#!/bin/bash
#SBATCH --job-name=t22_flank
#SBATCH --array=0-15
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --partition=TEST,NXFL,2204,2004
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/22_odgi_inject_cluster_pav/logs/flank_%a.log
# THEME 22 step 13 — LIFT-ARTIFACT confounder. A non-ref call is halLiftover FAILURE, which could be genuine absence
# OR divergent-sequence alignment failure (the non-ref RATE tracks divergence rho=0.69). Test: lift the +/-5kb FLANKS
# of each non-ref cluster to GRCm39. If both flanks lift to the SAME GRCm39 chrom within 50kb => the locus EXISTS in
# GRCm39 and the cluster is a genuine INSERTION between homologous flanks (NOT a divergent-block artifact).
set -euo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
CP=$B/analysis/claude_biomni_analysis/unique_pirna/cluster_pav
T22=$B/figures/analysis_figures/22_odgi_inject_cluster_pav; D=$T22/data
SIF=$B/cactus_v2.9.3.sif; HAL=$B/results/pangenome/output/mouse_17strain_pangenome.full.hal
STR=(129S1_SvImJ A_J AKR_J BALB_cJ C3H_HeJ C57BL_6NJ CAST_EiJ CBA_J DBA_2J FVB_NJ LP_J NOD_ShiLtJ NZO_HlLtJ PWK_PhJ SPRET_EiJ WSB_EiJ)
X=${STR[$SLURM_ARRAY_TASK_ID]}
awk 'BEGIN{OFS="\t"}{L=$2-5000; if(L<0)L=0; print $1,L,$2,$4"|L"; print $1,$3,$3+5000,$4"|R"}' $D/nonref/$X.nonref.bed > $D/nonref/$X.flanks.bed
singularity exec --bind /mnt "$SIF" halLiftover "$HAL" "$X" "$D/nonref/$X.flanks.bed" GRCm39 "$D/nonref/$X.flanks.lifted.bed" 2>/dev/null
/usr/bin/python3 - "$X" "$D" <<'PY'
import sys,collections
X,D=sys.argv[1],sys.argv[2]
pos=collections.defaultdict(lambda:collections.defaultdict(list))
for l in open(f"{D}/nonref/{X}.flanks.lifted.bed"):
    f=l.rstrip("\n").split("\t"); cid,side=f[3].rsplit("|",1); pos[cid][side].append((f[0],(int(f[1])+int(f[2]))//2))
tot=sum(1 for _ in open(f"{D}/nonref/{X}.nonref.bed")); both=clean=0
for cid,d in pos.items():
    if "L" in d and "R" in d:
        both+=1
        if any(lc==rc and abs(lm-rm)<50000 for lc,lm in d["L"] for rc,rm in d["R"]): clean+=1
open(f"{D}/nonref/{X}.flank_summary","w").write(f"{X}\t{tot}\t{both}\t{clean}\t{100*clean/max(tot,1):.0f}\n")
print(f"{X}: nonref={tot} both-flanks-lift={both} clean-insertion={clean} ({100*clean/max(tot,1):.0f}% genuine insertions)")
PY
rm -f $D/nonref/$X.flanks.bed $D/nonref/$X.flanks.lifted.bed
