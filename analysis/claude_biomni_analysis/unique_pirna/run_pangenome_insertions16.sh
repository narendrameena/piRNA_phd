#!/bin/bash
# 16-strain private-insertion extraction from the pangenome VCF (-> pangenome_te/ins16/{X}.private_insertions.fasta).
#SBATCH --job-name=panins16
#SBATCH --partition=2204,1804,2004,NXFL
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te/extract16.log
set -euo pipefail
BCF=/mnt/home3/miska/nm667/miniconda3/envs/rna2prot_snv/bin/bcftools
PY=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/python
V=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/pangenome/output/mouse_17strain_pangenome.vcf.gz
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
# all 16 strains; -a trims to carried alleles, --min-ac 1 keeps variant sites, -V snps drops SNP-only
"$BCF" view -s 129S1_SvImJ,AKR_J,A_J,BALB_cJ,C3H_HeJ,C57BL_6NJ,CAST_EiJ,CBA_J,DBA_2J,FVB_NJ,LP_J,NOD_ShiLtJ,NZO_HlLtJ,PWK_PhJ,SPRET_EiJ,WSB_EiJ -a --min-ac 1 -V snps "$V" \
 | "$BCF" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
 | "$PY" "$U/parse_insertions16.py"
echo "[done] $(date)"
