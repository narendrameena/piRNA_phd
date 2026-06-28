#!/bin/bash
# THEME 22 step 14 — reference-free CO-LOCATION of the 1,393 non-reference clusters. They don't lift to GRCm39, so
# lift each strain's non-ref clusters into EVERY OTHER strain (cactus HAL has all 16). A cluster that lifts to strain Y
# => Y carries the homologous insertion (SHARED); lifts to no other strain => strain-PRIVATE. Run via srun on a big-mem
# node, 8 lifts in parallel (halLiftover uses UDC mmap, low RAM/instance).
set -uo pipefail
B=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
D=$B/figures/analysis_figures/22_odgi_inject_cluster_pav/data
SIF=$B/cactus_v2.9.3.sif; HAL=$B/results/pangenome/output/mouse_17strain_pangenome.full.hal
mkdir -p $D/colo
STR=(129S1_SvImJ A_J AKR_J BALB_cJ C3H_HeJ C57BL_6NJ CAST_EiJ CBA_J DBA_2J FVB_NJ LP_J NOD_ShiLtJ NZO_HlLtJ PWK_PhJ SPRET_EiJ WSB_EiJ)
do_lift(){ local X=$1 Y=$2
  singularity exec --bind /mnt "$SIF" halLiftover "$HAL" "$X" "$D/nonref/$X.nonref.bed" "$Y" "$D/colo/${X}__${Y}.lifted.bed" 2>/dev/null
}
export -f do_lift; export SIF HAL D
for X in "${STR[@]}"; do for Y in "${STR[@]}"; do [ "$X" != "$Y" ] && echo "$X $Y"; done; done | xargs -P 8 -n 2 bash -c 'do_lift "$@"' _
echo "ALL_COLO_LIFTS_DONE  ($(ls $D/colo/*__*.lifted.bed 2>/dev/null | grep -v test | wc -l) pair files)"
