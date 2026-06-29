#!/usr/bin/env bash
# Render the CREATION-only source-locus figures (confound-fix §2a) — true new strain-private loci (clusters_at
# breadth<=3, read-verified), NOT the mislabeled top-FPM propagation set (Fig_srcmulti_*). Uses the precomputed
# creation projection (build_source_projection.py source_loci_master_creation.tsv) so make_source_pav_multi skips
# per-locus HAL liftover. Output: Fig_srccreation_{strain}_{TE}_chr{N}_{start}.{png,pdf,svg}
set -uo pipefail
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
PY=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
PROJ=$U/pangenome_te/source_projection_creation.tsv
cd "$U"
n=0
# genome-PROBED subset: clusters_at-creation AND genome-PAV breadth<=3 (true insertion-creation, absent in >=13);
# excludes loci that are genome-conserved (flank-false / propagation) despite restricted expression.
tail -n +2 pangenome_te/source_loci_master_creation_private.tsv | while IFS=$'\t' read -r carrier chrom_own start end te strand fpm cls gbreadth; do
  [ -z "$carrier" ] && continue
  N=${chrom_own##*chr}; te_san=${te//\//-}
  out="Fig_srccreation_${carrier}_${te_san}_chr${N}_${start}"
  n=$((n+1))
  echo "[$n] $out (FPM=$fpm)"
  $PY make_source_pav_multi.py "$carrier" "$chrom_own" "$start" "$end" "$te" "$strand" "$out" "$PROJ" \
    2>&1 | grep -E "carrier=|wrote|Error|Traceback|Exception" | head -4
done
echo "DONE rendering creation figures"
