PY=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
n=0
tail -n +2 pangenome_te/divergence_loci_final.tsv | while IFS=$'\t' read g39c g39s g39e breadth carrier maxfpm expressing genomen; do
  [ -z "$g39c" ] && continue
  n=$((n+1))
  slug="chr${g39c}_$((g39s/1000000))p$(((g39s%1000000)/1000))_${carrier}"
  lab="Divergence locus chr${g39c}:${g39s} — genome-conserved (${genomen}/16), strain-restricted expression (${breadth}/16)"
  $PY make_pav_locus_multi.py "$g39c" "$g39s" "$g39e" "$lab" _ "Fig_divergence_pav_${slug}" 2>&1 | grep -iE "present=|wrote|Error|Traceback"
  echo "=== [$n] $slug (carrier $carrier, breadth $breadth, genome $genomen) ==="
done
echo DIVERGENCE_V2_DONE
