#!/usr/bin/env bash
# Re-render the full PICB-locus catalogue (20 loci × main/multi/single) with the updated Panel A (strand height-split)
# + Panel B (TE-DNA/TE-piRNA/AS→TE header, family×strand top sub-bar, antisense-by-family + silencing bottom sub-bar).
PY=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
n=0
while IFS=$'\t' read -r C S E L SLUG; do
  [ -z "$C" ] && continue
  n=$((n+1)); echo "=== [$n] $SLUG ($C:$S-$E) ==="
  $PY make_pav_locus.py        "$C" "$S" "$E" "$L" _ "Fig_pav_locus_${SLUG}"       2>&1 | grep -iE "wrote|Error|Traceback" | sed 's/^/  [main]   /'
  $PY make_pav_locus_multi.py  "$C" "$S" "$E" "$L" _ "Fig_pav_locus_${SLUG}_multi" 2>&1 | grep -iE "wrote|Error|Traceback" | sed 's/^/  [multi]  /'
  $PY make_pav_locus_single.py "$C" "$S" "$E" "$L" "Fig_pav_locus_${SLUG}_single"  2>&1 | grep -iE "wrote|Error|Traceback" | sed 's/^/  [single] /'
done < /tmp/allloci.tsv
echo CATALOGUE_DONE
