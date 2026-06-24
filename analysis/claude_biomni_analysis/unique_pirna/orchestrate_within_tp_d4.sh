#!/bin/bash
# Orchestrate the WITHIN-TIMEPOINT re-derivation WITH within-tp D4:
#   1. wait for all 48 per-tp pools (build_pools16_pertp.py, PID arg)
#   2. submit the within-tp D4 SNP-variant array (run_step416_pertp.sh)
#   3. wait for the array -> 16 step4_16/{X}.snp_withintp.csv.gz
#   4. aggregate -> unique16/snp_variant_refinement_withintp.csv
#   5. within_tp_rederive.py -> within-tp final_classified_clean_2read.csv.gz (cross-tp backed up)
set -uo pipefail
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
PY=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
BUILDER_PID=${1:-0}
LOG=$U/within_tp_d4_orchestrator.log
exec >>"$LOG" 2>&1
echo "================ orchestrator start $(date) (builder pid=$BUILDER_PID) ================"

# 1. wait for all 48 per-tp pools
echo "[1] waiting for 48 per-tp pools ..."
while true; do
  n=$(ls "$U"/pools_pertp/*.pool.txt.gz 2>/dev/null | wc -l)
  alive=$(ps -p "$BUILDER_PID" -o pid= 2>/dev/null | wc -l)
  echo "  pools=$n/48 builder_alive=$alive $(date)"
  [ "$n" -ge 48 ] && break
  if [ "$BUILDER_PID" != "0" ] && [ "$alive" -eq 0 ] && [ "$n" -lt 48 ]; then
    echo "ERROR: builder pid=$BUILDER_PID gone at $n/48 pools -> aborting"; exit 1
  fi
  sleep 120
done
echo "[1] all 48 per-tp pools present $(date)"

# 2. submit the within-tp D4 array
echo "[2] submitting within-tp D4 array (run_step416_pertp.sh) ..."
JID=$(sbatch --parsable "$U/run_step416_pertp.sh") || { echo "ERROR sbatch failed"; exit 1; }
echo "  job=$JID"

# 3. wait for the array to finish
echo "[3] waiting for array $JID ..."
while squeue -j "$JID" -h 2>/dev/null | grep -q .; do sleep 60; done
sleep 10
nout=$(ls "$U"/step4_16/*.snp_withintp.csv.gz 2>/dev/null | wc -l)
echo "[3] array $JID done $(date); outputs=$nout/16"
if [ "$nout" -lt 16 ]; then
  echo "ERROR: only $nout/16 within-tp snp outputs -> check step4_16/snp_withintp_*.log"; exit 1
fi

# 4. aggregate within-tp snp set
echo "[4] aggregating within-tp snp set ..."
"$PY" "$U/aggregate_within_tp_snp.py" || { echo "ERROR aggregate"; exit 1; }

# 5. within-tp re-derivation (with within-tp D4)
echo "[5] within-tp re-derivation (within-tp D4) ..."
"$PY" "$U/within_tp_rederive.py" || { echo "ERROR within_tp_rederive"; exit 1; }

echo "================ orchestrator DONE $(date) ================"
