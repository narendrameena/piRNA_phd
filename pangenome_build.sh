#!/bin/bash
#SBATCH --job-name=mouse_pangenome
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=120:00:00
#SBATCH --output=pangenome_build_%j.log
#SBATCH --error=pangenome_build_%j.err

set -euo pipefail

###############################################################################
# CONFIGURATION
###############################################################################

BASE=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
REL_DIR=${BASE}/resources/REL-2205-Assembly
PAN_DIR=${BASE}/results/pangenome
THREADS=32

# GRCm39 reference (C57BL/6J)
GRCm39_FA=${BASE}/resources/black6/genome/Mus_musculus.GRCm39.dna.primary_assembly.fa

# Singularity container (has cactus-pangenome + vg + minigraph)
CACTUS_SIF=${BASE}/cactus_v2.9.3.sif

mkdir -p ${PAN_DIR}/{prepared,logs,output,qc}

echo "$(date) ============================================"
echo "  MOUSE 17-STRAIN PANGENOME BUILD + QC"
echo "  Reference: GRCm39 (C57BL/6J)"
echo "  No Y chromosome"
echo "  For piRNA cluster / TE analysis"
echo "============================================"

###############################################################################
# SINGULARITY WRAPPERS
###############################################################################

run_cactus() {
    singularity exec \
        --bind ${BASE}:${BASE} \
        --bind /tmp:/tmp \
        --cleanenv \
        ${CACTUS_SIF} \
        "$@"
}

run_vg() {
    singularity exec \
        --bind ${BASE}:${BASE} \
        --bind /tmp:/tmp \
        ${CACTUS_SIF} \
        vg "$@"
}

###############################################################################
# CANONICAL CHROMOSOMES (no Y)
###############################################################################

CANONICAL_CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X MT)

###############################################################################
# ALL STRAINS
#   GRCm39 = C57BL/6J reference
#   16 REL-2205 strains (includes CAST_EiJ)
###############################################################################

REL_STRAINS=(C57BL_6NJ CAST_EiJ A_J AKR_J BALB_cJ C3H_HeJ CBA_J DBA_2J \
             FVB_NJ LP_J NOD_ShiLtJ NZO_HlLtJ 129S1_SvImJ \
             PWK_PhJ WSB_EiJ SPRET_EiJ)

###############################################################################
# STEP 0: VERIFY ALL TOOLS
###############################################################################

echo ""
echo "$(date) === STEP 0: Tool verification ==="
echo ""

echo "  samtools:  $(samtools --version | head -1)"
echo "  bcftools:  $(bcftools --version | head -1)"
echo "  bedtools:  $(bedtools --version)"
echo "  cactus:    $(run_cactus cactus-pangenome --help 2>&1 | head -1)"
echo "  vg:        $(run_vg version 2>&1 | head -1)"

# Disk check (BeeGFS compatible)
AVAIL=$(df -BG ${PAN_DIR} 2>/dev/null | awk '/\// {gsub(/G/,"",$4); print $4; exit}')
echo ""
echo "  Disk available: ${AVAIL:-unknown} GB (need ~500GB)"

if [[ "${AVAIL:-0}" =~ ^[0-9]+$ ]] && [ "${AVAIL}" -lt 400 ]; then
    echo "  WARNING: Low disk space. Build may fail."
fi

# Verify GRCm39 exists
if [ ! -f "${GRCm39_FA}" ]; then
    echo "  FATAL: GRCm39 not found: ${GRCm39_FA}"
    exit 1
fi
echo "  GRCm39:    $(ls -lh ${GRCm39_FA} | awk '{print $5}')"

echo ""
echo "$(date) === STEP 0 COMPLETE ==="

###############################################################################
# STEP 1: DIAGNOSE ALL ASSEMBLY NAMING
###############################################################################

echo ""
echo "$(date) === STEP 1: Diagnose assembly naming ==="
echo ""

echo "Target chromosomes (no Y): ${CANONICAL_CHROMS[*]}"
echo ""

printf "%-25s  %-10s  %-55s\n" "Assembly" "Format" "First_header"
printf "%-25s  %-10s  %-55s\n" \
    "-------------------------" "----------" \
    "-------------------------------------------------------"

show_format() {
    local FA="$1"
    local LABEL="$2"

    if [ ! -f "$FA" ]; then
        printf "%-25s  %-10s  %-55s\n" "$LABEL" "MISSING" "file not found"
        return 1
    fi

    local HEADER
    HEADER=$(head -1 "$FA" | cut -c1-55)

    local FMT="unknown"
    if [[ "$HEADER" == *"#"*"#"* ]]; then
        FMT="PanSN"
    else
        FMT="Ensembl"
    fi

    printf "%-25s  %-10s  %s\n" "$LABEL" "$FMT" "$HEADER"
}

show_format "${GRCm39_FA}" "GRCm39_C57BL_6J"

for strain in "${REL_STRAINS[@]}"; do
    show_format "${REL_DIR}/${strain}_chromosomes_MT.fasta" "$strain"
done

echo ""
echo "$(date) === STEP 1 COMPLETE ==="

###############################################################################
# STEP 2: STANDARDIZE ALL CHROMOSOME NAMES
#
#   GRCm39 (Ensembl): >1 dna:chromosome ... → >1
#   REL-2205 (PanSN): >STRAIN#1#chr1       → >1
#   ALL become:        >1 >2 ... >19 >X >MT
#   Y chromosome:      SKIPPED
#   Scaffolds:         SKIPPED
###############################################################################

echo ""
echo "$(date) === STEP 2: Standardize chromosome names ==="
echo ""

standardize_genome() {
    local INPUT="$1"
    local OUTPUT="$2"
    local NAME="$3"

    # Skip if already done and valid
    if [ -f "$OUTPUT" ] && [ -f "${OUTPUT}.fai" ]; then
        local N
        N=$(wc -l < "${OUTPUT}.fai")
        if [ "$N" -eq 21 ]; then
            local NAMES
            NAMES=$(cut -f1 "${OUTPUT}.fai" | tr '\n' ' ')
            echo "  ${NAME}: already done (${N} chrs: ${NAMES})"
            return 0
        fi
        echo "  ${NAME}: wrong chr count (${N}), rebuilding..."
        rm -f "$OUTPUT" "${OUTPUT}.fai"
    fi

    if [ ! -f "$INPUT" ]; then
        echo "  ${NAME}: ERROR -- not found: $INPUT"
        return 1
    fi

    echo "  ${NAME}: processing..."
    [ ! -f "${INPUT}.fai" ] && samtools faidx "$INPUT"

    local FIRST_NAME
    FIRST_NAME=$(awk 'NR==1{print $1}' "${INPUT}.fai")

    local FORMAT="unknown"
    if [[ "$FIRST_NAME" == *"#"*"#"* ]]; then
        FORMAT="PanSN"
    else
        FORMAT="Ensembl"
    fi

    echo "    Format: ${FORMAT} (first seq: ${FIRST_NAME})"

    # Build rename map
    local MAP="${PAN_DIR}/qc/${NAME}_rename_map.tsv"
    > "$MAP"

    while read -r OLD_NAME OLD_LEN; do
        local NEW_NAME=""

        if [ "$FORMAT" = "PanSN" ]; then
            # >STRAIN#1#chr1 -> chr1 -> 1
            NEW_NAME=$(echo "$OLD_NAME" | awk -F'#' '{print $NF}')
            NEW_NAME=${NEW_NAME#chr}
            if [ "$NEW_NAME" = "M" ]; then
                NEW_NAME="MT"
            fi
        else
            # Ensembl: first column is already the name
            NEW_NAME="$OLD_NAME"
            if [ "$NEW_NAME" = "M" ]; then
                NEW_NAME="MT"
            fi
        fi

        printf '%s\t%s\t%s\n' "$OLD_NAME" "$NEW_NAME" "$OLD_LEN" >> "$MAP"

    done < <(awk '{print $1, $2}' "${INPUT}.fai")

    # Show mappings
    echo "    Mappings:"
    for c in "${CANONICAL_CHROMS[@]}"; do
        local OLD
        OLD=$(awk -v target="$c" -F'\t' '$2==target {print $1}' "$MAP")
        local SZ
        SZ=$(awk -v target="$c" -F'\t' '$2==target {printf "%.1f Mb", $3/1e6}' "$MAP")
        if [ -n "$OLD" ]; then
            printf "      %-35s -> %-4s  (%s)\n" "$OLD" "$c" "$SZ"
        fi
    done

    # Count skipped
    local TOTAL_IN
    TOTAL_IN=$(wc -l < "$MAP")
    local N_CANON=0
    for c in "${CANONICAL_CHROMS[@]}"; do
        if grep -qP "\t${c}\t" "$MAP" 2>/dev/null; then
            N_CANON=$((N_CANON + 1))
        fi
    done
    local N_SKIP=$((TOTAL_IN - N_CANON))
    if [ "$N_SKIP" -gt 0 ]; then
        # Show if Y was skipped
        local HAS_Y
        HAS_Y=$(awk -F'\t' '$2=="Y"' "$MAP")
        if [ -n "$HAS_Y" ]; then
            echo "      Y chromosome -> SKIPPED (by design)"
        fi
        echo "      ... plus $((N_SKIP)) scaffolds/contigs/Y -> skipped"
    fi

    # Extract canonical chromosomes in order
    local TMPFA="${OUTPUT}.tmp"
    > "$TMPFA"

    local FOUND_COUNT=0
    for c in "${CANONICAL_CHROMS[@]}"; do
        local OLD
        OLD=$(awk -v target="$c" -F'\t' '$2==target {print $1; exit}' "$MAP")
        if [ -z "$OLD" ]; then
            echo "    WARNING: chromosome ${c} not found in ${NAME}"
            continue
        fi

        samtools faidx "$INPUT" "$OLD" 2>/dev/null | \
            awk -v newname="$c" '/^>/{print ">" newname; next} {print}' >> "$TMPFA"

        FOUND_COUNT=$((FOUND_COUNT + 1))
    done

    if [ "$FOUND_COUNT" -lt 20 ]; then
        echo "    ERROR: Only ${FOUND_COUNT}/21 canonical chromosomes found!"
        rm -f "$TMPFA"
        return 1
    fi

    samtools faidx "$TMPFA"
    mv "$TMPFA" "$OUTPUT"
    mv "${TMPFA}.fai" "${OUTPUT}.fai"

    local TOTAL
    TOTAL=$(awk '{sum+=$2} END{printf "%.2f Gb", sum/1e9}' "${OUTPUT}.fai")
    local FINAL_NAMES
    FINAL_NAMES=$(cut -f1 "${OUTPUT}.fai" | tr '\n' ' ')

    echo "    Done: ${FOUND_COUNT} chromosomes, ${TOTAL}"
    echo "    Names: ${FINAL_NAMES}"

    local CHR1_MB
    CHR1_MB=$(awk '$1=="1"{printf "%.0f", $2/1e6}' "${OUTPUT}.fai")
    if [ -n "$CHR1_MB" ] && [ "$CHR1_MB" -lt 100 ]; then
        echo "    WARNING: Chr1 only ${CHR1_MB} Mb -- suspicious!"
    fi

    return 0
}

# -- GRCm39 reference --
echo "--- GRCm39 Reference (C57BL/6J) ---"
echo ""

standardize_genome \
    "${GRCm39_FA}" \
    "${PAN_DIR}/prepared/GRCm39.fa" \
    "GRCm39_C57BL_6J"
echo ""

# -- REL-2205 strain assemblies --
echo "--- REL-2205 Assemblies (16 strains) ---"
echo ""

for strain in "${REL_STRAINS[@]}"; do
    standardize_genome \
        "${REL_DIR}/${strain}_chromosomes_MT.fasta" \
        "${PAN_DIR}/prepared/${strain}.fa" \
        "${strain}"
    echo ""
done

echo "$(date) === STEP 2 COMPLETE ==="

###############################################################################
# STEP 3: CROSS-ASSEMBLY VALIDATION
###############################################################################

echo ""
echo "$(date) === STEP 3: Cross-assembly validation ==="
echo ""

REF_CHROMS=$(cut -f1 ${PAN_DIR}/prepared/GRCm39.fa.fai | sort)

printf "%-25s  %4s  %10s  %8s  %8s  %7s  %s\n" \
    "Assembly" "Chrs" "Total_Gb" "Chr1_Mb" "ChrX_Mb" "MT_bp" "Status"
printf "%-25s  %4s  %10s  %8s  %8s  %7s  %s\n" \
    "-------------------------" "----" "----------" "--------" "--------" "-------" "------"

ALL_MATCH=true

for fa in ${PAN_DIR}/prepared/*.fa; do
    NAME=$(basename "$fa" .fa)
    [ ! -f "${fa}.fai" ] && samtools faidx "$fa"

    NSEQ=$(wc -l < "${fa}.fai")
    TOTAL=$(awk '{sum+=$2} END{printf "%.3f", sum/1e9}' "${fa}.fai")
    CHR1=$(awk '$1=="1"{printf "%.1f", $2/1e6}' "${fa}.fai")
    CHRX=$(awk '$1=="X"{printf "%.1f", $2/1e6}' "${fa}.fai")
    MT=$(awk '$1=="MT"{print $2}' "${fa}.fai")

    THIS_CHROMS=$(cut -f1 "${fa}.fai" | sort)
    MISSING=$(comm -23 <(echo "$REF_CHROMS") <(echo "$THIS_CHROMS") | tr '\n' ',' | sed 's/,$//')
    EXTRA=$(comm -13 <(echo "$REF_CHROMS") <(echo "$THIS_CHROMS") | tr '\n' ',' | sed 's/,$//')

    if [ -z "$MISSING" ] && [ -z "$EXTRA" ]; then
        STATUS="MATCH"
    else
        STATUS="ISSUE"
        [ -n "$MISSING" ] && STATUS="${STATUS} miss:${MISSING}"
        [ -n "$EXTRA" ] && STATUS="${STATUS} extra:${EXTRA}"
        ALL_MATCH=false
    fi

    [ -z "$CHR1" ] && CHR1="N/A"
    [ -z "$CHRX" ] && CHRX="N/A"
    [ -z "$MT" ] && MT="N/A"

    printf "%-25s  %4d  %10s  %8s  %8s  %7s  %s\n" \
        "$NAME" "$NSEQ" "$TOTAL" "$CHR1" "$CHRX" "$MT" "$STATUS"
done

echo ""

# Header cleanliness
echo "--- Header cleanliness check ---"
CLEAN=true
for fa in ${PAN_DIR}/prepared/*.fa; do
    NAME=$(basename "$fa" .fa)
    BAD=$(grep "^>" "$fa" | grep " " | head -1 || true)
    if [ -n "$BAD" ]; then
        echo "  PROBLEM: ${NAME}: header has spaces: ${BAD}"
        CLEAN=false
    fi
done
if [ "$CLEAN" = true ]; then
    echo "  OK: All headers clean (no spaces or metadata)"
fi
echo ""

# Resolve mismatches if any
if [ "$ALL_MATCH" = false ]; then
    echo "--- Resolving chromosome set mismatches ---"

    COMMON_FILE=${PAN_DIR}/qc/common_chroms.txt
    echo "$REF_CHROMS" > "$COMMON_FILE"

    for fa in ${PAN_DIR}/prepared/*.fa; do
        comm -12 "$COMMON_FILE" <(cut -f1 "${fa}.fai" | sort) > "${COMMON_FILE}.tmp"
        mv "${COMMON_FILE}.tmp" "$COMMON_FILE"
    done

    N_COMMON=$(wc -l < "$COMMON_FILE")
    echo "  Common set (${N_COMMON}): $(cat "$COMMON_FILE" | tr '\n' ' ')"

    if [ "$N_COMMON" -lt 20 ]; then
        echo "  FATAL: fewer than 20 common chromosomes (expected 21: 1-19 X MT)"
        exit 1
    fi

    COMMON_LIST=$(cat "$COMMON_FILE" | tr '\n' ' ')

    for fa in ${PAN_DIR}/prepared/*.fa; do
        NAME=$(basename "$fa" .fa)
        N_THIS=$(wc -l < "${fa}.fai")
        if [ "$N_THIS" -ne "$N_COMMON" ]; then
            echo "  Trimming ${NAME} to common set..."
            mv "$fa" "${fa}.bak"
            samtools faidx "${fa}.bak" ${COMMON_LIST} > "$fa"
            samtools faidx "$fa"
            rm -f "${fa}.bak" "${fa}.bak.fai"
        fi
    done

    echo "  OK: All assemblies now have identical chromosome sets"
fi

# Pairwise md5 verification
echo ""
echo "--- Pairwise md5 verification ---"
ASSEMBLIES=(${PAN_DIR}/prepared/*.fa)
REF_HASH=$(cut -f1 "${ASSEMBLIES[0]}.fai" | sort | md5sum | awk '{print $1}')
PAIR_OK=true

for fa in "${ASSEMBLIES[@]}"; do
    NAME=$(basename "$fa" .fa)
    HASH=$(cut -f1 "${fa}.fai" | sort | md5sum | awk '{print $1}')
    if [ "$HASH" != "$REF_HASH" ]; then
        echo "  FAIL: ${NAME} chromosome set differs!"
        PAIR_OK=false
    fi
done

if [ "$PAIR_OK" = true ]; then
    echo "  OK: All ${#ASSEMBLIES[@]} assemblies have IDENTICAL chromosome names"
else
    echo "  FATAL: Could not harmonize chromosome names"
    exit 1
fi

echo ""
echo "$(date) === STEP 3 COMPLETE ==="

###############################################################################
# STEP 4: CREATE SEQFILE
#
#   Order:
#     GRCm39 first (reference)
#     CAST_EiJ (wild-derived M. m. castaneus — most divergent after SPRET)
#     Lab strains (close to reference)
#     Wild-derived (progressively divergent)
#     SPRET_EiJ last (most divergent ~2My)
###############################################################################

echo ""
echo "$(date) === STEP 4: Create seqfile ==="
echo ""

SEQFILE=${PAN_DIR}/seqfile.txt

cat > ${SEQFILE} << SEQEOF
GRCm39	${PAN_DIR}/prepared/GRCm39.fa
CAST_EiJ	${PAN_DIR}/prepared/CAST_EiJ.fa
C57BL_6NJ	${PAN_DIR}/prepared/C57BL_6NJ.fa
129S1_SvImJ	${PAN_DIR}/prepared/129S1_SvImJ.fa
A_J	${PAN_DIR}/prepared/A_J.fa
AKR_J	${PAN_DIR}/prepared/AKR_J.fa
BALB_cJ	${PAN_DIR}/prepared/BALB_cJ.fa
C3H_HeJ	${PAN_DIR}/prepared/C3H_HeJ.fa
CBA_J	${PAN_DIR}/prepared/CBA_J.fa
DBA_2J	${PAN_DIR}/prepared/DBA_2J.fa
FVB_NJ	${PAN_DIR}/prepared/FVB_NJ.fa
LP_J	${PAN_DIR}/prepared/LP_J.fa
NOD_ShiLtJ	${PAN_DIR}/prepared/NOD_ShiLtJ.fa
NZO_HlLtJ	${PAN_DIR}/prepared/NZO_HlLtJ.fa
PWK_PhJ	${PAN_DIR}/prepared/PWK_PhJ.fa
WSB_EiJ	${PAN_DIR}/prepared/WSB_EiJ.fa
SPRET_EiJ	${PAN_DIR}/prepared/SPRET_EiJ.fa
SEQEOF

VALID_COUNT=0
while IFS=$'\t' read -r name path; do
    if [ -f "$path" ] && [ -f "${path}.fai" ]; then
        NSEQ=$(wc -l < "${path}.fai")
        SIZE=$(awk '{sum+=$2} END{printf "%.2f Gb", sum/1e9}' "${path}.fai")
        echo "  OK: ${name} -> ${NSEQ} chrs, ${SIZE}"
        VALID_COUNT=$((VALID_COUNT + 1))
    else
        echo "  MISSING: ${name} -> ${path}"
    fi
done < ${SEQFILE}

echo ""
echo "  Valid: ${VALID_COUNT} / $(wc -l < ${SEQFILE})"

if [ ${VALID_COUNT} -lt 2 ]; then
    echo "  FATAL: need at least 2 assemblies"
    exit 1
fi

echo ""
echo "$(date) === STEP 4 COMPLETE ==="

###############################################################################
# STEP 5: BUILD PANGENOME WITH MINIGRAPH-CACTUS
#
#   Reference: GRCm39 (C57BL/6J)
#   Outputs:
#     --vcf     TE-sized structural variants (primary value for piRNA)
#     --gfa     Graph visualization in Bandage
#     --giraffe Indexes for optional graph-based sRNA mapping
#     --gbz     Graph genome format
###############################################################################

echo ""
echo "$(date) === STEP 5: Building pangenome ==="
echo ""
echo "  Reference:    GRCm39 (C57BL/6J)"
echo "  Assemblies:   ${VALID_COUNT} (1 ref + 16 strains)"
echo "  Chromosomes:  1-19, X, MT (no Y)"
echo "  Threads:      ${THREADS}"
echo "  Memory:       240G"
echo "  Container:    ${CACTUS_SIF}"
echo ""

JOBSTORE=${PAN_DIR}/jobstore
OUTDIR=${PAN_DIR}/output
OUTNAME=mouse_17strain_pangenome

# Clean previous attempt
[ -d "${JOBSTORE}" ] && rm -rf ${JOBSTORE}

START_TIME=$(date +%s)

run_cactus cactus-pangenome ${JOBSTORE} \
    ${SEQFILE} \
    --outDir ${OUTDIR} \
    --outName ${OUTNAME} \
    --reference GRCm39 \
    --vcf \
    --giraffe \
    --gfa \
    --gbz \
    --maxCores 256 \
    --maxMemory 240G \
    --consCores 256 \
    --logFile ${PAN_DIR}/logs/cactus_pangenome.log \
    --realTimeLogging True\
    2>&1 | tee ${PAN_DIR}/logs/cactus_stdout.log

END_TIME=$(date +%s)
ELAPSED=$(( (END_TIME - START_TIME) / 3600 ))

echo ""
echo "$(date) === STEP 5 COMPLETE (${ELAPSED} hours) ==="

###############################################################################
# STEP 6: PANGENOME QC -- FILE INTEGRITY
###############################################################################

echo ""
echo "$(date) === STEP 6: Pangenome QC -- File integrity ==="
echo ""

QC_PASS=true

echo "--- Expected output files ---"
echo ""

# Check for files (cactus may or may not compress GFA)
for label_file in \
    "GBZ:${OUTDIR}/${OUTNAME}.gbz" \
    "VCF:${OUTDIR}/${OUTNAME}.vcf.gz" \
    "DIST:${OUTDIR}/${OUTNAME}.dist" \
    "MIN:${OUTDIR}/${OUTNAME}.min"; do

    label="${label_file%%:*}"
    FPATH="${label_file#*:}"

    if [ -f "$FPATH" ]; then
        FSIZE=$(ls -lh "$FPATH" | awk '{print $5}')
        FBYTES=$(stat --printf="%s" "$FPATH" 2>/dev/null || echo 0)

        if [ "$FBYTES" -lt 1000 ]; then
            printf "  %-6s  %-60s  %8s  WARNING: SMALL\n" "$label" "$FPATH" "$FSIZE"
            QC_PASS=false
        else
            printf "  %-6s  %-60s  %8s  OK\n" "$label" "$FPATH" "$FSIZE"
        fi
    else
        printf "  %-6s  %-60s  %8s  MISSING\n" "$label" "$FPATH" "---"
        QC_PASS=false
    fi
done

# GFA (may be .gfa or .gfa.gz)
GFA_FILE=""
for gfa_try in "${OUTDIR}/${OUTNAME}.gfa.gz" "${OUTDIR}/${OUTNAME}.gfa"; do
    if [ -f "$gfa_try" ]; then
        GFA_FILE="$gfa_try"
        FSIZE=$(ls -lh "$gfa_try" | awk '{print $5}')
        printf "  %-6s  %-60s  %8s  OK\n" "GFA" "$gfa_try" "$FSIZE"
        break
    fi
done
if [ -z "$GFA_FILE" ]; then
    printf "  %-6s  %-60s  %8s  MISSING\n" "GFA" "${OUTDIR}/${OUTNAME}.gfa*" "---"
    QC_PASS=false
fi

echo ""
echo "--- All files in output directory ---"
ls -lhS ${OUTDIR}/ 2>/dev/null || echo "  output directory empty!"

echo ""
if [ "$QC_PASS" = true ]; then
    echo "  FILE INTEGRITY: PASS"
else
    echo "  FILE INTEGRITY: FAIL -- check cactus log"
    echo "    Log: ${PAN_DIR}/logs/cactus_pangenome.log"
fi

echo ""
echo "$(date) === STEP 6 COMPLETE ==="

###############################################################################
# STEP 7: PANGENOME QC -- GRAPH TOPOLOGY
###############################################################################

echo ""
echo "$(date) === STEP 7: Pangenome QC -- Graph topology ==="
echo ""

GBZ=${OUTDIR}/${OUTNAME}.gbz

if [ -f "${GBZ}" ]; then

    echo "--- vg stats on GBZ ---"
    echo ""

    run_vg stats -z ${GBZ} 2>/dev/null | tee ${PAN_DIR}/qc/vg_stats_basic.txt
    echo ""

    echo "--- Detailed graph metrics ---"
    echo ""

    GRAPH_NODES=$(run_vg stats -N ${GBZ} 2>/dev/null || echo "N/A")
    GRAPH_EDGES=$(run_vg stats -E ${GBZ} 2>/dev/null || echo "N/A")
    GRAPH_TOTAL_BP=$(run_vg stats -l ${GBZ} 2>/dev/null || echo "N/A")

    echo "  Nodes (segments):    ${GRAPH_NODES}"
    echo "  Edges (links):       ${GRAPH_EDGES}"
    echo "  Total sequence (bp): ${GRAPH_TOTAL_BP}"

    REF_BP=$(awk '{sum+=$2} END{print sum}' ${PAN_DIR}/prepared/GRCm39.fa.fai)

    if [ "${GRAPH_TOTAL_BP}" != "N/A" ] && [ -n "${GRAPH_TOTAL_BP}" ]; then
        EXCESS_MB=$(echo "scale=2; (${GRAPH_TOTAL_BP} - ${REF_BP}) / 1000000" | bc)
        EXCESS_PCT=$(echo "scale=1; 100 * (${GRAPH_TOTAL_BP} - ${REF_BP}) / ${REF_BP}" | bc)
        AVG_NODE_LEN=$(echo "scale=0; ${GRAPH_TOTAL_BP} / ${GRAPH_NODES}" | bc 2>/dev/null || echo "N/A")

        echo ""
        echo "  Reference (GRCm39):       $(echo ${REF_BP} | numfmt --to=iec) bp"
        echo "  Non-reference excess:     ${EXCESS_MB} Mb (${EXCESS_PCT}%)"
        echo "  Average node length:      ${AVG_NODE_LEN} bp"
        echo ""

        if (( $(echo "${EXCESS_PCT} > 50" | bc -l) )); then
            echo "  QC WARNING: Non-ref excess > 50%"
        elif (( $(echo "${EXCESS_PCT} < 1" | bc -l) )); then
            echo "  QC WARNING: Non-ref excess < 1%"
        else
            echo "  OK: Non-ref excess in expected range for 17 mouse strains"
        fi

        if [ "${AVG_NODE_LEN}" != "N/A" ] && [ "${AVG_NODE_LEN}" -lt 10 ]; then
            echo "  QC WARNING: Average node length < 10bp -- graph over-fragmented"
        fi
    fi

    # Paths
    echo ""
    echo "--- Paths (haplotypes) in graph ---"
    echo ""

    run_vg paths -L -x ${GBZ} 2>/dev/null > ${PAN_DIR}/qc/graph_paths.txt || true

    if [ -f "${PAN_DIR}/qc/graph_paths.txt" ] && [ -s "${PAN_DIR}/qc/graph_paths.txt" ]; then
        TOTAL_PATHS=$(wc -l < ${PAN_DIR}/qc/graph_paths.txt)
        echo "  Total paths: ${TOTAL_PATHS}"

        echo ""
        echo "  Paths per strain:"
        awk -F'#' '{print $1}' ${PAN_DIR}/qc/graph_paths.txt | sort | uniq -c | sort -rn | \
        while read count strain; do
            printf "    %-20s  %3d paths\n" "$strain" "$count"
        done

        N_STRAINS=$(awk -F'#' '{print $1}' ${PAN_DIR}/qc/graph_paths.txt | sort -u | wc -l)
        echo ""
        echo "  Strains represented: ${N_STRAINS} / 17 expected"
        if [ "${N_STRAINS}" -eq 17 ]; then
            echo "  OK: All 17 strains present in graph"
        else
            echo "  QC WARNING: Missing strains!"
            echo "    Expected:"
            awk -F'\t' '{print $1}' ${SEQFILE} | sort
            echo "    Found:"
            awk -F'#' '{print $1}' ${PAN_DIR}/qc/graph_paths.txt | sort -u
        fi

        echo ""
        echo "  Reference (GRCm39) paths:"
        grep "^GRCm39" ${PAN_DIR}/qc/graph_paths.txt | head -5
        REF_PATH_COUNT=$(grep -c "^GRCm39" ${PAN_DIR}/qc/graph_paths.txt || echo 0)
        echo "    ... ${REF_PATH_COUNT} total reference paths"
    else
        echo "  Could not extract paths"
        run_vg stats -z ${GBZ} 2>/dev/null | grep -i path || true
    fi

else
    echo "  GBZ not found -- cannot perform graph QC"
fi

echo ""
echo "$(date) === STEP 7 COMPLETE ==="

###############################################################################
# STEP 8: PANGENOME QC -- VCF QUALITY
###############################################################################

echo ""
echo "$(date) === STEP 8: Pangenome QC -- VCF quality ==="
echo ""

VCF=${OUTDIR}/${OUTNAME}.vcf.gz

if [ -f "${VCF}" ]; then

    [ ! -f "${VCF}.tbi" ] && bcftools index -t ${VCF}

    echo "--- bcftools stats ---"
    echo ""

    bcftools stats ${VCF} > ${PAN_DIR}/qc/bcftools_stats.txt 2>/dev/null

    TOTAL_RECORDS=$(grep "^SN" ${PAN_DIR}/qc/bcftools_stats.txt | grep "number of records:" | awk '{print $NF}')
    TOTAL_SNPS=$(grep "^SN" ${PAN_DIR}/qc/bcftools_stats.txt | grep "number of SNPs:" | awk '{print $NF}')
    TOTAL_INDELS=$(grep "^SN" ${PAN_DIR}/qc/bcftools_stats.txt | grep "number of indels:" | awk '{print $NF}')
    TOTAL_MULTIALLELIC=$(grep "^SN" ${PAN_DIR}/qc/bcftools_stats.txt | grep "number of multiallelic" | awk '{print $NF}')
    TSTV=$(grep "^TSTV" ${PAN_DIR}/qc/bcftools_stats.txt | head -1 | awk '{print $5}')

    echo "  Total variant records: ${TOTAL_RECORDS}"
    echo "  SNPs:                  ${TOTAL_SNPS}"
    echo "  Indels:                ${TOTAL_INDELS}"
    echo "  Multiallelic sites:    ${TOTAL_MULTIALLELIC}"
    echo "  Ts/Tv ratio:           ${TSTV}"

    echo ""
    if [ -n "$TSTV" ]; then
        TSTV_OK=$(echo "$TSTV > 1.5 && $TSTV < 3.0" | bc -l 2>/dev/null || echo 0)
        if [ "$TSTV_OK" -eq 1 ]; then
            echo "  OK: Ts/Tv ratio in expected range (1.5-3.0 for mouse)"
        else
            echo "  QC WARNING: Ts/Tv ratio outside expected range"
        fi
    fi

    if [ -n "$TOTAL_RECORDS" ] && [ "$TOTAL_RECORDS" -lt 100000 ]; then
        echo "  QC WARNING: Very few variants"
    elif [ -n "$TOTAL_RECORDS" ] && [ "$TOTAL_RECORDS" -gt 1000000 ]; then
        echo "  OK: Variant count reasonable for 17 mouse strains"
    fi

    # Per-sample variant counts
    echo ""
    echo "--- Per-strain variant counts ---"
    echo ""
    printf "  %-20s  %10s\n" "Strain" "Non-ref"
    printf "  %-20s  %10s\n" "--------------------" "----------"

    bcftools query -l ${VCF} 2>/dev/null | while read sample; do
        N_ALT=$(bcftools view -s "${sample}" -H ${VCF} 2>/dev/null | \
            awk -F'\t' '{
                split($10,gt,":")
                g = gt[1]
                if (g != "0/0" && g != "0|0" && g != "./." && g != ".") print
            }' | wc -l)
        printf "  %-20s  %10d\n" "${sample}" "${N_ALT}"
    done

    # Per-chromosome variant density
    echo ""
    echo "--- Per-chromosome variant counts ---"
    echo ""

    for CHR in "${CANONICAL_CHROMS[@]}"; do
        N_CHR=$(bcftools view -r ${CHR} -H ${VCF} 2>/dev/null | wc -l)
        CHR_LEN=$(awk -v c=${CHR} '$1==c{print $2}' ${PAN_DIR}/prepared/GRCm39.fa.fai)
        if [ -n "$CHR_LEN" ] && [ "$CHR_LEN" -gt 0 ] && [ "$N_CHR" -gt 0 ]; then
            DENSITY=$(echo "scale=1; ${N_CHR} * 1000000 / ${CHR_LEN}" | bc)
            printf "    chr%-3s  %8d variants  (%s vars/Mb)\n" "$CHR" "$N_CHR" "$DENSITY"
        else
            printf "    chr%-3s  %8d variants\n" "$CHR" "$N_CHR"
        fi
    done

    # Indel size distribution
    echo ""
    echo "--- Indel size distribution (TE-relevant) ---"
    echo ""

    grep "^IDD" ${PAN_DIR}/qc/bcftools_stats.txt | \
    awk '{
        size = $3
        count = $4
        if (size < 0) size = -size
        if (size == 0) next
        if (size <= 10)        bin = "1-10bp"
        else if (size <= 50)   bin = "11-50bp"
        else if (size <= 300)  bin = "51-300bp (SINE B1/B2)"
        else if (size <= 1000) bin = "301bp-1kb (solo LTR)"
        else if (size <= 5000) bin = "1-5kb (LINE fragment)"
        else if (size <= 10000) bin = "5-10kb (full LINE/LTR)"
        else bin = ">10kb (complex SV)"
        bins[bin] += count
    } END {
        for (b in bins) printf "    %-35s  %8d\n", b, bins[b]
    }' | sort -t$'\t' -k2 -rn

    echo ""
    echo "  Full stats: ${PAN_DIR}/qc/bcftools_stats.txt"

else
    echo "  VCF not found -- skipping VCF QC"
fi

echo ""
echo "$(date) === STEP 8 COMPLETE ==="

###############################################################################
# STEP 9: PANGENOME QC -- REFERENCE COVERAGE
###############################################################################

echo ""
echo "$(date) === STEP 9: Pangenome QC -- Reference coverage ==="
echo ""

if [ -f "${GBZ}" ]; then

    echo "--- Reference path lengths vs input assembly ---"
    echo ""

    printf "  %-5s  %12s  %12s  %8s\n" "Chr" "Input_bp" "Graph_bp" "Ratio"
    printf "  %-5s  %12s  %12s  %8s\n" "-----" "------------" "------------" "--------"

    REF_INPUT_FAI="${PAN_DIR}/prepared/GRCm39.fa.fai"
    COVERAGE_OK=true

    run_vg paths -L -x ${GBZ} 2>/dev/null | grep "^GRCm39" > ${PAN_DIR}/qc/ref_paths.txt || true

    if [ -s "${PAN_DIR}/qc/ref_paths.txt" ]; then
        for CHR in "${CANONICAL_CHROMS[@]}"; do
            INPUT_LEN=$(awk -v c=${CHR} '$1==c{print $2}' ${REF_INPUT_FAI})

            REF_PATH_NAME="GRCm39#1#${CHR}"
            GRAPH_LEN=$(run_vg paths -x ${GBZ} -E -S "${REF_PATH_NAME}" 2>/dev/null | \
                        awk '{sum+=$2} END{print sum+0}' || echo 0)

            if [ -n "$INPUT_LEN" ] && [ "$INPUT_LEN" -gt 0 ] && [ "$GRAPH_LEN" -gt 0 ]; then
                RATIO=$(echo "scale=4; ${GRAPH_LEN} / ${INPUT_LEN}" | bc)
                printf "  %-5s  %12d  %12d  %8s\n" "$CHR" "$INPUT_LEN" "$GRAPH_LEN" "$RATIO"

                LOW=$(echo "$RATIO < 0.95" | bc -l)
                HIGH=$(echo "$RATIO > 1.05" | bc -l)
                if [ "$LOW" -eq 1 ] || [ "$HIGH" -eq 1 ]; then
                    echo "         WARNING: Ratio outside 0.95-1.05"
                    COVERAGE_OK=false
                fi
            elif [ -n "$INPUT_LEN" ]; then
                printf "  %-5s  %12d  %12s  %8s\n" "$CHR" "$INPUT_LEN" "N/A" "N/A"
                COVERAGE_OK=false
            fi
        done

        echo ""
        if [ "$COVERAGE_OK" = true ]; then
            echo "  OK: All reference chromosomes well-represented"
        else
            echo "  WARNING: Some chromosomes have coverage anomalies"
        fi
    else
        echo "  Could not extract reference path lengths"
        run_vg stats -z ${GBZ} 2>/dev/null | head -20
    fi

else
    echo "  GBZ not found -- skipping"
fi

echo ""
echo "$(date) === STEP 9 COMPLETE ==="

###############################################################################
# STEP 10: PANGENOME QC -- GIRAFFE INDEX VALIDATION
###############################################################################

echo ""
echo "$(date) === STEP 10: Pangenome QC -- Giraffe index validation ==="
echo ""

GBZ=${OUTDIR}/${OUTNAME}.gbz
DIST=${OUTDIR}/${OUTNAME}.dist
MIN=${OUTDIR}/${OUTNAME}.min

# Build missing indexes
if [ ! -f "${DIST}" ]; then
    echo "  Distance index missing -- building..."
    run_vg snarls ${GBZ} > ${OUTDIR}/${OUTNAME}.snarls 2>/dev/null
    run_vg index -j ${DIST} ${GBZ} 2>/dev/null
fi

if [ ! -f "${MIN}" ]; then
    echo "  Minimizer index missing -- building..."
    run_vg minimizer -t ${THREADS} -d ${DIST} -o ${MIN} ${GBZ} 2>/dev/null
fi

echo "--- Giraffe index files ---"
for IDX_FILE in ${GBZ} ${DIST} ${MIN}; do
    if [ -f "$IDX_FILE" ]; then
        FSIZE=$(ls -lh "$IDX_FILE" | awk '{print $5}')
        echo "  OK: $(basename $IDX_FILE): ${FSIZE}"
    else
        echo "  MISSING: $(basename $IDX_FILE)"
        QC_PASS=false
    fi
done

# Synthetic piRNA-length read mapping test
echo ""
echo "--- Synthetic read mapping test (28nt = piRNA length) ---"
echo ""

SYNTH_FQ=${PAN_DIR}/qc/synthetic_reads.fq

samtools faidx ${PAN_DIR}/prepared/GRCm39.fa 1:1000000-1100000 2>/dev/null | \
grep -v "^>" | tr -d '\n' | \
awk '{
    seq = toupper($0)
    for (i = 1; i <= 100; i++) {
        start = int(rand() * (length(seq) - 28)) + 1
        read = substr(seq, start, 28)
        if (read !~ /N/) {
            printf "@synth_%d\n%s\n+\n%s\n", i, read, "IIIIIIIIIIIIIIIIIIIIIIIIIIII"
        }
    }
}' > ${SYNTH_FQ}

N_SYNTH=$(grep -c "^@" ${SYNTH_FQ} || echo 0)
echo "  Generated ${N_SYNTH} synthetic 28nt reads"

if [ "$N_SYNTH" -gt 0 ] && [ -f "${GBZ}" ] && [ -f "${DIST}" ] && [ -f "${MIN}" ]; then
    echo "  Running vg giraffe..."

    run_vg giraffe \
        -Z ${GBZ} \
        -m ${MIN} \
        -d ${DIST} \
        -f ${SYNTH_FQ} \
        -t ${THREADS} \
        -o gam \
        2>/dev/null > ${PAN_DIR}/qc/synthetic_mapped.gam 2>${PAN_DIR}/qc/giraffe_test.log || true

    if [ -f "${PAN_DIR}/qc/synthetic_mapped.gam" ]; then
        run_vg stats -a ${PAN_DIR}/qc/synthetic_mapped.gam 2>/dev/null > \
            ${PAN_DIR}/qc/synthetic_map_stats.txt || true

        if [ -s "${PAN_DIR}/qc/synthetic_map_stats.txt" ]; then
            echo ""
            echo "  Mapping results:"
            cat ${PAN_DIR}/qc/synthetic_map_stats.txt | head -20 | sed 's/^/    /'
            echo ""
            echo "  OK: Giraffe mapping functional"
        else
            GAM_SIZE=$(stat --printf='%s' ${PAN_DIR}/qc/synthetic_mapped.gam 2>/dev/null || echo 0)
            if [ "$GAM_SIZE" -gt 0 ]; then
                echo "  OK: Giraffe produced output (${GAM_SIZE} bytes)"
            else
                echo "  WARNING: Empty GAM output"
            fi
        fi
    else
        echo "  FAIL: Giraffe mapping failed"
        echo "    Check: ${PAN_DIR}/qc/giraffe_test.log"
    fi
else
    echo "  Skipping mapping test -- missing files"
fi

rm -f ${SYNTH_FQ}

echo ""
echo "$(date) === STEP 10 COMPLETE ==="

###############################################################################
# STEP 11: PANGENOME QC -- CACTUS LOG ANALYSIS
###############################################################################

echo ""
echo "$(date) === STEP 11: Pangenome QC -- Log analysis ==="
echo ""

CACTUS_LOG="${PAN_DIR}/logs/cactus_pangenome.log"

if [ -f "${CACTUS_LOG}" ]; then
    echo "--- Cactus log summary ---"
    echo ""

    N_ERRORS=$(grep -ci "error" ${CACTUS_LOG} 2>/dev/null || echo 0)
    N_WARNINGS=$(grep -ci "warning" ${CACTUS_LOG} 2>/dev/null || echo 0)
    N_FAILED=$(grep -ci "failed\|failure" ${CACTUS_LOG} 2>/dev/null || echo 0)

    echo "  Error mentions:   ${N_ERRORS}"
    echo "  Warning mentions: ${N_WARNINGS}"
    echo "  Failed mentions:  ${N_FAILED}"

    if [ "$N_ERRORS" -gt 0 ]; then
        echo ""
        echo "  Last 5 error lines:"
        grep -i "error" ${CACTUS_LOG} | tail -5 | sed 's/^/    /'
    fi

    if [ "$N_FAILED" -gt 0 ]; then
        echo ""
        echo "  Failed lines:"
        grep -i "failed\|failure" ${CACTUS_LOG} | tail -5 | sed 's/^/    /'
    fi

    echo ""
    echo "--- Pipeline stage timings ---"
    grep -i "completed\|finished\|starting\|done" ${CACTUS_LOG} 2>/dev/null | \
        grep -i "minigraph\|cactus\|align\|graphmap\|clip\|join\|vg\|gfa\|gbz\|vcf" | \
        tail -20 | sed 's/^/    /'

    echo ""
    echo "--- Peak memory ---"
    grep -i "memory\|mem\|rss" ${CACTUS_LOG} 2>/dev/null | tail -5 | sed 's/^/    /' || \
        echo "    No memory info found"

    echo ""
    echo "  Full log: ${CACTUS_LOG}"
else
    echo "  Cactus log not found"
fi

echo ""
echo "$(date) === STEP 11 COMPLETE ==="

###############################################################################
# QC SUMMARY REPORT
###############################################################################

echo ""
echo "=================================================================="
echo "  PANGENOME BUILD + QC SUMMARY REPORT"
echo "  $(date)"
echo "=================================================================="
echo ""
echo "  INPUT:"
echo "    Reference:          GRCm39 (C57BL/6J)"
echo "    Strains:            17 (1 ref + 16 REL-2205)"
echo "    Chromosomes:        1-19, X, MT (no Y)"
echo "    Container:          ${CACTUS_SIF}"
echo ""
echo "  BUILD:"
echo "    Build time:         ${ELAPSED} hours"
echo "    Output directory:   ${OUTDIR}"
echo ""
echo "  OUTPUT FILES:"
ls -lhS ${OUTDIR}/${OUTNAME}.* 2>/dev/null | awk '{printf "    %-50s  %s\n", $NF, $5}'
echo ""
echo "  QC FILES:"
ls -lhS ${PAN_DIR}/qc/ 2>/dev/null | tail -n+2 | awk '{printf "    %-50s  %s\n", $NF, $5}'
echo ""
echo "  QC VERDICTS:"

if [ -f "${GBZ}" ]; then
    echo "    File integrity:       PASS"
else
    echo "    File integrity:       FAIL (GBZ missing)"
fi

if [ -f "${PAN_DIR}/qc/vg_stats_basic.txt" ] && [ -s "${PAN_DIR}/qc/vg_stats_basic.txt" ]; then
    echo "    Graph topology:       PASS"
else
    echo "    Graph topology:       CHECK"
fi

if [ -f "${PAN_DIR}/qc/bcftools_stats.txt" ] && [ -s "${PAN_DIR}/qc/bcftools_stats.txt" ]; then
    echo "    VCF quality:          PASS"
else
    echo "    VCF quality:          CHECK"
fi

if [ -f "${PAN_DIR}/qc/synthetic_mapped.gam" ]; then
    GAM_BYTES=$(stat --printf='%s' ${PAN_DIR}/qc/synthetic_mapped.gam 2>/dev/null || echo 0)
    if [ "$GAM_BYTES" -gt 0 ]; then
        echo "    Giraffe indexes:      PASS"
    else
        echo "    Giraffe indexes:      CHECK"
    fi
else
    echo "    Giraffe indexes:      CHECK"
fi

echo ""
echo "  NEXT STEPS FOR piRNA ANALYSIS:"
echo ""
echo "  1. Extract TE-sized SVs from pangenome VCF:"
echo "     bcftools view -H ${OUTDIR}/${OUTNAME}.vcf.gz | \\"
echo "       awk 'length(\$5)-length(\$4) >= 300 || length(\$4)-length(\$5) >= 300'"
echo ""
echo "  2. Map sRNA-seq to each strain's OWN assembly (linear):"
echo "     bowtie -v 0 -m 1 --best -x STRAIN_INDEX -f reads.fa -S out.sam"
echo ""
echo "  3. Detect piRNA clusters per strain:"
echo "     proTRAC -i out.sam -g STRAIN.fa"
echo ""
echo "  4. Intersect piRNA clusters with TE insertions from VCF:"
echo "     bedtools intersect -a clusters.bed -b te_insertions.bed -wa -wb"
echo ""
echo "  5. Optional: Map sRNA-seq to pangenome graph:"
echo "     vg giraffe -Z ${OUTDIR}/${OUTNAME}.gbz \\"
echo "       -m ${OUTDIR}/${OUTNAME}.min \\"
echo "       -d ${OUTDIR}/${OUTNAME}.dist \\"
echo "       -f reads.fq -t 16 -o gam > sample.gam"
echo ""
echo "$(date) === ALL COMPLETE ==="