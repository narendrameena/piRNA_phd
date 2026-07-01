cat > ~/scratch/inProgress/mice_PiRNA/pangenome_build.sh << 'ENDOFSCRIPT'
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
GENOME_DIR=${BASE}/resources/black6/genome
PAN_DIR=${BASE}/results/pangenome
THREADS=32

# Singularity container (has cactus-pangenome + vg + minigraph)
CACTUS_SIF=${BASE}/cactus_v2.9.3.sif

# Activate conda for samtools/bcftools/bedtools
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake

mkdir -p ${PAN_DIR}/{prepared,logs,output,qc}

echo "$(date) ============================================"
echo "  MOUSE 17-STRAIN PANGENOME BUILD"
echo "  For piRNA cluster / TE analysis"
echo "============================================"

###############################################################################
# SINGULARITY WRAPPERS
#
#   cactus-pangenome and vg live inside the container.
#   samtools, bcftools, bedtools run natively from conda.
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

AVAIL=$(df -BG ${PAN_DIR} | awk 'NR==2{print \$4}' | sed 's/G//')
echo ""
echo "  Disk available: ${AVAIL} GB (need ~500GB)"

if [ "${AVAIL}" -lt 400 ]; then
    echo "  WARNING: Low disk space. Build may fail."
fi

echo ""
echo "$(date) === STEP 0 COMPLETE ==="

###############################################################################
# STEP 1: DIAGNOSE ALL ASSEMBLY NAMING
#
#   Your assemblies have 3 naming conventions:
#
#   T2T files:        >1 unmasked:chromosome primary_assembly:...
#                     Ensembl format — first word = chromosome, rest = metadata
#
#   REL-2205 strains: >129S1_SvImJ#1#chr1
#                     PanSN format — STRAIN#HAPLOTYPE#chrN
#                     ⚠ cactus adds its own STRAIN#HAP# → double prefix → CRASH
#
#   ALL must become:  >1  >2  >3 ... >19  >X  >Y  >MT
###############################################################################

echo ""
echo "$(date) === STEP 1: Diagnose assembly naming ==="
echo ""

CANONICAL_CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT)
echo "Target names: ${CANONICAL_CHROMS[*]}"
echo ""

printf "%-25s  %-10s  %-50s\n" "Assembly" "Format" "First_header"
printf "%-25s  %-10s  %-50s\n" \
    "─────────────────────────" "──────────" \
    "──────────────────────────────────────────────────"

show_format() {
    local FA="\$1"
    local LABEL="\$2"

    if [ ! -f "$FA" ]; then
        printf "%-25s  %-10s  %-50s\n" "$LABEL" "MISSING" "file not found"
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

show_format "${GENOME_DIR}/C57BL_6J_T2T.fa" "C57BL_6J_T2T"
show_format "${GENOME_DIR}/CAST_EiJ_T2T.fa" "CAST_EiJ_T2T"

REL_STRAINS=(C57BL_6NJ A_J AKR_J BALB_cJ C3H_HeJ CBA_J DBA_2J \
             FVB_NJ LP_J NOD_ShiLtJ NZO_HlLtJ 129S1_SvImJ \
             PWK_PhJ WSB_EiJ SPRET_EiJ)

for strain in "${REL_STRAINS[@]}"; do
    show_format "${REL_DIR}/${strain}_chromosomes_MT.fasta" "$strain"
done

echo ""
echo "$(date) === STEP 1 COMPLETE ==="

###############################################################################
# STEP 2: STANDARDIZE ALL CHROMOSOME NAMES
#
#   PanSN:   >STRAIN#1#chr1  → extract after last # → strip chr → >1
#   Ensembl: >1 unmasked:... → take first word only → >1
#
#   Uses samtools faidx per-chromosome (safe for multi-GB files)
#   Outputs ONLY canonical chromosomes in canonical order
#
#   piRNA relevance:
#     - ALL chromosomes: piRNA clusters on every chr
#     - X: pachytene piRNA clusters enriched
#     - MT: mitochondrial piRNAs exist
#     - Strip scaffolds/contigs: noise for piRNA mapping
###############################################################################

echo ""
echo "$(date) === STEP 2: Standardize chromosome names ==="
echo ""

standardize_genome() {
    local INPUT="\$1"
    local OUTPUT="\$2"
    local NAME="\$3"

    # Skip if done and validated
    if [ -f "$OUTPUT" ] && [ -f "${OUTPUT}.fai" ]; then
        local N
        N=$(wc -l < "${OUTPUT}.fai")
        if [ "$N" -ge 19 ]; then
            local NAMES
            NAMES=$(cut -f1 "${OUTPUT}.fai" | tr '\n' ' ')
            echo "  ${NAME}: already done (${N} chrs: ${NAMES})"
            return 0
        fi
        echo "  ${NAME}: incomplete output, rebuilding..."
        rm -f "$OUTPUT" "${OUTPUT}.fai"
    fi

    if [ ! -f "$INPUT" ]; then
        echo "  ${NAME}: ERROR — not found: $INPUT"
        return 1
    fi

    echo "  ${NAME}: processing..."
    [ ! -f "${INPUT}.fai" ] && samtools faidx "$INPUT"

    # Detect format from first sequence name
    local FIRST_NAME
    FIRST_NAME=$(awk 'NR==1{print $1}' "${INPUT}.fai")

    local FORMAT="unknown"
    if [[ "$FIRST_NAME" == *"#"*"#"* ]]; then
        FORMAT="PanSN"
    else
        FORMAT="Ensembl"
    fi

    echo "    Format: ${FORMAT} (first seq: ${FIRST_NAME})"

    # Build rename map: old_name → new_name
    local MAP="${PAN_DIR}/qc/${NAME}_rename_map.tsv"
    > "$MAP"

    while read -r OLD_NAME OLD_LEN; do
        local NEW_NAME=""

        if [ "$FORMAT" = "PanSN" ]; then
            # >STRAIN#1#chr1 → extract after last # → chr1 → strip chr → 1
            NEW_NAME=$(echo "$OLD_NAME" | awk -F'#' '{print $NF}')
            NEW_NAME=${NEW_NAME#chr}
            if [ "$NEW_NAME" = "M" ]; then
                NEW_NAME="MT"
            fi
        else
            # >1 or >X or >MT — first column of .fai is already the name
            NEW_NAME="$OLD_NAME"
            if [ "$NEW_NAME" = "M" ]; then
                NEW_NAME="MT"
            fi
        fi

        printf '%s\t%s\t%s\n' "$OLD_NAME" "$NEW_NAME" "$OLD_LEN" >> "$MAP"

    done < <(awk '{print \$1, $2}' "${INPUT}.fai")

    # Show mappings for canonical chromosomes
    echo "    Mappings:"
    for c in "${CANONICAL_CHROMS[@]}"; do
        local OLD
        OLD=$(awk -v target="$c" -F'\t' '\$2==target {print $1}' "$MAP")
        local SZ
        SZ=$(awk -v target="$c" -F'\t' '\$2==target {printf "%.1f Mb", $3/1e6}' "$MAP")
        if [ -n "$OLD" ]; then
            printf "      %-35s → %-4s  (%s)\n" "$OLD" "$c" "$SZ"
        fi
    done

    # Count skipped non-canonical sequences
    local TOTAL_IN
    TOTAL_IN=$(wc -l < "$MAP")
    local N_CANON=0
    for c in "${CANONICAL_CHROMS[@]}"; do
        if grep -q "	${c}	" "$MAP" 2>/dev/null; then
            N_CANON=$((N_CANON + 1))
        fi
    done
    local N_SKIP=$((TOTAL_IN - N_CANON))
    if [ "$N_SKIP" -gt 0 ]; then
        echo "      ... plus ${N_SKIP} scaffolds/contigs → skipped"
    fi

    # Extract each canonical chromosome and rename header
    local TMPFA="${OUTPUT}.tmp"
    > "$TMPFA"

    local FOUND_COUNT=0
    for c in "${CANONICAL_CHROMS[@]}"; do
        local OLD
        OLD=$(awk -v target="$c" -F'\t' '\$2==target {print $1; exit}' "$MAP")
        if [ -z "$OLD" ]; then
            continue
        fi

        samtools faidx "$INPUT" "$OLD" 2>/dev/null | \
            awk -v newname="$c" '/^>/{print ">" newname; next} {print}' >> "$TMPFA"

        FOUND_COUNT=$((FOUND_COUNT + 1))
    done

    if [ "$FOUND_COUNT" -lt 19 ]; then
        echo "    ERROR: Only ${FOUND_COUNT} canonical chromosomes found!"
        rm -f "$TMPFA"
        return 1
    fi

    # Index and finalize
    samtools faidx "$TMPFA"
    mv "$TMPFA" "$OUTPUT"
    mv "${TMPFA}.fai" "${OUTPUT}.fai"

    local TOTAL
    TOTAL=$(awk '{sum+=$2} END{printf "%.2f Gb", sum/1e9}' "${OUTPUT}.fai")
    local FINAL_NAMES
    FINAL_NAMES=$(cut -f1 "${OUTPUT}.fai" | tr '\n' ' ')

    echo "    ✓ ${FOUND_COUNT} chromosomes, ${TOTAL}"
    echo "    ✓ Names: ${FINAL_NAMES}"

    # Sanity check chr1 size
    local CHR1_MB
    CHR1_MB=$(awk '\$1=="1"{printf "%.0f", $2/1e6}' "${OUTPUT}.fai")
    if [ -n "$CHR1_MB" ] && [ "$CHR1_MB" -lt 100 ]; then
        echo "    ⚠ Chr1 only ${CHR1_MB} Mb — suspicious!"
    fi

    return 0
}

# ── T2T assemblies ──
echo "--- T2T Assemblies ---"
echo "    Format: >1 unmasked:chromosome... → take first word → >1"
echo ""

standardize_genome \
    "${GENOME_DIR}/C57BL_6J_T2T.fa" \
    "${PAN_DIR}/prepared/C57BL_6J.fa" \
    "C57BL_6J_T2T"
echo ""

standardize_genome \
    "${GENOME_DIR}/CAST_EiJ_T2T.fa" \
    "${PAN_DIR}/prepared/CAST_EiJ.fa" \
    "CAST_EiJ_T2T"
echo ""

# ── REL-2205 strain assemblies ──
echo "--- REL-2205 Assemblies ---"
echo "    Format: >STRAIN#1#chr1 → strip STRAIN#1#chr → >1"
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
#
#   CRITICAL GATE: every assembly must have IDENTICAL chromosome names.
#   Pangenome build WILL fail if any mismatch exists.
###############################################################################

echo ""
echo "$(date) === STEP 3: Cross-assembly validation ==="
echo ""

# Reference chromosome set
REF_CHROMS=$(cut -f1 ${PAN_DIR}/prepared/C57BL_6J.fa.fai | sort)

printf "%-25s  %4s  %10s  %8s  %8s  %7s  %s\n" \
    "Assembly" "Chrs" "Total_Gb" "Chr1_Mb" "ChrX_Mb" "MT_bp" "Status"
printf "%-25s  %4s  %10s  %8s  %8s  %7s  %s\n" \
    "─────────────────────────" "────" "──────────" "────────" "────────" "───────" "──────"

ALL_MATCH=true

for fa in ${PAN_DIR}/prepared/*.fa; do
    NAME=$(basename "$fa" .fa)
    [ ! -f "${fa}.fai" ] && samtools faidx "$fa"

    NSEQ=$(wc -l < "${fa}.fai")
    TOTAL=$(awk '{sum+=$2} END{printf "%.3f", sum/1e9}' "${fa}.fai")
    CHR1=$(awk '\$1=="1"{printf "%.1f", $2/1e6}' "${fa}.fai")
    CHRX=$(awk '\$1=="X"{printf "%.1f", $2/1e6}' "${fa}.fai")
    MT=$(awk '\$1=="MT"{print $2}' "${fa}.fai")

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

# Header cleanliness check — no spaces allowed in FASTA headers
echo "--- Header cleanliness check ---"
CLEAN=true
for fa in ${PAN_DIR}/prepared/*.fa; do
    NAME=$(basename "$fa" .fa)
    BAD=$(grep "^>" "$fa" | grep " " | head -1)
    if [ -n "$BAD" ]; then
        echo "  PROBLEM: ${NAME}: header has spaces: ${BAD}"
        CLEAN=false
    fi
done
if [ "$CLEAN" = true ]; then
    echo "  OK: All headers clean (no spaces or metadata)"
fi
echo ""

# If mismatches found, find common chromosome set and trim
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

    if [ "$N_COMMON" -lt 19 ]; then
        echo "  FATAL: fewer than 19 common chromosomes"
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

# Pairwise verification via md5
echo ""
echo "--- Pairwise md5 verification ---"
ASSEMBLIES=(${PAN_DIR}/prepared/*.fa)
REF_HASH=$(cut -f1 "${ASSEMBLIES[0]}.fai" | sort | md5sum | awk '{print \$1}')
PAIR_OK=true

for fa in "${ASSEMBLIES[@]}"; do
    NAME=$(basename "$fa" .fa)
    HASH=$(cut -f1 "${fa}.fai" | sort | md5sum | awk '{print \$1}')
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
#   Order matters for cactus-pangenome:
#     C57BL_6J first (T2T reference — most complete)
#     CAST_EiJ (T2T wild-derived M. m. castaneus)
#     Lab strains (close to reference)
#     Wild-derived (progressively divergent)
#     SPRET_EiJ last (most divergent ~2My, different TE landscape)
#
#   piRNA relevance:
#     T2T reference captures centromeric/pericentromeric piRNA clusters
#     SPRET last because most divergent TE repertoire
###############################################################################

echo ""
echo "$(date) === STEP 4: Create seqfile ==="
echo ""

SEQFILE=${PAN_DIR}/seqfile.txt

cat > ${SEQFILE} << SEQEOF
C57BL_6J	${PAN_DIR}/prepared/C57BL_6J.fa
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
echo ""
echo "  cactus will create PanSN paths:"
echo "    C57BL_6J#1#1, C57BL_6J#1#2, ..., SPRET_EiJ#1#MT"

if [ ${VALID_COUNT} -lt 2 ]; then
    echo "  FATAL: need at least 2 assemblies"
    exit 1
fi

echo ""
echo "$(date) === STEP 4 COMPLETE ==="

###############################################################################
# STEP 5: BUILD PANGENOME WITH MINIGRAPH-CACTUS
#
#   PARAMETER CHOICES FOR piRNA / TE BIOLOGY:
#
#   --reference C57BL_6J
#     T2T genome = most complete assembly
#     Captures centromeric sequence (piRNA-producing)
#     Captures satellite arrays (piRNA targets)
#     Captures pericentromeric TEs (piRNA cluster-rich)
#
#   --vcf
#     Variants across 17 strains
#     TE insertions appear as large INDELs (300bp-10kb)
#     Intersect with piRNA clusters to find strain-specific
#     TE insertions creating new piRNA defense loci
#
#   --giraffe (produces .gbz + .dist + .min)
#     For mapping sRNA-seq reads to pangenome graph
#     vg giraffe handles multi-mapping in TE repeats better
#     than linear aligners — CRITICAL for piRNA reads (24-32nt)
#
#   --gfa
#     Graph viewable in Bandage
#     Visualize piRNA cluster loci as graph bubbles
#
#   minigraph captures SVs >= ~50bp:
#     SINEs (B1/B2):  150-300bp  captured
#     LINEs (L1):     1-7kb      captured
#     LTRs (IAP):     5-7kb      captured
#     LTRs (ETn):     5-9kb      captured
#     ALL TE families relevant to piRNA defense
#
#   cactus base-level alignment:
#     Resolves individual TE copies in tandem arrays
#     Critical for piRNA clusters which are TE-dense
###############################################################################

echo ""
echo "$(date) === STEP 5: Building pangenome ==="
echo ""
echo "  Reference:    C57BL_6J (T2T)"
echo "  Assemblies:   ${VALID_COUNT}"
echo "  Threads:      ${THREADS}"
echo "  Memory:       240G"
echo "  Container:    ${CACTUS_SIF}"
echo ""

JOBSTORE=${PAN_DIR}/jobstore
OUTDIR=${PAN_DIR}/output
OUTNAME=mouse_17strain_pirna

# Clean previous attempt
[ -d "${JOBSTORE}" ] && rm -rf ${JOBSTORE}

START_TIME=$(date +%s)

run_cactus cactus-pangenome ${JOBSTORE} \
    ${SEQFILE} \
    --outDir ${OUTDIR} \
    --outName ${OUTNAME} \
    --reference C57BL_6J \
    --vcf \
    --giraffe \
    --gfa \
    --gbz \
    --maxCores ${THREADS} \
    --maxMemory 240G \
    --consCores ${THREADS} \
    --logFile ${PAN_DIR}/logs/cactus_pangenome.log \
    --realTimeLogging \
    2>&1 | tee ${PAN_DIR}/logs/cactus_stdout.log

END_TIME=$(date +%s)
ELAPSED=$(( (END_TIME - START_TIME) / 3600 ))

echo ""
echo "$(date) === STEP 5 COMPLETE (${ELAPSED} hours) ==="

###############################################################################
# STEP 6: VERIFY / BUILD VG GIRAFFE INDEXES
#
#   Needed to map piRNA (sRNA-seq) reads to pangenome.
#   --giraffe flag above should produce them, but verify.
#
#   piRNA reads are 24-32nt — very short.
#   Graph-based mapping disambiguates TE multi-mappers
#   using haplotype context.
###############################################################################

echo ""
echo "$(date) === STEP 6: Verify vg indexes ==="
echo ""

GBZ=${OUTDIR}/${OUTNAME}.gbz

if [ ! -f "${GBZ}" ]; then
    echo "FATAL: GBZ not produced!"
    echo "Check: ${PAN_DIR}/logs/cactus_pangenome.log"
    ls -lh ${OUTDIR}/ 2>/dev/null || echo "  output dir empty"
    exit 1
fi

echo "  GBZ:  $(ls -lh ${GBZ} | awk '{print \$5}')"

# Distance index
DIST=${OUTDIR}/${OUTNAME}.dist
if [ ! -f "${DIST}" ]; then
    echo "  Building distance index..."
    run_vg snarls ${GBZ} > ${OUTDIR}/${OUTNAME}.snarls
    run_vg index -j ${DIST} ${GBZ}
fi
echo "  DIST: $(ls -lh ${DIST} | awk '{print \$5}')"

# Minimizer index
MIN=${OUTDIR}/${OUTNAME}.min
if [ ! -f "${MIN}" ]; then
    echo "  Building minimizer index..."
    run_vg minimizer -t ${THREADS} -d ${DIST} -o ${MIN} ${GBZ}
fi
echo "  MIN:  $(ls -lh ${MIN} | awk '{print \$5}')"

echo ""
echo "$(date) === STEP 6 COMPLETE ==="

###############################################################################
# STEP 7: PANGENOME STATISTICS
#
#   Key metrics for piRNA/TE biology:
#   - Graph excess = non-reference TE insertions
#   - SV size distribution matches TE family sizes
#   - Per-strain variants (divergent strains = more TEs)
###############################################################################

echo ""
echo "$(date) === STEP 7: Pangenome statistics ==="
echo ""

echo "=================================================================="
echo "  MOUSE 17-STRAIN PANGENOME — piRNA / TE STATISTICS"
echo "=================================================================="

echo ""
echo "--- Output files ---"
ls -lh ${OUTDIR}/${OUTNAME}.*

REF_BP=$(awk '{sum+=$2} END{print sum}' ${PAN_DIR}/prepared/C57BL_6J.fa.fai)

# GFA statistics
GFA=${OUTDIR}/${OUTNAME}.gfa
if [ -f "${GFA}" ]; then
    echo ""
    echo "--- Graph topology ---"

    N_SEGMENTS=$(grep -c '^S' ${GFA} || echo 0)
    N_LINKS=$(grep -c '^L' ${GFA} || echo 0)
    N_PATHS=$(grep -c '^[PW]' ${GFA} || echo 0)
    GRAPH_BP=$(awk '/^S/{sum += length($3)} END{print sum}' ${GFA})

    EXCESS_MB=$(echo "scale=2; (${GRAPH_BP} - ${REF_BP}) / 1000000" | bc)
    EXCESS_PCT=$(echo "scale=1; 100 * (${GRAPH_BP} - ${REF_BP}) / ${REF_BP}" | bc)

    echo "  Segments (nodes):     ${N_SEGMENTS}"
    echo "  Links (edges):        ${N_LINKS}"
    echo "  Paths (haplotypes):   ${N_PATHS}"
    echo ""
    echo "  Total graph sequence: $(echo ${GRAPH_BP} | numfmt --to=iec) bp"
    echo "  Reference (T2T):      $(echo ${REF_BP} | numfmt --to=iec) bp"
    echo "  Non-reference excess: ${EXCESS_MB} Mb (${EXCESS_PCT}%)"
    echo ""
    echo "  The excess contains strain-specific TE insertions"
    echo "  that may form novel piRNA clusters"
fi

# VCF statistics
VCF=${OUTDIR}/${OUTNAME}.vcf.gz
if [ -f "${VCF}" ]; then
    echo ""
    echo "--- Variant statistics ---"

    TOTAL_VAR=$(bcftools view -H ${VCF} 2>/dev/null | wc -l)
    N_SNP=$(bcftools view -v snps -H ${VCF} 2>/dev/null | wc -l)
    N_INDEL=$(bcftools view -v indels -H ${VCF} 2>/dev/null | wc -l)

    echo "  Total variants: ${TOTAL_VAR}"
    echo "  SNPs:           ${N_SNP}"
    echo "  Indels:         ${N_INDEL}"

    echo ""
    echo "  TE-relevant SV size distribution:"
    echo ""

    bcftools query -f '%INFO/SVLEN\n' ${VCF} 2>/dev/null | \
    awk '{
        v = (\$1 < 0) ? -\$1 : \$1
    }
    v >= 50 && v < 300    { sine++ }
    v >= 300 && v < 1000  { small++ }
    v >= 1000 && v < 3000 { trunc++ }
    v >= 3000 && v < 7000 { full++ }
    v >= 7000 && v < 15000{ ltr++ }
    v >= 15000            { mega++ }
    END {
        print "    50-300bp    (SINE B1/B2):           " sine+0
        print "    300bp-1kb   (small TE / solo LTR):  " small+0
        print "    1-3kb       (truncated LINE):       " trunc+0
        print "    3-7kb       (full-length LINE/LTR): " full+0
        print "    7-15kb      (IAP/ETn with LTRs):    " ltr+0
        print "    >15kb       (complex/tandem TE):    " mega+0
    }'

    # Per-strain variant counts
    echo ""
    echo "--- Variants per strain (vs C57BL_6J reference) ---"
    echo ""

    bcftools query -l ${VCF} 2>/dev/null | while read sample; do
        N_ALT=$(bcftools view -s "${sample}" -H ${VCF} 2>/dev/null | \
                awk -F'\t' '{
                    split(\$10,gt,":")
                    if(gt[1]!="0/0" && gt[1]!="0|0" && gt[1]!="./.") print
                }' | wc -l)
        printf "    %-20s %8d variants\n" "${sample}" "${N_ALT}"
    done
fi

echo ""
echo "$(date) === STEP 7 COMPLETE ==="

###############################################################################
# STEP 8: EXTRACT piRNA-RELEVANT TE VARIANTS
#
#   Filter pangenome VCF for TE-sized insertions/deletions:
#   - New TE insertions → may create strain-specific piRNA clusters
#   - TE deletions → loss of piRNA cluster material
#   - Size bins match mouse TE families:
#     B1/B2 SINEs:  ~150-300bp
#     L1 LINEs:     ~1-7kb (full ~6.5kb)
#     IAP LTRs:     ~5-7kb
#     ETn LTRs:     ~5-9kb
#     Solo LTRs:    ~300-500bp
###############################################################################

echo ""
echo "$(date) === STEP 8: Extract piRNA-relevant TE variants ==="
echo ""

if [ -f "${VCF}" ]; then

    # TE-sized variants (300bp - 15kb)
    echo "  Extracting TE-sized structural variants..."

    bcftools view -H ${VCF} 2>/dev/null | \
    awk -F'\t' '{
        ref_len = length(\$4)
        alt_len = length(\$5)
        sv_size = (alt_len > ref_len) ? alt_len - ref_len : ref_len - alt_len
        sv_type = (alt_len > ref_len) ? "INS" : "DEL"
        if (sv_size >= 300 && sv_size <= 15000) {
            print \$1"\t"\$2-1"\t"\$2+ref_len"\t"sv_type"_"sv_size"bp\t"sv_size"\t."
        }
    }' | sort -k1,1 -k2,2n > ${OUTDIR}/te_sized_variants.bed

    N_TE=$(wc -l < ${OUTDIR}/te_sized_variants.bed)
    echo "  Total TE-sized variants: ${N_TE}"

    # Separate by type
    grep "INS" ${OUTDIR}/te_sized_variants.bed > ${OUTDIR}/te_insertions.bed 2>/dev/null || true
    grep "DEL" ${OUTDIR}/te_sized_variants.bed > ${OUTDIR}/te_deletions.bed 2>/dev/null || true

    N_INS=$(wc -l < ${OUTDIR}/te_insertions.bed 2>/dev/null || echo 0)
    N_DEL=$(wc -l < ${OUTDIR}/te_deletions.bed 2>/dev/null || echo 0)

    echo "    Insertions: ${N_INS} (new TEs in some strains)"
    echo "    Deletions:  ${N_DEL} (TE losses in some strains)"

    echo ""
    echo "  Size distribution:"
    awk '{
        size = \$5
        if (size < 500)       bin = "300-500bp  (SINE/solo LTR)"
        else if (size < 1500) bin = "500-1.5kb  (truncated element)"
        else if (size < 4000) bin = "1.5-4kb    (LINE fragment)"
        else if (size < 7000) bin = "4-7kb      (full LINE/IAP)"
        else                  bin = "7-15kb     (large LTR/tandem)"
        print bin
    }' ${OUTDIR}/te_sized_variants.bed | sort | uniq -c | sort -rn | \
    awk '{printf "    %6d  %s\n", \$1, substr(\$0, index(\$0,\$2))}'

    # Per-chromosome distribution
    echo ""
    echo "  Per-chromosome TE variant counts:"
    cut -f1 ${OUTDIR}/te_sized_variants.bed | sort -V | uniq -c | \
    awk '{printf "    chr%-3s  %5d TE-sized variants\n", \$2, \$1}'

else
    echo "  VCF not found — skipping"
fi

echo ""
echo "$(date) === STEP 8 COMPLETE ==="

###############################################################################
# STEP 9: VARIANT DENSITY ALONG CHROMOSOMES
#
#   High-density regions = structural variation hotspots
#   Many overlap piRNA cluster loci (piRNA clusters are TE traps)
###############################################################################

echo ""
echo "$(date) === STEP 9: Variant density ==="
echo ""

if [ -f "${VCF}" ]; then
    echo "chr	position_Mb	variant_count" > ${OUTDIR}/variant_density_per_Mb.tsv

    for CHR in $(seq 1 19) X; do
        bcftools view -r ${CHR} -H ${VCF} 2>/dev/null | \
        awk -v chr=${CHR} '{print chr"\t"int(\$2/1000000)}' | \
        sort | uniq -c | \
        awk '{print \$2"\t"\$3"\t"$1}' >> ${OUTDIR}/variant_density_per_Mb.tsv
    done

    echo "  Saved: ${OUTDIR}/variant_density_per_Mb.tsv"

    echo ""
    echo "  Top 10 variant-dense 1Mb windows:"
    tail -n+2 ${OUTDIR}/variant_density_per_Mb.tsv | \
    sort -t$'\t' -k3 -rn | head -10 | \
    awk -F'\t' '{printf "    chr%-3s  %3d-%3d Mb  %6d variants\n", \$1, \$2, \$2+1, \$3}'
fi

echo ""
echo "$(date) === STEP 9 COMPLETE ==="

###############################################################################
# STEP 10: piRNA ANALYSIS-READY FILES
#
#   Create files needed for downstream piRNA cluster detection
###############################################################################

echo ""
echo "$(date) === STEP 10: piRNA analysis prep files ==="
echo ""

# Chromosome sizes for bedtools
awk '{print \$1"\t"$2}' ${PAN_DIR}/prepared/C57BL_6J.fa.fai > ${OUTDIR}/chrom.sizes
echo "  chrom.sizes"

# Strain list
awk -F'\t' '{print $1}' ${SEQFILE} > ${OUTDIR}/strain_list.txt
echo "  strain_list.txt"

# Per-strain genome sizes
echo ""
echo "  Per-strain genome sizes:"
for fa in ${PAN_DIR}/prepared/*.fa; do
    NAME=$(basename "$fa" .fa)
    SIZE=$(awk '{sum+=$2} END{printf "%.3f Gb", sum/1e9}' "${fa}.fai")
    DIFF=$(awk -v ref=${REF_BP} '{sum+=$2} END{printf "%+.1f Mb", (sum-ref)/1e6}' "${fa}.fai")
    printf "    %-20s  %s  (%s vs ref)\n" "$NAME" "$SIZE" "$DIFF"
done

echo ""
echo "$(date) === STEP 10 COMPLETE ==="

###############################################################################
# FINAL SUMMARY
###############################################################################

echo ""
echo "=================================================================="
echo "  PANGENOME BUILD COMPLETE"
echo "  17 mouse strains — optimized for piRNA / TE analysis"
echo "=================================================================="
echo ""
echo "  Output: ${OUTDIR}"
echo ""
echo "  PANGENOME FILES:"
echo "    ${OUTNAME}.gbz              Graph genome (vg giraffe mapping)"
echo "    ${OUTNAME}.gfa              GFA graph (view in Bandage)"
echo "    ${OUTNAME}.vcf.gz           All variants across 17 strains"
echo "    ${OUTNAME}.dist             Distance index"
echo "    ${OUTNAME}.min              Minimizer index"
echo ""
echo "  piRNA / TE FILES:"
echo "    te_sized_variants.bed       TE-sized SVs (300bp-15kb)"
echo "    te_insertions.bed           New TE insertions"
echo "    te_deletions.bed            TE deletions"
echo "    variant_density_per_Mb.tsv  Hotspot detection"
echo "    chrom.sizes                 For bedtools"
echo ""
echo "  QC FILES:"
echo "    ${PAN_DIR}/qc/*_rename_map.tsv    Chromosome rename maps"
echo ""
echo "  NEXT STEPS FOR piRNA ANALYSIS:"
echo ""
echo "  1. MAP sRNA-seq reads to pangenome:"
echo "     singularity exec ${CACTUS_SIF} vg giraffe \\"
echo "       -Z ${OUTDIR}/${OUTNAME}.gbz \\"
echo "       -m ${OUTDIR}/${OUTNAME}.min \\"
echo "       -d ${OUTDIR}/${OUTNAME}.dist \\"
echo "       -f reads.fq.gz -t 16 -o gam > sample.gam"
echo ""
echo "  2. SURJECT to strain coordinates:"
echo "     singularity exec ${CACTUS_SIF} vg surject \\"
echo "       -x ${OUTDIR}/${OUTNAME}.gbz \\"
echo "       -b -R 'C57BL_6J' sample.gam > sample.bam"
echo ""
echo "  3. DETECT piRNA clusters:"
echo "     Run proTRAC or PICB on surjected BAMs per strain"
echo ""
echo "  4. INTERSECT clusters with TE variants:"
echo "     bedtools intersect -a pirna_clusters.bed \\"
echo "       -b ${OUTDIR}/te_insertions.bed -wa -wb"
echo ""
echo "  5. COMPARE with liftover-based piRNA cluster analysis"
echo ""
echo "  File sizes:"
ls -lh ${OUTDIR}/${OUTNAME}.* 2>/dev/null
echo ""
echo "$(date) === ALL COMPLETE ==="
ENDOFSCRIPT

echo "Script saved to: ~/scratch/inProgress/mice_PiRNA/pangenome_build.sh"
echo "Lines: $(wc -l < ~/scratch/inProgress/mice_PiRNA/pangenome_build.sh)"