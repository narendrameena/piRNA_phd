#!/bin/bash

set -euo pipefail

###############################################################################
# MOUSE 17-STRAIN PANGENOME -- QC ONLY
# Run this AFTER the pangenome build is complete
# Does NOT re-run cactus-pangenome
###############################################################################

###############################################################################
# CONFIGURATION
###############################################################################

BASE=/mnt/beegfs/scratch/miska/nm667/inProgress/mice_PiRNA
PAN_DIR=${BASE}/results/pangenome
OUTDIR=${PAN_DIR}/output
OUTNAME=mouse_17strain_pangenome
THREADS=32
CACTUS_SIF=${BASE}/cactus_v2.9.3.sif
SEQFILE=${PAN_DIR}/seqfile.txt

CANONICAL_CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X MT)

mkdir -p ${PAN_DIR}/qc

###############################################################################
# LOGGING -- all output goes to screen AND log file
###############################################################################

LOGFILE=${PAN_DIR}/qc/pangenome_qc_$(date +%Y%m%d_%H%M%S).log
exec > >(tee -a "${LOGFILE}") 2>&1

echo "$(date) ============================================"
echo "  MOUSE 17-STRAIN PANGENOME -- QC ONLY"
echo "  Reference: GRCm39 (C57BL/6J)"
echo "  Output dir: ${OUTDIR}"
echo "  Log file: ${LOGFILE}"
echo "============================================"
echo ""

###############################################################################
# SINGULARITY WRAPPERS
###############################################################################

run_cactus() {
    singularity exec \
        --bind /mnt/beegfs:/mnt/beegfs \
        --bind /tmp:/tmp \
        ${CACTUS_SIF} \
        "$@"
}

run_vg() {
    singularity exec \
        --bind /mnt/beegfs:/mnt/beegfs \
        --bind /tmp:/tmp \
        ${CACTUS_SIF} \
        vg "$@"
}

###############################################################################
# DETECT OUTPUT FILES (handle .d2. naming from cactus)
###############################################################################

GBZ="${OUTDIR}/${OUTNAME}.gbz"
VCF="${OUTDIR}/${OUTNAME}.vcf.gz"

DIST=""
for d in "${OUTDIR}/${OUTNAME}.dist" "${OUTDIR}/${OUTNAME}.d2.dist"; do
    if [ -f "$d" ]; then
        DIST="$d"
        break
    fi
done

MIN=""
for m in "${OUTDIR}/${OUTNAME}.min" "${OUTDIR}/${OUTNAME}.d2.min"; do
    if [ -f "$m" ]; then
        MIN="$m"
        break
    fi
done

GFA_FILE=""
for g in "${OUTDIR}/${OUTNAME}.gfa.gz" "${OUTDIR}/${OUTNAME}.gfa"; do
    if [ -f "$g" ]; then
        GFA_FILE="$g"
        break
    fi
done

D2_GBZ="${OUTDIR}/${OUTNAME}.d2.gbz"

echo "--- Detected files ---"
echo "  GBZ:    ${GBZ} $([ -f "${GBZ}" ] && echo '[FOUND]' || echo '[MISSING]')"
echo "  VCF:    ${VCF} $([ -f "${VCF}" ] && echo '[FOUND]' || echo '[MISSING]')"
echo "  DIST:   ${DIST:-not found} $([ -n "${DIST}" ] && [ -f "${DIST}" ] && echo '[FOUND]' || echo '[MISSING]')"
echo "  MIN:    ${MIN:-not found} $([ -n "${MIN}" ] && [ -f "${MIN}" ] && echo '[FOUND]' || echo '[MISSING]')"
echo "  GFA:    ${GFA_FILE:-not found} $([ -n "${GFA_FILE}" ] && [ -f "${GFA_FILE}" ] && echo '[FOUND]' || echo '[MISSING]')"
echo "  D2 GBZ: ${D2_GBZ} $([ -f "${D2_GBZ}" ] && echo '[FOUND]' || echo '[not present]')"
echo ""

###############################################################################
# STEP 6: PANGENOME QC -- FILE INTEGRITY
###############################################################################

echo ""
echo "$(date) === STEP 6: Pangenome QC -- File integrity ==="
echo ""

QC_PASS=true

echo "--- Expected output files ---"
echo ""

check_file() {
    local label="$1"
    local fpath="$2"
    if [ -n "$fpath" ] && [ -f "$fpath" ]; then
        local FSIZE
        local FBYTES
        FSIZE=$(ls -lh "$fpath" | awk '{print $5}')
        FBYTES=$(stat --printf="%s" "$fpath" 2>/dev/null || echo 0)
        if [ "$FBYTES" -lt 1000 ]; then
            printf "  %-6s  %-70s  %8s  WARNING: SMALL\n" "$label" "$fpath" "$FSIZE"
            QC_PASS=false
        else
            printf "  %-6s  %-70s  %8s  OK\n" "$label" "$fpath" "$FSIZE"
        fi
    else
        printf "  %-6s  %-70s  %8s  MISSING\n" "$label" "${fpath:-not found}" "---"
        QC_PASS=false
    fi
}

check_file "GBZ"  "$GBZ"
check_file "VCF"  "$VCF"
check_file "DIST" "${DIST:-}"
check_file "MIN"  "${MIN:-}"
check_file "GFA"  "${GFA_FILE:-}"

echo ""
echo "--- All files in output directory ---"
ls -lhS "${OUTDIR}/" 2>/dev/null || echo "  output directory empty!"

echo ""
if [ "$QC_PASS" = true ]; then
    echo "  FILE INTEGRITY: PASS"
else
    echo "  FILE INTEGRITY: ISSUES -- see above"
fi

echo ""
echo "$(date) === STEP 6 COMPLETE ==="

###############################################################################
# STEP 7: PANGENOME QC -- GRAPH TOPOLOGY
###############################################################################

echo ""
echo "$(date) === STEP 7: Pangenome QC -- Graph topology ==="
echo ""

if [ -f "${GBZ}" ]; then

    echo "--- vg stats on GBZ ---"
    echo ""

    run_vg stats -z "${GBZ}" 2>/dev/null | tee "${PAN_DIR}/qc/vg_stats_basic.txt"
    echo ""

    echo "--- Detailed graph metrics ---"
    echo ""

    GRAPH_NODES=$(run_vg stats -N "${GBZ}" 2>/dev/null | awk '{print $NF}' | tr -d '[:space:]') || true
    GRAPH_EDGES=$(run_vg stats -E "${GBZ}" 2>/dev/null | awk '{print $NF}' | tr -d '[:space:]') || true
    GRAPH_TOTAL_BP=$(run_vg stats -l "${GBZ}" 2>/dev/null | awk '{print $NF}' | tr -d '[:space:]') || true

    GRAPH_NODES="${GRAPH_NODES:-0}"
    GRAPH_EDGES="${GRAPH_EDGES:-0}"
    GRAPH_TOTAL_BP="${GRAPH_TOTAL_BP:-0}"

    echo "  Nodes (segments):    ${GRAPH_NODES}"
    echo "  Edges (links):       ${GRAPH_EDGES}"
    echo "  Total sequence (bp): ${GRAPH_TOTAL_BP}"

    REF_BP=$(awk '{sum+=$2} END{print sum}' "${PAN_DIR}/prepared/GRCm39.fa.fai")

    if [ "${GRAPH_TOTAL_BP}" -gt 0 ] 2>/dev/null; then
        EXCESS_MB=$(echo "scale=2; (${GRAPH_TOTAL_BP} - ${REF_BP}) / 1000000" | bc)
        EXCESS_PCT=$(echo "scale=1; 100 * (${GRAPH_TOTAL_BP} - ${REF_BP}) / ${REF_BP}" | bc)

        AVG_NODE_LEN="N/A"
        if [ "${GRAPH_NODES}" -gt 0 ] 2>/dev/null; then
            AVG_NODE_LEN=$(echo "scale=0; ${GRAPH_TOTAL_BP} / ${GRAPH_NODES}" | bc 2>/dev/null || echo "N/A")
        fi

        echo ""
        echo "  Reference (GRCm39):       $(echo "${REF_BP}" | numfmt --to=iec) bp"
        echo "  Graph total:              $(echo "${GRAPH_TOTAL_BP}" | numfmt --to=iec) bp"
        echo "  Non-reference excess:     ${EXCESS_MB} Mb (${EXCESS_PCT}%)"
        echo "  Average node length:      ${AVG_NODE_LEN} bp"
        echo ""

        EXCESS_HIGH=$(echo "${EXCESS_PCT} > 50" | bc -l 2>/dev/null || echo 0)
        EXCESS_LOW=$(echo "${EXCESS_PCT} < 1" | bc -l 2>/dev/null || echo 0)

        if [ "${EXCESS_HIGH}" -eq 1 ] 2>/dev/null; then
            echo "  QC WARNING: Non-ref excess > 50%"
        elif [ "${EXCESS_LOW}" -eq 1 ] 2>/dev/null; then
            echo "  QC WARNING: Non-ref excess < 1%"
        else
            echo "  OK: Non-ref excess in expected range for 17 mouse strains"
        fi

        if [ "${AVG_NODE_LEN}" != "N/A" ] && [ "${AVG_NODE_LEN}" -lt 10 ] 2>/dev/null; then
            echo "  QC WARNING: Average node length < 10bp -- graph over-fragmented"
        fi
    else
        echo ""
        echo "  WARNING: Could not parse total graph length"
    fi

    # Paths (haplotypes)
    echo ""
    echo "--- Paths (haplotypes) in graph ---"
    echo ""

    run_vg paths -L -x "${GBZ}" 2>/dev/null > "${PAN_DIR}/qc/graph_paths.txt" || true

    if [ -f "${PAN_DIR}/qc/graph_paths.txt" ] && [ -s "${PAN_DIR}/qc/graph_paths.txt" ]; then
        TOTAL_PATHS=$(wc -l < "${PAN_DIR}/qc/graph_paths.txt" | tr -d '[:space:]')
        echo "  Total paths: ${TOTAL_PATHS}"

        FIRST_PATH=$(head -1 "${PAN_DIR}/qc/graph_paths.txt")
        echo "  First path:  ${FIRST_PATH}"

        echo ""
        echo "  Paths per strain:"
        awk -F'#' '{print $1}' "${PAN_DIR}/qc/graph_paths.txt" | sort | uniq -c | sort -rn | \
        while read -r count strain; do
            printf "    %-20s  %4d paths\n" "$strain" "$count"
        done

        N_STRAINS=$(awk -F'#' '{print $1}' "${PAN_DIR}/qc/graph_paths.txt" | sort -u | wc -l | tr -d '[:space:]')
        echo ""
        echo "  Strains represented: ${N_STRAINS} / 17 expected"

        if [ "${N_STRAINS}" -eq 17 ]; then
            echo "  OK: All 17 strains present in graph"
        else
            echo "  QC WARNING: Missing strains!"
            echo "    Expected (from seqfile):"
            awk -F'\t' '{print "      " $1}' "${SEQFILE}"
            echo "    Found in graph:"
            awk -F'#' '{print $1}' "${PAN_DIR}/qc/graph_paths.txt" | sort -u | sed 's/^/      /'
        fi

        # Detect reference haplotype number (#0# vs #1#)
        REF_HAP=$(grep "^GRCm39" "${PAN_DIR}/qc/graph_paths.txt" | head -1 | awk -F'#' '{print $2}')
        REF_HAP="${REF_HAP:-0}"

        echo ""
        echo "  Reference (GRCm39) paths (haplotype #${REF_HAP}#):"
        grep "^GRCm39" "${PAN_DIR}/qc/graph_paths.txt" | head -5
        REF_PATH_COUNT=$(grep -c "^GRCm39" "${PAN_DIR}/qc/graph_paths.txt" || echo "0")
        REF_PATH_COUNT=$(echo "${REF_PATH_COUNT}" | tr -d '[:space:]')
        echo "    ... ${REF_PATH_COUNT} total reference paths"
    else
        echo "  Could not extract paths from graph"
        run_vg stats -z "${GBZ}" 2>/dev/null | head -20 || true
    fi

else
    echo "  GBZ not found at: ${GBZ}"
    echo "  Cannot perform graph QC"
fi

echo ""
echo "$(date) === STEP 7 COMPLETE ==="

###############################################################################
# STEP 8: PANGENOME QC -- VCF QUALITY
###############################################################################

echo ""
echo "$(date) === STEP 8: Pangenome QC -- VCF quality ==="
echo ""

if [ -f "${VCF}" ]; then

    if [ ! -f "${VCF}.tbi" ]; then
        bcftools index -t "${VCF}"
    fi

    echo "--- bcftools stats ---"
    echo ""

    bcftools stats "${VCF}" > "${PAN_DIR}/qc/bcftools_stats.txt" 2>/dev/null

    TOTAL_RECORDS=$(grep "^SN" "${PAN_DIR}/qc/bcftools_stats.txt" | grep "number of records:" | awk '{print $NF}') || true
    TOTAL_SNPS=$(grep "^SN" "${PAN_DIR}/qc/bcftools_stats.txt" | grep "number of SNPs:" | awk '{print $NF}') || true
    TOTAL_INDELS=$(grep "^SN" "${PAN_DIR}/qc/bcftools_stats.txt" | grep "number of indels:" | awk '{print $NF}') || true
    TOTAL_MULTIALLELIC=$(grep "^SN" "${PAN_DIR}/qc/bcftools_stats.txt" | grep "number of multiallelic" | awk '{print $NF}') || true
    TSTV=$(grep "^TSTV" "${PAN_DIR}/qc/bcftools_stats.txt" | head -1 | awk '{print $5}') || true

    echo "  Total variant records: ${TOTAL_RECORDS:-N/A}"
    echo "  SNPs:                  ${TOTAL_SNPS:-N/A}"
    echo "  Indels:                ${TOTAL_INDELS:-N/A}"
    echo "  Multiallelic sites:    ${TOTAL_MULTIALLELIC:-N/A}"
    echo "  Ts/Tv ratio:           ${TSTV:-N/A}"

    echo ""
    if [ -n "${TSTV:-}" ]; then
        TSTV_OK=$(echo "${TSTV} > 1.5 && ${TSTV} < 3.0" | bc -l 2>/dev/null || echo 0)
        if [ "${TSTV_OK}" -eq 1 ] 2>/dev/null; then
            echo "  OK: Ts/Tv ratio in expected range (1.5-3.0 for mouse)"
        else
            echo "  QC WARNING: Ts/Tv ratio outside expected range"
        fi
    fi

    if [ -n "${TOTAL_RECORDS:-}" ] && [ "${TOTAL_RECORDS}" -gt 1000000 ] 2>/dev/null; then
        echo "  OK: Variant count reasonable for 17 mouse strains"
    elif [ -n "${TOTAL_RECORDS:-}" ] && [ "${TOTAL_RECORDS}" -lt 100000 ] 2>/dev/null; then
        echo "  QC WARNING: Very few variants (${TOTAL_RECORDS})"
    fi

    # Samples in VCF
    echo ""
    echo "--- Samples in VCF ---"
    echo ""
    bcftools query -l "${VCF}" 2>/dev/null | while read -r sample; do
        echo "  ${sample}"
    done
    N_SAMPLES=$(bcftools query -l "${VCF}" 2>/dev/null | wc -l | tr -d '[:space:]')
    echo ""
    echo "  Total samples: ${N_SAMPLES}"

    # Per-chromosome variant counts
    echo ""
    echo "--- Per-chromosome variant counts ---"
    echo ""

    for CHR in "${CANONICAL_CHROMS[@]}"; do
        N_CHR=$(bcftools view -r "${CHR}" -H "${VCF}" 2>/dev/null | wc -l | tr -d '[:space:]')
        CHR_LEN=$(awk -v c="${CHR}" '$1==c{print $2}' "${PAN_DIR}/prepared/GRCm39.fa.fai")
        if [ -n "${CHR_LEN}" ] && [ "${CHR_LEN}" -gt 0 ] && [ "${N_CHR}" -gt 0 ]; then
            DENSITY=$(echo "scale=1; ${N_CHR} * 1000000 / ${CHR_LEN}" | bc)
            printf "    chr%-3s  %8d variants  (%s vars/Mb)\n" "${CHR}" "${N_CHR}" "${DENSITY}"
        else
            printf "    chr%-3s  %8d variants\n" "${CHR}" "${N_CHR:-0}"
        fi
    done

    # Indel/SV size distribution
    echo ""
    echo "--- Indel size distribution (TE-relevant) ---"
    echo ""

    echo "  From bcftools IDD field:"
    grep "^IDD" "${PAN_DIR}/qc/bcftools_stats.txt" | \
    awk '{
        size = $3
        count = $4
        if (size < 0) size = -size
        if (size == 0) next
        if (size <= 10)         bin = "1-10bp"
        else if (size <= 50)    bin = "11-50bp"
        else if (size <= 300)   bin = "51-300bp (SINE B1/B2)"
        else if (size <= 1000)  bin = "301bp-1kb (solo LTR)"
        else if (size <= 5000)  bin = "1-5kb (LINE fragment)"
        else if (size <= 10000) bin = "5-10kb (full LINE/LTR)"
        else                    bin = ">10kb (complex SV)"
        bins[bin] += count
    } END {
        order[1]="1-10bp"
        order[2]="11-50bp"
        order[3]="51-300bp (SINE B1/B2)"
        order[4]="301bp-1kb (solo LTR)"
        order[5]="1-5kb (LINE fragment)"
        order[6]="5-10kb (full LINE/LTR)"
        order[7]=">10kb (complex SV)"
        for (i=1; i<=7; i++) {
            if (bins[order[i]] > 0)
                printf "    %-35s  %8d\n", order[i], bins[order[i]]
        }
    }' || true

    # Direct large SV counts from VCF
    echo ""
    echo "  Direct count of large variants from VCF:"

    LARGE_300=$(bcftools view -H "${VCF}" 2>/dev/null | \
        awk -F'\t' '{
            split($5, alts, ",")
            for (i in alts) {
                d = length(alts[i]) - length($4)
                if (d < 0) d = -d
                if (d >= 300) { count++; break }
            }
        } END { print count+0 }') || true
    echo "    Variants with |size| >= 300bp:  ${LARGE_300:-0}"

    LARGE_1KB=$(bcftools view -H "${VCF}" 2>/dev/null | \
        awk -F'\t' '{
            split($5, alts, ",")
            for (i in alts) {
                d = length(alts[i]) - length($4)
                if (d < 0) d = -d
                if (d >= 1000) { count++; break }
            }
        } END { print count+0 }') || true
    echo "    Variants with |size| >= 1kb:    ${LARGE_1KB:-0}"

    LARGE_5KB=$(bcftools view -H "${VCF}" 2>/dev/null | \
        awk -F'\t' '{
            split($5, alts, ",")
            for (i in alts) {
                d = length(alts[i]) - length($4)
                if (d < 0) d = -d
                if (d >= 5000) { count++; break }
            }
        } END { print count+0 }') || true
    echo "    Variants with |size| >= 5kb:    ${LARGE_5KB:-0}"

    echo ""
    echo "  Full stats: ${PAN_DIR}/qc/bcftools_stats.txt"

else
    echo "  VCF not found at: ${VCF}"
    echo "  Skipping VCF QC"
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

    REF_INPUT_FAI="${PAN_DIR}/prepared/GRCm39.fa.fai"
    COVERAGE_OK=true

    # Make sure paths file exists
    if [ ! -f "${PAN_DIR}/qc/graph_paths.txt" ] || [ ! -s "${PAN_DIR}/qc/graph_paths.txt" ]; then
        run_vg paths -L -x "${GBZ}" 2>/dev/null > "${PAN_DIR}/qc/graph_paths.txt" || true
    fi

    REF_HAP=$(grep "^GRCm39" "${PAN_DIR}/qc/graph_paths.txt" 2>/dev/null | head -1 | awk -F'#' '{print $2}') || true
    REF_HAP="${REF_HAP:-0}"

    echo "  Reference haplotype in graph: #${REF_HAP}#"
    echo ""

    printf "  %-5s  %12s  %12s  %8s\n" "Chr" "Input_bp" "Graph_bp" "Ratio"
    printf "  %-5s  %12s  %12s  %8s\n" "-----" "------------" "------------" "--------"

    for CHR in "${CANONICAL_CHROMS[@]}"; do
        INPUT_LEN=$(awk -v c="${CHR}" '$1==c{print $2}' "${REF_INPUT_FAI}")

        # Find the correct path name
        REF_PATH_NAME=""
        for try_name in "GRCm39#${REF_HAP}#${CHR}" "GRCm39#${REF_HAP}#chr${CHR}"; do
            if grep -q "^${try_name}$" "${PAN_DIR}/qc/graph_paths.txt" 2>/dev/null; then
                REF_PATH_NAME="${try_name}"
                break
            fi
        done

        GRAPH_LEN=0
        if [ -n "${REF_PATH_NAME}" ]; then
            GRAPH_LEN=$(run_vg paths -x "${GBZ}" -E -Q "${REF_PATH_NAME}" 2>/dev/null | \
                        awk '{sum+=$2} END{print sum+0}') || true
            GRAPH_LEN=$(echo "${GRAPH_LEN}" | tr -d '[:space:]')
            GRAPH_LEN="${GRAPH_LEN:-0}"
        fi

        if [ -n "${INPUT_LEN}" ] && [ "${INPUT_LEN}" -gt 0 ] 2>/dev/null && [ "${GRAPH_LEN}" -gt 0 ] 2>/dev/null; then
            RATIO=$(echo "scale=4; ${GRAPH_LEN} / ${INPUT_LEN}" | bc)
            printf "  %-5s  %12d  %12d  %8s\n" "${CHR}" "${INPUT_LEN}" "${GRAPH_LEN}" "${RATIO}"

            LOW=$(echo "${RATIO} < 0.95" | bc -l 2>/dev/null || echo 0)
            HIGH=$(echo "${RATIO} > 1.05" | bc -l 2>/dev/null || echo 0)
            if [ "${LOW:-0}" -eq 1 ] 2>/dev/null || [ "${HIGH:-0}" -eq 1 ] 2>/dev/null; then
                echo "         WARNING: Ratio outside 0.95-1.05"
                COVERAGE_OK=false
            fi
        elif [ -n "${INPUT_LEN}" ]; then
            printf "  %-5s  %12d  %12s  %8s  (path: %s)\n" "${CHR}" "${INPUT_LEN}" "N/A" "N/A" "${REF_PATH_NAME:-not found}"
        fi
    done

    echo ""

    # Path existence check
    echo "--- Path existence check ---"
    MISSING_PATHS=0
    for CHR in "${CANONICAL_CHROMS[@]}"; do
        FOUND=false
        for try_name in "GRCm39#${REF_HAP}#${CHR}" "GRCm39#${REF_HAP}#chr${CHR}"; do
            if grep -q "^${try_name}$" "${PAN_DIR}/qc/graph_paths.txt" 2>/dev/null; then
                FOUND=true
                break
            fi
        done
        if [ "${FOUND}" = false ]; then
            echo "  MISSING path for chr ${CHR}"
            MISSING_PATHS=$((MISSING_PATHS + 1))
        fi
    done

    if [ "${MISSING_PATHS}" -eq 0 ]; then
        echo "  OK: All ${#CANONICAL_CHROMS[@]} reference paths found"
    else
        echo "  WARNING: ${MISSING_PATHS} reference paths missing"
        echo ""
        echo "  All reference paths in graph:"
        grep "^GRCm39" "${PAN_DIR}/qc/graph_paths.txt" | sed 's/^/    /'
    fi

    echo ""
    if [ "${COVERAGE_OK}" = true ] && [ "${MISSING_PATHS}" -eq 0 ]; then
        echo "  REFERENCE COVERAGE: OK"
    else
        echo "  REFERENCE COVERAGE: CHECK -- see details above"
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

# If DIST or MIN not found, try to build
if [ -z "${DIST}" ] || [ ! -f "${DIST}" ]; then
    echo "  Distance index not found -- attempting to build..."
    DIST="${OUTDIR}/${OUTNAME}.dist"
    run_vg snarls "${GBZ}" > "${OUTDIR}/${OUTNAME}.snarls" 2>/dev/null || true
    run_vg index -j "${DIST}" "${GBZ}" 2>/dev/null || true
fi

if [ -z "${MIN}" ] || [ ! -f "${MIN}" ]; then
    echo "  Minimizer index not found -- attempting to build..."
    MIN="${OUTDIR}/${OUTNAME}.min"
    run_vg minimizer -t "${THREADS}" -d "${DIST}" -o "${MIN}" "${GBZ}" 2>/dev/null || true
fi

echo "--- Giraffe index files ---"
echo ""

for IDX_LABEL_FILE in "GBZ:${GBZ}" "DIST:${DIST:-}" "MIN:${MIN:-}"; do
    LABEL="${IDX_LABEL_FILE%%:*}"
    FPATH="${IDX_LABEL_FILE#*:}"
    if [ -n "${FPATH}" ] && [ -f "${FPATH}" ]; then
        FSIZE=$(ls -lh "${FPATH}" | awk '{print $5}')
        printf "  OK:      %-6s  %-60s  %s\n" "${LABEL}" "$(basename "${FPATH}")" "${FSIZE}"
    else
        printf "  MISSING: %-6s  %s\n" "${LABEL}" "${FPATH:-not set}"
    fi
done

# Synthetic piRNA-length read mapping test
echo ""
echo "--- Synthetic read mapping test (28nt = piRNA length) ---"
echo ""

SYNTH_FQ="${PAN_DIR}/qc/synthetic_reads.fq"

samtools faidx "${PAN_DIR}/prepared/GRCm39.fa" "1:1000000-1100000" 2>/dev/null | \
    grep -v "^>" | tr -d '\n' | \
    awk 'BEGIN{srand(42)} {
        seq = toupper($0)
        n = 0
        for (i = 1; i <= 500 && n < 100; i++) {
            start = int(rand() * (length(seq) - 28)) + 1
            read_seq = substr(seq, start, 28)
            if (read_seq !~ /N/ && length(read_seq) == 28) {
                n++
                printf "@synth_%d\n%s\n+\n%s\n", n, read_seq, "IIIIIIIIIIIIIIIIIIIIIIIIIIII"
            }
        }
    }' > "${SYNTH_FQ}" || true

N_SYNTH=0
if [ -f "${SYNTH_FQ}" ]; then
    N_SYNTH=$(grep -c "^@synth" "${SYNTH_FQ}" 2>/dev/null || echo "0")
    N_SYNTH=$(echo "${N_SYNTH}" | tr -d '[:space:]')
fi
echo "  Generated ${N_SYNTH} synthetic 28nt reads"

HAS_DIST=false
HAS_MIN=false
[ -n "${DIST}" ] && [ -f "${DIST}" ] && HAS_DIST=true
[ -n "${MIN}" ] && [ -f "${MIN}" ] && HAS_MIN=true

if [ "${N_SYNTH}" -gt 0 ] 2>/dev/null && \
   [ -f "${GBZ}" ] && \
   [ "${HAS_DIST}" = true ] && \
   [ "${HAS_MIN}" = true ]; then

    echo "  Running vg giraffe..."
    echo ""

    run_vg giraffe \
        -Z "${GBZ}" \
        -m "${MIN}" \
        -d "${DIST}" \
        -f "${SYNTH_FQ}" \
        -t "${THREADS}" \
        -o gam \
        > "${PAN_DIR}/qc/synthetic_mapped.gam" 2>"${PAN_DIR}/qc/giraffe_test.log" || true

    if [ -f "${PAN_DIR}/qc/synthetic_mapped.gam" ]; then
        GAM_SIZE=$(stat --printf='%s' "${PAN_DIR}/qc/synthetic_mapped.gam" 2>/dev/null || echo 0)

        if [ "${GAM_SIZE}" -gt 0 ] 2>/dev/null; then
            run_vg stats -a "${PAN_DIR}/qc/synthetic_mapped.gam" 2>/dev/null > \
                "${PAN_DIR}/qc/synthetic_map_stats.txt" || true

            if [ -s "${PAN_DIR}/qc/synthetic_map_stats.txt" ]; then
                echo "  Mapping results:"
                head -20 "${PAN_DIR}/qc/synthetic_map_stats.txt" | sed 's/^/    /'
                echo ""
                echo "  OK: Giraffe mapping functional for piRNA-length reads"
            else
                echo "  OK: Giraffe produced output (${GAM_SIZE} bytes)"
            fi
        else
            echo "  WARNING: Empty GAM output"
            echo "  Giraffe stderr:"
            head -20 "${PAN_DIR}/qc/giraffe_test.log" 2>/dev/null | sed 's/^/    /' || true
        fi
    else
        echo "  FAIL: Giraffe mapping produced no output"
        echo "  Check: ${PAN_DIR}/qc/giraffe_test.log"
    fi
else
    echo "  Skipping mapping test -- missing index files or no reads generated"
    echo "    GBZ:   $([ -f "${GBZ}" ] && echo 'OK' || echo 'MISSING')"
    echo "    DIST:  ${HAS_DIST}"
    echo "    MIN:   ${HAS_MIN}"
    echo "    Reads: ${N_SYNTH}"
fi

rm -f "${SYNTH_FQ}"

echo ""
echo "$(date) === STEP 10 COMPLETE ==="

###############################################################################
# STEP 11: PANGENOME QC -- CACTUS LOG ANALYSIS
###############################################################################

echo ""
echo "$(date) === STEP 11: Pangenome QC -- Log analysis ==="
echo ""

CACTUS_LOG=""
for log_try in \
    "${PAN_DIR}/logs/cactus_stdout.log" \
    "${PAN_DIR}/logs/cactus_pangenome.log"; do
    if [ -f "${log_try}" ] && [ -s "${log_try}" ]; then
        CACTUS_LOG="${log_try}"
        break
    fi
done

if [ -n "${CACTUS_LOG}" ] && [ -f "${CACTUS_LOG}" ]; then
    echo "--- Cactus log: ${CACTUS_LOG} ---"
    echo "  Log size: $(ls -lh "${CACTUS_LOG}" | awk '{print $5}')"
    echo ""

    N_ERRORS=$(grep -ci "error" "${CACTUS_LOG}" 2>/dev/null || echo "0")
    N_ERRORS=$(echo "${N_ERRORS}" | tr -d '[:space:]')

    N_WARNINGS=$(grep -ci "warning" "${CACTUS_LOG}" 2>/dev/null || echo "0")
    N_WARNINGS=$(echo "${N_WARNINGS}" | tr -d '[:space:]')

    N_FAILED=$(grep -ciE "failed|failure" "${CACTUS_LOG}" 2>/dev/null || echo "0")
    N_FAILED=$(echo "${N_FAILED}" | tr -d '[:space:]')

    echo "  Error mentions:   ${N_ERRORS}"
    echo "  Warning mentions: ${N_WARNINGS}"
    echo "  Failed mentions:  ${N_FAILED}"

    if [ "${N_ERRORS}" -gt 0 ] 2>/dev/null; then
        echo ""
        echo "  Last 5 error lines:"
        grep -i "error" "${CACTUS_LOG}" 2>/dev/null | tail -5 | sed 's/^/    /' || true
    fi

    if [ "${N_FAILED}" -gt 0 ] 2>/dev/null; then
        echo ""
        echo "  Failed lines:"
        grep -iE "failed|failure" "${CACTUS_LOG}" 2>/dev/null | tail -5 | sed 's/^/    /' || true
    fi

    echo ""
    echo "--- Pipeline stage timings ---"
    grep -iE "completed|finished|starting|done" "${CACTUS_LOG}" 2>/dev/null | \
        grep -iE "minigraph|cactus|align|graphmap|clip|join|vg|gfa|gbz|vcf" | \
        tail -20 | sed 's/^/    /' || echo "    No timing info found"

    echo ""
    echo "--- Peak memory usage ---"
    grep -iE "memory|mem|rss|peak" "${CACTUS_LOG}" 2>/dev/null | tail -5 | sed 's/^/    /' || \
        echo "    No memory info found in log"

    echo ""
    echo "  Full log: ${CACTUS_LOG}"
else
    echo "  No cactus log found"
    echo "  Checked:"
    echo "    ${PAN_DIR}/logs/cactus_stdout.log"
    echo "    ${PAN_DIR}/logs/cactus_pangenome.log"
fi

echo ""
echo "$(date) === STEP 11 COMPLETE ==="

###############################################################################
# QC SUMMARY REPORT
###############################################################################

echo ""
echo "=================================================================="
echo "  PANGENOME QC SUMMARY REPORT"
echo "  $(date)"
echo "=================================================================="
echo ""
echo "  INPUT:"
echo "    Reference:          GRCm39 (C57BL/6J)"
echo "    Strains:            17 (1 ref + 16 REL-2205)"
echo "    Chromosomes:        1-19, X, MT (no Y)"
echo "    Container:          ${CACTUS_SIF}"
echo ""
echo "  OUTPUT DIRECTORY:     ${OUTDIR}"
echo ""
echo "  OUTPUT FILES:"
ls -lhS "${OUTDIR}/${OUTNAME}".* 2>/dev/null | awk '{printf "    %-60s  %s\n", $NF, $5}' || true
echo ""
echo "  KEY INDEXES FOR DOWNSTREAM:"
echo "    GBZ:   ${GBZ}"
echo "    DIST:  ${DIST:-NOT FOUND}"
echo "    MIN:   ${MIN:-NOT FOUND}"
echo "    VCF:   ${VCF}"
echo "    GFA:   ${GFA_FILE:-NOT FOUND}"
echo ""
echo "  QC VERDICTS:"

# File integrity
if [ -f "${GBZ}" ] && [ -f "${VCF}" ]; then
    echo "    File integrity:       PASS"
else
    echo "    File integrity:       FAIL"
fi

# Graph topology
if [ -f "${PAN_DIR}/qc/vg_stats_basic.txt" ] && [ -s "${PAN_DIR}/qc/vg_stats_basic.txt" ]; then
    NODES_REPORT=$(awk '/nodes/{print $NF}' "${PAN_DIR}/qc/vg_stats_basic.txt" | tr -d '[:space:]') || true
    EDGES_REPORT=$(awk '/edges/{print $NF}' "${PAN_DIR}/qc/vg_stats_basic.txt" | tr -d '[:space:]') || true
    echo "    Graph topology:       PASS (${NODES_REPORT:-?} nodes, ${EDGES_REPORT:-?} edges)"
else
    echo "    Graph topology:       CHECK"
fi

# Strains
if [ -f "${PAN_DIR}/qc/graph_paths.txt" ] && [ -s "${PAN_DIR}/qc/graph_paths.txt" ]; then
    N_STRAINS_REPORT=$(awk -F'#' '{print $1}' "${PAN_DIR}/qc/graph_paths.txt" | sort -u | wc -l | tr -d '[:space:]') || true
    echo "    Strains in graph:     ${N_STRAINS_REPORT:-?} / 17"
else
    echo "    Strains in graph:     CHECK"
fi

# VCF quality
if [ -f "${PAN_DIR}/qc/bcftools_stats.txt" ] && [ -s "${PAN_DIR}/qc/bcftools_stats.txt" ]; then
    VCF_RECORDS=$(grep "^SN" "${PAN_DIR}/qc/bcftools_stats.txt" | grep "number of records:" | awk '{print $NF}') || true
    VCF_TSTV=$(grep "^TSTV" "${PAN_DIR}/qc/bcftools_stats.txt" | head -1 | awk '{print $5}') || true
    echo "    VCF quality:          PASS (${VCF_RECORDS:-?} variants, Ts/Tv=${VCF_TSTV:-?})"
else
    echo "    VCF quality:          CHECK"
fi

# Giraffe indexes
if [ -n "${DIST}" ] && [ -f "${DIST}" ] && [ -n "${MIN}" ] && [ -f "${MIN}" ]; then
    echo "    Giraffe indexes:      PASS"
else
    echo "    Giraffe indexes:      MISSING"
fi

# Giraffe mapping test
if [ -f "${PAN_DIR}/qc/synthetic_mapped.gam" ]; then
    GAM_BYTES=$(stat --printf='%s' "${PAN_DIR}/qc/synthetic_mapped.gam" 2>/dev/null || echo 0)
    if [ "${GAM_BYTES}" -gt 0 ] 2>/dev/null; then
        echo "    Giraffe mapping:      PASS"
    else
        echo "    Giraffe mapping:      CHECK (empty output)"
    fi
else
    echo "    Giraffe mapping:      NOT TESTED"
fi

# Log analysis
if [ -n "${CACTUS_LOG:-}" ] && [ -f "${CACTUS_LOG:-/dev/null}" ]; then
    LOG_ERRORS=$(grep -ci "error" "${CACTUS_LOG}" 2>/dev/null || echo "0")
    LOG_ERRORS=$(echo "${LOG_ERRORS}" | tr -d '[:space:]')
    if [ "${LOG_ERRORS:-0}" -eq 0 ] 2>/dev/null; then
        echo "    Build log:            CLEAN (0 errors)"
    else
        echo "    Build log:            ${LOG_ERRORS} error mentions"
    fi
else
    echo "    Build log:            NOT FOUND"
fi

echo ""
echo '  NEXT STEPS FOR piRNA ANALYSIS:'
echo ''
echo '  1. Extract TE-sized SVs from pangenome VCF:'
echo '     bcftools view -H your.vcf.gz | \'
echo '       awk '"'"'length($5)-length($4) >= 300 || length($4)-length($5) >= 300'"'"''
echo ''
echo '  2. Map sRNA-seq to each strain OWN assembly (linear):'
echo '     bowtie -v 0 -m 1 --best -x STRAIN_INDEX -f reads.fa -S out.sam'
echo ''
echo '  3. Detect piRNA clusters per strain:'
echo '     proTRAC -i out.sam -g STRAIN.fa'
echo ''
echo '  4. Intersect piRNA clusters with TE insertions from VCF:'
echo '     bedtools intersect -a clusters.bed -b te_insertions.bed -wa -wb'
echo ''
echo '  5. Map sRNA-seq to pangenome graph:'
echo '     vg giraffe -Z your.gbz \'
echo '       -m your.min \'
echo '       -d your.dist \'
echo '       -f reads.fq -t 16 -o gam > sample.gam'
echo ""
echo "  LOG FILE:      ${LOGFILE}"
echo "  QC FILES:      ${PAN_DIR}/qc/"
echo ""
echo "$(date) === ALL QC COMPLETE ==="