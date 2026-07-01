#!/bin/bash

BASE=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA
REL_DIR=${BASE}/resources/REL-2205-Assembly
GENOME_DIR=${BASE}/resources/black6/genome

echo "=== T2T Assemblies ==="
for f in C57BL_6J_T2T.fa CAST_EiJ_T2T.fa; do
    FULL="${GENOME_DIR}/${f}"
    if [ -f "$FULL" ]; then
        SIZE=$(ls -lh "$FULL" | awk '{print \$5}')
        FIRST=$(head -1 "$FULL")
        echo "  ✓ ${f}  (${SIZE})  header: ${FIRST:0:60}"
    else
        echo "  ✗ ${f}  NOT FOUND at: ${FULL}"
    fi
done

echo ""
echo "=== REL-2205 Assemblies ==="
for strain in C57BL_6NJ A_J AKR_J BALB_cJ C3H_HeJ CBA_J DBA_2J \
              FVB_NJ LP_J NOD_ShiLtJ NZO_HlLtJ 129S1_SvImJ \
              PWK_PhJ WSB_EiJ SPRET_EiJ; do
    FULL="${REL_DIR}/${strain}_chromosomes_MT.fasta"
    if [ -f "$FULL" ]; then
        SIZE=$(ls -lh "$FULL" | awk '{print \$5}')
        FIRST=$(head -1 "$FULL")
        echo "  ✓ ${strain}  (${SIZE})  header: ${FIRST:0:60}"
    else
        echo "  ✗ ${strain}  NOT FOUND at: ${FULL}"
    fi
done

echo ""
echo "=== ALL files in REL-2205 directory ==="
ls -lh ${REL_DIR}/*chromosomes_MT.fasta 2>/dev/null | awk '{print "  ", $NF, \$5}'

echo ""
echo "=== ALL files in genome directory ==="
ls -lh ${GENOME_DIR}/*.fa 2>/dev/null | awk '{print "  ", $NF, \$5}'

echo ""
echo "=== Installed tools ==="
for tool in samtools bcftools bedtools vg cactus-pangenome minigraph singularity apptainer; do
    if command -v $tool &>/dev/null; then
        VER=$($tool --version 2>&1 | head -1 || echo "found")
        echo "  ✓ $tool: $VER"
    else
        echo "  ✗ $tool: NOT INSTALLED"
    fi
done

echo ""
echo "=== Conda environments ==="
conda env list 2>/dev/null || echo "  conda not available"

echo ""
echo "=== Available HPC modules (cactus/vg) ==="
module avail 2>&1 | grep -iE "cactus|vg|singularity|apptainer" || echo "  none found"

echo ""
echo "=== Disk space ==="
df -h ${BASE} | awk 'NR==2{print "  Available:", \$4, "of", \$2}'
