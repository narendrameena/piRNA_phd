#!/bin/bash
#SBATCH --job-name=panins
#SBATCH --partition=2204,1804,2004
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te/extract.log
set -euo pipefail
BCF=/mnt/home3/miska/nm667/miniconda3/envs/rna2prot_snv/bin/bcftools
PY=/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/python
V=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/pangenome/output/mouse_17strain_pangenome.vcf.gz
U=/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna
# subset to 3 pilot strains, trim alleles not carried by them, keep sites variant among the 3, drop
# SNP-only sites; stream the (CHROM POS REF ALT[multi] GT...) to the private-insertion extractor.
"$BCF" view -s C57BL_6NJ,CAST_EiJ,SPRET_EiJ -a --min-ac 1 -V snps "$V" \
 | "$BCF" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
 | "$PY" "$U/parse_insertions.py"
echo "[done] $(date)"
