#!/usr/bin/env python3
"""
Re-run only the SV extraction step with corrected bcftools norm -m - to properly
handle multi-allelic sites. Then regenerate Fig3.
Bug: original awk used length($5) on the entire comma-joined multi-allelic ALT
field, inflating apparent INS count and suppressing DEL classification.
Fix: bcftools norm -m - decomposes multi-allelic sites into biallelic records
so each REF/ALT comparison is for the specific allele C57BL_6NJ carries.
"""
import os, subprocess
import pandas as pd

OUT_DIR = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome"
PAN_VCF = "/mnt/beegfs/scratch/miska/nm667/inProgress/mice_PiRNA/results/pangenome/output/mouse_17strain_pangenome.vcf.gz"
vcf_private_bed = os.path.join(OUT_DIR, "C57BL_6NJ_private_sv_GRCm39coords.bed")

print("=== Re-running SV step with bcftools norm -m - fix ===")
print(f"VCF: {PAN_VCF}")

sv_cmd = f"""
export PATH="/mnt/home3/miska/nm667/miniconda3/bin:$PATH"
bcftools view -s C57BL_6NJ \\
  "{PAN_VCF}" 2>/dev/null | \\
bcftools norm -m - 2>/dev/null | \\
bcftools view -c 1 -e 'GT[0]="0/0" || GT[0]="0|0" || GT[0]="./."' 2>/dev/null | \\
awk '/^#/{{next}} {{
  ref=length($4); alt=length($5)
  diff=(alt>ref) ? alt-ref : ref-alt
  if(diff>=300) {{
    sv_type=(alt>ref)?"INS":"DEL"
    printf "%s\\t%d\\t%d\\t%s\\t%d\\t%s\\n",$1,$2-1,$2+length($4)-1,$3,diff,sv_type
  }}
}}' > {vcf_private_bed}
echo "Done"
"""

print("Extracting C57BL_6NJ SVs (this may take several minutes)...")
result = subprocess.run(sv_cmd, shell=True, capture_output=True, text=True, timeout=900)
print(f"  {result.stdout.strip()}")
if result.returncode != 0:
    print(f"  STDERR: {result.stderr[:500]}")

if not os.path.exists(vcf_private_bed) or os.path.getsize(vcf_private_bed) == 0:
    print("ERROR: no output produced")
    exit(1)

sv_df = pd.read_csv(vcf_private_bed, sep="\t", header=None,
                    names=["chr","start","end","id","sv_size_bp","sv_type"])

print(f"\nC57BL_6NJ SVs (≥300bp) with non-ref genotype: {len(sv_df)}")
print(f"Insertions (INS): {(sv_df['sv_type']=='INS').sum()}")
print(f"Deletions  (DEL): {(sv_df['sv_type']=='DEL').sum()}")
print(f"Median SV size:   {sv_df['sv_size_bp'].median():.0f} bp")
print(f"Max SV size:      {sv_df['sv_size_bp'].max():,} bp")

print("\nSize distribution:")
bins = [(300,500,"300-499bp"), (500,1000,"0.5-1kb"), (1000,5000,"1-5kb"),
        (5000,10000,"5-10kb"), (10000,100000,"10-100kb"), (100000,10**9,">100kb")]
for lo, hi, label in bins:
    n = ((sv_df["sv_size_bp"] >= lo) & (sv_df["sv_size_bp"] < hi)).sum()
    print(f"  {label:15s}: {n:6d}")

sv_df.to_csv(os.path.join(OUT_DIR, "C57BL_6NJ_TE_sized_SVs_GRCm39.csv"), index=False)
print(f"\nSaved updated CSV.")
print("=== DONE — now run pangenome_figures.py to regenerate Fig3 ===")
