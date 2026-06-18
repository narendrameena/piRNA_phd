# Fig_pav_locus_newC2cd4b_chr9 — C2cd4b region (chr9, 7/16)

**Locus:** C2cd4b region (chr9, 7/16)  ·  GRCm39 chr9:67680073–67700208  ·  slug `newC2cd4b_chr9`
**Variant:** main — top-expressing strain shown in detail (Panel B) + base-resolution Panel C
**Shows:** one PICB piRNA cluster at this locus across the 16-strain pangenome × 3 timepoints (E16.5→P12.5→P20.5).
 Panel A = per-strain × timepoint cluster FPM (log) + genome PAV (● present / ○ absent); Panel B = per-strain
 per-timepoint primary-sRNA coverage + the TE-family×strand top bar and the antisense-to-TE **silencing** second bar
 (length = %on-TE); Panel C = base resolution (5′-U = 1U; red 5′ arrow = antisense-to-TE silencing).

**How generated** — full raw-FASTQ→figure pipeline in [`../PIPELINE.md`](../PIPELINE.md). In brief:
strain sRNA-seq FASTQ → cutadapt (20–36 nt, adapter TGGAATTCTCGGGTGCCAAGG) → STAR (exact, multimap-permissive,
EndToEnd) on each strain's UNMASKED REL-2205 genome → PICB cluster calling → pangenome projection (GRCm39) →
RepeatMasker TE families + halLiftover genome PAV → `fetch_primary` (primary 24–32 nt reads in the cluster:
strand, sense/antisense-to-TE, 1U, per-family) → this figure.

- **Code:** `code/make_pav_locus.py` + `code/pav_clusters.py`  (driver `code/render_catalogue_updated.sh`)
- **Command:** `make_pav_locus.py 9 67680073 67700208 "C2cd4b region (chr9, 7/16)" _ Fig_pav_locus_newC2cd4b_chr9`
- **Data tables:** `data/picb_pangenome_clusters.tsv` (cluster projection, all strains), `data/locus_genome_pav.tsv`
  (genome presence/absence), `data/catalogue_loci.tsv` (this locus)
- **Raw input:** `results/STAR_srna_strain_wise/{strain}/{strain}-{tp}.{rep}/Aligned.sortedByCoord.out.bam`
- **Formats:** PDF + SVG + PNG.
