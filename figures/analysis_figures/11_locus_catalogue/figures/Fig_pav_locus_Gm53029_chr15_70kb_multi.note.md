# Fig_pav_locus_Gm53029_chr15_70kb_multi — chr15 70kb (Gm53029)

**Locus:** chr15 70kb (Gm53029)  ·  GRCm39 chr15:74231581–74302138  ·  slug `Gm53029_chr15_70kb`
**Variant:** multi — every present strain stacked in Panel B
**Shows:** one PICB piRNA cluster at this locus across the 16-strain pangenome × 3 timepoints (E16.5→P12.5→P20.5).
 Panel A = per-strain × timepoint cluster FPM (log) + genome PAV (● present / ○ absent); Panel B = per-strain
 per-timepoint primary-sRNA coverage + the TE-family×strand top bar and the antisense-to-TE **silencing** second bar
 (length = %on-TE); Panel C = base resolution (5′-U = 1U; red 5′ arrow = antisense-to-TE silencing).

**How generated** — full raw-FASTQ→figure pipeline in [`../PIPELINE.md`](../PIPELINE.md). In brief:
strain sRNA-seq FASTQ → cutadapt (20–36 nt, adapter TGGAATTCTCGGGTGCCAAGG) → STAR (exact, multimap-permissive,
EndToEnd) on each strain's UNMASKED REL-2205 genome → PICB cluster calling → pangenome projection (GRCm39) →
RepeatMasker TE families + halLiftover genome PAV → `fetch_primary` (primary 24–32 nt reads in the cluster:
strand, sense/antisense-to-TE, 1U, per-family) → this figure.

- **Code:** `code/make_pav_locus_multi.py` + `code/pav_clusters.py`  (driver `code/render_catalogue_updated.sh`)
- **Command:** `make_pav_locus_multi.py 15 74231581 74302138 "chr15 70kb (Gm53029)" _ Fig_pav_locus_Gm53029_chr15_70kb_multi`
- **Data tables:** `data/picb_pangenome_clusters.tsv` (cluster projection, all strains), `data/locus_genome_pav.tsv`
  (genome presence/absence), `data/catalogue_loci.tsv` (this locus)
- **Raw input:** `results/STAR_srna_strain_wise/{strain}/{strain}-{tp}.{rep}/Aligned.sortedByCoord.out.bam`
- **Formats:** PDF + SVG + PNG.
