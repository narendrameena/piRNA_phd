# Fig_snp_allele_test

**TEST — the exact-'unique' SNP-alleles are standing variation (shared ±1–3 SNP with ~11 strains), not innovation**

- **Shows:** for the 4,394 piRNAs that the EXACT definition calls unique but SNP-aware rejects: A how many OTHER strains express a 1–3 mm allele (same stage); B mismatch distance to the nearest expressed allele; C exact-unique composition (clean innovation vs SNP-allele standing variation).
- **Result:** these SNP-alleles are shared (as a 1–3 mm allele) with a **mean of 11.1 of the 15 other strains** — the distribution **peaks at all 15** — and **88 % differ by a single SNP**. So they are SNPs in widely-shared / conserved piRNAs (STANDING GENETIC VARIATION), not strain-specific INNOVATION. They make up **29 %** of the exact-sequence "unique" set; the remaining 71 % is clean innovation (= the SNP-aware set).
- **Why trustworthy:** uses the per-tp 1–3 mm SNP-allele records (`snp_variant_refinement_withintp.csv`, classify_step416_pertp.py) restricted to the 4,394; sharing = distinct `variant_strain` per candidate. **Interpretation BioMNI 3/3-verified**: SNP-aware definition preferred (separates innovation from standing variation).
- **How:** `code/Fig_snp_allele_test.py`.
- **Data:** `data/SourceData_Fig_snp_allele_test.csv.gz`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
