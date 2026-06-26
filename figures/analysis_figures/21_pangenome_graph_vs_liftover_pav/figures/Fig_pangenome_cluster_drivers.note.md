# Fig_pangenome_cluster_drivers

**Drivers of strain-variable piRNA clusters: TE insertions (genetic novelty), gene context, antisense pseudogene fragments, and piC-DoG readthrough**

- **Shows:** A TE families at the genetically-variable (liftover-dispensable/private) loci; B developmental class × gene context; C antisense vs sense pseudogene fragments in pachytene clusters; D piC-DoG readthrough per strain (whole-testis bulk RNA).
- **Result:** the genetically-novel minority (clusters not conserved across strains) arises at **young TE insertions** — **L1MdTf/L1MdA LINE-1 + IAP LTR** dominate. **Pachytene** clusters are **intergenic** (A-MYB-driven); **pre-pachytene/hybrid** clusters are **3′UTR-enriched** (mRNA-derived). **1,285 antisense** pseudogene fragments sit in pachytene clusters (candidate *trans*-acting mRNA-silencing piRNA sources, Loubalova 2025). **~65 %** of downstream clusters show whole-testis-RNA **readthrough** from the upstream gene's 3′ end (piC-DoG; Konstantinidou 2024) — consistent across wild and classical strains.
- **Why trustworthy:** SV/TE from the **precomputed deconstructed pangenome VCF** (≥300 bp ref/alt length difference) typed by blastn vs the mouse TE consensus; dev-class from per-tp halLiftover clusters (all 3 timepoints incl. fetal E16.5); gene-context vs GRCm39 Ensembl gene models; PGF = blastn of intergenic pachytene-cluster sequence vs GRCm39 CDS; readthrough = `samtools bedcov` of pooled whole-testis bulk RNA (gap coverage ≥10 % of last-exon, cluster expressed).
- **How:** `code/05_sv_te_variable_loci.sh`, `04_devclass_genecontext.sh`, `07_pgf_pachytene.sh`, `06_readthrough.slurm`; `Fig_pangenome_cluster_drivers.py`.
- **Data:** `data/sv_te_hits.tsv`, `master_devclass.tsv`, `master_genecontext.tsv`, `pgf_hits.tsv`, `readthrough_*.tsv`.

Full pipeline: [`PIPELINE.md`](../PIPELINE.md).
