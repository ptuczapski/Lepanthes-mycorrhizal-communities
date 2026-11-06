# Fungal communities associated with sympatric *Lepanthes* orchids

Reproducible analysis scripts for a study of root-associated fungal communities in four sympatric epiphytic *Lepanthes* orchids sampled in the Monteverde cloud forest, Costa Rica.

**Repository status:** Analysis scripts accompanying a manuscript intended for *Ecology and Evolution*.

## Study overview

The study examines fungal community composition and diversity in roots of *Lepanthes cribbii*, *L. falx-bellica*, *L. mentosa*, and *L. monteverdensis*. Full-length fungal ITS amplicons were sequenced using PacBio technology. The workflow includes ASV inference and taxonomic assignment, fungal guild annotation, construction of two analysis datasets, phylogenetic diversity, ordination, multivariate statistical testing, ASV accumulation curves, taxon prevalence, and an phylogenetic sensitivity analysis.

Two datasets are used throughout:

- **GF — General Fungi:** The broader root-associated fungal ASV dataset.
- **CMF — Classified/putative mycorrhizal fungi:** A subset selected using fungal guild annotations.

The scripts document the analyses originally performed. They are not a single-command, fully automated pipeline: some steps depend on intermediate files, cluster software, and manual table processing.

## Repository layout

```text
.
├── README.md
├── CITATION.cff
├── LICENSE                 # Add the selected license
├── .gitignore              # Add before uploading data
└── scripts/
    ├── 1.dada2.R
    ├── 2.funguild.and.CMFdataset.R
    ├── 3.dada2toqiime2.sh
    ├── 4.remsingletons.sh
    ├── 5.ASVtableGF_nosingletons.R
    ├── 6.ASVtableCMF_nosingletons.R
    ├── 7.phylogeny.sh
    ├── 8.alphabetadiversity.sh
    ├── 9.GF_analysis_nosingletons.R
    ├── 10.CMF_analysis_nosingletons.R
    ├── 11.Accumulation_curves_calculation.R
    ├── 12.Fungal_taxa_prevalence.R
    └── 13.Gblocks.sh
```

Files with `.sh` extensions in this archive include shell/SLURM commands; consult their contents rather than assuming they are plain-text data files. The script numbers indicate the broad workflow order, but some later scripts are independent or optional analyses.

## Script guide

| Script | Role | Principal inputs / dependencies |
| --- | --- | --- |
| `1.dada2.R` | Process amplicon reads, infer ASVs, remove chimeras, and assign taxonomy | Demultiplexed FASTQ files; DADA2; UNITE reference |
| `2.funguild.and.CMFdataset.R` | Annotate fungal guilds and derive the mycorrhizal subset | ASV table, taxonomy, representative sequences; FUNGuildR |
| `3.dada2toqiime2.sh` | Convert and import DADA2-derived data into QIIME 2 | ASV/BIOM tables, representative sequences, metadata |
| `4.remsingletons.sh` | Filter features by prevalence within each orchid-species subset | QIIME 2 feature tables and sample metadata |
| `5.ASVtableGF_nosingletons.R` | Assemble the GF table after species-level feature filtering | Exported GF feature tables; manual intermediate-table steps |
| `6.ASVtableCMF_nosingletons.R` | Assemble the CMF table after species-level feature filtering | Exported CMF feature tables; manual intermediate-table steps |
| `7.phylogeny.sh` | Align ASV sequences and reconstruct/filter a rooted phylogeny | Representative sequences and QIIME 2 feature table |
| `8.alphabetadiversity.sh` | Alpha rarefaction, phylogenetic diversity, and QIIME 2 beta-diversity analyses | GF/CMF tables, rooted tree, metadata |
| `9.GF_analysis_nosingletons.R` | GF ordination, PERMDISP, and PERMANOVA | GF table, taxonomy, tree, metadata |
| `10.CMF_analysis_nosingletons.R` | CMF ordination, PERMDISP, and PERMANOVA | CMF table, taxonomy, tree, metadata |
| `11.Accumulation_curves_calculation.R` | Species-specific ASV accumulation with randomized sample orders | GF or CMF table and metadata |
| `12.Fungal_taxa_prevalence.R` | Fungal group prevalence and ASV richness | GF table, taxonomy, metadata, FUNGuild annotations, tree |
| `13.Gblocks.sh` | Sensitivity analysis using Gblocks-filtered ITS alignment and phylogenetic diversity | Representative sequences, GF feature table, metadata; Gblocks and QIIME 2 |

## Important analysis choices

- **Feature filtering:** Script 4 filters ASVs by occurrence across samples within orchid-species subsets; Consult the actual QIIME 2 commands for the precise rule.
- **Rarefaction:** The downstream GF and CMF analyses use sampling depths of **679** and **445** reads per sample, respectively, where specified.
- **Community dissimilarity:** Analyses include presence/absence Jaccard and unweighted UniFrac distances.
- **Ordination:** PCoA and NMDS are used to visualize community composition. NMDS dimensionality differs by dataset and distance metric; see Scripts 9 and 10.
- **Statistical tests:** PERMDISP and PERMANOVA evaluate dispersion and group effects. The R PERMANOVA models use **9,999 permutations** where specified.
- **Accumulation:** Script 11 uses **1,000 randomized sample orders** to summarize cumulative ASV richness; percentile envelopes describe sample-order variation.

## Software

The workflow uses R version 4.3.2 and the packages `dada2`, `phyloseq`, `Biostrings`, `vegan`, `ggplot2`, `tidyverse`, `data.table`, and `FUNGuildR`, as applicable; it also uses QIIME 2, MAFFT, FastTree, Gblocks, and cluster shell/SLURM commands.

## Input data and availability

The 13 scripts alone do not contain all data needed for full reproduction. Relevant inputs include:

1. Demultiplexed sequencing reads and the corresponding sequence archive accession: **BioProject PRJNA1345165, accession numbers: SRX31179501 - SRX31179909**.
2. Sample metadata (`qiime2metadata.tsv`) with the fields used in the analyses, feature tables and ASV FASTA files: **: https://doi.org/10.71927/uga.27592**.

## Running the analyses

1. Obtain the input data and install the software versions used in the original study.
2. Inspect each script's working-directory assumptions and update paths for the local or HPC environment.
3. Run the preprocessing and guild-classification scripts (1–2).
4. Import data and perform feature filtering in QIIME 2 (3–4).
5. Reconstruct the combined GF and CMF tables (5–6), including the documented manual intermediate-table operations where applicable.
6. Reconstruct the phylogeny and calculate diversity metrics (7–8).
7. Run the dataset-specific analyses (9–10) and the additional accumulation/prevalence analyses (11–12).
8. Run the Gblocks sensitivity analysis (13) if reproducing that supplementary comparison.

## Reproducibility and provenance

These scripts are intended to document the analyses performed for the associated manuscript. Annotations or path adjustments made for public release should not be interpreted as evidence that every script has been rerun in a clean environment. When reproducing the published results, consult the manuscript methods and original parameter settings, and record any deviations.

## Citation

If you use this code, cite the associated manuscript **[ADD MANUSCRIPT CITATION]** and the archived software release **[ADD DOI WHEN AVAILABLE]**. See (CITATION.cff) for machine-readable citation metadata.

## License

**[ADD LICENSE FILE AND UPDATE THIS SECTION]**. Until a license is added, public availability on GitHub does not automatically grant permission to reuse or redistribute the code.

## Contact

For questions about the analysis, open an issue in this repository **(https://github.com/ptuczapski/Lepanthes-mycorrhizal-communities/)**.
