# aquifer-defensome

Code accompanying **"The defensome of prokaryotes in aquifers"**
(Li, P. *et al.*, *Nature Communications* **16**, 6482 (2025)).

The repository contains every script used in the study: the upstream bioinformatic
pipeline that goes from raw metagenomic reads to the 27,578 dereplicated MAGs and their
mobile-genetic-element / defense / anti-defense annotation, and the R analysis scripts
that produce the statistics and the figures.

> 中文版说明见 **[README_zh.md](README_zh.md)**。

---

## What is in here

```
aquifer-defensome/
├── pipelines/          # upstream bioinformatics (Slurm + conda, shell/python)
│   ├── 01_reads_to_MAG_catalog/       QC, assembly, binning, dRep, GTDB-Tk, CheckM2
│   ├── 02_gene_prediction/            Prodigal ORFs, Prokka annotation
│   ├── 03_abundance_and_function/     DIAMOND screening, coverm abundance
│   ├── 04_mobile_genetic_elements/    geNomad, IntegronFinder, ICEfinder2, minced, BLAST
│   ├── 05_defense_systems/            DefenseFinder 2.0.0, DefensePredictor
│   ├── 06_phylogeny_and_traits/       anti-CRISPR phylogeny, OGT prediction
│   └── 07_table_postprocessing/       defense-island and integron tables (python/awk)
├── analyses/           # downstream statistics and figures (R)
│   ├── 01_main_figures/               Fig01 - Fig06
│   └── 02_supplementary_figures/      FigS02_B - FigS18
└── docs/               # workflow, layout, file mapping, data and software notes
```

Every R script is named after the figure it produces, keeps that mapping in its header,
and is numbered in manuscript order. A one-to-one table of script -> figure -> original
file name is in [docs/03_file_mapping.md](docs/03_file_mapping.md).

## Documentation

| Document | Content |
| --- | --- |
| [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md) | End-to-end workflow, per-stage inputs and outputs, run order and inter-script dependencies |
| [docs/02_directory_layout.md](docs/02_directory_layout.md) | What lives where, and the naming conventions |
| [docs/03_file_mapping.md](docs/03_file_mapping.md) | Every script: previous path -> new path, purpose, plus the one code correction |
| [docs/04_data_availability.md](docs/04_data_availability.md) | Working directory and input tables of every R script |
| [docs/05_software_and_environment.md](docs/05_software_and_environment.md) | Cluster setup, conda environments, R packages |

## Reproducing the analyses

**1. Upstream pipelines.** The scripts under `pipelines/` are the Slurm job scripts used
for the production runs (`#SBATCH` headers included). The site-specific absolute paths
have been replaced by variables: export `LUSTRE_HOME`, `HOME_ROOT`, `DB_ROOT` and, where
needed, `SHARED_HOME`, adapt the queue names, then submit in the order given in
[docs/01](docs/01_pipeline_overview.md).

**2. R analyses.** The scripts require R >= 4.2 and the packages listed in
[docs/05](docs/05_software_and_environment.md). Input tables are read relative to
`GTR_DATA`, which defaults to the current working directory:

```bash
cd /path/to/data && Rscript /path/to/aquifer-defensome/analyses/01_main_figures/Fig01.R
# or
Rscript -e 'Sys.setenv(GTR_DATA = "/path/to/data"); source("analyses/01_main_figures/Fig01.R")'
```

Four scripts additionally read from `GTR_REV`, the working copy of the peer-review
revision analyses.

### Run order matters for two scripts

Most scripts are independent, but two of them also build objects that other scripts
consume, so they must be run first:

| Run this | before | because it defines |
| --- | --- | --- |
| `analyses/02_supplementary_figures/FigS10.R` | `Fig02.R`, `Fig03.R`, `Fig05_A-D.R`, `FigS02_B.R`, `FigS07.R` | `ds_density_tax` (and `md_p_val1` for Fig02) |
| `analyses/02_supplementary_figures/FigS17-18.R` | `Fig06.R` | `prophage_anti_coloca_4` and the gene-map tables |

The full table is at the end of [docs/01](docs/01_pipeline_overview.md).

## Data availability

The repository is code-only. The MAG catalog and the annotation tables are provided as
Supplementary Data of the paper; the raw reads are available from the NCBI SRA under the
accession given in the paper. The directory layout the scripts expect is documented in
[docs/04_data_availability.md](docs/04_data_availability.md).

## Note on the previous layout

The previous version of this repository was two loose folders (`R codes/`,
`Shell scripts/`). It has been reorganized into the layered structure above: no script was
dropped, all 23 R scripts and 27 pipeline scripts are present, paths were made portable,
and one script was repaired (a full-width colon that prevented `Fig.S17-18.R` from
parsing). No absolute path, user home directory or account name remains in the
repository -- everything is handled by the `GTR_DATA`, `GTR_REV`, `LUSTRE_HOME`,
`HOME_ROOT`, `DB_ROOT`, `SHARED_HOME`, `GTR_DI_DIR` and `GTR_INTEGRON_DIR` variables.

## Citation

Li, P. *et al.* The defensome of prokaryotes in aquifers. *Nat. Commun.* **16**, 6482
(2025).
