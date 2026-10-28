# Workflow overview

The study is organized in two layers. The **upstream layer** (`pipelines/`) turns raw
metagenomic reads into the MAG catalog and the annotation tables with Slurm job scripts.
The **downstream layer** (`analyses/`) turns those tables into the statistics and figures
of the paper with R. The layers communicate only through files, so any step can be re-run
on its own.

```
   raw reads (NCBI SRA)
        |
   [1] reads -> MAG catalog -------------------------------------------------+
       QC, assembly, binning, dRep, GTDB-Tk, CheckM2                          |
        |                                                                      |
   [2] gene prediction ----+                                                  |
       Prodigal ORFs        |   [3] abundance and function ---+               |
       Prokka (1,626 MAGs)  |       DIAMOND ARG/MRG/VFG       |               |
        |                   |       coverm MAG abundance      |               |
        |                   |                                |               |
   [4] mobile genetic elements                    [5] defense systems          |
       geNomad (phage/plasmid)                        DefenseFinder 2.0.0     |
       IntegronFinder, ICEfinder2                     DefensePredictor        |
       minced + spacer BLAST                                                  |
        |                                                                      |
   [6] phylogeny and traits ---+                                             |
       anti-CRISPR tree         |   [7] table postprocessing ---+              |
       OGT prediction           |       defense islands         |              |
        |                       |       integron DG tables      |              |
        v                       v                              v              v
   [8] R analyses: analyses/01_main_figures, analyses/02_supplementary_figures
       statistics, figures, source-data tables
```

## Stage 1 — reads to MAG catalog

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 1.1 | `01_quality_control.sh` | metaWRAP (`read_qc`) | per-sample cleaned reads |
| 1.2 | `02_assembly.sh` | metaWRAP / megahit | contigs (1 kb and 0.5 kb sets) |
| 1.3 | `03_binning.sh` | metaWRAP (metabat2, maxbin2, concoct) | initial bins + bin refinement |
| 1.4 | `04_dereplicate_bins.sh` | dRep + CheckM | dereplicated MAGs (98% ANI, 30% coverage) |
| 1.5 | `05_gtdbtk_taxonomy.sh` | GTDB-Tk r207 | `gtdbtk.*.summary.tsv` |
| 1.6 | `06_checkm2_quality.sh` | CheckM2 | completeness / contamination |

The dereplicated catalog contains **27,578 MAGs**; the 1,626 MAGs that pass the quality
threshold are used where per-MAG annotations are needed (Prokka, DefensePredictor).

## Stage 2 — gene prediction

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 2.1 | `07_prodigal_orf.sh` | Prodigal 2.6.3 (meta) | ORF nucleotide and protein sets |
| 2.2 | `22_prokka_annotate.sh` | Prokka | GFF/FFN/FAA of the 1,626 high-quality MAGs |

## Stage 3 — abundance and functional screening

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 3.1 | `08_diamond_screen.sh` | DIAMOND | hits against ARG (CARD), MRG, mibPOP and VFG databases |
| 3.2 | `09_coverm_abundance.sh` | bowtie2 + coverm | per-sample MAG abundance (RPKM/TPM) |

## Stage 4 — mobile genetic elements

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 4.1 | `10_genomad.sh` | geNomad | plasmids, viruses/proviruses in the MAG contigs |
| 4.2 | `11_integron_finder.sh` | IntegronFinder | integron calls |
| 4.3 | `12_icefinder.sh` | ICEfinder2 | ICE/IME calls |
| 4.4 | `30_minced_crispr.sh` | minced | CRISPR arrays and spacers |
| 4.5 | `31_blast_spacers.sh` | BLAST+ 2.11 | spacer -> prophage hits (historical phage-host links) |

## Stage 5 — defense systems

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 5.1 | `20_defensefinder.sh` | DefenseFinder 2.0.0 (`--preserve-raw`) | `merged_mags_*_defense_finder_systems.tsv` and `..._genes.tsv`, plus the `DIs/` output |
| 5.2 | `21_defense_predictor.sh` | DefensePredictor | per-MAG defense-gene predictions |

`merged_mags_*_defense_finder_systems.tsv` is the table every defensome analysis reuses
(rows with `activity == "Defense"` / `"Antidefense"`).

## Stage 6 — phylogeny and traits

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 6.1 | `40_acr_phylogeny.sh` | seqkit, muscle, trimAl, FastTree | `all_Acr.faa`, `503acr_aln.trim.afa`, tree |
| 6.2 | `50_ogt_prediction.sh` | OGT_prediction | predicted growth temperature per MAG |

## Stage 7 — table post-processing

```
DIs/*.tsv                        per-batch DefenseFinder island output
   | 01_collect_DI_hits.py        keep Defense rows, standardise contig IDs
   v
defense_islands_summary4.csv
   | 02_DI_summary.py             per-island summary      (defense_islands.txt)
   | 03_DI_nondG_orf_ids.py       ORF IDs of each island
   | 04_DI_nondG_filter.sh        subtract DG IDs, seqkit extraction
   v
01all_DIorf.txt, 03filtered_01all_DIorf.txt, 04non_DG_orf.faa

integron DefenseFinder output
   | 10_integron_DG_extract.sh    pull the DS rows of integron contigs
   | 11_DS_filter.py              clean the hit_id -> genome mapping
   | 12_DS_bind.py                bind and count defense genes per genome
   v
all_integron_dss.txt, all_integron_dss2.txt, per-genome counts

DefensePredictor output
   | 20_DP_filter.sh              merge, keep log-odds >= 7.2, drop sequence columns
   v
1626mag_defensepredictor_processed.csv
```

## Stage 8 — R analyses

Main figures:

| Script | Manuscript | Content |
| --- | --- | --- |
| `Fig01.R` | Fig. 1 | DS prevalence in 27,578 MAGs, DG novelty, DS number/density, effector-type treemap |
| `Fig02.R` | Fig. 2 | DS number in CPR vs other phyla vs Pseudomonadota, phylum heatmap |
| `Fig03.R` | Fig. 3 | VHR/DHR, MGE load and DS number, DPANN and CPR clades |
| `Fig04.R` | Fig. 4 | ARG content vs DS number (incl. MazEF, pan-immune panel) |
| `Fig05_A-D.R` | Fig. 5a-d | MGE counts, phylum sankey, plasmid:prophage ratio, DS density on MGEs, O/E and OR |
| `Fig05_E-I.R` | Fig. 5e-i | defense islands: gene number, DG proportion, co-localization network, Fisher tests |
| `Fig06.R` | Fig. 6 | figure assembly fragment (gene maps drawn from `FigS17-18.R` objects) |

Supplementary figures:

| Script | Manuscript | Content |
| --- | --- | --- |
| `FigS02_B.R` | Fig. S2b | MAG abundance by niche breadth, stacked composition per zone |
| `FigS02_C.R` | Fig. S2c | DS density by geological setting, well type, confined/unconfined, depth |
| `FigS02_D.R` | Fig. S2d | upset plot of defense families across the seven zones |
| `FigS03.R` | Fig. S3 | DG number vs genome size, DG distribution, percentage, density |
| `FigS04.R` | Fig. S4 | archaeal defensome, DPANN comparisons, DG occurrence |
| `FigS06.R` | Fig. S6 | DefenseFinder vs DefensePredictor |
| `FigS07.R` | Fig. S7 | prophage number vs MAG size, VHR vs size, supercarriers |
| `FigS08.R` | Fig. S8 | same analysis as S2_C plus the source-data export |
| `FigS09.R` | Fig. S9 | plaque assay (wet-lab validation) |
| `FigS10.R` | Fig. S10 | phylum/genus DS heatmaps and CPR clades; **also the data-prep script** |
| `FigS12.R` | Fig. S12 | DS distribution in selected phyla |
| `FigS13.R` | Fig. S13 | MRG vs DS number |
| `FigS14.R` | Fig. S14 | VFG vs DS number |
| `FigS15.R` | Fig. S15 | prophage/plasmid counts per phylum, sankey |
| `FigS16.R` | Fig. S16 | COG composition of non-defense genes in defense islands |
| `FigS17-18.R` | Fig. S17, S18 | gene-level maps; **also the data-prep script for Fig06.R** |

No script is provided in the repository for Fig. S1, S2a, S5 and S11 (schematic,
sampling-map and review-figure panels were assembled outside R).

## Inter-script dependencies

R scripts are self-contained except for two data-prep scripts. Run them first:

| Data-prep script | Defines | Needed by |
| --- | --- | --- |
| `FigS10.R` | `ds_density_tax` | `Fig02.R`, `Fig03.R`, `Fig05_A-D.R`, `FigS02_B.R`, `FigS07.R` |
| `FigS10.R` | `md_p_val1` | `Fig02.R` |
| `FigS17-18.R` | `prophage_anti_coloca_4` and the gene-map tables | `Fig06.R` |

In practice: `source("analyses/02_supplementary_figures/FigS10.R")` and
`source("analyses/02_supplementary_figures/FigS17-18.R")` in the same session before
running the dependent scripts.
