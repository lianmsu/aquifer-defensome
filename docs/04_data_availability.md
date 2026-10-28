# Data availability and input files

The repository contains **code only**. Every table the scripts read is listed below
so that the analyses can be re-run once the data are in place.

Two roots are used:

* `GTR_DATA` — the data directory of the study (default: the current working directory;
  `Sys.setenv(GTR_DATA = "/path/to/defensome")`). Each script `setwd()`s into the
  sub-directory named in the *Working directory* column.
* `GTR_REV` — the working copy of the peer-review revision analyses, used by four
  scripts (default: the current working directory).

The sub-directory names are the ones used during the study (some are in Chinese, e.g.
`02去冗余信息` is the MAG dereplication folder, `01样品信息` the sample metadata folder).
A `<GTR_DATA>` tree with the directories referenced by the scripts:

```
<GTR_DATA>/
├── groundwater/                      main analysis directory of the 27,578 MAGs
│   ├── 01样品信息/                   sample metadata
│   ├── 02去冗余信息/                 dereplicated MAG table (Wdb_ANI99.csv), env.xlsx
│   ├── 02ANI99MAGs_defense/          DefenseFinder output, DIs/
│   ├── 03mag_info/                   per-MAG DS/MGE/taxonomy summary (ds_density_tax*.csv)
│   ├── 01defensefinder_info/         DefenseFinder system definitions
│   ├── 04mag_DI_gene/, 04mag_DI_net/ defense-island annotation and network input
│   ├── 05MGE/                        prophage/plasmid tables
│   ├── 10mag_arg/                    DIAMOND/MRG/VFG heat tables
│   ├── 31gene_circle/                gene maps of representative contigs
│   ├── 4VHR_DHR/                     VHR table
│   ├── s10mag_dg/, s1mag_dg/, s12mag_abundance/   DG counts and MAG abundance
│   ├── archaea/{01MAG_info,01mag_info,02defense}/ archaeal subsection
│   └── 88defense_finder_database/    DefenseFinder database version
└── RefSeq/                           NCBI RefSeq reference set used for comparison

<GTR_REV>/
├── 01defensome/NC_revision/1additional_analyses/{ICEs,defensepredictor,genomesize_phage}
├── 01defensome/NC_revision/2experiments/plaque_results_20250508/   plaque assay (Fig. S9)
└── 03eu_pre/11eukdetect/                                            eukdetect output
```

## Where the data come from

| Data | Source |
| --- | --- |
| Raw metagenomic reads | NCBI SRA, accession given in the paper |
| Assembly, binning, dereplication, taxonomy, quality | regenerate with `pipelines/01_reads_to_MAG_catalog/` |
| Prodigal ORFs, Prokka annotations | regenerate with `pipelines/02_gene_prediction/` |
| DIAMOND functional screening, coverm abundance | regenerate with `pipelines/03_abundance_and_function/` |
| geNomad / IntegronFinder / ICEfinder2 / minced / BLAST | regenerate with `pipelines/04_mobile_genetic_elements/` |
| DefenseFinder / DefensePredictor | regenerate with `pipelines/05_defense_systems/` |
| Acr phylogeny, OGT | regenerate with `pipelines/06_phylogeny_and_traits/` |
| Defense-island and integron tables | regenerate with `pipelines/07_table_postprocessing/` |
| `bac.gtdb.r207.xlsx`, `GW_DATA_lpw.xlsx`, `VHR.xlsx`, `env.xlsx` | curated metadata tables of the study |
| PaCRISPR reference Acr sequences | PaCRISPR database |

## Per-script inputs


### `analyses/01_main_figures/Fig01.R`

**Defensome overview** — Prevalence of defense systems in the 27,578 MAGs, novelty of defense genes, defense-gene distribution across phyla, DS number and DS density, and the effector-type treemap.

* Manuscript: Fig. 1
* Working directory: `groundwater/88defense_finder_database`, `groundwater/02去冗余信息`, `groundwater/02ANI99MAGs_defense`, `RefSeq`, `groundwater/01defensefinder_info`, `groundwater/03mag_info`
* Reads:
    * `uniq_DG.txt`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `df-RefSeq.csv`
    * `DS_compare2.txt`
    * `df_List_Systems_lpw.csv`
    * `genome_stats.csv`
    * `magsize_fenlei.txt`
* Writes: `1_DSs_distribution.pdf`, `DSs_percentage.pdf`, `DSs_density_magsize.pdf`, `DSs_magsize_number2.pdf`, `dg_novel.pdf`, `dg_phyla.pdf`

### `analyses/01_main_figures/Fig02.R`

**Defensome of CPR, Pseudomonadota and other phyla** — DS number across CPR / other phyla / Pseudomonadota, the phylum-level DS frequency heatmap and the per-phylum comparisons.

* Manuscript: Fig. 2
* Working directory: `groundwater/02ANI99MAGs_defense`, `groundwater/02去冗余信息`, `groundwater/01样品信息`
* Reads:
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `bac.gtdb.r207.xlsx`
    * `gtdbtk.bac120.summary.tsv`
* Writes: `3phyla_dss.pdf`

### `analyses/01_main_figures/Fig03.R`

**VHR/DHR and the defensome** — Relationship between virus-to-host ratio (VHR), prophage/plasmid load and DS number, plus DPANN and CPR clade characteristics.

* Manuscript: Fig. 3
* Working directory: `groundwater/4VHR_DHR`, `groundwater/03mag_info`, `groundwater/05MGE`, `groundwater/02去冗余信息`, `groundwater/02ANI99MAGs_defense`, `groundwater/archaea/02defense`
* Reads:
    * `VHR.xlsx`
    * `ds_density_tax_r220.csv`
    * `virus_taxonomy_count.tsv`
    * `plasmid_taxonomy_count.tsv`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
* Writes: `VHR_DHR.pdf`, `3CPR_clade_dss.pdf`, `3_CPR_clade_DSrichness.pdf`, `1_dpann_phyla_DSrichness_10.pdf`, `1_p11_DPANN_dss_liuxinglv.pdf`

### `analyses/01_main_figures/Fig04.R`

**Antibiotic resistance genes and the defensome** — Regression between ARG (sarg) content and DS number, including the MazEF case and the pan-immune ARG panel.

* Manuscript: Fig. 4
* Working directory: `groundwater/10mag_arg`, `groundwater/02ANI99MAGs_defense`, `groundwater/02去冗余信息`
* Reads:
    * `sarg_gene_out1_type.csv`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `dg_counts_matrix.csv`
    * `gtdbtk.bac120.summary.tsv`
* Writes: `1_arg_dss_regresion.pdf`, `1_MazEF_arg3.pdf`, `1_pan-immune_arg2.pdf`

### `analyses/01_main_figures/Fig05_A-D.R`

**Mobile genetic elements and the defensome** — Prophage/plasmid counts, phylum-level sankey, plasmid:prophage ratio, DS density on chromosomes vs MGEs, and the observed/expected (O/E) plus odds-ratio statistics.

* Manuscript: Fig. 5a-d
* Working directory: `GTR_REV/01defensome/NC_revision/1additional_analyses/ICEs`, `groundwater/02去冗余信息`, `groundwater/05MGE`, `groundwater/02ANI99MAGs_defense`, `groundwater/03mag_info`
* Reads:
    * `virus_taxonomy_count.tsv`
    * `all_prophage_contig.tsv`
    * `plasmid_taxonomy_count.tsv`
    * `complete_integron_contig.id`
    * `all_ICE_contigid.txt`
    * `gtdbtk.bac120.summary.tsv`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `genome_stats.csv`
    * `ds_density_tax_r220.csv`
    * `all_integron_dss2.txt`
    * `all_ICE_dss2.txt`
    * `all_prophage.tsv`
    * `all_prophage_dss2.txt`
    * `all_plasmid_contig.tsv`
    * `all_virus_dss_counts_matrix.csv`
    * `all_plasmid_dss_counts_matrix.csv`
    * `all_prophage_dss_counts_matrix.csv`
    * `all_integron_dss_counts_matrix.csv`
    * `all_ICE_dss_counts_matrix.csv`
* Writes: `1_mge_phyla_sankey_new.pdf`, `1mge_plasmid_prophage_bili_new.pdf`, `1_mge_ds_density_new.pdf`, `1_mag_ds_oe.pdf`

### `analyses/01_main_figures/Fig05_E-I.R`

**Defense islands** — Defense-island gene number for CPR vs other phyla, by geological environment, the proportion of defense genes in DIs, and DS co-localization network plus Fisher tests.

* Manuscript: Fig. 5e-i
* Working directory: `groundwater/02ANI99MAGs_defense`, `groundwater/02去冗余信息`, `groundwater`, `groundwater/04mag_DI_net`
* Reads:
    * `defense_islands_summary4.csv`
    * `gtdbtk.bac120.summary.tsv`
    * `GW_DATA_lpw.xlsx`
    * `ds_counts_matrix.csv`
* Writes: `8CPR_DIs_genenumber.pdf`, `8DI-length-(genes).pdf`, `8dizhihuanjing_DIs_genenumber.pdf`, `8DGs-proportion-in-DIs.pdf`, `8_di_ds_oe.pdf`

### `analyses/01_main_figures/Fig06.R`

**Arms-race figure assembly** — Assembly fragment of Fig. 6: draws the prophage/anti-defense gene maps from the objects built by FigS17-18.R.

* Manuscript: Fig. 6
* Working directory: `groundwater/31gene_circle`, `groundwater/31gene_circle/acr_cas_two_contig`
* Reads: -
* Writes: `p_prophage_anti_colo4.pdf`, `p_acrcascontig.pdf (objects built in FigS17-18.R)`

### `analyses/02_supplementary_figures/FigS02_B.R`

**MAG abundance by niche breadth** — MAG abundance grouped by niche breadth and the stacked abundance composition per zone.

* Manuscript: Fig. S2b
* Working directory: `groundwater/s12mag_abundance`, `groundwater`
* Reads:
    * `bac.rpkm.tsv`
    * `ds_density_tax_r220.csv`
    * `GW_DATA_lpw.xlsx`
* Writes: `abs_stack.pdf`, `abs_stack_zone2.pdf`

### `analyses/02_supplementary_figures/FigS02_C.R`

**DS density across geological settings** — DS density by geological environment zone, well type, confined/unconfined aquifer and well depth; upset plot of defense families. Superseded by FigS08.R (which adds the source-data export).

* Manuscript: Fig. S2c
* Working directory: `groundwater/02去冗余信息`, `groundwater`, `groundwater/02ANI99MAGs_defense`, `groundwater/03mag_info`, `groundwater/01样品信息`
* Reads:
    * `env.xlsx`
    * `GW_DATA_lpw.xlsx`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `genome_stats.csv`
    * `GW_DATA.xlsx`
* Writes: `4_geographic_analyses.pdf`, `4_log10DSs_density_jing-type2.pdf`, `4_p_DSs_log10density_dizhihuanjing2.pdf`, `DSs_density_chengya2.pdf`

### `analyses/02_supplementary_figures/FigS02_D.R`

**Defense families across zones** — Upset plot of the defense families found in the seven aquifer zones (archaeal MAGs).

* Manuscript: Fig. S2d
* Working directory: `groundwater/02去冗余信息`, `groundwater`, `groundwater/archaea/02defense`
* Reads:
    * `env.xlsx`
    * `GW_DATA_lpw.xlsx`
* Writes: `4_upset_zone_ds_families.pdf`

### `analyses/02_supplementary_figures/FigS03.R`

**Defense genes vs genome size** — DG number vs MAG size, DG distribution per MAG, DG percentage and DG density.

* Manuscript: Fig. S3
* Working directory: `groundwater/s10mag_dg`, `groundwater/02去冗余信息`, `groundwater/03mag_info`, `groundwater/s1mag_dg`
* Reads:
    * `dg_counts_matrix.csv`
    * `Wdb_ANI99.csv`
    * `genome_stats.csv`
* Writes: `s1_DGnumber-density_magsize.pdf`, `s1_DGs_distribution.pdf`, `s1_DGs_percentage.pdf`, `s1_p3_DGs_density_magsize.pdf`, `s1_p4_DGs_number_magsize.pdf`

### `analyses/02_supplementary_figures/FigS04.R`

**Defensome of archaea** — DS distribution in archaea, defense-family composition and treemap, DPANN vs other phyla comparisons, DS prevalence and DG occurrence.

* Manuscript: Fig. S4
* Working directory: `groundwater/archaea/02defense`, `groundwater/archaea/01MAG_info`, `groundwater/01defensefinder_info`, `groundwater/archaea/01mag_info`
* Reads:
    * `archaea_ds_counts_matrix.csv`
    * `2269_metainfo.xlsx`
    * `zhong_70_10.xlsx`
    * `quality_report_checkm2.tsv`
    * `df_List_Systems_lpw.csv`
    * `archaea_dg_counts_matrix.csv`
    * `md_archaea.csv`
* Writes: `1_p1_DSs_distribution.pdf`, `1_p2_DSs_percentage.pdf`, `1_DPANN_phyla_all2.pdf`, `1_g1_superphylum_compare.pdf`, `1_DGs_distribution.pdf`, `1_DG_chuxian.pdf`

### `analyses/02_supplementary_figures/FigS06.R`

**DefenseFinder vs DefensePredictor** — Compares the defense genes called by DefenseFinder and DefensePredictor, including the genes identified by both tools.

* Manuscript: Fig. S6
* Working directory: `GTR_REV/01defensome/NC_revision/1additional_analyses/defensepredictor`, `groundwater/4VHR_DHR`
* Reads:
    * `defensepredictor_genecount.csv`
    * `1626dg_counts_matrix.csv`
    * `dgs_mergemags_with_coords.txt`
    * `1626mag_defensepredictor_processed_filtered.csv`
* Writes: `1_DP_DF_compare.pdf`

### `analyses/02_supplementary_figures/FigS07.R`

**VHR and genome size** — Prophage number per MAG vs MAG size, VHR vs genome size, and the features of 'supercarrier' MAGs.

* Manuscript: Fig. S7
* Working directory: `groundwater/4VHR_DHR`, `groundwater/03mag_info`, `groundwater/05MGE`, `GTR_REV/01defensome/NC_revision/1additional_analyses/genomesize_phage`, `groundwater/02去冗余信息`, `groundwater/02ANI99MAGs_defense`
* Reads:
    * `VHR.xlsx`
    * `ds_density_tax_r220.csv`
    * `virus_taxonomy_count.tsv`
    * `plasmid_taxonomy_count.tsv`
    * `VHR-DHR.csv`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
* Writes: `1_prophage_MAGsize.pdf`

### `analyses/02_supplementary_figures/FigS08.R`

**DS density across geological settings (with source data)** — Same analysis as FigS02_C.R plus the export of the source-data tables (figs8*.csv).

* Manuscript: Fig. S8
* Working directory: `groundwater/02去冗余信息`, `groundwater`, `groundwater/02ANI99MAGs_defense`, `groundwater/03mag_info`, `groundwater/01样品信息`
* Reads:
    * `env.xlsx`
    * `GW_DATA_lpw.xlsx`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `genome_stats.csv`
    * `GW_DATA.xlsx`
* Writes: `4_geographic_analyses.pdf`, `figs8ab.csv`, `figs8c.csv`, `figs8d.csv`, `figs8e.csv`

### `analyses/02_supplementary_figures/FigS09.R`

**Plaque assay experiment** — Quantification and plotting of the three plaque-assay replicates (wet-lab validation).

* Manuscript: Fig. S9
* Working directory: `GTR_REV/01defensome/NC_revision/2experiments/plaque_results_20250508`
* Reads:
    * `plaque_1.xlsx`
    * `plaque_2.xlsx`
    * `plaque_3.xlsx`
* Writes: `plaque_count_plot_2.pdf`

### `analyses/02_supplementary_figures/FigS10.R`

**DS distribution across phyla and genera (data-prep script)** — Phylum- and genus-level DS frequency heatmaps, DS richness of the major phyla, and the CPR clade analysis. Also builds `ds_density_tax` (reused by Fig02, Fig03, Fig05_A-D, FigS02_B and FigS07) and `md_p_val1` (reused by Fig02) - run this script first.

* Manuscript: Fig. S10
* Working directory: `groundwater/02ANI99MAGs_defense`, `groundwater/02去冗余信息`, `groundwater/01样品信息`, `groundwater/03mag_info`
* Reads:
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `bac.gtdb.r207.xlsx`
    * `gtdbtk.bac120.summary.tsv`
    * `genome_stats.csv`
    * `ds_density_tax_r220.csv`
    * `gtdbtk.bac120.summary_ncbi2.txt`
* Writes: `3phyla_dss.pdf`, `3phyla_dss_richness.pdf`, `3genera_dss.pdf`, `3genera_dss_richness.pdf`, `3CPR_DSs_number.pdf`, `3CPR_clade_dss.pdf`, `3_CPR_clade_DSrichness.pdf`

### `analyses/02_supplementary_figures/FigS12.R`

**DS distribution in selected phyla** — DS distribution for Chlamydiota, Patescibacteria, Campylobacterota, Nitrospirota, Myxococcota and other phyla of interest.

* Manuscript: Fig. S12
* Working directory: `groundwater/02去冗余信息`, `groundwater/02ANI99MAGs_defense`
* Reads:
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `gtdbtk.bac120.summary.tsv`
* Writes: `11_DSs_distribution_p__Campylobacterota.pdf`, `11_DSs_distribution_p__Myxococcota.pdf`, `11_DSs_distribution_p__Patescibacteria.pdf`

### `analyses/02_supplementary_figures/FigS13.R`

**Mobile resistance genes and the defensome** — Regression between MRG content and DS number, with the MazEF panel.

* Manuscript: Fig. S13
* Working directory: `groundwater/10mag_arg`, `groundwater/02ANI99MAGs_defense`, `groundwater/02去冗余信息`
* Reads:
    * `mrg_gene_out1_type.csv`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `dg_counts_matrix.csv`
    * `gtdbtk.bac120.summary.tsv`
* Writes: `1_arg_dss_regresion.pdf`, `1_MazEF_mrg3.pdf`, `2_MazEF_mrg2.pdf`

### `analyses/02_supplementary_figures/FigS14.R`

**Virulence factor genes and the defensome** — Regression between VFG content and DS number.

* Manuscript: Fig. S14
* Working directory: `groundwater/10mag_arg`, `groundwater/02ANI99MAGs_defense`, `groundwater/02去冗余信息`
* Reads:
    * `merged_VFG_output_file_num.csv`
    * `Wdb_ANI99.csv`
    * `ds_counts_matrix.csv`
    * `dg_counts_matrix.csv`
    * `gtdbtk.bac120.summary.tsv`
* Writes: `3_VFG_dss_regresion.pdf`, `03_ds_VFG2.pdf`

### `analyses/02_supplementary_figures/FigS15.R`

**Prophages and plasmids across phyla** — Counts of prophages and plasmids per phylum and the corresponding sankey diagram.

* Manuscript: Fig. S15
* Working directory: `GTR_REV/01defensome/NC_revision/1additional_analyses/ICEs`, `groundwater/02去冗余信息`, `GTR_REV/03eu_pre/11eukdetect`
* Reads:
    * `virus_taxonomy_count.tsv`
    * `all_prophage_contig.tsv`
    * `plasmid_taxonomy_count.tsv`
    * `complete_integron_contig.id`
    * `all_ICE_contigid.txt`
    * `gtdbtk.bac120.summary.tsv`
* Writes: `1_mge_phyla_sankey_new.pdf`

### `analyses/02_supplementary_figures/FigS16.R`

**Function of non-defense genes in defense islands** — COG composition of the eggNOG annotation of the non-defense genes found inside defense islands.

* Manuscript: Fig. S16
* Working directory: `groundwater/04mag_DI_gene`
* Reads:
    * `out.emapper.annotations.xlsx`
* Writes: `COG composition figure`

### `analyses/02_supplementary_figures/FigS17-18.R`

**Gene-level maps of representative elements (data-prep script)** — Gene-circle/gene-arrow maps of representative prophages, defense-antidefense co-localization, Acr-Cas contig pairs, and the T4SS/TxSS defense islands. Also builds the objects that Fig06.R draws - run it before Fig06.R.

* Manuscript: Fig. S17, Fig. S18
* Working directory: `groundwater/31gene_circle`, `groundwater/31gene_circle/eggnog_10prophage_anticolocation`, `groundwater/31gene_circle/virus_ds_anti_both`, `groundwater/31gene_circle/acr_cas_two_contig`, `groundwater/31gene_circle/DIs_with-T4SS`
* Reads:
    * `virus_3anti.gff`
    * `Ji15_bin.49_virus_anti.txt`
    * `Ji15_bin.49_defense.txt`
    * `coloca3_anti.gff`
    * `coloca3.txt`
    * `out.emapper.annotations`
    * `65virus.gff`
    * `65virus_anti_annotation.txt`
    * `65virus_DG_annotation.txt`
    * `two_special_contig.gff`
    * `twocontig_anti_annotation.txt`
    * `twocontig_DG_annotation.txt`
    * `di119.gff`
    * `DIs_nonDG_gene_CDD.txt`
    * `di119_DG_annotation.txt`
    * `allDIs_t.gff`
    * `DIs_nonDG_gene_CDD_TxSS.txt`
    * `allDIs_DG_annotation.txt`
* Writes: `p_virus_3contig.pdf`, `p_mag_virus.pdf`, `p_prophage_anti_colo2.pdf`, `p_prophage_anti_colo4.pdf`, `p_acrcascontig.pdf`, `p_crispr.pdf`, `p_DI-with-T4SS.pdf`, `p_allDI-with-TxSS.pdf`, `p_allDI-with-TxSS_2.pdf`

