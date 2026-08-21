# Bacterial gene-position analysis pipeline

This directory contains the R pipeline used to convert published bacterial proteomics workbooks into identifier-mapped, chromosome-position-aware data and to test how gene position relates to mean expression, growth-dependent expression, gene-set membership, and promoter strength.

The pipeline analyzes eight datasets from *Escherichia coli*, *Bacillus subtilis*, and *Vibrio natriegens*. It runs five gene universes: all eligible genes, PanKB core genes, PanKB accessory genes (including rare genes), transcription/translation genes, and all other genes.

## Repository contents and data availability

This public package contains the analysis scripts, raw-workbook extraction manifest, compact reference tables, QC reports, and analysis results. It intentionally does **not** redistribute third-party supplementary workbooks, database exports, reference GFF files, the local R package library, regenerable intermediate datasets, or the separate manuscript-figure collection.

See [`THIRD_PARTY_DATA.md`](THIRD_PARTY_DATA.md) for the required filenames, official source links, and retrieval notes for the source files used in this analysis. After placing those inputs under `data/raw_data/`, the pipeline recreates the omitted intermediate directories without changing the source files.

## Quick start

1.  Follow [`THIRD_PARTY_DATA.md`](THIRD_PARTY_DATA.md) to obtain the required inputs, keep the specified filenames, and place them under `Revision_1/data/raw_data/` using the structure below.
2.  Install the required R packages. Optionally, place a project-local library at `Revision_1/R_libs/`; that directory is not included in this submission.
3.  From the directory containing `Revision_1/`, run:

``` powershell
Rscript .\Revision_1\scripts\run_pipeline.R
```

The driver creates intermediate data under `Revision_1/data/`, QC tables under `Revision_1/qc/`, and results under:

``` text
Revision_1/gene_set_analysis/
|-- all/
|-- core/
|-- accessory/
|-- transcription_translation/
`-- other/
```

Each result directory contains `figs/` and `tables/`. `run_pipeline.R` locates the project from its own path; `GENE_POSITION_ROOT` may be set explicitly as an override. The optional UniProt download step is not called by the driver, so either generate the three TSVs first or supply previously downloaded copies as described in `THIRD_PARTY_DATA.md`.

## Required input layout

Third-party files shown below are required at runtime but are not stored in this public package. Generated directories such as `raw_extracted/`, `raw_mapped/`, `molecular_weight/`, `normalized_proteomics/`, `position_mapped_proteomics/`, and `analysis_ready_proteomics/` may be absent before the first run. The included `qc/` and `gene_set_analysis/` directories are the retained results from the reported run and may be regenerated.

``` text
<project-root>/
`-- Revision_1/
    |-- README.md
    |-- scripts/
    |-- R_libs/                         # optional local R library
    `-- data/
        `-- raw_data/
            |-- proteomics/
            |   |-- manifest.csv
            |   |-- Schmidt_raw.xlsx
            |   |-- Peebo_raw.xlsx
            |   |-- Valgepea_raw.xlsx
            |   |-- Li_raw.xlsx
            |   |-- Goelzer_raw.xlsx
            |   `-- Zhu_raw.xlsx
            |-- reference_gff/
            |   |-- GCF_000005845.2.gff
            |   |-- GCF_000009045.1.gff
            |   |-- GCF_001456255.1.gff
            |   `-- U00096.2.gff3
            |-- promoter/
            |   |-- LB_promoter_regions.csv
            |   `-- M9_promoter_regions.csv
            |-- PanKB/
            |   |-- Ecoli_focal_gene_info.csv
            |   |-- Ecoli_gene_annotations.csv
            |   |-- Bsub_focal_gene_info.csv
            |   |-- Bsub_Ecoli_gene_annotations.csv
            |   |-- Vnat_focal_gene_info.csv
            |   `-- Vnat_Ecoli_gene_annotations.csv
            `-- molecular_weight/
                |-- zhu_E_coli_uniprot_molecular_weight_raw.tsv
                |-- zhu_B_subtilis_uniprot_molecular_weight_raw.tsv
                `-- zhu_V_natriegens_uniprot_molecular_weight_raw.tsv
```

## Obtaining the inputs

### Proteomics workbooks

Download each publication's supplementary workbook, rename it as shown, and place it in `data/raw_data/proteomics/`. The supplied `manifest.csv` specifies the exact sheets and cell ranges used; retain it unchanged unless the source workbooks change.

| File | Source | Data used |
|----|----|----|
| `Schmidt_raw.xlsx` | [Download Schmidt et al. (2016) Supplementary Tables S1–S29](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.3418/MediaObjects/41587_2016_BFnbt3418_MOESM18_ESM.xlsx) | Table S5 abundance and Table S23 condition metadata |
| `Peebo_raw.xlsx` | [Peebo et al. (2015)](https://doi.org/10.1039/C4MB00721B) Supplementary table | Protein concentrations and growth rates |
| `Valgepea_raw.xlsx` | [Valgepea et al. (2013)](https://doi.org/10.1039/C3MB70119K) Supplementary table | Table S2 protein abundance; RNA columns are provenance only |
| `Li_raw.xlsx` | [Download Li et al. (2014) Supplementary table](https://www.cell.com/cms/10.1016/j.cell.2014.02.033/attachment/de4645b4-4bea-4197-abd2-d45249e51aaa/mmc1.xlsx) | Gene-level protein synthesis |
| `Goelzer_raw.xlsx` | [Download Goelzer et al. (2015) Supplementary table](https://ars.els-cdn.com/content/image/1-s2.0-S1096717615001317-mmc3.xlsx) | PYR, S, TS, CH, and CHG sheets |
| `Zhu_raw.xlsx` | [Download Zhu et al. (2025) Supplementary table](https://www.pnas.org/doi/suppl/10.1073/pnas.2427091122/suppl_file/pnas.2427091122.sd01.xlsx) | Dataset S01 for all three species |

### Reference GFF files

Download the GFF3 annotations for the exact RefSeq assemblies from NCBI and place them in `data/raw_data/reference_gff/`:

| File | Reference and analyzed sequence |
|----|----|
| `GCF_000005845.2.gff` | *E. coli* K-12 MG1655; `NC_000913.3` |
| `GCF_000009045.1.gff` | *B. subtilis* 168; `NC_000964.3` |
| `GCF_001456255.1.gff` | *V. natriegens* ATCC 14048 chromosome 1; `NZ_CP009977.1` |

The promoter mapper separately requires historical `U00096.2.gff3`, because the promoter tables use that MG1655 coordinate system. Do not substitute the newer `GCF_000005845.2.gff`.

### Promoter tables

Download the [LB promoter regions](https://ecolipromoterdb.com/data_tables/LB_promoter_regions.csv) and [M9 promoter regions](https://ecolipromoterdb.com/data_tables/M9_promoter_regions.csv) from E. coli Promoter Database "PromoterDB" and place them in `data/raw_data/promoter/`.

### PanKB tables

The driver does not download PanKB data. For each species, obtain the complete gene-annotation CSV and the focal-genome gene CSV from these endpoints:

``` text
https://pankb.org/pangenome_analyses/gene_annotation/csv/?species=SPECIES_KEY
https://pankb.org/gene_function/genome_info/genes/csv/?species=SPECIES_KEY&genome_id=ASSEMBLY
```

| Species key | Assembly | Focal file | Annotation file |
|----|----|----|----|
| `Escherichia_coli` | `GCF_000005845.2` | `Ecoli_focal_gene_info.csv` | `Ecoli_gene_annotations.csv` |
| `Bacillus_subtilis` | `GCF_000009045.1` | `Bsub_focal_gene_info.csv` | `Bsub_Ecoli_gene_annotations.csv` |
| `Vibrio_natriegens` | `GCF_001456255.1` | `Vnat_focal_gene_info.csv` | `Vnat_Ecoli_gene_annotations.csv` |

The annotation stage joins focal genes to PanKB families, maps `original_locus_tag` to the reference identifier, and combines PanKB `accessory` and `rare` as the binary accessory class. PanKB data retrieval date: 2026-08-06.

### Zhu molecular-weight TSVs

The three raw UniProt TSVs may be supplied as a fixed snapshot or regenerated after workbook extraction:

``` powershell
Rscript .\Revision_1\scripts\extract_raw_proteomics.R
Rscript .\Revision_1\scripts\download_molecular_weight_api.R
```

The downloader uses `curl.exe`, UniProt ID mapping, and UniParc recovery and writes to `data/raw_data/molecular_weight/`. The preparation stage writes derived molecular-weight CSVs to `data/molecular_weight/`.

## Software

The pipeline was checked with R 4.4.2 on Windows and uses:

``` text
betareg, car, dplyr, ggbreak, ggh4x, ggplot2, gridExtra, jsonlite,
lmtest, patchwork, performance, purrr, readr, readxl, scales,
stringr, tibble, tidyr, viridis
```

Scripts prepend `Revision_1/R_libs/` to `.libPaths()`. The project-local directory is optional and is not included here; install missing packages in that directory or in a normal R library. The optional downloader also needs Internet access and `curl.exe`; PDF output uses Cairo where requested.

## Workflow

``` mermaid
flowchart TD
    RAW["raw_data/proteomics<br/>six workbooks + manifest"] --> EX["1. Extract declared sheets/ranges"]
    EX --> RE["data/raw_extracted"]
    GFF["raw_data/reference_gff<br/>three RefSeq GFFs"] --> MAP["2. Map stable reference identifiers"]
    RE --> MAP
    MAP --> RM["data/raw_mapped + mapping QC"]
    RE --> DL["Optional UniProt download"]
    DL --> MWR["raw_data/molecular_weight"]
    MWR --> MW["3. Prepare molecular-weight tables"]
    RM --> N["4. Normalize eight datasets"]
    MW --> N
    N --> POS["5. Attach chromosome positions + eligibility"]
    GFF --> POS
    GFF --> TT["6a. Build transcription/translation set"]
    PK["raw_data/PanKB<br/>six CSVs"] --> JOIN["6b. Add gene-set annotations"]
    TT --> JOIN
    POS --> JOIN
    JOIN --> AR["data/analysis_ready_proteomics"]
    PROM["raw promoters + U00096.2 GFF3"] --> PM["7. Map promoters to b-numbers"]
    AR --> ANALYZE["8. Analyze five gene universes"]
    PM --> ANALYZE
    ANALYZE --> OUT["gene_set_analysis/{all, core, accessory,<br/>transcription_translation, other}"]
```

## Scripts

| Script | Role |
|----|----|
| `run_pipeline.R` | Windows entry point; runs all active stages and all five gene universes. |
| `extract_raw_proteomics.R` | Extracts only manifest-declared workbook ranges and writes extraction QC. |
| `mapping_helpers.R` | Shared GFF parsing, identifier normalization, alias, and position helpers; sourced by other scripts. |
| `build_cross_strain_mapping.R` | Maps raw rows to stable identifiers, applies only unique common-name fallbacks, and writes mapping QC. |
| `download_molecular_weight_api.R` | Optional UniProt/UniParc downloader; not called by the driver. |
| `prepare_molecular_weight_tables.R` | Builds one UniProt molecular-weight record per Zhu source protein. |
| `proteomics_normalization.R` | Normalizes Schmidt, Peebo, Valgepea, Li, and Goelzer abundance without using position. |
| `zhu_normalization.R` | Converts Zhu mass fractions to copy fractions using matched molecular weights. |
| `attach_reference_positions.R` | Attaches reference positions and assigns analysis eligibility. |
| `build_transcription_translation_gene_set.R` | Builds RNA-polymerase/ribosomal-protein membership from GFF CDS products. |
| `build_raw_rebuild_joined_proteomics.R` | Adds PanKB and transcription/translation classes to eligible rows. |
| `promoter.R` | Maps promoter peaks to intragenic or nearest downstream genes within 500 bp. |
| `final_gene_position_gene_sets.R` | Runs the main statistics and figures once per selected gene universe. |
| `run_gene_set_analyses.R` | During the all-gene run, creates supplementary, transcription/translation (TT) attenuation, and interaction outputs. |

## Analysis rules

- Primary identifiers are UniProt accessions or supplied reference locus tags. A normalized common name is used only when all matching aliases identify one reference locus. Ambiguous and unmatched records remain unresolved.
- Li compound records containing `+` are excluded before normalization because their abundance cannot be assigned to one gene.
- Normalization divides Schmidt, Peebo, Valgepea, and retained Li values by the condition total, then adds `1e-8`. Goelzer normalizes technical replicates, averages technical and biological replicates, replaces missing conditions with zero, then adds `1e-8`. Zhu divides mass fraction by matched molecular mass, renormalizes to copy fraction, averages biological replicates, then adds `1e-8`.
- A row is analysis-eligible only if its identifier mapping is one-to-one, the reference identifier is unique within that dataset, and its reference position is finite. Position is joined only by species and stable identifier.
- Normalized position is midpoint distance from oriC, scaled from zero at oriC to one at the terminus.
- Consensus analyses generally retain a gene observed in any source dataset. Only consensus beta regression requires at least three E. coli source datasets or both B. subtilis datasets. Optional downsampling uses that same beta-regression universe.
- The two retained `supp_*` PDFs each have five pages. Panel display limits use the 1st-99th response percentiles, but the OLS line, equation, zero-slope P-value, and `n` use all finite observations. Annotations are bottom-left; point styling matches the promoter scatter (`alpha = 0.05`, `stroke = 0`).
- Core/accessory forest plots show the `norm_pos`-by-accessory interaction, meaning accessory positional slope minus core positional slope. All-gene and TT-excluded models use circles distinguished by color.
- TT attenuation uses 10,000 equal-size random removals. It reports `delta_TT = beta_all - beta_nonTT`, its null-standardized Z score, and an empirical P-value. Separate outputs summarize individual datasets and consensus datasets; the consensus table also includes the E. coli promoter-strength comparison.
- Condition downsampling is disabled by default (`ENABLE_DOWNSAMPLING <- FALSE`). Set it to `TRUE` in `final_gene_position_gene_sets.R` to run 100 replicates per condition count. It runs only for the `all` gene universe and uses the same minimum dataset-presence rules as the consensus beta regressions.

## Outputs and QC

Each gene-set directory contains the main correlation, position-region, beta-regression, consensus, promoter, and stability outputs. The cross-gene-set supplementary module runs only during the `all` analysis and writes these additional files under `gene_set_analysis/all/`:

``` text
figs/
|-- supp_lnE_by_position_all_gene_sets.pdf       # five pages, one per gene universe
|-- supp_k_by_position_all_gene_sets.pdf         # five pages, one per gene universe
|-- tt_attenuation_individual.{pdf,png}
|-- tt_attenuation_consensus.{pdf,png}
|-- fig_interaction_coefficients_forest_individual.{pdf,png}
`-- fig_interaction_coefficients_forest_consensus.{pdf,png}

tables/
|-- gene_set_dataset_counts.csv
|-- tt_attenuation_individual.csv
|-- tt_attenuation_consensus.csv
|-- interaction_models_group_slopes_summary.csv
`-- interaction_models_interaction_term_summary.csv
```

When downsampling is enabled, the `all/` directory also receives `figs/downsampling_effect_size_evolution.pdf`, `tables/downsampling_simulation_raw.csv`, and `tables/downsampling_simulation_summary.csv`. The driver creates or overwrites expected outputs but does not clean old results, so remove obsolete output files before a clean reproducibility run.
