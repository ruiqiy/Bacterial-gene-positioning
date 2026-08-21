# Gene Position

This repository contains the R scripts and intermediate data used for the manuscript "Bacterial chromosomal gene 
positioning is shaped by selection on both mean and growth-dependent expression. The analysis focuses on converting 
raw proteomics data to protein copy fraction data, calculating promoter strengths and coordinates, and performing 
comprehensive statistical analyses. The scripts were written with assistance of generative artificial intelligence.

## Workflow overview

The data processing and analysis pipeline consists of four main scripts:

1. **`proteomics_normalization.R`**
   - Converts raw data to protein copy fraction data for 5 of the datasets.
   - Outputs intermediate `.rds` files.

2. **`zhu_normalization.R`**
   - Converts raw data to protein copy fraction data for the 3 datasets published by Zhu et al.
   - Outputs intermediate `.rds` files.

3. **`promoter.R`**
   - Calculates promoter strengths and genomic coordinates.
   - Outputs intermediate `.rds` files (`LB_promoter.rds`, `M9_promoter.rds`, etc.).

4. **`final_gene_position.R`**
   - Takes the `.rds` files generated in the previous steps.
   - Performs the main statistical analyses and generates the tables and figures used in the manuscript.

## Data Availability

The raw data used in this study was sourced from the following publications/repositories. GFF files were obtained 
from NCBI. Molecular weight data files (e.g., `All-polypeptides-of-E.-coli-K-12-substr.-MG1655.txt`, `Ecoli_MW.txt`) included
columns "Proteins", "Genes", and "Molecular-Weight-From-Sequence", and were obtained from Biocyc.

*   **Dataset 1 (Schmidt et al.)**
    *   **Source/Link:** https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.3418/MediaObjects/41587_2016_BFnbt3418_MOESM18_ESM.xlsx
    Table S23 Growth condition and Growth rate column for `schmidt_meta.csv`. Table S5 Protein copies/cell part for `schmidt_data_copy.csv`.
    *   **Files used:** `schmidt_meta.csv`, `schmidt_data_copy.csv`, `All-polypeptides-of-E.-coli-K-12-substr.-MG1655.txt`

*   **Dataset 2 (Peebo et al.)**
    *   **Source/Link:** https://doi.org/10.1039/C4MB00721B Supplementary table
    *   **Files used:** `peebo_copy.csv`

*   **Dataset 3 (Valgepea et al.)**
    *   **Source/Link:** https://doi.org/10.1039/C3MB70119K Table S2 molecules/cell data.
    *   **Files used:** `Valgepea_copy.csv`

*   **Dataset 4 (Li et al.)**
    *   **Source/Link:** https://www.cell.com/cms/10.1016/j.cell.2014.02.033/attachment/de4645b4-4bea-4197-abd2-d45249e51aaa/mmc1.xlsx
    *   **Files used:** `Li_copy.csv`

*   **Dataset 5 (Goelzer et al.)**
    *   **Source/Link:** https://ars.els-cdn.com/content/image/1-s2.0-S1096717615001317-mmc3.xlsx Sheets PYR, S, TS, CH, CHG.
    *   **Files used:** `Goelzer_data_copy.xlsx`, `All-polypeptides-of-B.-subtilis-subtilis-168.txt`

*   **Datasets 6-8 (Zhu et al.)**
    *   **Source/Link:** https://www.pnas.org/doi/suppl/10.1073/pnas.2427091122/suppl_file/pnas.2427091122.sapp.pdf for metadata.
    https://www.pnas.org/doi/suppl/10.1073/pnas.2427091122/suppl_file/pnas.2427091122.sd01.xlsx for `*Ecoli*.csv`, `*Vnat*.csv`, `*Bsub*.csv`
    *   **Files used:** `*Ecoli*.csv`, `*Vnat*.csv`, `*Bsub*.csv`, `Ecoli_MW.txt`, `Vnat_MW.txt`, `Bsub_MW.txt` and `GCF_001043215.1.gff`, 
    `GCF_001456255.1.gff`, `GCF_000009045.1.gff`

*   **Promoter & Genomic Reference Data**
    *   **Source/Link:** https://ecolipromoterdb.com/data_tables/LB_promoter_regions.csv for `LB_promoter_regions.csv`
    https://ecolipromoterdb.com/data_tables/M9_promoter_regions.csv for `M9_promoter_regions.csv`
    *   **Files used:** `U00096.2.gff3`, `LB_promoter_regions.csv`, `M9_promoter_regions.csv`
