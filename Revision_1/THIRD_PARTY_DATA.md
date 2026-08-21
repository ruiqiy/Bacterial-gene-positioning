# Obtaining third-party input data

The pipeline requires source files from published supplementary materials and public biological databases. These files are intentionally not redistributed in this repository. Download them from the official sources below, retain the required local filenames, and place them under `Revision_1/data/raw_data/`.

## Published proteomics workbooks

Place all six files in `data/raw_data/proteomics/`. The included `manifest.csv` records the exact worksheets and cell ranges read by `scripts/extract_raw_proteomics.R`.

| Required filename | Official source |
|------------------------------------|------------------------------------|
| `Schmidt_raw.xlsx` | [Schmidt et al. (2016), Supplementary Tables S1-S29](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.3418/MediaObjects/41587_2016_BFnbt3418_MOESM18_ESM.xlsx) |
| `Peebo_raw.xlsx` | [Peebo et al. (2015), supplementary files](https://doi.org/10.1039/C4MB00721B) |
| `Valgepea_raw.xlsx` | [Valgepea et al. (2013), supplementary files](https://doi.org/10.1039/C3MB70119K) |
| `Li_raw.xlsx` | [Li et al. (2014), supplementary table](https://www.cell.com/cms/10.1016/j.cell.2014.02.033/attachment/de4645b4-4bea-4197-abd2-d45249e51aaa/mmc1.xlsx) |
| `Goelzer_raw.xlsx` | [Goelzer et al. (2015), supplementary table](https://ars.els-cdn.com/content/image/1-s2.0-S1096717615001317-mmc3.xlsx) |
| `Zhu_raw.xlsx` | [Zhu et al. (2025), Dataset S01](https://www.pnas.org/doi/suppl/10.1073/pnas.2427091122/suppl_file/pnas.2427091122.sd01.xlsx) |

## Reference genome annotations

Download GFF3 annotations for the exact accessions below from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/), rename them as required, and place them in `data/raw_data/reference_gff/`. Do not substitute a different assembly version.

| Required filename | NCBI record | Analyzed sequence |
|------------------------|------------------------|------------------------|
| `GCF_000005845.2.gff` | [*E. coli* K-12 MG1655](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/) | `NC_000913.3` |
| `GCF_000009045.1.gff` | [*B. subtilis* 168](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000009045.1/) | `NC_000964.3` |
| `GCF_001456255.1.gff` | [*V. natriegens* ATCC 14048](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001456255.1/) | `NZ_CP009977.1` |
| `U00096.2.gff3` | [historical MG1655 sequence `U00096.2`](https://www.ncbi.nlm.nih.gov/nuccore/U00096.2) | `U00096.2` |

The historical `U00096.2.gff3` is specifically required for promoter-coordinate mapping; it is not interchangeable with `GCF_000005845.2.gff`.

## PromoterDB tables

Place these files in `data/raw_data/promoter/`:

| Required filename | Official source |
|------------------------------------|------------------------------------|
| `LB_promoter_regions.csv` | [PromoterDB LB table](https://ecolipromoterdb.com/data_tables/LB_promoter_regions.csv) |
| `M9_promoter_regions.csv` | [PromoterDB M9 table](https://ecolipromoterdb.com/data_tables/M9_promoter_regions.csv) |

## PanKB exports

The snapshots used here were retrieved on 2026-08-06. Place all six CSV files in `data/raw_data/PanKB/`. Use the species key and assembly in the query URLs, then rename each download as shown.

``` text
https://pankb.org/pangenome_analyses/gene_annotation/csv/?species=SPECIES_KEY
https://pankb.org/gene_function/genome_info/genes/csv/?species=SPECIES_KEY&genome_id=ASSEMBLY
```

| Species key         | Assembly          | Required filename                 |
|------------------------|------------------------|------------------------|
| `Escherichia_coli`  | `GCF_000005845.2` | `Ecoli_gene_annotations.csv`      |
| `Escherichia_coli`  | `GCF_000005845.2` | `Ecoli_focal_gene_info.csv`       |
| `Bacillus_subtilis` | `GCF_000009045.1` | `Bsub_Ecoli_gene_annotations.csv` |
| `Bacillus_subtilis` | `GCF_000009045.1` | `Bsub_focal_gene_info.csv`        |
| `Vibrio_natriegens` | `GCF_001456255.1` | `Vnat_Ecoli_gene_annotations.csv` |
| `Vibrio_natriegens` | `GCF_001456255.1` | `Vnat_focal_gene_info.csv`        |

PanKB requests citation of its current publication; see the [PanKB About page](https://pankb.org/about/) for citation and access information.

## UniProt molecular-weight snapshots

These inputs can be regenerated from `Zhu_raw.xlsx` after workbook extraction:

``` powershell
Rscript .\Revision_1\scripts\extract_raw_proteomics.R
Rscript .\Revision_1\scripts\download_molecular_weight_api.R
```

The downloader writes the following files to `data/raw_data/molecular_weight/`. The reported run used UniProt release `2026_02`; API results may change with later releases.

| Required filename                                   |
|-----------------------------------------------------|
| `zhu_E_coli_uniprot_molecular_weight_raw.tsv`       |
| `zhu_B_subtilis_uniprot_molecular_weight_raw.tsv`   |
| `zhu_V_natriegens_uniprot_molecular_weight_raw.tsv` |

UniProt applies the [CC BY 4.0 license](https://www.uniprot.org/help/license/) to copyrightable database content; retain attribution when sharing derived records.
