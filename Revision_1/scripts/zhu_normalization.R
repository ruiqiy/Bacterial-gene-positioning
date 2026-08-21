#!/usr/bin/env Rscript

# Normalize Zhu mass-fraction measurements after stable identifiers have been
# attached to the raw rows. Molecular weights come from the identifier-based
# tables under data/molecular_weight. No GFF or genomic position is used here.

project_root <- normalizePath(
  Sys.getenv(
    "GENE_POSITION_ROOT",
    unset = getwd()
  ),
  winslash = "/"
)
revision_root <- file.path(project_root, "Revision_1")
input_root <- file.path(revision_root, "data", "raw_mapped")
output_root <- file.path(revision_root, "data", "normalized_proteomics")
mw_file <- file.path(
  revision_root,
  "data",
  "molecular_weight",
  "zhu_molecular_weight_all.csv"
)
qc_root <- file.path(revision_root, "qc", "normalization")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(
  file.path(revision_root, "R_libs"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

metadata <- tribble(
  ~ProjectID,   ~Batch, ~SampleID, ~OriginalID, ~Species,         ~Condition_Raw,                    ~Temp_C, ~GrowthRate_hr,
  # Vibrio natriegens (Batch C - 37°C)
  "WD039LQ",    "C",    "C1",      "V1",        "V_natriegens",   "LB3 broth",                        37,      2.78,
  "WD039LQ",    "C",    "C2",      "V2",        "V_natriegens",   "0.4% glucose + 0.4% CAA",          37,      2.05,
  "WD039LQ",    "C",    "C3",      "V3",        "V_natriegens",   "0.4% glucose",                     37,      1.65,
  "WD039LQ",    "C",    "C4",      "V4",        "V_natriegens",   "0.4% glycerol",                    37,      0.86,
  "WD039LQ",    "C",    "C5",      "V5",        "V_natriegens",   "0.4% fructose",                    37,      1.29,
  "WD039LQ",    "C",    "C6",      "V6",        "V_natriegens",   "50 mM acetate",                    37,      0.26,
  
  # Vibrio natriegens (Batch L - 30°C + control)
  "XB02027LQ",  "L",    "L1",      "V1",        "V_natriegens",   "LB3 broth",                        30,      2.10,
  "XB02027LQ",  "L",    "L2",      "V2",        "V_natriegens",   "0.4% glucose + 0.4% CAA",          30,      1.63,
  "XB02027LQ",  "L",    "L3",      "V3",        "V_natriegens",   "0.4% glucose",                     30,      1.30,
  "XB02027LQ",  "L",    "L4",      "V4",        "V_natriegens",   "0.4% fructose",                    30,      1.11,
  "XB02027LQ",  "L",    "L5",      "V5",        "V_natriegens",   "0.4% glycerol",                    30,      0.70,
  "XB02027LQ",  "L",    "L6",      "V6",        "V_natriegens",   "0.4% glucose",                     37,      1.64, # Replicate of C3
  
  # Bacillus subtilis (Batch F - 37°C)
  "AD066LQ",    "F",    "F6",      "B6",        "B_subtilis",     "LB broth",                         37,      1.84,
  "AD066LQ",    "F",    "F7",      "B7",        "B_subtilis",     "0.4% glucose + 0.4% CAA",          37,      1.38,
  "AD066LQ",    "F",    "F8",      "B8",        "B_subtilis",     "0.4% glucose + 0.4% CAA",          37,      1.38, # Replicate of F7
  "AD066LQ",    "F",    "F9",      "B9",        "B_subtilis",     "0.4% glucose + 1 mM glutamate",    37,      0.80,
  "AD066LQ",    "F",    "F10",     "B10",       "B_subtilis",     "0.4% mannose + 1 mM glutamate",    37,      0.64,
  "AD066LQ",    "F",    "F11",     "B11",       "B_subtilis",     "0.4% ribose + 1 mM glutamate",     37,      0.44,
  
  # Escherichia coli (Batch F - 37°C)
  "AD067LQ",    "F",    "F12",     "E12",       "E_coli",         "0.2% glucose",                     37,      0.92,
  "AD067LQ",    "F",    "F13",     "E13",       "E_coli",         "0.2% glucose + 1 mM glutamate",    37,      1.00,
  "AD067LQ",    "F",    "F13-1",   "E13-1",     "E_coli",         "LB broth",                         37,      1.90,
  "AD067LQ",    "F",    "F13-2",   "E13-2",     "E_coli",         "0.2% glucose + 0.2% CAA",          37,      1.20,
  "AD067LQ",    "F",    "F13-3",   "E13-3",     "E_coli",         "0.2% glycerol",                    37,      0.67,
  "AD067LQ",    "F",    "F13-4",   "E13-4",     "E_coli",         "0.2% mannose",                     37,      0.42,
  
  # Escherichia coli (Batch K - 30°C + control)
  "XB02028LQ",  "K",    "K1",      "E1",        "E_coli",         "LB broth",                         30,      1.205,
  "XB02028LQ",  "K",    "K2",      "E2",        "E_coli",         "0.2% glucose + 0.2% CAA",          30,      0.915,
  "XB02028LQ",  "K",    "K3",      "E3",        "E_coli",         "0.2% glucose",                     30,      0.644,
  "XB02028LQ",  "K",    "K4",      "E4",        "E_coli",         "0.2% glycerol",                    30,      0.418,
  "XB02028LQ",  "K",    "K5",      "E5",        "E_coli",         "0.2% mannose",                     30,      0.395,
  "XB02028LQ",  "K",    "K6",      "E6",        "E_coli",         "0.2% glucose",                     37,      0.949, # Replicate of F12
  
  # Bacillus subtilis (Batch G - 43°C + control)
  "XA01797LQ",  "G",    "G1",      "B1",        "B_subtilis",     "LB broth",                         43,      2.57,
  "XA01797LQ",  "G",    "G2",      "B2",        "B_subtilis",     "0.4% glucose + 0.4% CAA",          43,      1.91,
  "XA01797LQ",  "G",    "G3",      "B3",        "B_subtilis",     "0.4% glucose + 1 mM glutamate",    43,      1.17,
  "XA01797LQ",  "G",    "G4",      "B4",        "B_subtilis",     "0.4% mannose + 1 mM glutamate",    43,      1.01,
  "XA01797LQ",  "G",    "G5",      "B5",        "B_subtilis",     "0.4% ribose + 1 mM glutamate",     43,      0.79,
  "XA01797LQ",  "G",    "G6",      "B6",        "B_subtilis",     "0.4% glucose + 1 mM glutamate",    37,      0.85 # Replicate of F9
)

species_config <- tribble(
  ~species_code, ~dataset, ~files,
  "E_coli", "zhu_E_coli", list(c("Ecoli_batch_F.csv", "Ecoli_batch_K.csv")),
  "B_subtilis", "zhu_B_subtilis", list(c("Bsub_batch_F.csv", "Bsub_batch_G.csv")),
  "V_natriegens", "zhu_V_natriegens", list(c("Vnat_batch_C.csv", "Vnat_batch_L.csv"))
)

molecular_weight <- read_csv(mw_file, show_col_types = FALSE) |>
  filter(mw_status == "mapped", is.finite(molecular_weight_Da)) |>
  select(dataset, source_protein_id, molecular_weight_Da, mw_source)
if (anyDuplicated(molecular_weight[c("dataset", "source_protein_id")])) {
  stop("Updated Zhu molecular-weight table has duplicated dataset/protein keys")
}

mapping_fields <- c(
  "source_identifier_raw",
  "source_identifier_normalized",
  "source_gene_name",
  "reference_identifier",
  "candidate_reference_identifiers",
  "mapping_method",
  "mapping_confidence",
  "mapping_status",
  "mapping_cardinality",
  "dataset",
  "species"
)

normalization_qc <- list()

for (config_index in seq_len(nrow(species_config))) {
  config <- species_config[config_index, ]
  species_metadata <- filter(metadata, Species == config$species_code)

  batches <- map(
    unlist(config$files[[1L]], use.names = FALSE),
    function(filename) {
      data <- read_csv(
        file.path(input_root, filename),
        show_col_types = FALSE,
        name_repair = "minimal"
      ) |>
        rename(gene = `gene name`, locus_tag = `gene locus`, source_protein_id = `protein ID`) |>
        mutate(
          gene = as.character(gene),
          locus_tag = as.character(locus_tag),
          source_protein_id = as.character(source_protein_id)
        )

      sample_columns <- intersect(names(data), metadata$SampleID)
      group_columns <- c(
        "gene",
        "locus_tag",
        "source_protein_id",
        mapping_fields
      )
      data |>
        group_by(across(all_of(group_columns))) |>
        summarise(
          across(all_of(sample_columns), ~ mean(.x, na.rm = TRUE)),
          .groups = "drop"
        )
    }
  )

  join_columns <- c("gene", "locus_tag", "source_protein_id", mapping_fields)
  abundance <- reduce(
    batches,
    full_join,
    by = join_columns,
    relationship = "one-to-one"
  )
  rows_before_mw <- nrow(abundance)

  abundance <- abundance |>
    left_join(
      filter(molecular_weight, dataset == config$dataset) |>
        select(-dataset),
      by = "source_protein_id",
      relationship = "many-to-one"
    )
  missing_mw <- sum(is.na(abundance$molecular_weight_Da))
  abundance <- filter(abundance, !is.na(molecular_weight_Da))

  sample_columns <- intersect(names(abundance), species_metadata$SampleID)
  for (sample_column in sample_columns) {
    relative_copies <- abundance[[sample_column]] /
      abundance$molecular_weight_Da
    denominator <- sum(relative_copies, na.rm = TRUE)
    abundance[[sample_column]] <- relative_copies / denominator
  }

  nonmissing_expression <- as.matrix(abundance[sample_columns])
  storage.mode(nonmissing_expression) <- "numeric"
  keep_expressed <- rowSums(
    apply(nonmissing_expression, 2, function(x) replace(x, is.na(x), 0))
  ) > 0
  abundance <- abundance[keep_expressed, , drop = FALSE]
  abundance <- abundance |>
    mutate(across(all_of(sample_columns), ~ replace_na(.x, 0) + 1e-8))

  long <- abundance |>
    select(
      gene,
      locus_tag,
      source_protein_id,
      all_of(mapping_fields),
      all_of(sample_columns)
    ) |>
    pivot_longer(
      cols = all_of(sample_columns),
      names_to = "SampleID",
      values_to = "Expression"
    ) |>
    left_join(species_metadata, by = "SampleID", relationship = "many-to-one")

  output <- long |>
    group_by(
      gene,
      locus_tag,
      source_protein_id,
      across(all_of(mapping_fields)),
      Condition_Raw,
      Temp_C
    ) |>
    summarise(
      Mean_Expr = mean(Expression, na.rm = TRUE),
      Mean_GR = mean(GrowthRate_hr, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      abundance_column = paste0(
        config$dataset,
        "__mu_",
        gsub("\\.", "_", sprintf("%.6f", Mean_GR))
      )
    ) |>
    select(
      gene,
      locus_tag,
      source_protein_id,
      all_of(mapping_fields),
      abundance_column,
      Mean_Expr
    ) |>
    pivot_wider(names_from = abundance_column, values_from = Mean_Expr)

  abundance_columns <- names(output)[str_detect(names(output), fixed("__mu_"))]
  output$n_cond <- rowSums(!is.na(output[abundance_columns]))
  saveRDS(output, file.path(output_root, paste0(config$dataset, ".rds")))

  normalization_qc[[config$dataset]] <- tibble(
    dataset = config$dataset,
    source_rows_after_batch_merge = rows_before_mw,
    molecular_weight_mapped = rows_before_mw - missing_mw,
    molecular_weight_missing = missing_mw,
    molecular_weight_coverage = (rows_before_mw - missing_mw) / rows_before_mw,
    expressed_rows_retained = nrow(output),
    mapped_one_to_one_retained = sum(output$mapping_status == "mapped_one_to_one"),
    unresolved_rows_retained = sum(output$mapping_status != "mapped_one_to_one"),
    abundance_columns = length(abundance_columns)
  )
}

write_csv(
  bind_rows(normalization_qc),
  file.path(qc_root, "zhu_molecular_weight_and_normalization_qc.csv"),
  na = ""
)
message("Wrote three position-independent Zhu datasets using updated molecular weights")
