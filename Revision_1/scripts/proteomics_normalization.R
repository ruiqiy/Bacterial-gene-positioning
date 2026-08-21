#!/usr/bin/env Rscript

# Normalize the five non-Zhu protein-abundance datasets. Stable reference
# identifiers have already been attached to every source row by
# build_cross_strain_mapping.R. This stage deliberately does not read a GFF,
# calculate a genomic position, or discard a row because its mapping failed.

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
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

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

mapping_columns <- c(
  "source_file",
  "source_row_number",
  "source_row_id",
  "source_gene_name",
  "source_identifier_raw",
  "source_identifier_normalized",
  "source_gene_name_normalized",
  "reference_identifier",
  "candidate_reference_identifiers",
  "ambiguity_n",
  "had_ambiguous_evidence",
  "mapping_method",
  "mapping_confidence",
  "mapping_status",
  "mapping_cardinality",
  "identifier_well_formed",
  "dataset",
  "species"
)

format_mu <- function(dataset, growth_rate, suffix = NULL) {
  label <- gsub("\\.", "_", sprintf("%.6f", round(growth_rate, 6)))
  if (!is.null(suffix)) label <- paste0(label, suffix)
  paste0(dataset, "__mu_", label)
}

copy_fraction <- function(values, replace_nonfinite = FALSE, pseudo = 1e-8) {
  matrix <- as.matrix(values)
  storage.mode(matrix) <- "numeric"
  if (replace_nonfinite) matrix[!is.finite(matrix)] <- 0
  totals <- colSums(matrix, na.rm = TRUE)
  totals[totals == 0] <- NA_real_
  sweep(matrix, 2, totals, "/") + pseudo
}

# Schmidt 2016: protein copies per cell -> copy fraction.
schmidt <- read_csv(
  file.path(input_root, "schmidt_data_copy.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
)
schmidt_meta_raw <- read_csv(
  file.path(revision_root, "data", "raw_extracted", "schmidt_meta.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
)
schmidt_meta <- schmidt_meta_raw |>
  transmute(
    condition_raw = as.character(.data[[names(schmidt_meta_raw)[1L]]]),
    growth_rate = as.numeric(.data[[names(schmidt_meta_raw)[2L]]]),
    condition_key = str_to_lower(str_squish(condition_raw))
  )

schmidt_conditions <- setdiff(
  names(schmidt),
  c("Uniprot Accession", "Gene", mapping_columns)
)
schmidt_condition_info <- tibble(
  condition = schmidt_conditions,
  condition_key = str_to_lower(str_squish(schmidt_conditions))
) |>
  left_join(schmidt_meta, by = "condition_key")
if (anyNA(schmidt_condition_info$growth_rate)) {
  stop("Schmidt abundance columns did not all match extracted condition metadata")
}

schmidt_fraction <- copy_fraction(schmidt[schmidt_conditions])
colnames(schmidt_fraction) <- format_mu(
  "schmidt_2016",
  schmidt_condition_info$growth_rate
)
schmidt_out <- schmidt |>
  transmute(
    gene = as.character(Gene),
    uniprot_accession = as.character(`Uniprot Accession`),
    across(any_of(mapping_columns))
  ) |>
  mutate(n_cond = rowSums(!is.na(as.matrix(schmidt[schmidt_conditions])))) |>
  bind_cols(as.data.frame(schmidt_fraction, check.names = FALSE))
saveRDS(
  schmidt_out,
  file.path(output_root, "schmidt_2016_pseudo.rds")
)

# Peebo 2015: protein concentration -> copy fraction. Repeated growth rates
# remain separate conditions with deterministic replicate suffixes.
peebo <- read_csv(
  file.path(input_root, "peebo_copy.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
)
peebo_conditions <- names(peebo)[str_detect(
  names(peebo),
  "^[0-9]+(?:\\.[0-9]+)?(?:\\.\\.\\.[0-9]+)?$"
)]
peebo_condition_info <- tibble(
  condition = peebo_conditions,
  growth_rate = as.numeric(sub("\\.\\.\\..*$", "", peebo_conditions))
) |>
  mutate(mu_label = sprintf("%.6f", round(growth_rate, 6))) |>
  group_by(mu_label) |>
  mutate(
    replicate_index = row_number(),
    replicate_count = n(),
    suffix = if_else(
      replicate_count > 1L,
      paste0("_rep", replicate_index),
      ""
    )
  ) |>
  ungroup()
peebo_fraction <- copy_fraction(
  peebo[peebo_conditions],
  replace_nonfinite = TRUE
)
colnames(peebo_fraction) <- format_mu(
  "peebo_2015",
  peebo_condition_info$growth_rate,
  peebo_condition_info$suffix
)
peebo_out <- peebo |>
  transmute(
    gene = as.character(Gene),
    locus_tag = as.character(BID),
    across(any_of(mapping_columns))
  ) |>
  mutate(n_cond = rowSums(!is.na(peebo_fraction))) |>
  bind_cols(as.data.frame(peebo_fraction, check.names = FALSE))
saveRDS(peebo_out, file.path(output_root, "peebo_2015_pseudo.rds"))

# Valgepea 2013: protein molecules per cell -> copy fraction. RNA columns are
# retained in the extracted CSV for provenance but are not an analysis input.
valgepea <- read_csv(
  file.path(input_root, "Valgepea_copy.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
)
valgepea_protein_columns <- names(valgepea)[str_starts(names(valgepea), "Protein")]
valgepea_matrix <- as.matrix(valgepea[valgepea_protein_columns])
storage.mode(valgepea_matrix) <- "numeric"
valgepea_matrix[!is.finite(valgepea_matrix)] <- 0
valgepea_keep <- rowSums(valgepea_matrix) > 0
valgepea <- valgepea[valgepea_keep, , drop = FALSE]
valgepea_matrix <- valgepea_matrix[valgepea_keep, , drop = FALSE]
valgepea_fraction <- copy_fraction(
  as.data.frame(valgepea_matrix),
  replace_nonfinite = TRUE
)
valgepea_growth_rates <- as.numeric(sub("^Protein", "", valgepea_protein_columns))
colnames(valgepea_fraction) <- format_mu(
  "valgepea_2013_prot",
  valgepea_growth_rates
)
valgepea_out <- valgepea |>
  transmute(
    gene = str_split_fixed(as.character(gene), ",", 2)[, 1],
    locus_tag = as.character(KEGG_ID),
    across(any_of(mapping_columns))
  ) |>
  mutate(n_cond = rowSums(!is.na(valgepea_matrix))) |>
  bind_cols(as.data.frame(valgepea_fraction, check.names = FALSE))
saveRDS(
  valgepea_out,
  file.path(output_root, "Proteinvalgepea_2013_pseudo.rds")
)

# Li 2014: compound records were excluded during identifier mapping.
# Low-confidence bracketed measurements are admitted to the published "all"
# dataset; n_cond continues to count high-confidence values.
li <- read_csv(
  file.path(input_root, "Li_copy.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
)
li_conditions <- c(
  "MOPS complete",
  "MOPS minimal",
  "MOPS complete without methionine"
)
li_long <- li |>
  pivot_longer(
    cols = all_of(li_conditions),
    names_to = "condition",
    values_to = "value_raw"
  ) |>
  mutate(
    value_text = as.character(value_raw),
    low_confidence = str_detect(value_text, "\\["),
    value_all = suppressWarnings(
      as.numeric(str_remove_all(value_text, "\\[|\\]"))
    ),
    value_high_confidence = if_else(low_confidence, NA_real_, value_all)
  )

li_row_status <- li_long |>
  group_by(source_row_id, Gene) |>
  summarise(
    keep_high_confidence = !all(is.na(value_high_confidence)),
    keep_all = !all(is.na(value_all) | value_all == 0),
    n_cond = sum(!is.na(value_high_confidence)),
    .groups = "drop"
  ) |>
  filter(keep_high_confidence | keep_all)

li_base <- li |>
  semi_join(li_row_status, by = c("source_row_id", "Gene")) |>
  left_join(
    select(li_row_status, source_row_id, Gene, n_cond),
    by = c("source_row_id", "Gene"),
    relationship = "one-to-one"
  )

li_all_wide <- li_long |>
  semi_join(li_row_status, by = c("source_row_id", "Gene")) |>
  select(source_row_id, Gene, condition, value_all) |>
  pivot_wider(names_from = condition, values_from = value_all) |>
  right_join(
    select(li_base, source_row_id, Gene),
    by = c("source_row_id", "Gene"),
    relationship = "one-to-one"
  )
li_expression <- as.matrix(li_all_wide[li_conditions])
storage.mode(li_expression) <- "numeric"
li_fraction <- copy_fraction(as.data.frame(li_expression))
li_growth_rates <- log(2) / (c(21.5, 56.3, 26.5) / 60)
colnames(li_fraction) <- format_mu("li_2014_all", li_growth_rates)

li_out <- li_base |>
  arrange(match(
    paste(source_row_id, Gene, sep = "\r"),
    paste(li_all_wide$source_row_id, li_all_wide$Gene, sep = "\r")
  )) |>
  transmute(
    gene = as.character(Gene),
    locus_tag = reference_identifier,
    across(any_of(mapping_columns)),
    n_cond
  ) |>
  bind_cols(as.data.frame(li_fraction, check.names = FALSE))
saveRDS(li_out, file.path(output_root, "li_2014_all_pseudo.rds"))

# Goelzer 2015: normalize each technical replicate to a copy fraction, average
# technical replicates within biological replicate, then average biological
# replicates. This is the inherited abundance formula; molecular weight is not
# part of it because the source values are already copies per cell.
goelzer_conditions <- tibble(
  condition = c("PYR", "S", "TS", "CH", "CHG"),
  growth_rate = c(0.3, 0.6, 0.9, 1.1, 1.5)
)

calculate_goelzer_condition <- function(condition) {
  data <- read_csv(
    file.path(input_root, paste0("Goelzer_", condition, ".csv")),
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  technical_columns <- names(data)[
    str_detect(names(data), "^Bio[0-9]+ Tech[0-9]+\\(")
  ]
  copies <- as.matrix(data[technical_columns])
  storage.mode(copies) <- "numeric"
  totals <- colSums(copies, na.rm = TRUE)
  totals[totals == 0] <- NA_real_
  replicate_fractions <- sweep(copies, 2, totals, "/")
  biological_id <- str_match(technical_columns, "^Bio([0-9]+)")[, 2]
  biological_fractions <- sapply(
    split(seq_along(technical_columns), biological_id),
    function(columns) rowMeans(
      replicate_fractions[, columns, drop = FALSE],
      na.rm = TRUE
    )
  )
  if (!is.matrix(biological_fractions)) {
    biological_fractions <- matrix(biological_fractions, ncol = 1L)
  }
  condition_fraction <- rowMeans(biological_fractions, na.rm = TRUE)
  condition_fraction[is.nan(condition_fraction)] <- NA_real_

  data |>
    transmute(
      gene = as.character(`% Gene_Name`),
      locus_tag = as.character(`% BSU Number`),
      source_identifier_raw,
      source_identifier_normalized,
      source_gene_name,
      reference_identifier,
      candidate_reference_identifiers,
      mapping_method,
      mapping_confidence,
      mapping_status,
      mapping_cardinality,
      dataset,
      species,
      !!condition := condition_fraction
    )
}

goelzer_tables <- map(
  goelzer_conditions$condition,
  calculate_goelzer_condition
)
goelzer_join_keys <- setdiff(
  names(goelzer_tables[[1L]]),
  goelzer_conditions$condition
)
goelzer_wide <- reduce(
  goelzer_tables,
  full_join,
  by = goelzer_join_keys,
  relationship = "one-to-one"
)
goelzer_matrix <- as.matrix(goelzer_wide[goelzer_conditions$condition])
storage.mode(goelzer_matrix) <- "numeric"
goelzer_keep <- apply(
  goelzer_matrix,
  1,
  function(values) any(!is.na(values) & values != 0)
)
goelzer_wide <- goelzer_wide[goelzer_keep, , drop = FALSE]
goelzer_original <- goelzer_matrix[goelzer_keep, , drop = FALSE]
goelzer_matrix <- goelzer_original
goelzer_matrix[is.na(goelzer_matrix)] <- 0
goelzer_matrix <- goelzer_matrix + 1e-8
colnames(goelzer_matrix) <- format_mu(
  "Goelzer_2015",
  goelzer_conditions$growth_rate
)
goelzer_out <- goelzer_wide |>
  select(-all_of(goelzer_conditions$condition)) |>
  mutate(n_cond = rowSums(!is.na(goelzer_original) & goelzer_original != 0)) |>
  bind_cols(as.data.frame(goelzer_matrix, check.names = FALSE))
saveRDS(goelzer_out, file.path(output_root, "Goelzer_2015_pseudo.rds"))

message("Wrote five position-independent normalized proteomics datasets to ", output_root)
