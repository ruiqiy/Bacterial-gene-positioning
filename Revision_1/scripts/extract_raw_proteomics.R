#!/usr/bin/env Rscript

# Extract the exact manifest-selected regions from the immutable supplementary
# workbooks. The generated CSV files are machine-readable inputs for the copied
# normalization scripts; the original workbooks are never modified.

project_root <- normalizePath(
  Sys.getenv("GENE_POSITION_ROOT", unset = getwd()),
  winslash = "/"
)
revision_root <- file.path(project_root, "Revision_1")
raw_root <- file.path(revision_root, "data", "raw_data", "proteomics")
out_root <- file.path(revision_root, "data", "raw_extracted")
qc_root <- file.path(revision_root, "qc", "raw_rebuild")

.libPaths(c(
  file.path(revision_root, "R_libs"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
})

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

manifest <- read_csv(
  file.path(raw_root, "manifest.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
)
names(manifest)[1] <- sub("^\ufeff", "", names(manifest)[1])
manifest <- manifest |>
  mutate(across(everything(), ~ trimws(as.character(.x))))

read_manifest_blocks <- function(manifest_row) {
  ranges <- trimws(strsplit(manifest_row$range, ",", fixed = TRUE)[[1]])
  blocks <- map(
    ranges,
    ~ read_excel(
      file.path(raw_root, manifest_row$source_file),
      sheet = manifest_row$sheet,
      range = .x,
      col_names = FALSE,
      .name_repair = "minimal"
    ) |>
      as.data.frame(check.names = FALSE)
  )
  row_counts <- vapply(blocks, nrow, integer(1))
  if (length(unique(row_counts)) != 1L) {
    stop(
      manifest_row$dataset, " / ", manifest_row$sheet,
      ": discontiguous manifest ranges have unequal row counts"
    )
  }
  setNames(
    bind_cols(blocks, .name_repair = "unique"),
    paste0("V", seq_len(sum(vapply(blocks, ncol, integer(1)))))
  )
}

drop_blank_rows <- function(data) {
  keep <- apply(
    data,
    1,
    function(row) any(!is.na(row) & trimws(as.character(row)) != "")
  )
  data[keep, , drop = FALSE]
}

write_extract <- function(data, filename) {
  write_csv(data, file.path(out_root, filename), na = "")
  file.path(out_root, filename)
}

qc_rows <- list()
record_qc <- function(manifest_row, data, output_file, id_column) {
  id <- as.character(data[[id_column]])
  qc_rows[[length(qc_rows) + 1L]] <<- tibble(
    dataset = manifest_row$dataset,
    source_file = manifest_row$source_file,
    sheet = manifest_row$sheet,
    selected_range = manifest_row$range,
    output_file = basename(output_file),
    extracted_rows = nrow(data),
    extracted_columns = ncol(data),
    missing_identifier_rows = sum(is.na(id) | trimws(id) == ""),
    duplicated_nonmissing_identifiers = sum(
      duplicated(id) & !is.na(id) & trimws(id) != ""
    )
  )
}

# Schmidt: A and C contain identifiers; G:AB contains protein copies/cell.
# Condition metadata is extracted independently from Table S23. The manifest
# header_row values are worksheet row numbers; because the selected metadata
# range begins on its header row, that header is row 1 of the extracted block.
schmidt_cfg <- manifest |> filter(dataset == "Schmidt")
stopifnot(nrow(schmidt_cfg) == 1L)
schmidt_block <- read_manifest_blocks(schmidt_cfg)
schmidt_headers <- trimws(as.character(schmidt_block[1, ]))
schmidt_selected <- schmidt_block[-1, , drop = FALSE]
names(schmidt_selected) <- schmidt_headers
schmidt_selected <- drop_blank_rows(schmidt_selected) |>
  transmute(
    `Uniprot Accession` = as.character(`Uniprot Accession`),
    Gene = as.character(Gene),
    across(-c(`Uniprot Accession`, Gene), as.character)
  )
schmidt_selected_file <- write_extract(
  schmidt_selected,
  "schmidt_raw_selected.csv"
)

schmidt_meta_cfg <- manifest |> filter(dataset == "Schmidt_metadata")
stopifnot(nrow(schmidt_meta_cfg) == 1L)
schmidt_meta_block <- read_manifest_blocks(schmidt_meta_cfg)
schmidt_meta_headers <- trimws(as.character(schmidt_meta_block[1, ]))
schmidt_meta_raw <- schmidt_meta_block[-1, , drop = FALSE]
names(schmidt_meta_raw) <- schmidt_meta_headers

required_schmidt_meta <- c("Growth condition", "Strain", "Growth rate (h-1)")
if (!all(required_schmidt_meta %in% names(schmidt_meta_raw))) {
  stop(
    "Schmidt metadata extract is missing required columns: ",
    paste(setdiff(required_schmidt_meta, names(schmidt_meta_raw)), collapse = ", ")
  )
}

schmidt_meta <- drop_blank_rows(schmidt_meta_raw) |>
  transmute(
    `Growth condition` = str_squish(as.character(`Growth condition`)),
    Strain = str_squish(as.character(Strain)),
    `Growth rate (h-1)` = suppressWarnings(
      as.numeric(as.character(`Growth rate (h-1)`))
    )
  ) |>
  filter(Strain == "BW25113", `Growth rate (h-1)` > 0) |>
  mutate(
    `Growth condition` = if_else(
      `Growth condition` == "Osmotic-stress glucose3",
      "Osmotic-stress glucose",
      `Growth condition`
    )
  ) |>
  select(`Growth condition`, `Growth rate (h-1)`)

if (nrow(schmidt_meta) != 20L) {
  stop("Expected 20 positive-growth BW25113 Schmidt conditions; found ", nrow(schmidt_meta))
}
if (anyDuplicated(schmidt_meta$`Growth condition`)) {
  stop("Filtered Schmidt metadata contains duplicated growth conditions")
}

schmidt_conditions <- as.character(schmidt_meta[[1]])
normalize_condition <- function(x) {
  str_to_lower(str_squish(as.character(x)))
}
condition_lookup <- setNames(
  schmidt_conditions,
  normalize_condition(schmidt_conditions)
)
available_conditions <- setdiff(
  names(schmidt_selected),
  c("Uniprot Accession", "Gene")
)
selected_condition_names <- unname(
  condition_lookup[normalize_condition(available_conditions)]
)
keep_condition <- !is.na(selected_condition_names)
schmidt_input <- schmidt_selected |>
  select(
    `Uniprot Accession`,
    Gene,
    all_of(available_conditions[keep_condition])
  )
names(schmidt_input)[-(1:2)] <- selected_condition_names[keep_condition]
schmidt_input_file <- write_extract(
  schmidt_input,
  "schmidt_data_copy.csv"
)
schmidt_meta_file <- write_extract(schmidt_meta, "schmidt_meta.csv")
record_qc(schmidt_cfg, schmidt_input, schmidt_input_file, "Uniprot Accession")
record_qc(
  schmidt_meta_cfg,
  schmidt_meta,
  schmidt_meta_file,
  "Growth condition"
)

# Peebo: the selected G:AZ region contains two parallel 23-column blocks.
# Retain the protein-concentration block for normalization and write the full
# selected region separately for provenance.
peebo_cfg <- manifest |> filter(dataset == "Peebo")
stopifnot(nrow(peebo_cfg) == 1L)
peebo_block <- read_manifest_blocks(peebo_cfg)
peebo_growth_headers <- trimws(as.character(peebo_block[1, ]))
peebo_data <- drop_blank_rows(peebo_block[-1, , drop = FALSE])
concentration_count <- (ncol(peebo_data) - 2L) / 2L
if (concentration_count != 23L) {
  stop("Unexpected Peebo selected-region width")
}
concentration_cols <- 3:(2 + concentration_count)
cost_cols <- (3 + concentration_count):ncol(peebo_data)
if (!identical(
  peebo_growth_headers[concentration_cols],
  peebo_growth_headers[cost_cols]
)) {
  stop("Peebo concentration and expression-cost growth-rate headers differ")
}
peebo_selected <- as_tibble(peebo_data)
names(peebo_selected) <- c(
  "Gene",
  "BID",
  paste0(
    "ProteinConcentration_",
    peebo_growth_headers[concentration_cols],
    "_",
    seq_len(concentration_count)
  ),
  paste0(
    "ProteinExpressionCost_",
    peebo_growth_headers[cost_cols],
    "_",
    seq_len(concentration_count)
  )
)
write_extract(peebo_selected, "peebo_raw_selected.csv")
peebo_input <- peebo_selected |>
  select(Gene, BID, starts_with("ProteinConcentration_"))
names(peebo_input)[-(1:2)] <- peebo_growth_headers[concentration_cols]
peebo_input_file <- write_extract(peebo_input, "peebo_copy.csv")
record_qc(peebo_cfg, peebo_input, peebo_input_file, "BID")

# Valgepea: combine rows 10 (growth rate) and 12 (measurement type) into a
# single header. Column A/B names come directly from row 12.
valgepea_cfg <- manifest |> filter(dataset == "Valgepea")
stopifnot(nrow(valgepea_cfg) == 1L)
valgepea_block <- read_manifest_blocks(valgepea_cfg)
if (nrow(valgepea_block) < 4L) stop("Valgepea range contains no data rows")
growth_header <- trimws(as.character(valgepea_block[1, ]))
measurement_header <- trimws(as.character(valgepea_block[3, ]))
valgepea_names <- character(ncol(valgepea_block))
valgepea_names[1:2] <- c("KEGG_ID", "gene")
for (j in 3:ncol(valgepea_block)) {
  prefix <- case_when(
    str_detect(measurement_header[j], regex("^mRNA abun\\.$", TRUE)) ~ "RNA",
    str_detect(measurement_header[j], regex("^Protein abun\\.$", TRUE)) ~ "Protein",
    TRUE ~ NA_character_
  )
  if (is.na(prefix)) {
    stop("Unexpected Valgepea measurement header: ", measurement_header[j])
  }
  valgepea_names[j] <- paste0(prefix, growth_header[j])
}
valgepea_input <- valgepea_block[-(1:3), , drop = FALSE]
names(valgepea_input) <- valgepea_names
valgepea_input <- drop_blank_rows(valgepea_input) |>
  as_tibble()
valgepea_input_file <- write_extract(
  valgepea_input,
  "Valgepea_copy.csv"
)
record_qc(valgepea_cfg, valgepea_input, valgepea_input_file, "KEGG_ID")

# Li: the raw selected range is already a one-row-header table.
li_cfg <- manifest |> filter(dataset == "Li")
stopifnot(nrow(li_cfg) == 1L)
li_block <- read_manifest_blocks(li_cfg)
li_headers <- trimws(as.character(li_block[1, ]))
li_input <- li_block[-1, , drop = FALSE]
names(li_input) <- li_headers
li_input <- drop_blank_rows(li_input) |> as_tibble()
li_input_file <- write_extract(li_input, "Li_copy.csv")
record_qc(li_cfg, li_input, li_input_file, "Gene")

# Goelzer: write one manifest-trimmed CSV for each condition sheet. The
# normalization script reads the technical-replicate columns directly.
goelzer_cfg <- manifest |> filter(dataset == "Goelzer")
for (i in seq_len(nrow(goelzer_cfg))) {
  cfg <- goelzer_cfg[i, ]
  block <- read_manifest_blocks(cfg)
  headers <- trimws(as.character(block[1, ]))
  data <- block[-1, , drop = FALSE]
  names(data) <- make.unique(headers, sep = "_")
  data <- drop_blank_rows(data) |> as_tibble()
  output_file <- write_extract(
    data,
    paste0("Goelzer_", cfg$sheet, ".csv")
  )
  record_qc(cfg, data, output_file, "% BSU Number")
}

# Zhu: each selected sheet is already a one-row-header table. Preserve both
# protein ID and locus tag in the extracted CSV.
zhu_cfg <- manifest |> filter(str_starts(dataset, "Zhu_"))
for (i in seq_len(nrow(zhu_cfg))) {
  cfg <- zhu_cfg[i, ]
  block <- read_manifest_blocks(cfg)
  headers <- trimws(as.character(block[1, ]))
  data <- block[-1, , drop = FALSE]
  names(data) <- headers
  data <- drop_blank_rows(data) |> as_tibble()
  output_file <- write_extract(
    data,
    paste0(str_replace_all(cfg$sheet, " ", "_"), ".csv")
  )
  record_qc(cfg, data, output_file, "gene locus")
}

extraction_qc <- bind_rows(qc_rows)
write_csv(
  extraction_qc,
  file.path(qc_root, "raw_extraction_qc.csv"),
  na = ""
)

message(
  "Extracted ", nrow(extraction_qc),
  " manifest rows into ", out_root
)
