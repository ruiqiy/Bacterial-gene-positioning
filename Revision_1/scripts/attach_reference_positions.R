#!/usr/bin/env Rscript

# Attach reference-genome coordinates to normalized rows by stable reference
# identifier. Mapping status is already present and abundance values are never
# recalculated here. Ambiguous mappings, duplicated reference identifiers, and
# missing reference positions are explicitly marked ineligible.

project_root <- normalizePath(
  Sys.getenv(
    "GENE_POSITION_ROOT",
    unset = getwd()
  ),
  winslash = "/"
)
revision_root <- file.path(project_root, "Revision_1")
gff_root <- file.path(revision_root, "data", "raw_data", "reference_gff")
input_root <- file.path(revision_root, "data", "normalized_proteomics")
output_root <- file.path(revision_root, "data", "position_mapped_proteomics")
qc_root <- file.path(revision_root, "qc", "position_mapping")
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
  library(tibble)
})

source(file.path(revision_root, "scripts", "mapping_helpers.R"))

dataset_config <- tribble(
  ~dataset, ~filename, ~species,
  "schmidt_2016", "schmidt_2016_pseudo.rds", "Escherichia coli",
  "peebo_2015", "peebo_2015_pseudo.rds", "Escherichia coli",
  "valgepea_2013_prot", "Proteinvalgepea_2013_pseudo.rds", "Escherichia coli",
  "li_2014_all", "li_2014_all_pseudo.rds", "Escherichia coli",
  "Goelzer_2015", "Goelzer_2015_pseudo.rds", "Bacillus subtilis",
  "zhu_E_coli", "zhu_E_coli.rds", "Escherichia coli",
  "zhu_B_subtilis", "zhu_B_subtilis.rds", "Bacillus subtilis",
  "zhu_V_natriegens", "zhu_V_natriegens.rds", "Vibrio natriegens"
)

reference_config <- tribble(
  ~species, ~reference_strain, ~reference_accession, ~gff_file,
  ~main_seqname, ~chromosome_length, ~oric_start, ~oric_end,
  "Escherichia coli", "K-12 MG1655", "GCF_000005845.2",
  "GCF_000005845.2.gff", "NC_000913.3", 4641652, 3925634, 3926011,
  "Bacillus subtilis", "subsp. subtilis 168", "GCF_000009045.1",
  "GCF_000009045.1.gff", "NC_000964.3", 4215606, 1751, 1938,
  "Vibrio natriegens", "ATCC 14048", "GCF_001456255.1",
  "GCF_001456255.1.gff", "NZ_CP009977.1", 3248023, 2920325, 2920800
)

reference_positions <- pmap_dfr(
  reference_config,
  function(
      species,
      reference_strain,
      reference_accession,
      gff_file,
      main_seqname,
      chromosome_length,
      oric_start,
      oric_end) {
    build_reference_position_map(
      file.path(gff_root, gff_file),
      species,
      reference_strain,
      reference_accession,
      main_seqname,
      chromosome_length,
      oric_start,
      oric_end
    )
  }
)
if (anyDuplicated(reference_positions[c("species", "reference_identifier")])) {
  stop("Reference-position map contains duplicated species/identifier keys")
}
write_csv(
  reference_positions,
  file.path(revision_root, "data", "reference_gene_positions.csv"),
  na = ""
)

qc <- list()
excluded <- list()

for (config_index in seq_len(nrow(dataset_config))) {
  config <- dataset_config[config_index, ]
  normalized <- readRDS(file.path(input_root, config$filename))
  normalized$.normalized_row_number <- seq_len(nrow(normalized))
  abundance_columns <- names(normalized)[grepl("__mu_", names(normalized), fixed = TRUE)]

  required <- c(
    "reference_identifier",
    "mapping_status",
    "mapping_method",
    "mapping_confidence",
    "mapping_cardinality"
  )
  if (!all(required %in% names(normalized))) {
    stop(
      config$dataset,
      ": normalized data lack mapping columns: ",
      paste(setdiff(required, names(normalized)), collapse = ", ")
    )
  }

  mapped <- normalized |>
    mutate(species = config$species) |>
    left_join(
      reference_positions,
      by = c("species", "reference_identifier"),
      relationship = "many-to-one"
    ) |>
    arrange(.normalized_row_number)

  duplicate_counts <- mapped |>
    filter(
      mapping_status == "mapped_one_to_one",
      !is.na(reference_identifier)
    ) |>
    count(reference_identifier, name = "dataset_reference_identifier_rows")

  mapped <- mapped |>
    left_join(
      duplicate_counts,
      by = "reference_identifier",
      relationship = "many-to-one"
    ) |>
    mutate(
      reference_identifier_duplicated = coalesce(
        dataset_reference_identifier_rows > 1L,
        FALSE
      ),
      reference_position_status = if_else(
        is.finite(norm_pos),
        "mapped_reference_position",
        "missing_reference_position"
      ),
      analysis_eligibility_status = case_when(
        mapping_status != "mapped_one_to_one" ~ mapping_status,
        is.na(reference_identifier) ~ "missing_reference_identifier",
        reference_identifier_duplicated ~ "duplicated_reference_identifier",
        reference_position_status != "mapped_reference_position" ~
          "missing_reference_position",
        TRUE ~ "eligible"
      ),
      analysis_eligible = analysis_eligibility_status == "eligible",
      position_mapping_method = if_else(
        reference_position_status == "mapped_reference_position",
        "reference_identifier_to_reference_GFF",
        NA_character_
      )
    ) |>
    select(-.normalized_row_number)

  if (nrow(mapped) != nrow(normalized)) {
    stop(config$dataset, ": position join changed normalized row count")
  }
  if (!all(vapply(
    abundance_columns,
    function(column) identical(mapped[[column]], normalized[[column]]),
    logical(1)
  ))) {
    stop(config$dataset, ": abundance changed while positions were attached")
  }

  excluded[[config$dataset]] <- mapped |>
    mutate(dataset = config$dataset, normalized_row_number = row_number()) |>
    filter(!analysis_eligible) |>
    select(
      dataset,
      normalized_row_number,
      any_of(c("gene", "uniprot_accession", "locus_tag", "source_protein_id")),
      reference_identifier,
      mapping_status,
      reference_identifier_duplicated,
      reference_position_status,
      analysis_eligibility_status
    )

  qc[[config$dataset]] <- tibble(
    dataset = config$dataset,
    normalized_rows = nrow(mapped),
    mapped_one_to_one = sum(mapped$mapping_status == "mapped_one_to_one"),
    unmatched = sum(mapped$mapping_status == "unmatched"),
    ambiguous_one_to_many = sum(
      mapped$mapping_status == "ambiguous_one_to_many"
    ),
    duplicated_reference_identifier_rows = sum(
      mapped$reference_identifier_duplicated
    ),
    missing_reference_position_rows = sum(
      mapped$mapping_status == "mapped_one_to_one" &
        mapped$reference_position_status == "missing_reference_position"
    ),
    analysis_eligible = sum(mapped$analysis_eligible),
    analysis_eligible_fraction = mean(mapped$analysis_eligible),
    abundance_columns = length(abundance_columns),
    abundance_unchanged = TRUE
  )

  saveRDS(mapped, file.path(output_root, config$filename))
}

write_csv(
  bind_rows(qc),
  file.path(qc_root, "reference_position_mapping_qc.csv"),
  na = ""
)
write_csv(
  bind_rows(excluded),
  file.path(qc_root, "position_analysis_exclusions.csv"),
  na = ""
)
message("Attached reference positions by stable identifier without changing abundance")
