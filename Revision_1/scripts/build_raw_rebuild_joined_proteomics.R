#!/usr/bin/env Rscript

project_root <- normalizePath(
  Sys.getenv("GENE_POSITION_ROOT", unset = getwd()),
  winslash = "/"
)
revision_root <- file.path(project_root, "Revision_1")
input_root <- file.path(revision_root, "data", "position_mapped_proteomics")
joined_root <- file.path(revision_root, "data", "analysis_ready_proteomics")
qc_root <- file.path(revision_root, "qc", "analysis_ready")

.libPaths(c(file.path(revision_root, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

source(file.path(revision_root, "scripts", "mapping_helpers.R"))

dir.create(joined_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

dataset_config <- tribble(
  ~dataset,             ~input_file,                      ~species,            ~source_strain,       ~source_identifier_column,
  "schmidt_2016",       "schmidt_2016_pseudo.rds",        "Escherichia coli",  "K-12 BW25113",       "uniprot_accession",
  "peebo_2015",         "peebo_2015_pseudo.rds",          "Escherichia coli",  "K-12 BW25113",       "locus_tag",
  "li_2014_all",        "li_2014_all_pseudo.rds",         "Escherichia coli",  "K-12 MG1655",        "gene",
  "valgepea_2013_prot", "Proteinvalgepea_2013_pseudo.rds", "Escherichia coli",  "K-12 MG1655",        "locus_tag",
  "Goelzer_2015",       "Goelzer_2015_pseudo.rds",        "Bacillus subtilis", "subsp. subtilis 168", "locus_tag",
  "zhu_B_subtilis",     "zhu_B_subtilis.rds",             "Bacillus subtilis", "subsp. subtilis 168", "locus_tag",
  "zhu_E_coli",         "zhu_E_coli.rds",                 "Escherichia coli",  "K-12 NCM3722",       "locus_tag",
  "zhu_V_natriegens",   "zhu_V_natriegens.rds",           "Vibrio natriegens", "ATCC 14048",        "locus_tag"
)

# Load raw PanKB data and construct core vs accessory categories
pankb_dir <- file.path(revision_root, "data", "raw_data", "PanKB")
species_pankb_files <- tribble(
  ~species,            ~focal_file,                 ~annot_file,
  "Escherichia coli",  "Ecoli_focal_gene_info.csv", "Ecoli_gene_annotations.csv",
  "Bacillus subtilis", "Bsub_focal_gene_info.csv",  "Bsub_Ecoli_gene_annotations.csv",
  "Vibrio natriegens", "Vnat_focal_gene_info.csv",  "Vnat_Ecoli_gene_annotations.csv"
)

pankb <- map_dfr(seq_len(nrow(species_pankb_files)), function(i) {
  sp <- species_pankb_files$species[i]
  focal <- read_csv(file.path(pankb_dir, species_pankb_files$focal_file[i]), show_col_types = FALSE)
  annot <- read_csv(file.path(pankb_dir, species_pankb_files$annot_file[i]), show_col_types = FALSE)

  focal |>
    filter(!is.na(original_locus_tag), original_locus_tag != "") |>
    inner_join(annot, by = "gene", suffix = c("_focal", "_annot")) |>
    mutate(
      species = sp,
      pankb_three_class = tolower(pangenomic_class),
      binary_core_accessory = if_else(pankb_three_class == "core", "core", "accessory"),
      reference_identifier = if_else(
        species == "Bacillus subtilis",
        str_replace(original_locus_tag, "^BSU_", "BSU"),
        original_locus_tag
      ),
      original_exact_match = normalize_logical(original_exact_match)
    )
}) |>
  transmute(
    species,
    reference_identifier,
    pankb_reference_genome_accession = genome_id,
    pankb_three_class,
    binary_core_accessory,
    prevalence_count = frequency,
    original_exact_match,
    pankb_family = gene,
    pankb_member_locus_tag = locus_tag,
    pankb_original_locus_tag = original_locus_tag
  ) |>
  arrange(species, reference_identifier, desc(original_exact_match), pankb_member_locus_tag) |>
  distinct(species, reference_identifier, .keep_all = TRUE)

# Load reference transcription and translation genes
reference_gene_set <- read_csv(
  file.path(revision_root, "data", "reference_transcription_translation_genes.csv"),
  show_col_types = FALSE
) |>
  select(
    species,
    reference_identifier,
    transcription_translation_category = category,
    transcription_translation_product = product_annotation
  )

join_audit <- list()

for (i in seq_len(nrow(dataset_config))) {
  cfg <- dataset_config[i, ]
  dataset_name <- cfg$dataset
  x <- readRDS(file.path(input_root, cfg$input_file))
  input_rows <- nrow(x)
  x$source_row_number <- seq_len(input_rows)

  x <- x |>
    filter(analysis_eligible) |>
    mutate(analysis_row_number = row_number())
  eligible_rows <- nrow(x)

  join_audit[[length(join_audit) + 1L]] <- tibble(
    dataset = dataset_name,
    stage = "identifier_eligibility_filter",
    rows_before = input_rows,
    rows_after = eligible_rows,
    row_count_change = eligible_rows - input_rows
  )

  normalized_identifier <- if (dataset_name == "li_2014_all") {
    normalize_gene_name(x$gene)
  } else {
    normalize_reference_identifier(cfg$species, x[[cfg$source_identifier_column]])
  }

  x <- x |>
    mutate(
      original_row_number = source_row_number,
      dataset = dataset_name,
      species = cfg$species,
      source_strain = cfg$source_strain,
      original_gene_name = as.character(gene),
      gene_common_name = as.character(gene),
      source_gene_name_normalized = normalize_gene_name(gene),
      original_gene_identifier = as.character(.data[[cfg$source_identifier_column]]),
      normalized_gene_identifier = normalized_identifier,
      cross_strain_mapping_method = mapping_method,
      cross_strain_mapping_confidence = mapping_confidence,
      cross_strain_mapping_status = mapping_status,
      cross_strain_mapping_cardinality = mapping_cardinality,
      cross_strain_evidence = candidate_reference_identifiers,
      gene = reference_identifier
    )

  # Join PanKB core/accessory annotations
  before_pankb <- nrow(x)
  x <- left_join(x, pankb, by = c("species", "reference_identifier")) |>
    mutate(
      pankb_match_status = case_when(
        is.na(reference_identifier) ~ "no_reference_identifier",
        is.na(pankb_three_class) ~ "reference_identifier_not_in_PanKB_crosswalk",
        TRUE ~ "matched"
      ),
      pankb_three_class_main = pankb_three_class,
      binary_core_accessory_main = binary_core_accessory
    )

  join_audit[[length(join_audit) + 1L]] <- tibble(
    dataset = dataset_name,
    stage = "PanKB_annotation",
    rows_before = before_pankb,
    rows_after = nrow(x),
    row_count_change = nrow(x) - before_pankb
  )

  # Join transcription/translation gene set annotations
  before_gene_set <- nrow(x)
  x <- left_join(x, reference_gene_set, by = c("species", "reference_identifier")) |>
    mutate(
      is_transcription_translation = case_when(
        is.na(reference_identifier) ~ NA,
        !is.na(transcription_translation_category) ~ TRUE,
        TRUE ~ FALSE
      ),
      is_other_gene = case_when(
        is.na(is_transcription_translation) ~ NA,
        TRUE ~ !is_transcription_translation
      ),
      transcription_translation_match_status = case_when(
        is.na(reference_identifier) ~ "unresolved_reference_mapping",
        is_transcription_translation ~ "member",
        TRUE ~ "mapped_nonmember"
      )
    )

  join_audit[[length(join_audit) + 1L]] <- tibble(
    dataset = dataset_name,
    stage = "transcription_translation_annotation",
    rows_before = before_gene_set,
    rows_after = nrow(x),
    row_count_change = nrow(x) - before_gene_set
  )

  saveRDS(x, file.path(joined_root, paste0(dataset_name, "_joined.rds")))
}

write_csv(bind_rows(join_audit), file.path(qc_root, "analysis_ready_join_row_count_audit.csv"))

message("Wrote annotated analysis-ready RDS files for ", nrow(dataset_config), " datasets.")
