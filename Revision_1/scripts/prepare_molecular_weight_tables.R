#!/usr/bin/env Rscript

# Build accession-keyed Zhu molecular-weight tables exclusively from the
# downloaded UniProtKB ID-mapping results. Proteins without a finite UniProt
# mass remain unresolved; this script does not read FASTA files or calculate
# molecular mass from amino-acid sequences.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) {
  normalizePath(sub("^--file=", "", script_arg), winslash = "/")
} else {
  normalizePath(
    "Revision_1/scripts/prepare_molecular_weight_tables.R",
    winslash = "/"
  )
}

default_project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/"
)
project_root <- normalizePath(
  Sys.getenv("GENE_POSITION_ROOT", unset = default_project_root),
  winslash = "/"
)
revision_root <- file.path(project_root, "Revision_1")
extracted_root <- file.path(revision_root, "data", "raw_extracted")
raw_mw_root <- file.path(revision_root, "data", "raw_data", "molecular_weight")
mw_root <- file.path(revision_root, "data", "molecular_weight")
dir.create(mw_root, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(
  file.path(revision_root, "R_libs"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

uniprot_release <- "2026_02"
uniprot_release_date <- "2026-06-10"
retrieval_date <- Sys.Date()

output_fields <- c(
  "dataset",
  "source_protein_id",
  "source_gene_name",
  "source_locus_identifier",
  "database_accession",
  "entry_name",
  "length_aa",
  "molecular_weight_Da",
  "molecular_weight_kDa",
  "mw_status",
  "mw_source",
  "database_release",
  "database_release_date",
  "retrieval_date"
)

prepare_uniprot <- function(
    dataset,
    raw_csv,
    preferred_taxon_id = NA_integer_,
    source_label = "UniProtKB_sequence_mass") {
  source_rows <- read_csv(
    file.path(extracted_root, raw_csv),
    show_col_types = FALSE,
    name_repair = "minimal"
  ) |>
    transmute(
      source_protein_id = str_trim(as.character(`protein ID`)),
      source_gene_name = as.character(`gene name`),
      source_locus_identifier = as.character(`gene locus`)
    )

  if (anyDuplicated(source_rows$source_protein_id)) {
    stop(dataset, ": duplicated source protein IDs in extracted input")
  }

  raw_results <- read_tsv(
    file.path(
      raw_mw_root,
      paste0(dataset, "_uniprot_molecular_weight_raw.tsv")
    ),
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  if (!"Organism (ID)" %in% names(raw_results)) {
    raw_results[["Organism (ID)"]] <- NA_integer_
  }
  if (!"Organism" %in% names(raw_results)) {
    raw_results[["Organism"]] <- NA_character_
  }

  result_rows <- raw_results |>
    transmute(
      source_protein_id = str_trim(as.character(From)),
      database_accession = as.character(Entry),
      entry_name = as.character(`Entry Name`),
      length_aa = suppressWarnings(as.integer(Length)),
      molecular_weight_Da = suppressWarnings(as.numeric(Mass)),
      organism_id = suppressWarnings(as.integer(`Organism (ID)`)),
      organism_name = as.character(Organism)
    ) |>
    filter(!is.na(source_protein_id), source_protein_id != "")

  conflicting_masses <- result_rows |>
    filter(is.finite(molecular_weight_Da)) |>
    distinct(source_protein_id, length_aa, molecular_weight_Da) |>
    count(source_protein_id, name = "distinct_length_mass_pairs") |>
    filter(distinct_length_mass_pairs > 1L)
  if (nrow(conflicting_masses) > 0L) {
    stop(
      dataset,
      ": conflicting UniProt length/mass values for source accessions: ",
      paste(head(conflicting_masses$source_protein_id, 20L), collapse = ", ")
    )
  }

  # RefSeq WP accessions can point to several organism-specific UniProt
  # entries. Prefer a finite mass, then the focal taxon when supplied, then a
  # deterministic accession. The conflict check above guarantees that this
  # choice cannot change the selected mass.
  preferred_results <- result_rows |>
    mutate(
      finite_mass_priority = if_else(is.finite(molecular_weight_Da), 1L, 2L),
      taxon_priority = if_else(
        !is.na(preferred_taxon_id) & organism_id == preferred_taxon_id,
        1L,
        2L,
        missing = 2L
      )
    ) |>
    arrange(
      source_protein_id,
      finite_mass_priority,
      taxon_priority,
      database_accession
    ) |>
    group_by(source_protein_id) |>
    slice(1L) |>
    ungroup() |>
    select(
      source_protein_id,
      database_accession,
      entry_name,
      length_aa,
      molecular_weight_Da
    )

  source_rows |>
    left_join(
      preferred_results,
      by = "source_protein_id",
      relationship = "one-to-one"
    ) |>
    mutate(
      dataset = dataset,
      molecular_weight_Da = round(molecular_weight_Da, 3),
      molecular_weight_kDa = round(molecular_weight_Da / 1000, 6),
      mw_status = case_when(
        is.finite(molecular_weight_Da) ~ "mapped",
        is.na(database_accession) ~ "not_mapped_to_UniProtKB",
        TRUE ~ "UniProtKB_entry_missing_mass"
      ),
      mw_source = if_else(
        mw_status == "mapped",
        source_label,
        NA_character_
      ),
      database_release = uniprot_release,
      database_release_date = uniprot_release_date,
      retrieval_date = retrieval_date
    ) |>
    select(all_of(output_fields))
}

tables <- list(
  zhu_E_coli = prepare_uniprot(
    "zhu_E_coli",
    "Ecoli_batch_F.csv"
  ),
  zhu_B_subtilis = prepare_uniprot(
    "zhu_B_subtilis",
    "Bsub_batch_F.csv"
  ),
  zhu_V_natriegens = prepare_uniprot(
    "zhu_V_natriegens",
    "Vnat_batch_C.csv",
    preferred_taxon_id = 1219067L,
    source_label = "UniProtKB_sequence_mass_via_RefSeq_cross_reference"
  )
)

for (dataset in names(tables)) {
  write_csv(
    tables[[dataset]],
    file.path(mw_root, paste0(dataset, "_molecular_weight.csv")),
    na = ""
  )
}

all_rows <- bind_rows(tables)
unresolved <- filter(all_rows, mw_status != "mapped")
source_qc <- all_rows |>
  mutate(mw_source = coalesce(mw_source, "unresolved_no_UniProt_mass")) |>
  count(dataset, mw_source, name = "rows") |>
  arrange(dataset, desc(rows), mw_source)
qc <- bind_rows(lapply(names(tables), function(dataset) {
  data <- tables[[dataset]]
  tibble(
    dataset = dataset,
    raw_rows = nrow(data),
    unique_source_protein_ids = n_distinct(data$source_protein_id),
    mw_mapped = sum(data$mw_status == "mapped"),
    mw_unmapped = sum(data$mw_status != "mapped"),
    mw_mapping_fraction = round(mean(data$mw_status == "mapped"), 8)
  )
}))

write_csv(all_rows, file.path(mw_root, "zhu_molecular_weight_all.csv"), na = "")
write_csv(qc, file.path(mw_root, "molecular_weight_download_qc.csv"), na = "")
write_csv(
  source_qc,
  file.path(mw_root, "molecular_weight_source_qc.csv"),
  na = ""
)
write_csv(
  unresolved,
  file.path(mw_root, "molecular_weight_unresolved.csv"),
  na = ""
)

message("Wrote ", nrow(all_rows), " UniProt-only molecular-weight rows to ", mw_root)
for (i in seq_len(nrow(qc))) {
  message(
    qc$dataset[i],
    " ", qc$mw_mapped[i], "/", qc$raw_rows[i],
    " (", sprintf("%.3f", 100 * qc$mw_mapping_fraction[i]), "%)"
  )
}
