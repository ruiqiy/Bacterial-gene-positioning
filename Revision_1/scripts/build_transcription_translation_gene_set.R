project_root <- normalizePath(
  Sys.getenv("GENE_POSITION_ROOT", unset = getwd()),
  winslash = "/"
)
revision_root <- file.path(project_root, "Revision_1")
gff_root <- file.path(revision_root, "data", "raw_data", "reference_gff")

.libPaths(c(
  file.path(revision_root, "R_libs"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

source(file.path(revision_root, "scripts", "mapping_helpers.R"))

reference_config <- tibble::tribble(
  ~species, ~strain, ~reference_genome_accession, ~gff_file,
  "Escherichia coli", "K-12 MG1655", "GCF_000005845.2", "GCF_000005845.2.gff",
  "Bacillus subtilis", "subsp. subtilis 168", "GCF_000009045.1", "GCF_000009045.1.gff",
  "Vibrio natriegens", "ATCC 14048", "GCF_001456255.1", "GCF_001456255.1.gff"
)

is_ribosomal_protein <- function(product) {
  annotated_as_subunit <- str_detect(
    product,
    regex(
      paste0(
        "^(?:",
        "(?:30S|50S|type [AB] 50S) ribosomal (?:subunit )?protein [LS][0-9]",
        "|ribosomal protein B?[LS][0-9]",
        "|K-turn RNA binding protein; alternative ribosomal protein",
        "|alternative ribosomal protein",
        "|RNA degradation presenting factor \\(ribosomal protein S1 homolog\\)",
        ")"
      ),
      ignore_case = TRUE
    )
  )
  modification_enzyme <- str_detect(
    product,
    regex(
      "methyl.*transferase|acetyltransferase|ligase|hydroxylase|modification protein|accessory factor|protease",
      ignore_case = TRUE
    )
  )
  annotated_as_subunit & !modification_enzyme
}

is_structural_rna_polymerase_subunit <- function(product) {
  str_detect(
    product,
    regex(
      paste0(
        "^(?:",
        "(?:DNA-directed )?RNA polymerase subunit (?:alpha|beta|beta'|omega)",
        "|RNA polymerase \\((?:alpha|beta|beta'|delta subunit)",
        "|omega(?: 1)? subunit of RNA polymerase",
        ")"
      ),
      ignore_case = TRUE
    )
  )
}

reference_gene_set <- purrr::pmap_dfr(
  reference_config,
  function(species, strain, reference_genome_accession, gff_file) {
    cds <- read_gff_features(file.path(gff_root, gff_file)) |>
      filter(feature == "CDS", !is.na(locus_tag), !is.na(product)) |>
      distinct(locus_tag, .keep_all = TRUE) |>
      mutate(
        category = case_when(
          is_structural_rna_polymerase_subunit(product) ~ "RNA polymerase",
          is_ribosomal_protein(product) ~ "ribosomal protein",
          TRUE ~ NA_character_
        )
      ) |>
      filter(!is.na(category))

    cds |>
      transmute(
        species,
        strain,
        reference_genome_accession,
        reference_identifier = normalize_reference_identifier(species, locus_tag),
        gene_name = gene,
        category,
        product_annotation = product,
        annotation_source = paste0("NCBI RefSeq ", reference_genome_accession, " GFF CDS product"),
        inclusion_rule = if_else(
          category == "RNA polymerase",
          "curated structural RNA-polymerase subunit product; sigma factors excluded",
          "curated product annotation identifying an actual small- or large-subunit ribosomal protein"
        ),
        mapping_method = "reference-locus annotation",
        mapping_confidence = "high"
      )
  }
) |>
  arrange(species, category, reference_identifier)

if (anyDuplicated(reference_gene_set[c("species", "reference_identifier")])) {
  stop("Reference transcription/translation gene-set identifiers are duplicated")
}

gene_set_qc <- reference_gene_set |>
  count(species, category, name = "reference_gene_count") |>
  tidyr::complete(
    species = reference_config$species,
    category = c("RNA polymerase", "ribosomal protein"),
    fill = list(reference_gene_count = 0)
  )

write_csv(
  reference_gene_set,
  file.path(revision_root, "data", "reference_transcription_translation_genes.csv")
)
write_csv(
  gene_set_qc,
  file.path(revision_root, "qc", "transcription_translation_reference_counts.csv")
)

message(
  "Wrote ",
  nrow(reference_gene_set),
  " reference transcription/translation genes."
)
