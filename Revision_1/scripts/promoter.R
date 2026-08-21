project_root <- normalizePath(
  Sys.getenv("GENE_POSITION_ROOT", unset = getwd()),
  winslash = "/"
)
revision_root <- file.path(project_root, "Revision_1")
promoter_root <- file.path(revision_root, "data", "promoter")
raw_root <- file.path(revision_root, "data", "raw_data", "promoter")
gff_root <- file.path(revision_root, "data", "raw_data", "reference_gff")
mapped_root <- file.path(promoter_root, "mapped")
qc_root <- file.path(revision_root, "qc")

.libPaths(c(
  file.path(revision_root, "R_libs"),
  .libPaths()
))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# The promoter coordinates and the original promoter.R were built against
# U00096.2. Keep its exact circular length so promoter-to-gene assignment is
# unchanged; only the value assigned to a promoter changes from gene name to
# the stable MG1655 b-number.
genome_len <- 4639675

parse_gff_identifiers <- function(gff_file) {
  gff <- read_tsv(
    gff_file,
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE
  )

  genes <- gff |>
    filter(X3 == "gene") |>
    select(start = X4, end = X5, strand = X7, attributes = X9) |>
    mutate(
      reference_identifier = str_extract(attributes, "locus_tag=[^;]+"),
      reference_identifier = str_remove(reference_identifier, "locus_tag="),
      reference_identifier = if_else(
        is.na(reference_identifier),
        str_remove(str_extract(attributes, "ID=[^;]+"), "ID=gene-"),
        reference_identifier
      )
    ) |>
    select(reference_identifier, start, end, gene_strand = strand)

  if (any(is.na(genes$reference_identifier))) {
    stop("U00096.2 GFF gene rows are missing stable reference identifiers")
  }
  if (anyDuplicated(genes$reference_identifier)) {
    stop("U00096.2 GFF contains duplicated gene reference identifiers")
  }
  genes
}

# This preserves the original promoter.R assignment logic:
#   1. expand to every gene whose body contains the promoter peak;
#   2. otherwise choose the first nearest downstream gene within 500 bp,
#      using promoter strand and circular-genome distance;
#   3. retain an unmapped row when neither rule finds a gene.
map_promoters_to_identifiers_expanded <- function(promoter_file, genes_df) {
  promoters <- read_csv(promoter_file, show_col_types = FALSE)
  required <- c(
    "name",
    "left_boundary_coordinate",
    "right_boundary_coordinate",
    "peak_coordinate",
    "peak_activity"
  )
  if (!all(required %in% names(promoters))) {
    stop(
      basename(promoter_file),
      " is missing required columns: ",
      paste(setdiff(required, names(promoters)), collapse = ", ")
    )
  }

  if (!"strand" %in% names(promoters)) {
    promoters <- promoters |>
      mutate(strand = str_sub(name, -1, -1))
  }
  if (any(!promoters$strand %in% c("+", "-"))) {
    stop(basename(promoter_file), " contains invalid promoter strands")
  }

  promoters$.source_row_index <- seq_len(nrow(promoters))
  results_list <- vector("list", nrow(promoters))

  for (i in seq_len(nrow(promoters))) {
    row <- promoters[i, ]
    peak <- row$peak_coordinate
    promoter_strand <- row$strand

    matched_genes <- data.frame(
      reference_identifier = character(),
      mapping_type = character(),
      distance = numeric(),
      stringsAsFactors = FALSE
    )

    intragenic_matches <- genes_df |>
      filter(start <= peak, end >= peak)

    if (nrow(intragenic_matches) > 0) {
      current_matches <- data.frame(
        reference_identifier = intragenic_matches$reference_identifier,
        mapping_type = "Intragenic",
        distance = 0
      )
      matched_genes <- bind_rows(matched_genes, current_matches)
    }

    if (nrow(matched_genes) == 0) {
      if (promoter_strand == "+") {
        candidates <- genes_df |>
          mutate(dist = (start - peak) %% genome_len) |>
          filter(dist > 0, dist <= 500) |>
          arrange(dist)
      } else {
        candidates <- genes_df |>
          mutate(dist = (peak - end) %% genome_len) |>
          filter(dist > 0, dist <= 500) |>
          arrange(dist)
      }

      if (nrow(candidates) > 0) {
        nearest_gene <- candidates[1, ]
        current_matches <- data.frame(
          reference_identifier = nearest_gene$reference_identifier,
          mapping_type = "Intergenic",
          distance = nearest_gene$dist
        )
        matched_genes <- bind_rows(matched_genes, current_matches)
      }
    }

    if (nrow(matched_genes) > 0) {
      expanded_rows <- row[rep(1, nrow(matched_genes)), ]
      expanded_rows$reference_identifier <- matched_genes$reference_identifier
      expanded_rows$mapping_type <- matched_genes$mapping_type
      expanded_rows$distance_to_gene <- matched_genes$distance
      results_list[[i]] <- expanded_rows
    } else {
      row$reference_identifier <- NA_character_
      row$mapping_type <- "Unmapped"
      row$distance_to_gene <- NA_real_
      results_list[[i]] <- row
    }
  }

  bind_rows(results_list)
}

dir.create(mapped_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

gff_file <- file.path(gff_root, "U00096.2.gff3")
if (!file.exists(gff_file)) {
  stop("Missing promoter reference GFF: ", gff_file)
}
gene_db <- parse_gff_identifiers(gff_file)

promoter_config <- tibble::tribble(
  ~medium, ~input_file, ~output_file,
  "LB", "LB_promoter_regions.csv", "LB_promoter_identifiers.csv",
  "M9", "M9_promoter_regions.csv", "M9_promoter_identifiers.csv"
)

mapped_data <- list()
qc_rows <- list()
unmapped_rows <- list()

for (i in seq_len(nrow(promoter_config))) {
  medium <- promoter_config$medium[[i]]
  input_path <- file.path(raw_root, promoter_config$input_file[[i]])
  output_path <- file.path(mapped_root, promoter_config$output_file[[i]])
  if (!file.exists(input_path)) {
    stop("Missing raw promoter file: ", input_path)
  }

  message("Mapping ", medium, " promoters directly to MG1655 b-numbers...")
  mapped <- map_promoters_to_identifiers_expanded(input_path, gene_db)
  source_rows <- n_distinct(mapped$.source_row_index)
  mapped_source_rows <- n_distinct(
    mapped$.source_row_index[!is.na(mapped$reference_identifier)]
  )

  qc_rows[[medium]] <- tibble::tibble(
    medium,
    input_promoters = source_rows,
    expanded_output_rows = nrow(mapped),
    mapped_promoters = mapped_source_rows,
    unmapped_promoters = source_rows - mapped_source_rows,
    promoter_mapping_fraction = mapped_source_rows / source_rows,
    intragenic_rows = sum(mapped$mapping_type == "Intragenic"),
    intergenic_rows = sum(mapped$mapping_type == "Intergenic"),
    unmapped_rows = sum(mapped$mapping_type == "Unmapped"),
    unique_reference_identifiers = n_distinct(
      mapped$reference_identifier,
      na.rm = TRUE
    )
  )
  unmapped_rows[[medium]] <- mapped |>
    filter(is.na(reference_identifier)) |>
    mutate(medium = medium, .before = 1)

  mapped_data[[medium]] <- mapped
  write_csv(
    mapped |> select(-.source_row_index),
    output_path,
    na = ""
  )
}

write_csv(
  bind_rows(qc_rows),
  file.path(qc_root, "promoter_identifier_mapping_qc.csv")
)
write_csv(
  bind_rows(unmapped_rows) |> select(-.source_row_index),
  file.path(qc_root, "promoter_identifier_mapping_unmapped.csv"),
  na = ""
)

message("Identifier-mapped promoter CSV files written to ", mapped_root)
