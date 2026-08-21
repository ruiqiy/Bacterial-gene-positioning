normalize_gene_name <- function(x) {
  tolower(trimws(as.character(x)))
}

normalize_reference_identifier <- function(species, x) {
  x <- trimws(as.character(x))
  x <- sub("\\..*", "", x)
  species <- rep_len(as.character(species), length(x))
  out <- rep(NA_character_, length(x))

  ecoli <- species == "Escherichia coli"
  bsub <- species == "Bacillus subtilis"
  vnat <- species == "Vibrio natriegens"

  out[ecoli] <- tolower(x[ecoli])
  out[bsub] <- sub("^BSU_", "BSU", toupper(x[bsub]))
  out[vnat] <- toupper(x[vnat])
  out[out == "" | is.na(x)] <- NA_character_
  out
}

is_well_formed_reference_identifier <- function(species, x) {
  x <- normalize_reference_identifier(species, x)
  species <- rep_len(as.character(species), length(x))
  ifelse(
    species == "Escherichia coli",
    grepl("^b[0-9]{4}$", x),
    ifelse(
      species == "Bacillus subtilis",
      grepl("^BSU[0-9]{5}$", x),
      grepl("^PN96_RS[0-9]{5}$", x)
    )
  )
}

normalize_logical <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

gff_attribute <- function(attributes, key) {
  pattern <- paste0("(^|;)", key, "=([^;]*)")
  value <- stringr::str_match(attributes, pattern)[, 3]
  vapply(
    value,
    function(z) if (is.na(z)) NA_character_ else utils::URLdecode(z),
    character(1)
  )
}

read_gff_features <- function(path) {
  gff <- utils::read.delim(
    path,
    comment.char = "#",
    header = FALSE,
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE
  )
  names(gff) <- c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "phase", "attributes"
  )

  dplyr::transmute(
    gff,
    seqname,
    source,
    feature,
    start = as.numeric(start),
    end = as.numeric(end),
    strand,
    attributes,
    gene = gff_attribute(attributes, "gene"),
    name = gff_attribute(attributes, "Name"),
    locus_tag = gff_attribute(attributes, "locus_tag"),
    old_locus_tag = gff_attribute(attributes, "old_locus_tag"),
    gene_biotype = gff_attribute(attributes, "gene_biotype"),
    gene_synonym = gff_attribute(attributes, "gene_synonym"),
    product = gff_attribute(attributes, "product"),
    protein_id = gff_attribute(attributes, "protein_id"),
    dbxref = gff_attribute(attributes, "Dbxref")
  )
}

build_reference_position_map <- function(
    path,
    species,
    reference_strain,
    reference_accession,
    main_seqname,
    chromosome_length,
    oric_start,
    oric_end) {
  features <- read_gff_features(path) |>
    dplyr::filter(
      seqname == main_seqname,
      feature %in% c("gene", "pseudogene"),
      !is.na(locus_tag)
    ) |>
    dplyr::mutate(
      gene_midpoint = (start + end) / 2,
      oric_midpoint = (oric_start + oric_end) / 2,
      distance_raw = abs(gene_midpoint - oric_midpoint),
      distance_from_oric = dplyr::if_else(
        distance_raw > chromosome_length / 2,
        chromosome_length - distance_raw,
        distance_raw
      ),
      norm_pos = distance_from_oric / (chromosome_length / 2),
      type = dplyr::case_when(
        norm_pos <= 1 / 3 ~ "oriC",
        norm_pos < 2 / 3 ~ "mid",
        TRUE ~ "ter"
      )
    )

  current_tags <- features |>
    dplyr::transmute(
      species,
      reference_strain,
      reference_genome_accession = reference_accession,
      reference_identifier = normalize_reference_identifier(species, locus_tag),
      reference_identifier_source = "locus_tag",
      seqname,
      start,
      end,
      reference_feature_type = feature,
      reference_gene_name = gene,
      norm_pos,
      type
    )

  old_tags <- features |>
    dplyr::filter(!is.na(old_locus_tag), old_locus_tag != "") |>
    dplyr::transmute(
      species,
      reference_strain,
      reference_genome_accession = reference_accession,
      reference_identifier = normalize_reference_identifier(
        species,
        old_locus_tag
      ),
      reference_identifier_source = "old_locus_tag",
      seqname,
      start,
      end,
      reference_feature_type = feature,
      reference_gene_name = gene,
      norm_pos,
      type
    )

  positions <- dplyr::bind_rows(current_tags, old_tags) |>
    dplyr::filter(!is.na(reference_identifier), reference_identifier != "") |>
    dplyr::distinct()

  conflicting <- positions |>
    dplyr::group_by(species, reference_identifier) |>
    dplyr::summarise(
      coordinate_count = dplyr::n_distinct(paste(seqname, start, end)),
      .groups = "drop"
    ) |>
    dplyr::filter(coordinate_count > 1L)
  if (nrow(conflicting) > 0L) {
    stop(
      species,
      ": reference identifiers map to multiple chromosome coordinates: ",
      paste(utils::head(conflicting$reference_identifier, 20), collapse = ", ")
    )
  }

  positions |>
    dplyr::arrange(
      species,
      reference_identifier,
      dplyr::desc(reference_identifier_source == "locus_tag")
    ) |>
    dplyr::group_by(species, reference_identifier) |>
    dplyr::slice(1L) |>
    dplyr::ungroup()
}

build_gff_alias_map <- function(gff_gene_rows) {
  canonical <- gff_gene_rows |>
    dplyr::filter(feature == "gene", !is.na(locus_tag), !is.na(gene)) |>
    dplyr::transmute(
      alias = normalize_gene_name(gene),
      reference_identifier = locus_tag,
      alias_type = "canonical_gene"
    )

  synonyms <- gff_gene_rows |>
    dplyr::filter(feature == "gene", !is.na(locus_tag), !is.na(gene_synonym)) |>
    dplyr::transmute(
      reference_identifier = locus_tag,
      alias = stringr::str_split(gene_synonym, ",")
    ) |>
    tidyr::unnest_longer(alias) |>
    dplyr::transmute(
      alias = normalize_gene_name(alias),
      reference_identifier,
      alias_type = "gene_synonym"
    )

  dplyr::bind_rows(canonical, synonyms) |>
    dplyr::filter(!is.na(alias), alias != "") |>
    dplyr::distinct()
}

summarize_alias_mapping <- function(query_names, alias_map) {
  query <- tibble::tibble(
    source_gene_name = as.character(query_names),
    source_gene_name_normalized = normalize_gene_name(query_names)
  ) |>
    dplyr::distinct()

  candidates <- query |>
    dplyr::left_join(
      alias_map,
      by = c("source_gene_name_normalized" = "alias"),
      relationship = "many-to-many"
    ) |>
    dplyr::group_by(source_gene_name, source_gene_name_normalized) |>
    dplyr::filter(
      if (any(alias_type == "canonical_gene", na.rm = TRUE)) {
        alias_type == "canonical_gene"
      } else {
        TRUE
      }
    ) |>
    dplyr::ungroup()

  candidates |>
    dplyr::group_by(source_gene_name, source_gene_name_normalized) |>
    dplyr::summarise(
      candidate_reference_identifiers = paste(
        sort(unique(stats::na.omit(reference_identifier))),
        collapse = ";"
      ),
      ambiguity_n = dplyr::n_distinct(reference_identifier, na.rm = TRUE),
      reference_identifier = if (
        ambiguity_n == 1
      ) {
        stats::na.omit(reference_identifier)[1]
      } else {
        NA_character_
      },
      mapping_evidence = paste(sort(unique(stats::na.omit(alias_type))), collapse = ";"),
      mapping_status = dplyr::case_when(
        ambiguity_n == 0 ~ "unmatched",
        ambiguity_n == 1 ~ "mapped_one_to_one",
        TRUE ~ "ambiguous_one_to_many"
      ),
      .groups = "drop"
    )
}

collapse_source_mapping <- function(
    data,
    dataset,
    species,
    source_strain,
    source_identifier_type,
    reference_strain,
    reference_accession,
    mapping_method,
    mapping_confidence,
    evidence_source) {
  data |>
    dplyr::mutate(
      dataset = dataset,
      species = species,
      source_strain = source_strain,
      source_identifier_type = source_identifier_type,
      source_identifier_raw = as.character(source_identifier_raw),
      source_identifier_normalized = as.character(source_identifier_normalized),
      source_gene_name = as.character(source_gene_name),
      source_gene_name_normalized = normalize_gene_name(source_gene_name),
      reference_strain = reference_strain,
      reference_genome_accession = reference_accession,
      reference_identifier = normalize_reference_identifier(
        species,
        reference_identifier
      ),
      mapping_method = mapping_method,
      mapping_confidence = mapping_confidence,
      evidence_source = evidence_source
    ) |>
    dplyr::group_by(
      dataset,
      species,
      source_strain,
      source_identifier_type,
      source_identifier_raw,
      source_identifier_normalized,
      source_gene_name,
      source_gene_name_normalized,
      reference_strain,
      reference_genome_accession,
      mapping_method,
      mapping_confidence,
      evidence_source
    ) |>
    dplyr::summarise(
      candidate_reference_identifiers = paste(
        sort(unique(stats::na.omit(reference_identifier))),
        collapse = ";"
      ),
      ambiguity_n = dplyr::n_distinct(reference_identifier, na.rm = TRUE),
      reference_identifier = if (
        ambiguity_n == 1
      ) {
        stats::na.omit(reference_identifier)[1]
      } else {
        NA_character_
      },
      mapping_status = dplyr::case_when(
        ambiguity_n == 0 ~ "unmatched",
        ambiguity_n == 1 ~ "mapped_one_to_one",
        TRUE ~ "ambiguous_one_to_many"
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      identifier_well_formed = is_well_formed_reference_identifier(
        species,
        reference_identifier
      ),
      mapping_status = dplyr::if_else(
        mapping_status == "mapped_one_to_one" & !identifier_well_formed,
        "malformed_identifier",
        mapping_status
      ),
      mapping_cardinality = dplyr::case_when(
        mapping_status == "mapped_one_to_one" ~ "one-to-one",
        mapping_status == "ambiguous_one_to_many" ~ "one-to-many",
        TRUE ~ "unresolved"
      )
    )
}

left_join_preserve_rows <- function(x, y, by, relationship = "many-to-one") {
  if (anyDuplicated(y[unname(by)])) {
    stop("Mapping table has duplicated join keys: ", paste(unname(by), collapse = ", "))
  }

  before <- nrow(x)
  x$.revision_row_id <- seq_len(before)
  out <- dplyr::left_join(x, y, by = by, relationship = relationship)

  if (nrow(out) != before || !identical(out$.revision_row_id, seq_len(before))) {
    stop("Join changed row count or row order")
  }

  out$.revision_row_id <- NULL
  out
}
