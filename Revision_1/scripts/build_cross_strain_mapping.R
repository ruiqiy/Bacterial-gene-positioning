curr_dir <- normalizePath(getwd(), winslash = "/")
revision_root <- if (basename(curr_dir) == "Revision_1") curr_dir else file.path(curr_dir, "Revision_1")

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
})

source(file.path(revision_root, "scripts", "mapping_helpers.R"))

extracted_root <- file.path(revision_root, "data", "raw_extracted")
mapped_raw_root <- file.path(revision_root, "data", "raw_mapped")
gff_root <- file.path(revision_root, "data", "raw_data", "reference_gff")
out_root <- file.path(revision_root, "data")
qc_root <- file.path(revision_root, "qc")

ecoli_gff <- read_gff_features(file.path(gff_root, "GCF_000005845.2.gff"))
ecoli_aliases <- build_gff_alias_map(ecoli_gff)
bsub_gff <- read_gff_features(file.path(gff_root, "GCF_000009045.1.gff"))
vnat_gff <- read_gff_features(file.path(gff_root, "GCF_001456255.1.gff"))

reference_alias_lookup <- bind_rows(
  mutate(ecoli_aliases, species = "Escherichia coli"),
  mutate(build_gff_alias_map(bsub_gff), species = "Bacillus subtilis"),
  mutate(build_gff_alias_map(vnat_gff), species = "Vibrio natriegens")
) |>
  mutate(
    reference_identifier = normalize_reference_identifier(
      species,
      reference_identifier
    )
  ) |>
  group_by(species, source_gene_name_normalized = alias) |>
  summarise(
    fallback_candidate_reference_identifiers = paste(
      sort(unique(na.omit(reference_identifier))),
      collapse = ";"
    ),
    fallback_ambiguity_n = n_distinct(reference_identifier, na.rm = TRUE),
    fallback_alias_type = paste(sort(unique(alias_type)), collapse = ";"),
    .groups = "drop"
  ) |>
  mutate(
    fallback_reference_identifier = if_else(
      fallback_ambiguity_n == 1L,
      fallback_candidate_reference_identifiers,
      NA_character_
    )
  )

apply_common_name_fallback <- function(data) {
  if (!"identifier_well_formed" %in% names(data)) {
    data$identifier_well_formed <- NA
  }

  data |>
    mutate(
      source_gene_name_normalized = normalize_gene_name(source_gene_name)
    ) |>
    left_join(
      reference_alias_lookup,
      by = c("species", "source_gene_name_normalized"),
      relationship = "many-to-one"
    ) |>
    mutate(
      .use_common_name_fallback =
        coalesce(mapping_status, "unmatched") != "mapped_one_to_one" &
        coalesce(fallback_ambiguity_n, 0L) == 1L,
      reference_identifier = if_else(
        .use_common_name_fallback,
        fallback_reference_identifier,
        reference_identifier
      ),
      candidate_reference_identifiers = if_else(
        .use_common_name_fallback,
        fallback_candidate_reference_identifiers,
        candidate_reference_identifiers
      ),
      ambiguity_n = if_else(
        .use_common_name_fallback,
        1L,
        coalesce(as.integer(ambiguity_n), 0L)
      ),
      mapping_method = if_else(
        .use_common_name_fallback,
        paste0("unique_RefSeq_common_name_fallback_", fallback_alias_type),
        mapping_method
      ),
      mapping_confidence = if_else(
        .use_common_name_fallback,
        "high",
        mapping_confidence
      ),
      mapping_status = if_else(
        .use_common_name_fallback,
        "mapped_one_to_one",
        mapping_status
      ),
      mapping_cardinality = if_else(
        .use_common_name_fallback,
        "one-to-one",
        mapping_cardinality
      ),
      identifier_well_formed = if_else(
        .use_common_name_fallback,
        is_well_formed_reference_identifier(species, reference_identifier),
        identifier_well_formed
      )
    ) |>
    select(-starts_with("fallback_"), -.use_common_name_fallback)
}

direct_mapping <- function(
    data,
    dataset,
    species,
    source_strain,
    identifier_type,
    reference_strain,
    reference_accession,
    evidence_source) {
  parsed <- data |>
    transmute(
      source_gene_name,
      source_identifier_raw = as.character(source_identifier_raw),
      source_identifier_normalized = normalize_reference_identifier(
        species,
        source_identifier_raw
      ),
      reference_identifier = if_else(
        is_well_formed_reference_identifier(species, source_identifier_raw),
        normalize_reference_identifier(species, source_identifier_raw),
        NA_character_
      )
    )

  collapse_source_mapping(
    parsed,
    dataset = dataset,
    species = species,
    source_strain = source_strain,
    source_identifier_type = identifier_type,
    reference_strain = reference_strain,
    reference_accession = reference_accession,
    mapping_method = "author_provided_reference_locus_identifier",
    mapping_confidence = "high",
    evidence_source = evidence_source
  )
}

# Schmidt: use the raw UniProt accession, not the common gene name, as the
# cross-strain key. The MG1655 RefSeq GFF provides the accession-to-b-number
# crosswalk in CDS Dbxref attributes.
schmidt_raw <- read_csv(
  file.path(extracted_root, "schmidt_data_copy.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  transmute(
    source_identifier_raw = as.character(`Uniprot Accession`),
    source_identifier_normalized = toupper(trimws(source_identifier_raw)),
    source_gene_name = as.character(`Gene`)
  ) |>
  filter(!is.na(source_identifier_raw), source_identifier_raw != "")

ecoli_uniprot <- ecoli_gff |>
  filter(feature == "CDS", !is.na(locus_tag), !is.na(dbxref)) |>
  transmute(
    source_identifier_normalized = str_match(
      dbxref,
      "(?:^|,)UniProtKB/Swiss-Prot:([^,;]+)"
    )[, 2],
    reference_identifier = locus_tag
  ) |>
  filter(!is.na(source_identifier_normalized)) |>
  distinct()

schmidt_map <- schmidt_raw |>
  left_join(
    ecoli_uniprot,
    by = "source_identifier_normalized",
    relationship = "many-to-many"
  ) |>
  collapse_source_mapping(
    dataset = "schmidt_2016",
    species = "Escherichia coli",
    source_strain = "K-12 BW25113",
    source_identifier_type = "UniProt accession",
    reference_strain = "K-12 MG1655",
    reference_accession = "GCF_000005845.2",
    mapping_method = "UniProt_accession_to_MG1655_GFF_CDS_Dbxref",
    mapping_confidence = "high",
    evidence_source = paste(
      "Revision_1/data/raw_extracted/schmidt_data_copy.csv +",
      "Revision_1/data/raw_data/reference_gff/GCF_000005845.2.gff"
    )
  )

# Peebo: b-numbers are supplied directly in the raw workbook.
peebo_raw <- read_csv(
  file.path(extracted_root, "peebo_copy.csv"),
  show_col_types = FALSE,
  name_repair = "unique"
) |>
  transmute(
    source_gene_name = as.character(Gene),
    source_identifier_raw = as.character(BID)
  ) |>
  filter(!is.na(source_identifier_raw), source_identifier_raw != "")

peebo_map <- direct_mapping(
  peebo_raw,
  dataset = "peebo_2015",
  species = "Escherichia coli",
  source_strain = "K-12 BW25113",
  identifier_type = "MG1655 b-number supplied by authors",
  reference_strain = "K-12 MG1655",
  reference_accession = "GCF_000005845.2",
  evidence_source = "Revision_1/data/raw_extracted/peebo_copy.csv"
)

# Valgepea: KEGG_ID is the MG1655 b-number.
valgepea_raw <- read_csv(
  file.path(extracted_root, "Valgepea_copy.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  transmute(
    source_identifier_raw = as.character(KEGG_ID),
    source_gene_name = str_split_fixed(as.character(gene), ",", 2)[, 1]
  ) |>
  filter(!is.na(source_identifier_raw), source_identifier_raw != "")

valgepea_map <- direct_mapping(
  valgepea_raw,
  dataset = "valgepea_2013_prot",
  species = "Escherichia coli",
  source_strain = "K-12 MG1655",
  identifier_type = "KEGG b-number",
  reference_strain = "K-12 MG1655",
  reference_accession = "GCF_000005845.2",
  evidence_source = "Revision_1/data/raw_extracted/Valgepea_copy.csv"
)

# Li: compound records do not identify gene-specific abundance and are excluded.
# Retain only single-gene names for the canonical/synonym lookup.
li_queries <- read_csv(
  file.path(extracted_root, "Li_copy.csv"),
  show_col_types = FALSE,
  name_repair = "minimal"
) |>
  filter(coalesce(!str_detect(as.character(Gene), fixed("+")), TRUE)) |>
  transmute(source_gene_name = trimws(as.character(Gene))) |>
  filter(!is.na(source_gene_name), source_gene_name != "") |>
  distinct()

li_alias_map <- summarize_alias_mapping(li_queries$source_gene_name, ecoli_aliases)

li_map <- li_alias_map |>
  transmute(
    source_gene_name,
    source_identifier_raw = source_gene_name,
    source_identifier_normalized = source_gene_name_normalized,
    reference_identifier,
    mapping_evidence,
    candidate_reference_identifiers,
    ambiguity_n,
    mapping_status
  ) |>
  mutate(
    dataset = "li_2014_all",
    species = "Escherichia coli",
    source_strain = "K-12 MG1655",
    source_identifier_type = "common gene name (unique RefSeq canonical/synonym match)",
    source_gene_name_normalized = normalize_gene_name(source_gene_name),
    reference_strain = "K-12 MG1655",
    reference_genome_accession = "GCF_000005845.2",
    reference_identifier = normalize_reference_identifier(
      species,
      reference_identifier
    ),
    mapping_method = case_when(
      mapping_status != "mapped_one_to_one" ~ "unresolved_common_name",
      str_detect(mapping_evidence, "canonical_gene") ~ "unique_RefSeq_canonical_gene_name",
      TRUE ~ "unique_RefSeq_gene_synonym"
    ),
    mapping_confidence = if_else(
      mapping_status == "mapped_one_to_one",
      "high",
      "unresolved"
    ),
    evidence_source = paste(
      "Revision_1/data/raw_extracted/Li_copy.csv +",
      "Revision_1/data/raw_data/reference_gff/GCF_000005845.2.gff"
    ),
    identifier_well_formed = is_well_formed_reference_identifier(
      species,
      reference_identifier
    ),
    mapping_cardinality = case_when(
      mapping_status == "mapped_one_to_one" ~ "one-to-one",
      mapping_status == "ambiguous_one_to_many" ~ "one-to-many",
      TRUE ~ "unresolved"
    )
  )

read_zhu_species <- function(files) {
  map_dfr(
    files,
    function(filename) {
      read_csv(
        file.path(extracted_root, filename),
        show_col_types = FALSE,
        name_repair = "minimal"
      ) |>
        transmute(
          source_gene_name = as.character(`gene name`),
          source_identifier_raw = as.character(`gene locus`)
        )
    }
  ) |>
    filter(!is.na(source_identifier_raw), source_identifier_raw != "") |>
    distinct()
}

zhu_ecoli_map <- direct_mapping(
  read_zhu_species(c("Ecoli_batch_F.csv", "Ecoli_batch_K.csv")),
  dataset = "zhu_E_coli",
  species = "Escherichia coli",
  source_strain = "K-12 NCM3722",
  identifier_type = "MG1655 b-number supplied by authors",
  reference_strain = "K-12 MG1655",
  reference_accession = "GCF_000005845.2",
  evidence_source = "Revision_1/data/raw_extracted/Ecoli_batch_F.csv;Ecoli_batch_K.csv"
)

zhu_bsub_map <- direct_mapping(
  read_zhu_species(c("Bsub_batch_F.csv", "Bsub_batch_G.csv")),
  dataset = "zhu_B_subtilis",
  species = "Bacillus subtilis",
  source_strain = "subsp. subtilis 168",
  identifier_type = "legacy BSU locus tag",
  reference_strain = "subsp. subtilis 168",
  reference_accession = "GCF_000009045.1",
  evidence_source = "Revision_1/data/raw_extracted/Bsub_batch_F.csv;Bsub_batch_G.csv"
)

zhu_vnat_map <- direct_mapping(
  read_zhu_species(c("Vnat_batch_C.csv", "Vnat_batch_L.csv")),
  dataset = "zhu_V_natriegens",
  species = "Vibrio natriegens",
  source_strain = "ATCC 14048",
  identifier_type = "PN96_RS locus tag",
  reference_strain = "ATCC 14048",
  reference_accession = "GCF_001456255.1",
  evidence_source = "Revision_1/data/raw_extracted/Vnat_batch_C.csv;Vnat_batch_L.csv"
)

goelzer_raw <- map_dfr(
  c("PYR", "S", "TS", "CH", "CHG"),
  function(condition) {
    read_csv(
      file.path(extracted_root, paste0("Goelzer_", condition, ".csv")),
      show_col_types = FALSE,
      name_repair = "minimal"
    ) |>
      transmute(
        source_gene_name = as.character(`% Gene_Name`),
        source_identifier_raw = as.character(`% BSU Number`)
      )
  }
) |>
  filter(!is.na(source_identifier_raw), source_identifier_raw != "") |>
  distinct()

goelzer_map <- direct_mapping(
  goelzer_raw,
  dataset = "Goelzer_2015",
  species = "Bacillus subtilis",
  source_strain = "subsp. subtilis 168",
  identifier_type = "legacy BSU locus tag",
  reference_strain = "subsp. subtilis 168",
  reference_accession = "GCF_000009045.1",
  evidence_source = "Revision_1/data/raw_extracted/Goelzer_{PYR,S,TS,CH,CHG}.csv"
)

mapping <- bind_rows(
  schmidt_map,
  peebo_map,
  valgepea_map,
  li_map,
  zhu_ecoli_map,
  zhu_bsub_map,
  zhu_vnat_map,
  goelzer_map
) |>
  select(
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
    reference_identifier,
    candidate_reference_identifiers,
    mapping_method,
    mapping_confidence,
    mapping_status,
    mapping_cardinality,
    ambiguity_n,
    identifier_well_formed,
    evidence_source,
    any_of("mapping_evidence")
  ) |>
  arrange(dataset, source_gene_name_normalized, source_identifier_normalized)

mapping <- apply_common_name_fallback(mapping)

collapse_reference_candidates <- function(reference_ids, candidate_text) {
  values <- c(reference_ids, candidate_text)
  values <- unlist(strsplit(values[!is.na(values)], ";", fixed = TRUE))
  values <- trimws(values)
  values <- values[values != "" & values != "NA"]
  paste(sort(unique(values)), collapse = ";")
}

unresolved <- mapping |>
  filter(mapping_status != "mapped_one_to_one")

dir.create(mapped_raw_root, recursive = TRUE, showWarnings = FALSE)

identifier_lookup <- mapping |>
  filter(dataset != "li_2014_all") |>
  filter(!is.na(source_identifier_normalized)) |>
  group_by(dataset, source_identifier_normalized) |>
  summarise(
    candidate_reference_identifiers = collapse_reference_candidates(
      reference_identifier,
      candidate_reference_identifiers
    ),
    ambiguity_n = if_else(
      candidate_reference_identifiers == "",
      0L,
      lengths(strsplit(candidate_reference_identifiers, ";", fixed = TRUE))
    ),
    had_ambiguous_evidence = any(mapping_status == "ambiguous_one_to_many"),
    reference_identifier = if (
      ambiguity_n == 1L && !had_ambiguous_evidence
    ) candidate_reference_identifiers else NA_character_,
    mapping_method = paste(sort(unique(mapping_method)), collapse = ";"),
    mapping_confidence = if_else(
      ambiguity_n == 1L && !had_ambiguous_evidence,
      "high",
      "unresolved"
    ),
    mapping_status = case_when(
      had_ambiguous_evidence | ambiguity_n > 1L ~ "ambiguous_one_to_many",
      ambiguity_n == 1L ~ "mapped_one_to_one",
      TRUE ~ "unmatched"
    ),
    mapping_cardinality = case_when(
      mapping_status == "mapped_one_to_one" ~ "one-to-one",
      mapping_status == "ambiguous_one_to_many" ~ "one-to-many",
      TRUE ~ "zero"
    ),
    .groups = "drop"
  )

li_lookup <- li_map |>
  select(
    source_gene_name_normalized,
    reference_identifier,
    candidate_reference_identifiers,
    ambiguity_n,
    mapping_method,
    mapping_confidence,
    mapping_status,
    mapping_cardinality
  )

if (anyDuplicated(li_lookup$source_gene_name_normalized)) {
  stop("Li lookup contains duplicated normalized source gene names")
}

raw_file_config <- tribble(
  ~filename, ~dataset, ~species, ~gene_column, ~identifier_column, ~join_mode,
  "schmidt_data_copy.csv", "schmidt_2016", "Escherichia coli",
  "Gene", "Uniprot Accession", "identifier_upper",
  "peebo_copy.csv", "peebo_2015", "Escherichia coli",
  "Gene", "BID", "identifier",
  "Valgepea_copy.csv", "valgepea_2013_prot", "Escherichia coli",
  "gene", "KEGG_ID", "identifier",
  "Li_copy.csv", "li_2014_all", "Escherichia coli",
  "Gene", NA_character_, "gene_name",
  "Goelzer_PYR.csv", "Goelzer_2015", "Bacillus subtilis",
  "% Gene_Name", "% BSU Number", "identifier",
  "Goelzer_S.csv", "Goelzer_2015", "Bacillus subtilis",
  "% Gene_Name", "% BSU Number", "identifier",
  "Goelzer_TS.csv", "Goelzer_2015", "Bacillus subtilis",
  "% Gene_Name", "% BSU Number", "identifier",
  "Goelzer_CH.csv", "Goelzer_2015", "Bacillus subtilis",
  "% Gene_Name", "% BSU Number", "identifier",
  "Goelzer_CHG.csv", "Goelzer_2015", "Bacillus subtilis",
  "% Gene_Name", "% BSU Number", "identifier",
  "Ecoli_batch_F.csv", "zhu_E_coli", "Escherichia coli",
  "gene name", "gene locus", "identifier",
  "Ecoli_batch_K.csv", "zhu_E_coli", "Escherichia coli",
  "gene name", "gene locus", "identifier",
  "Bsub_batch_F.csv", "zhu_B_subtilis", "Bacillus subtilis",
  "gene name", "gene locus", "identifier",
  "Bsub_batch_G.csv", "zhu_B_subtilis", "Bacillus subtilis",
  "gene name", "gene locus", "identifier",
  "Vnat_batch_C.csv", "zhu_V_natriegens", "Vibrio natriegens",
  "gene name", "gene locus", "identifier",
  "Vnat_batch_L.csv", "zhu_V_natriegens", "Vibrio natriegens",
  "gene name", "gene locus", "identifier"
)

map_extracted_file <- function(cfg) {
  raw <- read_csv(
    file.path(extracted_root, cfg$filename),
    show_col_types = FALSE,
    name_repair = "unique"
  ) |>
    mutate(
      source_file = cfg$filename,
      source_row_number = row_number(),
      source_row_id = paste0(cfg$filename, ":", source_row_number)
    )
  input_rows <- nrow(raw)

  if (is.na(cfg$identifier_column[[1L]])) {
    raw <- raw |>
      filter(coalesce(
        !str_detect(as.character(.data[[cfg$gene_column]]), fixed("+")),
        TRUE
      )) |>
      mutate(
        source_gene_name = str_trim(as.character(.data[[cfg$gene_column]])),
        source_identifier_raw = source_gene_name,
        source_gene_name_normalized = normalize_gene_name(source_gene_name)
      ) |>
      left_join(
        li_lookup,
        by = "source_gene_name_normalized",
        relationship = "many-to-one"
      )
  } else {
    raw <- raw |>
      mutate(
        source_gene_name = as.character(.data[[cfg$gene_column]]),
        source_identifier_raw = as.character(.data[[cfg$identifier_column]]),
        source_identifier_normalized = if (cfg$join_mode == "identifier_upper") {
          str_to_upper(str_trim(source_identifier_raw))
        } else {
          normalize_reference_identifier(cfg$species, source_identifier_raw)
        }
      ) |>
      left_join(
        filter(identifier_lookup, dataset == cfg$dataset) |>
          select(-dataset),
        by = "source_identifier_normalized",
        relationship = "many-to-one"
      )
  }

  raw <- raw |>
    mutate(species = cfg$species) |>
    apply_common_name_fallback()

  mapped_all <- raw |>
    mutate(
      dataset = cfg$dataset,
      species = cfg$species,
      mapping_status = coalesce(mapping_status, "unmatched"),
      mapping_confidence = coalesce(mapping_confidence, "unresolved"),
      mapping_cardinality = coalesce(mapping_cardinality, "zero"),
      candidate_reference_identifiers = coalesce(
        candidate_reference_identifiers,
        ""
      )
    )
  
  write_csv(mapped_all, file.path(mapped_raw_root, cfg$filename), na = "")
  tibble(
    dataset = cfg$dataset,
    source_file = cfg$filename,
    extracted_rows = input_rows,
    excluded_compound_rows =
      input_rows - n_distinct(mapped_all$source_row_id),
    rows_submitted_to_mapping = nrow(mapped_all),
    output_rows = nrow(mapped_all),
    mapped_one_to_one = sum(
      mapped_all$mapping_status == "mapped_one_to_one"
    ),
    unmatched = sum(
      mapped_all$mapping_status == "unmatched"
    ),
    ambiguous_one_to_many = sum(
      mapped_all$mapping_status == "ambiguous_one_to_many"
    ),
    non_one_to_one = sum(
      mapped_all$mapping_status != "mapped_one_to_one"
    )
  )
}

raw_mapping_qc <- map_dfr(
  seq_len(nrow(raw_file_config)),
  ~ map_extracted_file(raw_file_config[.x, ])
)

mapping_qc <- raw_mapping_qc |>
  group_by(dataset) |>
  summarise(
    extracted_row_count = sum(extracted_rows),
    excluded_compound_rows = sum(excluded_compound_rows),
    source_row_count = sum(rows_submitted_to_mapping),
    output_row_count = sum(output_rows),
    mapped_one_to_one = sum(mapped_one_to_one),
    unmatched = sum(unmatched),
    ambiguous_one_to_many = sum(ambiguous_one_to_many),
    non_one_to_one = sum(non_one_to_one),
    one_to_one_fraction = mapped_one_to_one / source_row_count,
    large_unresolved_fraction = one_to_one_fraction < 0.90,
    .groups = "drop"
  )

write_csv(mapping, file.path(out_root, "cross_strain_gene_mapping.csv"))
write_csv(mapping_qc, file.path(qc_root, "cross_strain_mapping_qc.csv"))
write_csv(unresolved, file.path(qc_root, "cross_strain_mapping_unresolved.csv"))
write_csv(raw_mapping_qc, file.path(qc_root, "raw_row_identifier_mapping_qc.csv"))

message(
  "Wrote cross-strain mapping and row-level mapped raw tables for ",
  n_distinct(mapping$dataset),
  " datasets."
)
