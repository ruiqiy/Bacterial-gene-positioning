#!/usr/bin/env Rscript

# Download UniProt molecular-weight raw TSVs for Zhu datasets (including UniParc recovery)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) {
  normalizePath(sub("^--file=", "", script_arg), winslash = "/")
} else {
  normalizePath("Revision_1/scripts/download_molecular_weight_api.R", winslash = "/")
}
revision_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
extracted_root <- file.path(revision_root, "data", "raw_extracted")
raw_mw_root <- file.path(revision_root, "data", "raw_data", "molecular_weight")
dir.create(raw_mw_root, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(revision_root, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(jsonlite)
  library(dplyr)
})

datasets <- list(
  zhu_E_coli = list(csv = "Ecoli_batch_F.csv", from = "UniProtKB_AC-ID", fields = "accession,id,length,mass", batch_size = 1000L),
  zhu_B_subtilis = list(csv = "Bsub_batch_F.csv", from = "UniProtKB_AC-ID", fields = "accession,id,length,mass", batch_size = 1000L),
  zhu_V_natriegens = list(csv = "Vnat_batch_C.csv", from = "RefSeq_Protein", fields = "accession,id,length,mass,organism_name,organism_id,xref_proteomes,xref_refseq", batch_size = 5000L)
)

base_curl <- c("-s", "--fail", "-L", "--connect-timeout", "30", "--max-time", "300", "--retry", "3", "--retry-delay", "2")

download_batch <- function(chunk_ids, cfg) {
  ids_file <- tempfile(fileext = ".txt")
  res_file <- tempfile(fileext = ".json")
  writeLines(chunk_ids, ids_file)
  
  system2("curl.exe", c(base_curl, "-F", paste0("from=", cfg$from), "-F", "to=UniProtKB", "-F", paste0("ids=<", ids_file), "https://rest.uniprot.org/idmapping/run"), stdout = res_file)
  job_id <- tryCatch(fromJSON(res_file)$jobId, error = function(e) NULL)
  unlink(c(ids_file, res_file))
  if (is.null(job_id)) return(NULL)
  
  status_file <- tempfile(fileext = ".txt")
  finished <- FALSE
  for (attempt in 1:60) {
    Sys.sleep(2)
    system2("curl.exe", c("-s", paste0("https://rest.uniprot.org/idmapping/status/", job_id)), stdout = status_file)
    st_raw <- paste(readLines(status_file, warn = FALSE), collapse = "")
    if (grepl("FINISHED", st_raw, ignore.case = TRUE) || grepl("results", st_raw, ignore.case = TRUE)) {
      finished <- TRUE
      break
    }
  }
  unlink(status_file)
  if (!finished) return(NULL)
  
  stream_file <- tempfile(fileext = ".tsv")
  stream_url <- paste0("https://rest.uniprot.org/idmapping/uniprotkb/results/stream/", job_id, "?compressed=false&format=tsv&fields=", cfg$fields)
  system2("curl.exe", c(base_curl, stream_url), stdout = stream_file)
  df_chunk <- tryCatch(read_tsv(stream_file, show_col_types = FALSE), error = function(e) NULL)
  unlink(stream_file)
  return(df_chunk)
}

recover_uniparc <- function(unmapped_ids) {
  if (length(unmapped_ids) == 0L) return(NULL)
  cat("  Recovering", length(unmapped_ids), "unmapped IDs via UniParc search API...\n")
  rows <- list()
  for (i in seq_along(unmapped_ids)) {
    id <- unmapped_ids[[i]]
    tmp_json <- tempfile(fileext = ".json")
    url <- paste0("https://rest.uniprot.org/uniparc/search?query=", URLencode(id, reserved = TRUE), "&format=json&size=1")
    system2("curl.exe", c("-s", "-L", url), stdout = tmp_json)
    res <- tryCatch(fromJSON(tmp_json, simplifyVector = FALSE)$results[[1]], error = function(e) NULL)
    unlink(tmp_json)
    if (!is.null(res)) {
      rows[[length(rows) + 1L]] <- tibble(
        From = as.character(id),
        Entry = as.character(res$uniParcId),
        `Entry Name` = as.character(res$uniParcId),
        Length = as.numeric(res$sequence$length),
        Mass = as.numeric(res$sequence$molWeight),
        Organism = "Vibrio natriegens",
        `Organism (ID)` = 1219067L,
        Proteomes = NA_character_,
        RefSeq = as.character(id)
      )
    }
  }
  bind_rows(rows)
}

for (name in names(datasets)) {
  cfg <- datasets[[name]]
  source_df <- read_csv(file.path(extracted_root, cfg$csv), show_col_types = FALSE)
  ids <- unique(trimws(as.character(source_df[["protein ID"]])))
  chunks <- split(ids, ceiling(seq_along(ids) / cfg$batch_size))
  all_chunks <- list()
  
  cat("Downloading", name, "(", length(ids), "IDs)...\n")
  for (i in seq_along(chunks)) {
    chunk_ids <- chunks[[i]]
    df_chunk <- NULL
    for (retry in 1:3) {
      df_chunk <- download_batch(chunk_ids, cfg)
      if (!is.null(df_chunk)) break
      Sys.sleep(3)
    }
    if (!is.null(df_chunk) && nrow(df_chunk) > 0) all_chunks[[i]] <- df_chunk
  }
  
  combined <- bind_rows(all_chunks) |> distinct()
  
  mapped_froms <- if (nrow(combined) > 0) unique(combined$From) else character(0)
  unmapped <- setdiff(ids, mapped_froms)
  if (length(unmapped) > 0) {
    uniparc_df <- recover_uniparc(unmapped)
    if (!is.null(uniparc_df) && nrow(uniparc_df) > 0) {
      combined <- bind_rows(combined, uniparc_df) |> distinct()
    }
  }
  
  out_path <- file.path(
    raw_mw_root,
    paste0(name, "_uniprot_molecular_weight_raw.tsv")
  )
  write_tsv(combined, out_path)
  cat("Saved", name, "->", nrow(combined), "rows,", length(unique(combined$From)), "unique source IDs mapped.\n\n")
}
