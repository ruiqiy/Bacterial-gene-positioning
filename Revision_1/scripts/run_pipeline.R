#!/usr/bin/env Rscript

# Locate the project from this script unless GENE_POSITION_ROOT is supplied.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) {
  normalizePath(sub("^--file=", "", script_arg), winslash = "/")
} else {
  normalizePath("Revision_1/scripts/run_pipeline.R", winslash = "/")
}
default_project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/"
)
project_root <- normalizePath(
  Sys.getenv("GENE_POSITION_ROOT", unset = default_project_root),
  winslash = "/"
)
Sys.setenv(GENE_POSITION_ROOT = project_root)
setwd(project_root)

revision_root <- file.path(project_root, "Revision_1")
script_root <- file.path(revision_root, "scripts")
analysis_root <- file.path(revision_root, "gene_set_analysis")
rscript_bin <- file.path(R.home("bin"), "Rscript.exe")

run_stage <- function(name, script_name) {
  cat(sprintf("\n[%s] %s\n", name, script_name))
  script_path <- file.path(script_root, script_name)
  status <- system2(rscript_bin, args = shQuote(script_path))
  if (status != 0) {
    stop(sprintf("Stage failed: %s", name))
  }
}

# 1. Extract raw proteomics workbooks to CSVs
run_stage("1/8 Extract raw proteomics", "extract_raw_proteomics.R")

# 2. Map raw rows to stable reference identifiers
run_stage("2/8 Map raw rows to stable identifiers", "build_cross_strain_mapping.R")

# 3. Prepare molecular-weight tables
run_stage("3/8 Prepare molecular-weight tables", "prepare_molecular_weight_tables.R")

# 4. Normalize proteomics datasets
run_stage("4a/8 Normalize non-Zhu proteomics", "proteomics_normalization.R")
run_stage("4b/8 Normalize Zhu proteomics", "zhu_normalization.R")

# 5. Attach chromosome position coordinates
run_stage("5/8 Attach reference positions", "attach_reference_positions.R")

# 6. Build gene set annotations and analysis-ready proteomics
run_stage("6a/8 Build transcription/translation membership", "build_transcription_translation_gene_set.R")
run_stage("6b/8 Build analysis-ready proteomics", "build_raw_rebuild_joined_proteomics.R")

# 7. Map promoter region positions
run_stage("7/8 Map raw promoter regions", "promoter.R")

# 8. Run gene set analysis across 5 universes
gene_sets <- c("all", "core", "accessory", "transcription_translation", "other")
analysis_script <- file.path(script_root, "final_gene_position_gene_sets.R")

for (gene_set in gene_sets) {
  cat(sprintf("\n[8/8 Analyze] %s\n", gene_set))
  out_dir <- file.path(analysis_root, gene_set)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  Sys.setenv(GENE_SET_FILTER = gene_set)
  Sys.setenv(GENE_SET_OUTPUT_ROOT = out_dir)

  status <- system2(rscript_bin, args = shQuote(analysis_script))
  if (status != 0) {
    stop(sprintf("Gene-set analysis failed for: %s", gene_set))
  }
}

Sys.unsetenv("GENE_SET_FILTER")
Sys.unsetenv("GENE_SET_OUTPUT_ROOT")

message("\nActive raw-workbook/reference-identifier pipeline completed.")
