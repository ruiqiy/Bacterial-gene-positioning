setwd("/home/usr/Desktop/research/gene_position")
rm(list = ls())

## ---- setup ----
# Load packages
pkgs <- c("readr", "dplyr", "tidyr", "stringr", "tibble", "rtracklayer", "purrr")
lapply(pkgs, library, character.only = TRUE)

## ---- user-configurable settings ----
dataset_tag <- "schmidt_2016"

meta_file <- "schmidt_meta.csv"    # condition + growth rate
data_file <- "schmidt_data_copy.csv"    # gene + per-condition protein copy number / cell

# GFF with gene annotations (MG1655, BW25113, etc.)
gff_path  <- "GCF_000750555.1.gff"
out_file  <- paste0(dataset_tag, "_pseudo.rds")

## ---- helper: normalize condition labels (case + space) ----
normalize_label <- function(x) {
  x %>%
    stringr::str_replace_all("\\s+", " ") %>%  # collapse multiple spaces
    stringr::str_trim() %>%                   # trim leading/trailing spaces
    stringr::str_to_lower()                   # lowercase
}

## ---- GFF helper: build gene -> coords map (no locus_tag) ----
build_gene_coords <- function(gff_path) {
  if (!file.exists(gff_path)) {
    stop("GFF not found: ", gff_path)
  }
  
  gr <- rtracklayer::import(gff_path, format = "gff3")
  
  # keep only gene features
  genes_gr <- gr[mcols(gr)$type == "gene"]
  
  # gene symbol: use 'gene' column, fall back to 'Name'
  gene_sym <- dplyr::coalesce(
    as.character(mcols(genes_gr)$gene),
    as.character(mcols(genes_gr)$Name)
  )
  
  coords_raw <- tibble(
    gene    = gene_sym,
    seqname = as.character(seqnames(genes_gr)),
    start   = as.integer(start(genes_gr)),
    end     = as.integer(end(genes_gr))
  ) %>%
    filter(!is.na(gene))  # drop unnamed features
  
  # Check for genes with multiple non-NA coordinate entries
  dup_genes <- coords_raw %>%
    filter(!is.na(start), !is.na(end)) %>%
    count(gene) %>%
    filter(n > 1)
  
  if (nrow(dup_genes) > 0) {
    warning(
      "Multiple non-NA coordinate entries found in GFF for some gene names. ",
      "Keeping the first non-NA start/end record per gene. Example genes: ",
      paste(head(dup_genes$gene, 20), collapse = ", "),
      if (nrow(dup_genes) > 20) " ..."
    )
  }
  
  # For each gene, keep the *first* row with non-NA start & end
  coords_gene <- coords_raw %>%
    filter(!is.na(start), !is.na(end)) %>%
    group_by(gene) %>%
    dplyr::slice(1) %>%              # first occurrence in original GFF order
    ungroup()
  
  coords_gene
}

## ---- read metadata: condition + growth rate ----
meta_raw <- read_csv(meta_file, col_types = cols())

if (ncol(meta_raw) < 2) {
  stop("schmidt_meta.csv must have at least two columns: condition and growth rate.")
}

meta <- meta_raw %>%
  transmute(
    condition_raw = .data[[ names(meta_raw)[1] ]],
    growth_rate   = as.numeric(.data[[ names(meta_raw)[2] ]])
  ) %>%
  mutate(
    cond_norm = normalize_label(condition_raw)
  )

if (any(is.na(meta$growth_rate))) {
  warning("Some growth_rate values are NA in metadata; this will cause errors later.")
}

## ---- read proteomics: mass per cell ----
prot_raw <- read_csv(data_file, col_types = cols())

if (ncol(prot_raw) < 2) {
  stop("schmidy_data.csv must have at least two columns: gene and at least one condition.")
}

names(prot_raw)[1] <- "gene"
condition_cols <- setdiff(names(prot_raw), "gene")

## ---- align conditions using normalized labels ----
prot_cond_df <- tibble(
  condition_col = condition_cols,
  cond_norm     = normalize_label(condition_col)
)

cond_join <- prot_cond_df %>%
  left_join(meta, by = "cond_norm")

unmatched <- cond_join %>%
  filter(is.na(growth_rate))

if (nrow(unmatched) > 0) {
  stop(
    "These conditions in proteomics data could not be matched (after case/space normalization) ",
    "to schmidt_meta.csv:\n",
    paste(unmatched$condition_col, collapse = ", ")
  )
}

meta_use <- cond_join
# meta_use has: condition_col, cond_norm, condition_raw, growth_rate

## ---- convert copy per cell -> copy fraction + pseudo ----
copy_mat <- as.matrix(prot_raw[condition_cols])
storage.mode(copy_mat) <- "numeric"

col_totals <- colSums(copy_mat, na.rm = TRUE)
col_totals[col_totals == 0] <- NA_real_

copy_frac <- sweep(copy_mat, 2, col_totals, "/")

eps <- 1e-8
copy_frac_pseudo <- copy_frac + eps

n_cond <- rowSums(!is.na(copy_frac))

## ---- construct dataset__mu_ column names ----
mu_num   <- round(meta_use$growth_rate, 6)
mu_label <- sprintf("%.6f", mu_num)      # exactly 6 decimal places
mu_label <- gsub("\\.", "_", mu_label)   # e.g. 0.470000 -> 0_470000

new_names <- paste0(dataset_tag, "__mu_", mu_label)
colnames(copy_frac_pseudo) <- new_names

## ---- assemble compiled table (no coords yet) ----
compiled <- tibble(
  gene   = prot_raw$gene,
  n_cond = n_cond
) %>%
  bind_cols(as.data.frame(copy_frac_pseudo, check.names = FALSE))

## ---- build gene coords from GFF and join by gene ----
gene_coords <- build_gene_coords(gff_path)

compiled2 <- compiled %>%
  # gene symbol -> genomic coordinates
  left_join(gene_coords, by = "gene") %>%
  # Drop genes with no coordinates
  relocate(gene, seqname, start, end, n_cond)

# Optional: sort by genomic coordinate
if ("start" %in% names(compiled2)) {
  compiled2 <- compiled2 %>% arrange(start, gene)
}

chrom_len <- 4631469	#BW25113
merged <- compiled2 %>%
  mutate(gene_mid = if_else(
    start <= end,
    (start + end)/2,
    # wrap-around case:
    ((start + (end + chrom_len))/2) %% chrom_len	
  ))

merged$OriC_start <- 3918994 #BW25113
merged$OriC_end <- 3919371

sal_df2 <-merged %>% 
  mutate(
    ori_mid = if_else(
      OriC_start <= OriC_end,
      (OriC_start + OriC_end)/2,
      # wrap-around case:
      ((OriC_start + (OriC_end + chrom_len))/2) %% chrom_len	
    )
  )

sal_df2 <- sal_df2 %>%
  
  mutate(chrom_len = chrom_len) %>%
  ungroup() %>%
  
  mutate(
    dist_raw = abs(gene_mid - ori_mid),
    dist      = if_else(dist_raw > chrom_len/2,
                        chrom_len - dist_raw,
                        dist_raw),
    # 5) normalize so 0 at oriC, 1 at ter
    norm_pos  = dist / (chrom_len/2)
  )

sal_df2$type <- ifelse(sal_df2$norm_pos <= 1/3,"oriC",
                       ifelse(sal_df2$norm_pos<2/3,"mid","ter"))

## ---- write output ----
saveRDS(sal_df2, out_file)
message("Done. Wrote compiled copy-fraction table to: ", out_file)















## ---- setup ----
pkgs <- c("readr", "dplyr", "tidyr", "stringr", "tibble", "rtracklayer")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}))
lapply(pkgs, library, character.only = TRUE)

## ---- user-configurable settings ----
dataset_tag <- "peebo_2015"

# Input Peebo file
data_file <- "peebo_copy.csv"   # path to your Peebo file

# GFF with gene annotations (for coordinates, matched by gene name)
gff_path  <- "GCF_000750555.1.gff"

# Output file
out_file  <- paste0(dataset_tag, "_pseudo.rds")

## ---- helper: normalize gene names for joining ----
normalize_gene <- function(x) {
  x %>%
    stringr::str_trim() %>%
    stringr::str_to_lower()
}

## ---- GFF helper: build gene -> coords map (match by normalized gene name) ----
build_gene_coords <- function(gff_path) {
  if (!file.exists(gff_path)) {
    stop("GFF not found: ", gff_path)
  }
  
  gr <- rtracklayer::import(gff_path, format = "gff3")
  
  # keep only gene features
  genes_gr <- gr[mcols(gr)$type == "gene"]
  
  # gene symbol: use 'gene' column, fall back to 'Name'
  gene_sym <- dplyr::coalesce(
    as.character(mcols(genes_gr)$gene),
    as.character(mcols(genes_gr)$Name)
  )
  
  coords_raw <- tibble(
    gene_raw = gene_sym,
    gene     = normalize_gene(gene_sym),
    seqname  = as.character(seqnames(genes_gr)),
    start    = as.integer(start(genes_gr)),
    end      = as.integer(end(genes_gr))
  ) %>%
    filter(!is.na(gene))  # drop unnamed / un-normalizable features
  
  # Check for genes with multiple non-NA coordinate entries (after normalization)
  dup_genes <- coords_raw %>%
    filter(!is.na(start), !is.na(end)) %>%
    count(gene) %>%
    filter(n > 1)
  
  if (nrow(dup_genes) > 0) {
    warning(
      "Multiple non-NA coordinate entries found in GFF for some normalized gene names. ",
      "Keeping the first non-NA start/end record per normalized gene. Example genes: ",
      paste(head(dup_genes$gene, 20), collapse = ", "),
      if (nrow(dup_genes) > 20) " ..."
    )
  }
  
  # For each normalized gene, keep the first row with non-NA start & end
  coords_gene <- coords_raw %>%
    filter(!is.na(start), !is.na(end)) %>%
    group_by(gene) %>%
    summarise(
      seqname = dplyr::first(seqname),
      start   = dplyr::first(start),
      end     = dplyr::first(end),
      .groups = "drop"
    )
  
  coords_gene    # columns: gene (normalized), seqname, start, end
}

## ---- read Peebo data ----
peebo_raw <- read_csv(data_file, col_types = cols())

# Expect columns: "Gene name", "B number identifier", "Molecular weight [kDa]", then conditions
peebo <- peebo_raw %>%
  dplyr::rename(
    gene      = Gene,
    locus_tag = BID
  )

# Conditions = all columns except gene, locus_tag
cond_cols <- setdiff(names(peebo), c("gene", "locus_tag"))

if (length(cond_cols) == 0) {
  stop("No condition columns detected. Check that the Peebo file has concentration columns beyond the first three.")
}

## ---- build condition metadata: growth rates & replicate labels ----
cond_info <- tibble(cond_col = cond_cols) %>%
  mutate(
    # Extract numeric growth rate from column name:
    # "0.22"   -> "0.22"
    # "0.22.1" -> "0.22"
    mu_str = sub("^([0-9]+(?:\\.[0-9]+)?).*", "\\1", cond_col),
    mu_num = as.numeric(mu_str)
  )

if (any(is.na(cond_info$mu_num))) {
  warning("Some condition column names could not be parsed into growth rates; mu_num is NA for those.")
}

cond_info <- cond_info %>%
  mutate(
    mu_num_round = round(mu_num, 6),
    mu_label     = sprintf("%.6f", mu_num_round),   # 6 decimal places
    mu_label     = gsub("\\.", "_", mu_label)       # 0.470000 -> 0_470000
  ) %>%
  group_by(mu_label) %>%
  mutate(
    rep_index    = dplyr::row_number(),
    n_rep        = dplyr::n(),
    mu_label_full = if_else(
      n_rep > 1L,
      paste0(mu_label, "_rep", rep_index),          # e.g. 0_220000_rep1, 0_220000_rep2
      mu_label
    )
  ) %>%
  ungroup()

dup_mu <- cond_info %>% count(mu_label) %>% filter(n > 1)
if (nrow(dup_mu) > 0) {
  warning(
    "Some growth rates appear in multiple condition columns (replicates). ",
    "They are kept as separate conditions with suffix _rep#. Growth rates: ",
    paste(dup_mu$mu_label, collapse = ", ")
  )
}

new_names <- paste0(dataset_tag, "__mu_", cond_info$mu_label_full)

## ---- convert concentration -> copy fraction + pseudo ----
# Concentration matrix (molecules / fL)
conc_mat <- as.matrix(peebo[cond_cols])
storage.mode(conc_mat) <- "numeric"

# Fill 0 for empty / missing protein expression values
conc_mat[!is.finite(conc_mat)] <- 0  # NA / NaN / Inf -> 0

# Total concentration per condition (column)
# Note: Because volume is constant per condition, concentration fraction = copy fraction
col_totals <- colSums(conc_mat, na.rm = TRUE)
col_totals[col_totals == 0] <- NA_real_

# Copy fraction per gene per condition
copy_frac <- sweep(conc_mat, 2, col_totals, "/")

# Pseudo copy fraction (added after filling zeros)
eps <- 1e-8
copy_frac_pseudo <- copy_frac + eps

# Number of conditions with non-NA copy fraction per gene
n_cond <- rowSums(!is.na(copy_frac))

# Rename columns with dataset_tag + growth rate (with replicate suffix if needed)
colnames(copy_frac_pseudo) <- new_names

## ---- assemble compiled Peebo table (no coords yet) ----
compiled <- tibble(
  gene      = peebo$gene,       # original case from Peebo
  locus_tag = peebo$locus_tag,
  n_cond    = n_cond
) %>%
  bind_cols(as.data.frame(copy_frac_pseudo, check.names = FALSE))

## ---- add genomic coordinates from GFF via normalized gene name ----
gene_coords <- build_gene_coords(gff_path)   # gene in this df is already normalized

compiled2 <- compiled %>%
  mutate(gene_norm = normalize_gene(gene)) %>%                 # normalize Peebo gene
  left_join(gene_coords, by = c("gene_norm" = "gene")) %>%     # join on normalized names
  select(-gene_norm) %>%                                       # drop helper column
  relocate(gene, locus_tag, seqname, start, end, n_cond)

# Optional: sort by genomic position
if ("start" %in% names(compiled2)) {
  compiled2 <- compiled2 %>% arrange(start, gene)
}

chrom_len <- 4631469	#BW25113
merged <- compiled2 %>%
  mutate(gene_mid = if_else(
    start <= end,
    (start + end)/2,
    # wrap-around case:
    ((start + (end + chrom_len))/2) %% chrom_len	
  ))

merged$OriC_start <- 3918994 #BW25113
merged$OriC_end <- 3919371

sal_df2 <-merged %>% 
  mutate(
    ori_mid = if_else(
      OriC_start <= OriC_end,
      (OriC_start + OriC_end)/2,
      # wrap-around case:
      ((OriC_start + (OriC_end + chrom_len))/2) %% chrom_len	
    )
  )

sal_df2 <- sal_df2 %>%
  
  mutate(chrom_len = chrom_len) %>%
  ungroup() %>%
  
  mutate(
    dist_raw = abs(gene_mid - ori_mid),
    dist      = if_else(dist_raw > chrom_len/2,
                        chrom_len - dist_raw,
                        dist_raw),
    # 5) normalize so 0 at oriC, 1 at ter
    norm_pos  = dist / (chrom_len/2)
  )

sal_df2$type <- ifelse(sal_df2$norm_pos <= 1/3,"oriC",
                       ifelse(sal_df2$norm_pos<2/3,"mid","ter"))


## ---- write output ----
saveRDS(sal_df2, out_file)
message("Done. Wrote compiled Peebo pseudo copy-fraction table to: ", out_file)






















## ---- setup ----
pkgs <- c("readr", "dplyr", "tibble", "stringr", "rtracklayer")
lapply(pkgs, library, character.only = TRUE)

## ---- user-configurable settings ----
dataset_tag <- "valgepea_2013"

data_file <- "Valgepea_copy.csv"            # your simplified input
gff_path  <- "GCF_000005845.2.gff"     # GFF with locus_tag = KEGG/b-number
out_file_rna  <- paste0("RNA",dataset_tag, "_pseudo.rds")
out_file_prot  <- paste0("Protein",dataset_tag, "_pseudo.rds")

## ---- GFF helper: locus_tag (KEGG/b-number) -> coords + length + gene name ----
build_locus_coords <- function(gff_path) {
  if (!file.exists(gff_path)) {
    stop("GFF not found: ", gff_path)
  }
  
  gr <- rtracklayer::import(gff_path, format = "gff3")
  
  genes_gr <- gr[mcols(gr)$type == "gene"]
  
  # locus_tag usually holds b-numbers (b0002, b0003, ...)
  locus <- mcols(genes_gr)$locus_tag
  if (is.null(locus)) {
    attrs <- mcols(genes_gr)$attributes
    locus <- sub(".*locus_tag=([^;]+).*", "\\1", attrs)
  }
  
  # common gene name from GFF: prefer 'gene', fallback to 'Name'
  gene_sym <- as.character(mcols(genes_gr)$gene)
  
  coords_raw <- tibble(
    locus_tag = as.character(locus),
    gene_name = gene_sym,
    seqname   = as.character(seqnames(genes_gr)),
    start     = as.integer(start(genes_gr)),
    end       = as.integer(end(genes_gr))
  ) %>%
    dplyr::filter(!is.na(locus_tag))
  
  # Warn if multiple entries for the same locus_tag
  dup_locus <- coords_raw %>%
    dplyr::count(locus_tag) %>%
    dplyr::filter(n > 1)
  
  if (nrow(dup_locus) > 0) {
    warning(
      "Multiple GFF entries found for some locus_tags. ",
      "Using the first non-NA coordinates and gene name for each. ",
      "Example locus_tags: ",
      paste(head(dup_locus$locus_tag, 20), collapse = ", "),
      if (nrow(dup_locus) > 20) " ..."
    )
  }
  
  coords <- coords_raw %>%
    dplyr::group_by(locus_tag) %>%
    dplyr::summarise(
      gene_name = dplyr::first(gene_name[!is.na(gene_name)]),
      seqname   = dplyr::first(seqname[!is.na(seqname)]),
      start     = dplyr::first(start[!is.na(start)]),
      end       = dplyr::first(end[!is.na(end)]),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(length_bp = end - start + 1L)
  
  # Check for suspicious coordinates: start > end
  bad_coords <- coords %>%
    dplyr::filter(!is.na(start), !is.na(end), start > end)
  
  if (nrow(bad_coords) > 0) {
    warning(
      "Detected genes with start > end in GFF. Printing their coordinates."
    )
    print(
      bad_coords %>%
        dplyr::select(locus_tag, gene_name, seqname, start, end)
    )
  }
  
  coords
}

## ---- read simplified Valgepea data ----
val_raw <- readr::read_csv(data_file, col_types = cols())

# Header: KEGG_ID, gene, RNA0.11, ..., Protein0.49, MW
val_data <- val_raw   # data start at row 2; header is row 1

# Rename key columns
val_data <- val_data %>%
  dplyr::rename(
    locus_tag = KEGG_ID,
    gene      = gene
  )

## ---- identify RNA and protein columns ----
all_cols  <- names(val_data)
expr_cols <- setdiff(all_cols, c("locus_tag", "gene"))

mrna_cols <- expr_cols[stringr::str_starts(expr_cols, "RNA")]
prot_cols <- expr_cols[stringr::str_starts(expr_cols, "Protein")]

if (length(mrna_cols) == 0L || length(prot_cols) == 0L) {
  stop("Could not find RNA* and Protein* columns in Valgepea.csv.")
}


## ---- attach GFF coordinates via locus_tag and replace gene names ----
locus_coords <- build_locus_coords(gff_path)

gene_info <- val_data %>%
  dplyr::select(locus_tag, gene) %>%
  dplyr::left_join(locus_coords, by = "locus_tag") %>%
  # replace original gene column with common name from GFF when available
  dplyr::mutate(gene = dplyr::coalesce(gene_name, gene)) %>%
  dplyr::select(-gene_name)

## ---- build mRNA and protein matrices (molecules/cell) ----
# mRNA
mrna_mat <- as.matrix(val_data[mrna_cols])
storage.mode(mrna_mat) <- "numeric"
mrna_mat[!is.finite(mrna_mat)] <- 0   # NA / NaN / Inf -> 0
n_cond_mrna <- rowSums(!is.na(mrna_mat))

# Protein
prot_mat <- as.matrix(val_data[prot_cols])
storage.mode(prot_mat) <- "numeric"
prot_mat[!is.finite(prot_mat)] <- 0   # NA / NaN / Inf -> 0
n_cond_prot <- rowSums(!is.na(prot_mat))

## ---- filter genes with all-zero expression per assay ----
# mRNA: keep genes with at least one nonzero RNA value
keep_mrna <- rowSums(mrna_mat) > 0
gene_info_mrna <- gene_info[keep_mrna, , drop = FALSE]
mrna_mat       <- mrna_mat[keep_mrna, , drop = FALSE]
n_cond_mrna    <- n_cond_mrna[keep_mrna]

# Protein: keep genes with at least one nonzero Protein value
keep_prot <- rowSums(prot_mat) > 0
gene_info_prot <- gene_info[keep_prot, , drop = FALSE]
prot_mat        <- prot_mat[keep_prot, , drop = FALSE]
n_cond_prot     <- n_cond_prot[keep_prot]

## ---- condition metadata: growth rates from colnames ----
# RNA0.11 -> 0.11, Protein0.11 -> 0.11
cond_info_mrna <- tibble(cond_col = mrna_cols) %>%
  dplyr::mutate(
    mu_str = sub("^RNA", "", cond_col),
    mu_num = as.numeric(mu_str),
    mu_num_round = round(mu_num, 6),
    mu_label = sprintf("%.6f", mu_num_round),
    mu_label = gsub("\\.", "_", mu_label)
  )

cond_info_prot <- tibble(cond_col = prot_cols) %>%
  dplyr::mutate(
    mu_str = sub("^Protein", "", cond_col),
    mu_num = as.numeric(mu_str),
    mu_num_round = round(mu_num, 6),
    mu_label = sprintf("%.6f", mu_num_round),
    mu_label = gsub("\\.", "_", mu_label)
  )

new_names_mrna <- paste0(dataset_tag, "_rna__mu_",  cond_info_mrna$mu_label)
new_names_prot <- paste0(dataset_tag, "_prot__mu_", cond_info_prot$mu_label)

## ---- mRNA TPM computation (from molecules/cell + gene length) ----
length_kb <- gene_info_mrna$length_bp / 1000

# RPK = molecules/cell / length_kb
rpk_mat <- sweep(mrna_mat, 1, length_kb, "/")

# Sum of RPK per condition
rpk_sums <- colSums(rpk_mat, na.rm = TRUE)
rpk_sums[rpk_sums == 0] <- NA_real_

# TPM = RPK / sum(RPK) * 1e6
tpm_mat <- sweep(rpk_mat, 2, rpk_sums / 1e6, "/")

# Pseudo count for TPM
eps_tpm <- 1e-3
tpm_pseudo <- tpm_mat + eps_tpm

colnames(tpm_pseudo) <- new_names_mrna

## ---- protein copy fraction computation (from molecules/cell) ----
# Calculate column totals directly from the molecules/cell matrix
col_totals <- colSums(prot_mat, na.rm = TRUE)
col_totals[col_totals == 0] <- NA_real_

# Copy fraction per gene per condition
copy_frac <- sweep(prot_mat, 2, col_totals, "/")

# Pseudo count for copy fraction
eps_copy <- 1e-8
copy_frac_pseudo <- copy_frac + eps_copy

colnames(copy_frac_pseudo) <- new_names_prot

## ---- assemble compiled tables ----
compiled_protein <- gene_info_prot %>%
  dplyr::mutate(n_cond = n_cond_prot) %>%
  dplyr::bind_cols(as.data.frame(copy_frac_pseudo, check.names = FALSE)) %>%
  #dplyr::filter(!is.na(start), !is.na(end)) %>%
  dplyr::relocate(gene, locus_tag, seqname, start, end, length_bp, n_cond)

compiled_mrna <- gene_info_mrna %>%
  dplyr::mutate(n_cond = n_cond_mrna) %>%
  dplyr::bind_cols(as.data.frame(tpm_pseudo, check.names = FALSE)) %>%
  #dplyr::filter(!is.na(start), !is.na(end)) %>%
  dplyr::relocate(gene, locus_tag, seqname, start, end, length_bp, n_cond)



chrom_len <- 4641652 #MG1655

merged <- compiled_protein %>%
  mutate(gene_mid = if_else(
    start <= end,
    (start + end)/2,
    # wrap-around case:
    ((start + (end + chrom_len))/2) %% chrom_len	
  ))

merged$OriC_start <- 3925634 #MG1655
merged$OriC_end <- 3926011

sal_df2 <-merged %>% 
  mutate(
    ori_mid = if_else(
      OriC_start <= OriC_end,
      (OriC_start + OriC_end)/2,
      # wrap-around case:
      ((OriC_start + (OriC_end + chrom_len))/2) %% chrom_len	
    )
  )

sal_df2 <- sal_df2 %>%
  
  mutate(chrom_len = chrom_len) %>%
  ungroup() %>%
  
  mutate(
    dist_raw = abs(gene_mid - ori_mid),
    dist      = if_else(dist_raw > chrom_len/2,
                        chrom_len - dist_raw,
                        dist_raw),
    # 5) normalize so 0 at oriC, 1 at ter
    norm_pos  = dist / (chrom_len/2)
  )

sal_df2$type <- ifelse(sal_df2$norm_pos <= 1/3,"oriC",
                       ifelse(sal_df2$norm_pos<2/3,"mid","ter"))


saveRDS(sal_df2, out_file_prot)

chrom_len <- 4641652 #MG1655

merged <- compiled_mrna %>%
  mutate(gene_mid = if_else(
    start <= end,
    (start + end)/2,
    # wrap-around case:
    ((start + (end + chrom_len))/2) %% chrom_len	
  ))

merged$OriC_start <- 3925634 #MG1655
merged$OriC_end <- 3926011

sal_df2 <-merged %>% 
  mutate(
    ori_mid = if_else(
      OriC_start <= OriC_end,
      (OriC_start + OriC_end)/2,
      # wrap-around case:
      ((OriC_start + (OriC_end + chrom_len))/2) %% chrom_len	
    )
  )

sal_df2 <- sal_df2 %>%
  
  mutate(chrom_len = chrom_len) %>%
  ungroup() %>%
  
  mutate(
    dist_raw = abs(gene_mid - ori_mid),
    dist      = if_else(dist_raw > chrom_len/2,
                        chrom_len - dist_raw,
                        dist_raw),
    # 5) normalize so 0 at oriC, 1 at ter
    norm_pos  = dist / (chrom_len/2)
  )

sal_df2$type <- ifelse(sal_df2$norm_pos <= 1/3,"oriC",
                       ifelse(sal_df2$norm_pos<2/3,"mid","ter"))
saveRDS(sal_df2, out_file_rna)











## ---- setup ----
pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "tibble", "readr", "rtracklayer")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}))
lapply(pkgs, library, character.only = TRUE)

## ---- user-configurable settings ----
dataset_tag <- "li_2014"

li_file   <- "Li_copy.csv"   # Li Table S1

mw_file   <- "All-polypeptides-of-E.-coli-K-12-substr.-MG1655.txt"  # MW table
gff_path  <- "GCF_000005845.2.gff"                                  # MG1655 GFF

out_file  <- paste0(dataset_tag, "_all_pseudo.rds")

## ---- helper: build gene -> coords map (by normalized gene name) ----
build_gene_coords <- function(gff_path) {
  if (!file.exists(gff_path)) {
    stop("GFF not found: ", gff_path)
  }
  
  gr <- rtracklayer::import(gff_path, format = "gff3")
  genes_gr <- gr[mcols(gr)$type == "gene"]
  
  gene_sym <-as.character(mcols(genes_gr)$gene)
  
  locus <- mcols(genes_gr)$locus_tag
  if (is.null(locus)) {
    attrs <- mcols(genes_gr)$attributes
    locus <- sub(".*locus_tag=([^;]+).*", "\\1", attrs)
  }
  
  coords_raw <- tibble(
    gene_raw  = gene_sym,
    locus_tag = as.character(locus),
    seqname   = as.character(seqnames(genes_gr)),
    start     = as.integer(start(genes_gr)),
    end       = as.integer(end(genes_gr))
  ) %>%
    dplyr::filter(!is.na(gene_raw)) %>%
    dplyr::mutate(
      gene_norm = stringr::str_to_lower(stringr::str_trim(gene_raw))
    )
  
  dup_genes <- coords_raw %>%
    dplyr::filter(!is.na(start), !is.na(end)) %>%
    dplyr::count(gene_norm) %>%
    dplyr::filter(n > 1)
  
  if (nrow(dup_genes) > 0) {
    warning(
      "Multiple GFF entries for some genes; keeping first coordinates per gene. Example genes: ",
      paste(head(dup_genes$gene_norm, 20), collapse = ", "),
      if (nrow(dup_genes) > 20) " ..."
    )
  }
  
  coords_gene <- coords_raw %>%
    dplyr::filter(!is.na(start), !is.na(end)) %>%
    dplyr::group_by(gene_norm) %>%
    dplyr::summarise(
      gene_gff  = dplyr::first(gene_raw),
      locus_tag = dplyr::first(locus_tag),
      seqname   = dplyr::first(seqname),
      start     = dplyr::first(start),
      end       = dplyr::first(end),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(length_bp = end - start + 1L)
  
  coords_gene   # gene_norm, gene_gff, locus_tag, seqname, start, end, length_bp
}

## ---- helper: read MW table and build gene -> MW map ----
build_mw_map <- function(mw_file) {
  if (!file.exists(mw_file)) {
    stop("MW file not found: ", mw_file)
  }
  
  mw_raw <- readr::read_tsv(mw_file, col_types = cols())
  
  # Assumes Biocyc-style export:
  # columns like "Proteins", "Genes", "Molecular-Weight-From-Sequence"
  mw_tbl <- mw_raw %>%
    dplyr::rename(
      protein_name = Proteins,
      gene_mw      = Genes,
      mw_kda       = `Molecular-Weight-From-Sequence`
    ) %>%
    dplyr::filter(!is.na(gene_mw), !is.na(mw_kda)) %>%
    dplyr::mutate(
      gene_norm = stringr::str_to_lower(stringr::str_trim(gene_mw))
    )
  
  mw_map <- mw_tbl %>%
    dplyr::group_by(gene_norm) %>%
    dplyr::summarise(
      gene_mw      = dplyr::first(gene_mw[!is.na(gene_mw)]),
      protein_name = dplyr::first(protein_name[!is.na(protein_name)]),
      mw_kda       = dplyr::first(mw_kda[!is.na(mw_kda)]),
      n_matches    = dplyr::n(),
      .groups      = "drop"
    )
  
  dup_mw <- mw_map %>% dplyr::filter(n_matches > 1)
  if (nrow(dup_mw) > 0) {
    warning(
      "Multiple MW entries for some genes; keeping first per gene. Example genes: ",
      paste(head(dup_mw$gene_norm, 20), collapse = ", "),
      if (nrow(dup_mw) > 20) " ..."
    )
  }
  
  mw_map
}

## ---- read Li Table S1 ----
li_raw <- read_csv(li_file)

stopifnot("Gene" %in% names(li_raw))

## ---- long format + parse low-confidence vs all ----
li_long <- li_raw %>%
  dplyr::rename(gene = Gene) %>%
  tidyr::pivot_longer(
    cols      = -gene,
    names_to  = "condition",
    values_to = "value_raw"
  ) %>%
  dplyr::mutate(
    value_chr   = as.character(value_raw),
    low_conf    = stringr::str_detect(value_chr, "\\["),
    numeric_str = stringr::str_replace_all(value_chr, "\\[|\\]", ""),
    val_all     = suppressWarnings(as.numeric(numeric_str)),      # version 2: accept low-confidence
    val_hiConf  = dplyr::if_else(low_conf, NA_real_, val_all),    # version 1: low-confidence -> NA
    n_genes     = stringr::str_count(gene, stringr::fixed("+")) + 1L
  )

## ---- split geneA+geneB rows and divide synthesis rate (two versions) ----
li_split <- li_long %>%
  tidyr::separate_rows(gene, sep = "\\+") %>%
  dplyr::mutate(gene = stringr::str_trim(gene))

# Version 1: hiConf (low-confidence -> NA)
li_split_v1 <- li_split %>%
  dplyr::mutate(
    synth = dplyr::if_else(is.na(val_hiConf), NA_real_, val_hiConf / n_genes)
  ) %>%
  dplyr::select(gene, condition, synth)

# Version 2: all (low-confidence treated as normal values)
li_split_v2 <- li_split %>%
  dplyr::mutate(
    synth = dplyr::if_else(is.na(val_all), NA_real_, val_all / n_genes)
  ) %>%
  dplyr::select(gene, condition, synth)

## ---- drop genes with all NA across the 3 conditions (per version) ----
li_split_v1 <- li_split_v1 %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(!all(is.na(synth))) %>%
  dplyr::ungroup()

li_split_v2 <- li_split_v2 %>%
  dplyr::group_by(gene) %>%
  # keep genes that have at least one non-zero synthesis value
  # (NAs are ignored; they don't count as expression)
  dplyr::filter(!all(is.na(synth) | synth == 0)) %>%
  dplyr::ungroup()

## ---- wide matrices: gene x condition (two versions) ----
li_wide_v1 <- li_split_v1 %>%
  tidyr::pivot_wider(
    names_from  = condition,
    values_from = synth
  )

li_wide_v2 <- li_split_v2 %>%
  tidyr::pivot_wider(
    names_from  = condition,
    values_from = synth
  )

cond_cols <- setdiff(names(li_wide_v1), "gene")
if (length(cond_cols) == 0L) {
  stop("No condition columns detected in Li data after pivot_wider.")
}

## ---- align gene set across both versions ----
genes_union <- sort(unique(c(li_wide_v1$gene, li_wide_v2$gene)))

li_wide_v1_full <- tibble(gene = genes_union) %>%
  dplyr::left_join(li_wide_v1, by = "gene")

li_wide_v2_full <- tibble(gene = genes_union) %>%
  dplyr::left_join(li_wide_v2, by = "gene")

## ---- growth rates and output column names ----
doubling_times_min <- c(
  "MOPS complete"                    = 21.5,
  "MOPS minimal"                     = 56.3,
  "MOPS complete without methionine" = 26.5
)

cond_info <- tibble(condition = cond_cols) %>%
  dplyr::mutate(
    doubling_min = doubling_times_min[condition],
    mu_num       = log(2) / (doubling_min / 60),   # hr^-1
    mu_num_round = round(mu_num, 6),
    mu_label     = sprintf("%.6f", mu_num_round),
    mu_label     = gsub("\\.", "_", mu_label)
  )

if (any(is.na(cond_info$doubling_min))) {
  warning("Some conditions have no defined doubling time; mu_label will be NA there.")
}

new_names_v1 <- paste0(dataset_tag, "_hiConf__mu_", cond_info$mu_label)
new_names_v2 <- paste0(dataset_tag, "_all__mu_",    cond_info$mu_label)

## ---- MW map and GFF coords ----
mw_map      <- build_mw_map(mw_file)
gene_coords <- build_gene_coords(gff_path)

## ---- assemble gene_info aligned to union of genes ----
gene_info <- tibble(gene = genes_union) %>%
  dplyr::mutate(
    gene_norm = stringr::str_to_lower(stringr::str_trim(gene))
  ) %>%
  dplyr::left_join(mw_map,      by = "gene_norm") %>%
  dplyr::left_join(gene_coords, by = "gene_norm") %>%
  dplyr::mutate(
    length_bp = dplyr::coalesce(length_bp, end - start + 1L)
  )

## ---- compute copy fractions (two versions) ----
expr_mat_v1 <- as.matrix(li_wide_v1_full[cond_cols])
expr_mat_v2 <- as.matrix(li_wide_v2_full[cond_cols])

storage.mode(expr_mat_v1) <- "numeric"
storage.mode(expr_mat_v2) <- "numeric"

# n_cond for hiConf version (how many non-NA per gene in v1)
n_cond_v1 <- rowSums(!is.na(expr_mat_v1))

# Note: We skip MW multiplication to calculate copy fractions directly

# --- Version 1: hiConf ---
col_totals_v1 <- colSums(expr_mat_v1, na.rm = TRUE)
col_totals_v1[col_totals_v1 == 0] <- NA_real_
copy_frac_v1 <- sweep(expr_mat_v1, 2, col_totals_v1, "/")
eps <- 1e-8
copy_frac_pseudo_v1 <- copy_frac_v1 + eps
colnames(copy_frac_pseudo_v1) <- new_names_v1

# --- Version 2: all data (including low-confidence) ---
col_totals_v2 <- colSums(expr_mat_v2, na.rm = TRUE)
col_totals_v2[col_totals_v2 == 0] <- NA_real_
copy_frac_v2 <- sweep(expr_mat_v2, 2, col_totals_v2, "/")
copy_frac_pseudo_v2 <- copy_frac_v2 + eps
colnames(copy_frac_pseudo_v2) <- new_names_v2

## ---- assemble final compiled Li table with 6 copy-fraction columns ----
compiled_li <- gene_info %>%
  dplyr::mutate(n_cond = n_cond_v1) %>%   # n_cond from hiConf version
  dplyr::bind_cols(
    as.data.frame(copy_frac_pseudo_v1, check.names = FALSE),
    as.data.frame(copy_frac_pseudo_v2, check.names = FALSE)
  ) %>%
  dplyr::relocate(gene, locus_tag, seqname, start, end, length_bp, n_cond)

chrom_len <- 4641652 #MG1655

merged <- compiled_li %>%
  mutate(gene_mid = if_else(
    start <= end,
    (start + end)/2,
    # wrap-around case:
    ((start + (end + chrom_len))/2) %% chrom_len	
  ))

merged$OriC_start <- 3925634 #MG1655
merged$OriC_end <- 3926011

sal_df2 <-merged %>% 
  mutate(
    ori_mid = if_else(
      OriC_start <= OriC_end,
      (OriC_start + OriC_end)/2,
      # wrap-around case:
      ((OriC_start + (OriC_end + chrom_len))/2) %% chrom_len	
    )
  )

sal_df2 <- sal_df2 %>%
  
  mutate(chrom_len = chrom_len) %>%
  ungroup() %>%
  
  mutate(
    dist_raw = abs(gene_mid - ori_mid),
    dist      = if_else(dist_raw > chrom_len/2,
                        chrom_len - dist_raw,
                        dist_raw),
    # 5) normalize so 0 at oriC, 1 at ter
    norm_pos  = dist / (chrom_len/2)
  )

sal_df2$type <- ifelse(sal_df2$norm_pos <= 1/3,"oriC",
                       ifelse(sal_df2$norm_pos<2/3,"mid","ter"))

sal_df2<-sal_df2%>%dplyr::select(-"li_2014_hiConf__mu_1_934364",-"li_2014_hiConf__mu_0_738700",
                        -"li_2014_hiConf__mu_1_569390")
## ---- save ----
saveRDS(sal_df2, out_file)



















##############################################
## Goelzer B. subtilis: pseudo mass fraction
## - Input: copies / cell for 5 conditions (PYR, S, TS, CH, CHG)
## - Uses gene-specific molecular weights
## - For each condition:
##   copies × MW -> mass per gene per replicate
##   -> mass fraction per replicate
##   -> average over tech replicates (within each bio replicate)
##   -> average over 3 bio replicates
## - Across conditions:
##   drop genes with 0/NA in all conditions
##   add pseudo count 1e-8
##   add genomic coords from GFF
##   rename columns to dataset__mu_X
##############################################

## =========================
## 0) SETUP
## =========================
rm(list = ls())

# ---- libraries ----
pkgs <- c(
  "readr","readxl","dplyr","tidyr","stringr",
  "purrr","tibble","rtracklayer"
)
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}))
lapply(pkgs, library, character.only = TRUE)

## ---- user-editable paths ----
dataset_tag <- "Goelzer_2015"

xlsx_file <- "Goelzer_data_copy.xlsx"                      # your Excel file
mw_file   <- "All-polypeptides-of-B.-subtilis-subtilis-168.txt"
gff_path  <- "GCF_000009045.1.gff"                     # <-- put your B. subtilis GFF here
out_file  <- paste0(dataset_tag, "_pseudo.rds")

## ---- growth rates for each condition ----
cond_info <- tibble::tibble(
  cond_label   = c("PYR","S","TS","CH","CHG"),
  growth_rate  = c(0.3, 0.6, 0.9, 1.1, 1.5)
)

## =========================
## 1) HELPERS
## =========================

## ---- build gene -> coords (B. subtilis GFF) ----
build_gene_coords <- function(gff_path) {
  if (!file.exists(gff_path)) {
    stop("GFF not found: ", gff_path)
  }
  
  gr <- rtracklayer::import(gff_path, format = "gff3")
  
  # keep only gene features
  genes_gr <- gr[mcols(gr)$type == "gene"]
  
  gene_sym <- as.character(mcols(genes_gr)$gene)
  
  coords_raw <- tibble::tibble(
    gene_gff  = gene_sym,
    gene_norm = tolower(gene_sym),
    seqname   = as.character(seqnames(genes_gr)),
    start     = as.integer(start(genes_gr)),
    end       = as.integer(end(genes_gr))
  ) %>%
    dplyr::filter(!is.na(gene_norm))
  
  dup_genes <- coords_raw %>%
    dplyr::filter(!is.na(start), !is.na(end)) %>%
    dplyr::count(gene_norm) %>%
    dplyr::filter(n > 1)
  
  if (nrow(dup_genes) > 0) {
    warning(
      "Multiple coordinate entries in GFF for some genes; ",
      "keeping first per gene_norm. Example gene_norm: ",
      paste(head(dup_genes$gene_norm, 20), collapse = ", "),
      if (nrow(dup_genes) > 20) " ..."
    )
  }
  
  coords_gene <- coords_raw %>%
    dplyr::filter(!is.na(start), !is.na(end)) %>%
    dplyr::group_by(gene_norm) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  coords_gene
}

## ---- calculate condition-level copy fraction from replicate matrix ----
calc_copyfrac_for_condition <- function(df, tech_cols) {
  # df: data frame with replicate columns
  # tech_cols: names of replicate columns (Bio1 Tech1, etc.)
  
  # Copies per cell matrix
  copies_mat <- as.matrix(df[tech_cols])
  storage.mode(copies_mat) <- "numeric"
  
  # Copy fraction per replicate
  col_totals <- colSums(copies_mat, na.rm = TRUE)
  col_totals[col_totals == 0] <- NA_real_
  copy_frac_rep <- sweep(copies_mat, 2, col_totals, `/`)
  
  # Average over technical replicates within each biological replicate
  tech_info <- tibble::tibble(
    tech_col = tech_cols,
    bio_id   = stringr::str_match(tech_cols, "^Bio(\\d+)")[,2]
  )
  
  grp <- split(seq_along(tech_cols), tech_info$bio_id)
  
  copy_frac_bio <- sapply(grp, function(idx) {
    rowMeans(copy_frac_rep[, idx, drop = FALSE], na.rm = TRUE)
  })
  if (!is.matrix(copy_frac_bio)) {
    copy_frac_bio <- matrix(copy_frac_bio, ncol = length(grp))
  }
  
  # Average over biological replicates (3 bios)
  copy_frac_cond <- rowMeans(copy_frac_bio, na.rm = TRUE)
  copy_frac_cond[is.nan(copy_frac_cond)] <- NA_real_
  
  copy_frac_cond
}

## =========================
## 2) MOLECULAR WEIGHT TABLE
## =========================

mw_tbl_raw <- readr::read_tsv(mw_file, col_types = cols())

mw_tbl <- mw_tbl_raw %>%
  dplyr::transmute(
    gene_mw_raw = Genes,
    gene_norm   = tolower(as.character(Genes)),
    mw_seq      = as.numeric(`Molecular-Weight-From-Sequence`)
  ) %>%
  dplyr::filter(!is.na(gene_norm)) %>%
  dplyr::group_by(gene_norm) %>%
  dplyr::summarise(
    mw_seq = dplyr::first(mw_seq),
    .groups = "drop"
  )

## =========================
## 3) PER-CONDITION MASS FRACTIONS
## =========================

sheet_names <- cond_info$cond_label

cond_tables <- purrr::map(sheet_names, function(sheet_nm) {
  message("Processing sheet: ", sheet_nm)
  
  df <- readxl::read_excel(xlsx_file, sheet = sheet_nm)
  
  # Gene name column (e.g. "% Gene_Name ")
  gene_col <- grep("^% Gene_Name", names(df), value = TRUE)
  if (length(gene_col) != 1L) {
    stop("Could not uniquely identify '% Gene_Name' column in sheet ", sheet_nm)
  }
  
  df2 <- df %>%
    dplyr::rename(gene = !!gene_col) %>%
    dplyr::mutate(
      gene = as.character(gene),
      gene_norm = tolower(gene)
    )
  
  # Technical replicate columns: "Bio1 Tech1(XXX)" etc.
  tech_cols <- grep("^Bio[0-9]+ Tech[0-9]+\\(", names(df2), value = TRUE)
  if (length(tech_cols) == 0L) {
    stop("No Bio*_Tech* columns found in sheet ", sheet_nm)
  }
  
  # We skip the MW join and missing-value filtering completely!
  
  # Copy fraction for this condition
  copy_frac_cond <- calc_copyfrac_for_condition(df2, tech_cols)
  
  tibble::tibble(
    gene      = df2$gene,
    gene_norm = df2$gene_norm,
    !!sheet_nm := copy_frac_cond
  )
})

## =========================
## 4) COMBINE CONDITIONS, FILTER GENES, ADD PSEUDO
## =========================

# Wide table: gene + one column per condition (PYR, S, TS, CH, CHG)
cond_wide <- purrr::reduce(cond_tables, dplyr::full_join, by = c("gene","gene_norm"))

cond_cols <- sheet_names

# Matrix of condition-level copy fractions (as computed from replicates)
copy_mat <- as.matrix(cond_wide[cond_cols])
storage.mode(copy_mat) <- "numeric"

# 1) Remove genes that are always NA or 0 across all conditions (before pseudo)
keep_idx <- apply(copy_mat, 1, function(x) {
  any(!is.na(x) & x != 0)
})

cond_filtered <- cond_wide[keep_idx, , drop = FALSE]

# 2) Fill 0 for NA values in expression data, then add pseudo-count
eps <- 1e-8
copy_mat2 <- as.matrix(cond_filtered[cond_cols])
storage.mode(copy_mat2) <- "numeric"

# fill NA with 0
copy_mat2[is.na(copy_mat2)] <- 0

# add pseudo-count
copy_mat2 <- copy_mat2 + eps

cond_filtered[cond_cols] <- copy_mat2

# 3) Number of conditions with non-zero expression *before* pseudo-count
copy_mat_orig <- copy_mat[keep_idx, , drop = FALSE]
n_cond_vec <- rowSums(!is.na(copy_mat_orig) & copy_mat_orig != 0)

cond_filtered$n_cond <- n_cond_vec

## =========================
## 5) ADD COORDS FROM GFF
## =========================

gene_coords <- build_gene_coords(gff_path) %>%
  dplyr::select(gene_norm, seqname, start, end)

compiled <- cond_filtered %>%
  dplyr::left_join(gene_coords, by = "gene_norm") %>%
  dplyr::relocate(gene, seqname, start, end, n_cond) %>%
  dplyr::select(-gene_norm)

# Optional: sort by genomic coordinate when available
if ("start" %in% names(compiled)) {
  compiled <- compiled %>%
    dplyr::arrange(dplyr::desc(is.na(start)), start, gene)
}

## =========================
## 6) RENAME COLUMNS TO dataset__mu_*
## =========================

mu_labels <- cond_info %>%
  dplyr::mutate(
    mu_num   = round(growth_rate, 6),
    mu_label = sprintf("%.6f", mu_num),
    mu_label = gsub("\\.", "_", mu_label),
    new_name = paste0(dataset_tag, "__mu_", mu_label)
  )

# names = NEW, values = OLD  (new = old in rename())
rename_map <- setNames(mu_labels$cond_label, mu_labels$new_name)
# This creates: c("PYR" = "Goelzer_2015__mu_0_300000", ...)


chrom_len <- 4215606	#Bacillus subtilis subsp. subtilis str. 168
merged <- compiled %>%
  mutate(gene_mid = if_else(
    start <= end,
    (start + end)/2,
    # wrap-around case:
    ((start + (end + chrom_len))/2) %% chrom_len	
  ))

merged$OriC_start <- 1751 #Bacillus subtilis	NC_000964.3
merged$OriC_end <- 1938

sal_df2 <-merged %>% 
  mutate(
    ori_mid = if_else(
      OriC_start <= OriC_end,
      (OriC_start + OriC_end)/2,
      # wrap-around case:
      ((OriC_start + (OriC_end + chrom_len))/2) %% chrom_len	
    )
  )

sal_df2 <- sal_df2 %>%
  
  mutate(chrom_len = chrom_len) %>%
  ungroup() %>%
  
  mutate(
    dist_raw = abs(gene_mid - ori_mid),
    dist      = if_else(dist_raw > chrom_len/2,
                        chrom_len - dist_raw,
                        dist_raw),
    # 5) normalize so 0 at oriC, 1 at ter
    norm_pos  = dist / (chrom_len/2)
  )

sal_df2$type <- ifelse(sal_df2$norm_pos <= 1/3,"oriC",
                       ifelse(sal_df2$norm_pos<2/3,"mid","ter"))


compiled_final <- sal_df2 %>%
  dplyr::rename(!!!rename_map)
## =========================
## 7) SAVE
## =========================

saveRDS(compiled_final, file = out_file)
message("Saved compiled dataset to: ", out_file)

