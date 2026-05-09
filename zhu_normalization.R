setwd("/home/usr/Desktop/research/gene_position")
rm(list = ls())
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)

# =========================================================================
# 1. USER CONFIGURATION (MANUAL INPUT REQUIRED)
# =========================================================================
# csv_pattern: A unique string to identify which CSV files belong to which species.
# Chrom_SeqID: The RefSeq ID found in the GFF 'seqname' column.

species_config <- tibble::tribble(
  ~Species,        ~csv_pattern,   ~GFF_File,                ~Chrom_SeqID,       ~Chrom_Len, ~OriC_Start, ~OriC_End, ~MW_File,
  
  # E. coli NCM3722 
  # RefSeq ID for CP011495 is usually NZ_CP011495.1
  "E_coli",        "Ecoli",        "GCF_001043215.1.gff",    "NZ_CP011495.1",    4678046,    501975,     502352,   "Ecoli_MW.txt",
  
  # V. natriegens (Chromosome 1)
  # RefSeq ID for CP009977 is usually NZ_CP009977.1
  "V_natriegens",  "Vnat",  "GCF_001456255.1.gff",    "NZ_CP009977.1",    3248023,    2920325,      2920800,   "Vnat_MW.txt",       
  
  # B. subtilis 168
  # Standard RefSeq ID
  "B_subtilis",    "Bsub",    "GCF_000009045.1.gff",    "NC_000964.3",      4215606,    1751,          1938, "Bsub_MW.txt"         
)

# -------------------------------------------------------------------------
# 1. COMPLETE TABLE S3 METADATA (Extracted from PNAS Supplementary)
# -------------------------------------------------------------------------

metadata_s3 <- tibble::tribble(
  ~ProjectID,   ~Batch, ~SampleID, ~OriginalID, ~Species,         ~Condition_Raw,                    ~Temp_C, ~GrowthRate_hr,
  # Vibrio natriegens (Batch C - 37°C)
  "WD039LQ",    "C",    "C1",      "V1",        "V_natriegens",   "LB3 broth",                        37,      2.78,
  "WD039LQ",    "C",    "C2",      "V2",        "V_natriegens",   "0.4% glucose + 0.4% CAA",          37,      2.05,
  "WD039LQ",    "C",    "C3",      "V3",        "V_natriegens",   "0.4% glucose",                     37,      1.65,
  "WD039LQ",    "C",    "C4",      "V4",        "V_natriegens",   "0.4% glycerol",                    37,      0.86,
  "WD039LQ",    "C",    "C5",      "V5",        "V_natriegens",   "0.4% fructose",                    37,      1.29,
  "WD039LQ",    "C",    "C6",      "V6",        "V_natriegens",   "50 mM acetate",                    37,      0.26,
  
  # Vibrio natriegens (Batch L - 30°C + control)
  "XB02027LQ",  "L",    "L1",      "V1",        "V_natriegens",   "LB3 broth",                        30,      2.10,
  "XB02027LQ",  "L",    "L2",      "V2",        "V_natriegens",   "0.4% glucose + 0.4% CAA",          30,      1.63,
  "XB02027LQ",  "L",    "L3",      "V3",        "V_natriegens",   "0.4% glucose",                     30,      1.30,
  "XB02027LQ",  "L",    "L4",      "V4",        "V_natriegens",   "0.4% fructose",                    30,      1.11,
  "XB02027LQ",  "L",    "L5",      "V5",        "V_natriegens",   "0.4% glycerol",                    30,      0.70,
  "XB02027LQ",  "L",    "L6",      "V6",        "V_natriegens",   "0.4% glucose",                     37,      1.64, # Replicate of C3
  
  # Bacillus subtilis (Batch F - 37°C)
  "AD066LQ",    "F",    "F6",      "B6",        "B_subtilis",     "LB broth",                         37,      1.84,
  "AD066LQ",    "F",    "F7",      "B7",        "B_subtilis",     "0.4% glucose + 0.4% CAA",          37,      1.38,
  "AD066LQ",    "F",    "F8",      "B8",        "B_subtilis",     "0.4% glucose + 0.4% CAA",          37,      1.38, # Replicate of F7
  "AD066LQ",    "F",    "F9",      "B9",        "B_subtilis",     "0.4% glucose + 1 mM glutamate",    37,      0.80,
  "AD066LQ",    "F",    "F10",     "B10",       "B_subtilis",     "0.4% mannose + 1 mM glutamate",    37,      0.64,
  "AD066LQ",    "F",    "F11",     "B11",       "B_subtilis",     "0.4% ribose + 1 mM glutamate",     37,      0.44,
  
  # Escherichia coli (Batch F - 37°C)
  "AD067LQ",    "F",    "F12",     "E12",       "E_coli",         "0.2% glucose",                     37,      0.92,
  "AD067LQ",    "F",    "F13",     "E13",       "E_coli",         "0.2% glucose + 1 mM glutamate",    37,      1.00,
  "AD067LQ",    "F",    "F13-1",   "E13-1",     "E_coli",         "LB broth",                         37,      1.90,
  "AD067LQ",    "F",    "F13-2",   "E13-2",     "E_coli",         "0.2% glucose + 0.2% CAA",          37,      1.20,
  "AD067LQ",    "F",    "F13-3",   "E13-3",     "E_coli",         "0.2% glycerol",                    37,      0.67,
  "AD067LQ",    "F",    "F13-4",   "E13-4",     "E_coli",         "0.2% mannose",                     37,      0.42,
  
  # Escherichia coli (Batch K - 30°C + control)
  "XB02028LQ",  "K",    "K1",      "E1",        "E_coli",         "LB broth",                         30,      1.205,
  "XB02028LQ",  "K",    "K2",      "E2",        "E_coli",         "0.2% glucose + 0.2% CAA",          30,      0.915,
  "XB02028LQ",  "K",    "K3",      "E3",        "E_coli",         "0.2% glucose",                     30,      0.644,
  "XB02028LQ",  "K",    "K4",      "E4",        "E_coli",         "0.2% glycerol",                    30,      0.418,
  "XB02028LQ",  "K",    "K5",      "E5",        "E_coli",         "0.2% mannose",                     30,      0.395,
  "XB02028LQ",  "K",    "K6",      "E6",        "E_coli",         "0.2% glucose",                     37,      0.949, # Replicate of F12
  
  # Bacillus subtilis (Batch G - 43°C + control)
  "XA01797LQ",  "G",    "G1",      "B1",        "B_subtilis",     "LB broth",                         43,      2.57,
  "XA01797LQ",  "G",    "G2",      "B2",        "B_subtilis",     "0.4% glucose + 0.4% CAA",          43,      1.91,
  "XA01797LQ",  "G",    "G3",      "B3",        "B_subtilis",     "0.4% glucose + 1 mM glutamate",    43,      1.17,
  "XA01797LQ",  "G",    "G4",      "B4",        "B_subtilis",     "0.4% mannose + 1 mM glutamate",    43,      1.01,
  "XA01797LQ",  "G",    "G5",      "B5",        "B_subtilis",     "0.4% ribose + 1 mM glutamate",     43,      0.79,
  "XA01797LQ",  "G",    "G6",      "B6",        "B_subtilis",     "0.4% glucose + 1 mM glutamate",    37,      0.85 # Replicate of F9
)

# ==============================================================================
# 3. LINEAR PROCESSING LOOP
# ==============================================================================

all_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

for(i in 1:nrow(species_config)) {
  
  # --- SETUP ---
  target_sp  <- species_config$Species[i]
  csv_pat    <- species_config$csv_pattern[i]
  gff_path   <- species_config$GFF_File[i]
  chrom_id   <- species_config$Chrom_SeqID[i]
  chrom_len  <- species_config$Chrom_Len[i]
  ori_start  <- species_config$OriC_Start[i]
  ori_end    <- species_config$OriC_End[i]
  mw_path <- species_config$MW_File[i]
  
  message("\n========================================================")
  message(" PROCESSING: ", target_sp)
  message("========================================================")
  
  # --------------------------------------------------------------------------
  # STEP 1: LOAD GFF & CLEAN COORDINATES
  # --------------------------------------------------------------------------
  message("1. Loading GFF: ", gff_path)
  
  if(!file.exists(gff_path)) { message("!! ERROR: GFF missing."); next }
  
  gff_raw <- read_tsv(gff_path, comment = "#", col_names = FALSE, show_col_types = FALSE)
  
  # Extract Name and Coordinates
  # - Filter for genes on the main chromosome
  # - Extract 'Name' (Common Name) strictly
  gff_genes <- gff_raw %>%
    select(seqname = X1, start = X4, end = X5, attributes = X9) %>%
    mutate(
      gene_name_raw = str_extract(attributes, "(?<=Name=)[^;]+"),
      # If Name is missing, use locus_tag as fallback to ensure we don't drop genes just because they lack a common name
      locus_tag_fallback = str_extract(attributes, "(?<=locus_tag=)[^;]+"),
      gene_name = coalesce(gene_name_raw, locus_tag_fallback)
    ) %>%
    filter(!is.na(gene_name)) %>%
    filter(seqname == chrom_id) 
  
  message("   -> Found ", nrow(gff_genes), " genes on main chromosome.")
  
  # --- HANDLE COORDINATE DUPLICATES ---
  # Keep the first coordinate entry for genes that appear multiple times
  dup_genes <- gff_genes %>%
    count(gene_name) %>%
    filter(n > 1)
  
  if(nrow(dup_genes) > 0) {
    message("   -> Warning: Found ", nrow(dup_genes), " genes with multiple GFF entries. Keeping the first occurrence.")
  }
  
  # Final Unique Coordinate Map (keep first occurrence per gene)
  gene_map <- gff_genes %>%
    group_by(gene_name) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    select(gene_name, start, end)
  
  message("   -> Final Coordinate Map size: ", nrow(gene_map), " unique gene names.")
  
  
  # --------------------------------------------------------------------------
  # STEP 2: LOAD CSVs & AVERAGE DUPLICATES
  # --------------------------------------------------------------------------
  sp_files <- all_files[str_detect(all_files, csv_pat)]
  if(length(sp_files) == 0) { message("!! No CSVs found."); next }
  
  message("2. Processing ", length(sp_files), " CSV files...")
  
  raw_list <- lapply(sp_files, function(f) {
    d <- read_csv(f, show_col_types = FALSE)
    if(!"gene" %in% names(d)) return(NULL)
    
    d <- d %>% mutate(gene = as.character(gene))
    
    # Only preserve 'gene locus' if the current species is B. subtilis
    grp_cols <- "gene"
    if(target_sp == "B_subtilis" && "gene locus" %in% names(d)) {
      grp_cols <- c("gene", "gene locus")
    }
    
    # Average duplicates within file (same gene name in CSV rows)
    if(n_distinct(d$gene) < nrow(d)) {
      d <- d %>%
        group_by(across(all_of(grp_cols))) %>%
        summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups="drop")
    }
    return(d)
  })
  
  raw_list <- raw_list[!sapply(raw_list, is.null)]
  
  # Determine join columns for the CSV merge based on the species
  join_cols <- "gene"
  if (target_sp == "B_subtilis" && "gene locus" %in% names(raw_list[[1]])) {
    join_cols <- c("gene", "gene locus")
  }
  
  # Merge all CSVs
  expression_df <- raw_list %>% reduce(full_join, by = join_cols)
  message("   -> Combined expression data has ", nrow(expression_df), " unique gene names.")
  
  # --------------------------------------------------------------------------
  # STEP 3: MERGE DATA + COORDINATES
  # --------------------------------------------------------------------------
  merged_df <- expression_df %>%
    inner_join(gene_map, by = c("gene" = "gene_name"))
  
  if(nrow(merged_df) == 0) {
    message("!! ERROR: 0 matches between CSV 'gene' and GFF 'Name'.")
    next
  }
  message("   -> Successfully mapped ", nrow(merged_df), " genes to coordinates.")
  
  # --------------------------------------------------------------------------
  # STEP 3.5: LOAD MW, JOIN, AND CONVERT MASS FRACTION TO COPY FRACTION
  # --------------------------------------------------------------------------
  message("3.5. Loading MW data and converting to Copy Fraction...")
  
  # Load MW data (Assumes BioCyc format: 'Gene' and 'Molecular-Weight-From-Sequence')
  mw_df <- read_tsv(mw_path, show_col_types = FALSE) %>%
    select(Gene, MW = `Molecular-Weight-From-Sequence`) %>%
    filter(!is.na(MW)) %>%
    mutate(Gene = as.character(Gene)) %>%
    group_by(Gene) %>% 
    summarise(MW = mean(MW, na.rm = TRUE), .groups = "drop") # Average duplicates
  
  # Conditional Join based explicitly on Species
  if (target_sp == "B_subtilis") {
    merged_df <- merged_df %>% left_join(mw_df, by = c("gene locus" = "Gene"))
  } else {
    merged_df <- merged_df %>% left_join(mw_df, by = c("gene" = "Gene"))
  }
  
  # Count missing MWs and report severity
  n_missing_mw <- sum(is.na(merged_df$MW))
  pct_missing <- (n_missing_mw / nrow(merged_df)) * 100
  message("   -> Missing MW data for ", n_missing_mw, " proteins (", round(pct_missing, 1), "%).")
  
  if (pct_missing > 10) {
    message("   !! WARNING: High percentage of missing MW data (>10%). Because these proteins are dropped, your denominator sum will be artificially low, which will inflate the estimated copy fraction of your remaining genes.")
  } else {
    message("   -> MW missing rate is low. The denominator will be largely unaffected, meaning the copy fraction estimation for the remaining genes will remain highly accurate.")
  }
  
  # Drop genes lacking MW (we cannot calculate copy fraction without it)
  merged_df <- merged_df %>% filter(!is.na(MW))
  
  # Identify data columns to apply the math
  sp_meta <- metadata_s3 %>% filter(Species == target_sp)
  data_cols <- intersect(names(merged_df), sp_meta$SampleID)
  
  # Math: Convert Mass Fraction to Copy Fraction
  for (col in data_cols) {
    # 1. Calculate relative copies (Mass Fraction / MW)
    rel_copies <- merged_df[[col]] / merged_df$MW
    
    # 2. Normalize to a fraction (Relative Copies / Sum of Relative Copies)
    merged_df[[col]] <- rel_copies / sum(rel_copies, na.rm = TRUE)
  }
  
  # Drop the MW column as it's no longer needed
  merged_df <- merged_df %>% select(-MW)
  
  
  # --------------------------------------------------------------------------
  # STEP 4: FILTER ZERO EXPRESSION (Before Pseudo Count)
  # --------------------------------------------------------------------------
  # Identify data columns for this species
  sp_meta <- metadata_s3 %>% filter(Species == target_sp)
  data_cols <- intersect(names(merged_df), sp_meta$SampleID)
  
  if(length(data_cols) == 0) {
    message("!! ERROR: No matching data columns found."); next 
  }
  
  n_before_filter <- nrow(merged_df)
  
  # Filter: Sum of expression across all conditions must be > 0
  # We convert NAs to 0 for this check (assuming NA = not detected)
  merged_df_clean <- merged_df %>%
    mutate(temp_sum = rowSums(across(all_of(data_cols), ~ replace_na(., 0)))) %>%
    filter(temp_sum > 0) %>%
    select(-temp_sum)
  
  n_removed <- n_before_filter - nrow(merged_df_clean)
  message("   -> Removed ", n_removed, " genes with 0 expression in all conditions.")
  message("   -> Remaining genes: ", nrow(merged_df_clean))
  
  
  # --------------------------------------------------------------------------
  # STEP 5: CALCULATE POSITIONS
  # --------------------------------------------------------------------------
  if (ori_start <= ori_end) {
    ori_mid <- (ori_start + ori_end) / 2
  } else {
    ori_mid <- ((ori_start + (ori_end + chrom_len)) / 2) %% chrom_len
  }
  
  pos_df <- merged_df_clean %>%
    mutate(
      gene_mid = if_else(start <= end, (start + end)/2, ((start + (end + chrom_len))/2) %% chrom_len),
      dist_raw = abs(gene_mid - ori_mid),
      dist     = if_else(dist_raw > chrom_len/2, chrom_len - dist_raw, dist_raw),
      norm_pos = dist / (chrom_len/2),
      type = case_when(
        norm_pos <= 1/3 ~ "oriC",
        norm_pos < 2/3 ~ "mid",
        TRUE ~ "ter"
      )
    )
  
  
  # --------------------------------------------------------------------------
  # STEP 6: AVERAGE BIOLOGICAL REPLICATES
  # --------------------------------------------------------------------------
  
  # ADD PSEUDO COUNT HERE (After zero filter, before log/averaging)
  pos_df_ready <- pos_df %>% 
    mutate(across(all_of(data_cols), ~ replace_na(., 0) + 1e-8))
  
  long_df <- pos_df_ready %>%
    select(gene = gene, norm_pos, type, all_of(data_cols)) %>%
    pivot_longer(cols = all_of(data_cols), names_to = "SampleID", values_to = "Expression") %>%
    left_join(sp_meta, by = "SampleID")
  
  # Average Replicates
  final_df <- long_df %>%
    group_by(gene, norm_pos, type, Condition_Raw, Temp_C) %>%
    summarise(
      Mean_Expr = mean(Expression, na.rm = TRUE),
      Mean_GR   = mean(GrowthRate_hr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(ColName = paste0("zhu_",target_sp, "__mu_", gsub("\\.", "_", sprintf("%.6f", Mean_GR)))) %>%
    select(gene, norm_pos, type, ColName, Mean_Expr) %>%
    pivot_wider(names_from = ColName, values_from = Mean_Expr) %>%
    mutate(n_cond = ncol(.) - 3)
  
  # --------------------------------------------------------------------------
  # STEP 7: SAVE
  # --------------------------------------------------------------------------
  out_name <- paste0("zhu_",target_sp, ".rds")
  saveRDS(final_df, out_name)
  message("   -> SAVED: ", out_name, " (", nrow(final_df), " rows)")
}