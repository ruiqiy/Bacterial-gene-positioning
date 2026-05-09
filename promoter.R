setwd("/home/usr/Desktop/research/gene_position")
rm(list = ls())

genome_len <- 4639675

# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# ==========================================
# 1. Parse the GFF3 File
# ==========================================
parse_gff <- function(gff_file) {
  # Read GFF3, skipping comment lines
  # U00096.2 GFF3 usually has standard 9 columns
  gff <- read_tsv(gff_file, comment = "#", col_names = FALSE, show_col_types = FALSE)
  
  # Filter for 'gene' features only to avoid redundancy with CDS/mRNA
  genes <- gff %>%
    filter(X3 == "gene") %>%
    select(start = X4, end = X5, strand = X7, attributes = X9) %>%
    mutate(
      # Extract common gene name (prefer 'gene=' tag)
      gene_name = str_extract(attributes, "gene=[^;]+"),
      gene_name = str_remove(gene_name, "gene="),
      
      # Fallback to Name if gene tag is missing
      gene_name = if_else(is.na(gene_name), 
                          str_remove(str_extract(attributes, "Name=[^;]+"), "Name="), 
                          gene_name),
      
      # Fallback to Locus Tag if name is missing (e.g., b0001)
      gene_name = if_else(is.na(gene_name), 
                          str_remove(str_extract(attributes, "ID=[^;]+"), "ID=gene-"), 
                          gene_name)
    ) %>%
    select(gene_name, start, end, gene_strand = strand)
  
  return(genes)
}

# Load gene annotations
cat("Parsing GFF3 annotation...\n")
# Ensure the file name matches your upload
gene_db <- parse_gff("U00096.2.gff3")

# ==========================================
# 2. Mapping Function (One-to-Many Logic)
# ==========================================
map_promoters_to_genes_expanded <- function(promoter_file, genes_df) {
  # Read promoter data
  promoters <- read_csv(promoter_file, show_col_types = FALSE)
  
  if (!"strand" %in% names(promoters)) {
    promoters <- promoters %>%
      mutate(strand = str_sub(name, -1, -1))
  }
  
  # Initialize a list to store result dataframes
  results_list <- list()
  
  # Iterate through each promoter
  for (i in 1:nrow(promoters)) {
    # Extract current promoter info
    row <- promoters[i, ]
    peak <- row$peak_coordinate
    promoter_strand <- row$strand
    
    # Placeholder for matched genes
    matched_genes <- data.frame(
      mapped_gene = character(),
      mapping_type = character(),
      distance = numeric(),
      stringsAsFactors = FALSE
    )
    
    # ---------------------------------------------------
    # Step A: Check for Intragenic (Overlap)
    # ---------------------------------------------------
    # Find genes where peak is strictly inside the gene body
    intragenic_matches <- genes_df %>%
      filter(start <= peak & end >= peak)
    
    if (nrow(intragenic_matches) > 0) {
      # If overlaps found, add all of them
      current_matches <- data.frame(
        mapped_gene = intragenic_matches$gene_name,
        mapping_type = "Intragenic",
        distance = 0
      )
      matched_genes <- bind_rows(matched_genes, current_matches)
    } 
    
    # ---------------------------------------------------
    # Step B: Check for Intergenic (Nearest Downstream)
    # ---------------------------------------------------
    # Only look for intergenic if no intragenic match was found
    
    if (nrow(matched_genes) == 0) {
      
      # Determine search direction based on PROMOTER strand
      if (promoter_strand == "+") {
        # Calculate distance from Peak -> Gene Start in the forward direction
        # The %% handles the wrap-around automatically
        candidates <- genes_df %>%
          mutate(
            # Distance: How far "forward" do I have to walk from Peak to hit Gene Start?
            dist = (start - peak) %% genome_len
          ) %>%
          filter(dist > 0 & dist <= 500) %>%
          arrange(dist)
        
      } else { # Promoter Strand is "-"
        # Calculate distance from Peak -> Gene End in the reverse direction
        candidates <- genes_df %>%
          mutate(
            # Distance: How far "backward" do I have to walk from Peak to hit Gene End?
            # (Peak - End) %% Genome_Length
            dist = (peak - end) %% genome_len
          ) %>%
          filter(dist > 0 & dist <= 500) %>%
          arrange(dist)
      }
      
      if (nrow(candidates) > 0) {
        # Select the closest gene (row 1). 
        # If there's a tie for distance (rare), take top one.
        nearest_gene <- candidates[1, ]
        
        current_matches <- data.frame(
          mapped_gene = nearest_gene$gene_name,
          mapping_type = "Intergenic",
          distance = nearest_gene$dist
        )
        matched_genes <- bind_rows(matched_genes, current_matches)
      }
    }
    
    # ---------------------------------------------------
    # Step C: Combine with Original Row
    # ---------------------------------------------------
    if (nrow(matched_genes) > 0) {
      # Replicate the original row for each match
      expanded_rows <- row[rep(1, nrow(matched_genes)), ]
      expanded_rows$mapped_gene <- matched_genes$mapped_gene
      expanded_rows$mapping_type <- matched_genes$mapping_type
      expanded_rows$distance_to_gene <- matched_genes$distance
      
      results_list[[i]] <- expanded_rows
    } else {
      # Keep the row even if no map found (fill NA)
      row$mapped_gene <- NA
      row$mapping_type <- "Unmapped"
      row$distance_to_gene <- NA
      results_list[[i]] <- row
    }
  }
  
  # Combine all results into one dataframe
  final_df <- bind_rows(results_list)
  return(final_df)
}

# ==========================================
# 3. Process LB and M9 Datasets
# ==========================================

cat("Mapping LB promoters...\n")
lb_mapped <- map_promoters_to_genes_expanded("LB_promoter_regions.csv", gene_db)
saveRDS(lb_mapped, "LB_promoter.rds")

cat("Mapping M9 promoters...\n")
m9_mapped <- map_promoters_to_genes_expanded("M9_promoter_regions.csv", gene_db)
saveRDS(m9_mapped, "M9_promoter.rds")

cat("Done! Files saved with expanded rows for multiple gene mappings.\n")