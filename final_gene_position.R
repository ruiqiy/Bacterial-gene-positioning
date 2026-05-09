# ===============================================
# Setup & Libraries
# ===============================================
setwd("/home/usr/Desktop/research/gene_position")
rm(list = ls())

ENABLE_DOWNSAMPLING           <- TRUE # Set to FALSE to skip the time-consuming down-sampling

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(car)
  library(betareg) 
  library(performance)
  library(patchwork)
  library(scales)
  library(ggh4x)
  library(gridExtra)
  library(viridis)
  library(lmtest)
})

# Create output directories
dir.create("figs", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# ===============================================
# Data Loading & Organization
# ===============================================
# Define file mapping
file_list <- c(
  schmidt_2016       = "schmidt_2016_pseudo.rds",
  peebo_2015         = "peebo_2015_pseudo.rds",
  li_2014_all        = "li_2014_all_pseudo.rds",
  valgepea_2013_prot = "Proteinvalgepea_2013_pseudo.rds",
  Goelzer_2015       = "Goelzer_2015_pseudo.rds",
  zhu_B_subtilis     = "zhu_B_subtilis.rds",
  zhu_E_coli         = "zhu_E_coli.rds",
  zhu_V_natriegens   = "zhu_V_natriegens.rds"
)

# Load all datasets into a named list
DATASETS <- lapply(file_list, readRDS)

if (!is.null(DATASETS$li_2014_all)) {
  DATASETS$li_2014_all$n_cond <- 3
}

# ===============================================
# Metadata Maps
# ===============================================
assay_map <- setNames(rep("proteomics", length(DATASETS)), names(DATASETS))

name_map <- c(
  schmidt_2016       = "Schmidt dataset",
  peebo_2015         = "Peebo dataset",
  li_2014_all        = "Li dataset",
  valgepea_2013_prot = "Valgepea dataset",
  Goelzer_2015       = "Goelzer dataset",
  zhu_B_subtilis     = "Zhu B. subtilis dataset",
  zhu_E_coli         = "Zhu E. coli dataset",
  zhu_V_natriegens   = "Zhu V. natriegens dataset"
)

# Species Grouping for Consensus Analysis & Plot Ordering
ecoli_ds <- c("schmidt_2016", "peebo_2015", "li_2014_all", "valgepea_2013_prot", "zhu_E_coli")
bsub_ds  <- c("Goelzer_2015", "zhu_B_subtilis")
vnat_ds  <- c("zhu_V_natriegens")

get_species_from_ds <- function(ds) {
  case_when(
    ds %in% ecoli_ds ~ "E. coli",
    ds %in% bsub_ds  ~ "B. subtilis",
    ds %in% vnat_ds  ~ "V. natriegens",
    TRUE ~ "Other"
  )
}

# ===============================================
# QC
# ===============================================
message("Datasets detected: ", paste(names(DATASETS), collapse = ", "))

# 1. Log Transform Expression Columns
DATASETS <- lapply(DATASETS, function(df) {
  mu_idx <- which(str_detect(names(df), "__mu_"))
  if (length(mu_idx) > 0) {
    df[mu_idx] <- lapply(df[mu_idx], log)
  }
  df
})

DATASETS_FULL <- DATASETS

# 2. Gene Filtering
floor_proteomics <- log(1e-8)
message("\n[QC] Starting Prevalence Filtering (Fixed Floors)...")

DATASETS <- map2(DATASETS, names(DATASETS), function(df, nm) {
  mu_cols <- grep("__mu_", names(df), value = TRUE)
  if (length(mu_cols) == 0) return(df)
  
  expr_mat   <- as.matrix(df[, mu_cols])
  is_present <- expr_mat > floor_proteomics
  n_cond     <- ncol(expr_mat)
  keep_mask  <- rowSums(is_present, na.rm = TRUE) >= (0.5 * n_cond)
  
  df_filtered <- df[keep_mask, ]
  
  cat(sprintf("Dataset: %-20s | Kept: %5d / %5d\n", nm, nrow(df_filtered), nrow(df)))
  return(df_filtered)
})


# ===============================================
# Precomputation: Env Meta & Gene Stats
# ===============================================
# Helpers
calc_M <- function(mu) { (0.619) * mu }

template_slope <- function(mu, y) {
  M <- (0.619) * mu
  ok <- is.finite(M) & is.finite(y)
  if (sum(ok) < 2) return(NA_real_)
  cov(M[ok], y[ok]) / var(M[ok])
}

safe_z <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  if (is.finite(s) && s > 0) (x - m) / s else rep(NA_real_, length(x))
}

p_to_stars_ns <- function(p) {
  case_when(is.na(p) ~ "NA", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "N.S.")
}

# Generate Stats Lists
env_meta_list <- list()
ds_gene_stats_list <- list()

for (nm in names(DATASETS)) {
  df <- DATASETS[[nm]]
  
  # Handle duplicates if any
  if (any(duplicated(names(df)))) names(df) <- make.unique(names(df), sep = "_dup")
  
  env_cols <- names(df)[str_detect(names(df), "__mu_")]
  if (!length(env_cols)) next
  
  # Env Meta
  env_meta <- tibble(env_col = env_cols) %>%
    mutate(
      mu_str     = str_replace(env_col, "^.*__mu_", ""),
      mu_num_str = str_replace(mu_str, "^(\\d+)_?(\\d+)?(?:.*)$", "\\1.\\2"),
      mu_num_str = if_else(str_detect(mu_num_str, "\\.$"), str_replace(mu_num_str, "\\.$", ""), mu_num_str),
      mu_num     = suppressWarnings(as.numeric(mu_num_str))
    ) %>%
    filter(!is.na(mu_num)) %>%
    distinct(env_col, .keep_all = TRUE)
  
  if (nrow(env_meta) < 2) next
  env_meta_list[[nm]] <- env_meta
  
  # Gene Stats
  long <- df %>%
    select(gene, all_of(env_meta$env_col)) %>%
    pivot_longer(-gene, names_to = "env_col", values_to = "ln_expr") %>%
    left_join(env_meta %>% select(env_col, mu_num), by = "env_col") %>%
    filter(is.finite(ln_expr), is.finite(mu_num))
  
  if (nrow(long) > 0) {
    ds_gene_stats_list[[nm]] <- long %>%
      group_by(gene) %>%
      summarise(
        n_env = n(),
        mean_ln = mean(ln_expr, na.rm = TRUE),
        slope = template_slope(mu_num, ln_expr),
        .groups = "drop"
      ) %>%
      mutate(slope = ifelse(is.finite(slope), slope, NA_real_))
  }
}

# Precompute Stats for UNFILTERED data 
ds_gene_stats_list_full <- list()

for (nm in names(DATASETS_FULL)) {
  df <- DATASETS_FULL[[nm]]
  
  # Handle duplicates if any
  if (any(duplicated(names(df)))) names(df) <- make.unique(names(df), sep = "_dup")
  
  # Use the same env_meta determined previously
  env_meta <- env_meta_list[[nm]]
  if (is.null(env_meta)) next
  
  long <- df %>%
    select(gene, all_of(env_meta$env_col)) %>%
    pivot_longer(-gene, names_to = "env_col", values_to = "ln_expr") %>%
    left_join(env_meta %>% select(env_col, mu_num), by = "env_col") %>%
    filter(is.finite(ln_expr), is.finite(mu_num))
  
  if (nrow(long) > 0) {
    ds_gene_stats_list_full[[nm]] <- long %>%
      group_by(gene) %>%
      summarise(
        n_env = n(),
        mean_ln = mean(ln_expr, na.rm = TRUE),
        slope = template_slope(mu_num, ln_expr),
        .groups = "drop"
      ) %>%
      mutate(slope = ifelse(is.finite(slope), slope, NA_real_))
  }
}

# ===============================================
# PLOTTING UTILITIES
# ===============================================
FONT_FAMILY <- "sans" 
BASE_SIZE   <- 14 # Enlarged for publication

strip_labeller_dataset <- function(x) {
  x_chr <- as.character(x)
  ifelse(!is.na(name_map[x_chr]), name_map[x_chr], x_chr)
}

theme_pub <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = FONT_FAMILY),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = base_size * 0.9),
      plot.title = element_text(face = "bold", size = base_size * 1.1),
      panel.spacing = unit(1.2, "lines")
    )
}


# ===============================================
# Expression-Position Correlation
# ===============================================
message("\n[Correlation Check] Calculating correlations with continuous gene position...")
cor_results <- list()
spearman_results <- list()

for (nm in names(DATASETS)) {
  stats <- ds_gene_stats_list[[nm]]
  if (is.null(stats)) next
  
  df <- DATASETS[[nm]]
  if (!"norm_pos" %in% names(df)) next
  
  dat <- stats %>%
    left_join(df %>% select(gene, norm_pos) %>% distinct(), by = "gene") %>%
    filter(is.finite(norm_pos))
  
  # Mean Expression vs norm_pos
  dat_mean <- dat %>% filter(is.finite(mean_ln))
  if (nrow(dat_mean) >= 10) {
    ct_mean <- cor.test(dat_mean$norm_pos, dat_mean$mean_ln, method="pearson")
    cor_results[[length(cor_results) + 1]] <- tibble(
      dataset = nm,
      metric = "Mean Expression",
      estimate = ct_mean$estimate,
      conf_low = ct_mean$conf.int[1],
      conf_high = ct_mean$conf.int[2],
      p_value = ct_mean$p.value
    )
  }
  
  # Slope vs norm_pos
  dat_slope <- dat %>% filter(is.finite(slope))
  if (nrow(dat_slope) >= 10) {
    ct_slope <- cor.test(dat_slope$norm_pos, dat_slope$slope, method="pearson")
    cor_results[[length(cor_results) + 1]] <- tibble(
      dataset = nm,
      metric = "Growth-Dependent Slope",
      estimate = ct_slope$estimate,
      conf_low = ct_slope$conf.int[1],
      conf_high = ct_slope$conf.int[2],
      p_value = ct_slope$p.value
    )
  }
  
  # Spearman Mean Expression vs norm_pos
  if (nrow(dat_mean) >= 10) {
    ct_mean_sp <- cor.test(dat_mean$norm_pos, dat_mean$mean_ln, method="spearman", exact=FALSE)
    spearman_results[[length(spearman_results) + 1]] <- tibble(
      dataset = nm, metric = "Mean Expression",
      estimate = ct_mean_sp$estimate, p_value = ct_mean_sp$p.value
    )
  }
  
  # Spearman Slope vs norm_pos
  if (nrow(dat_slope) >= 10) {
    ct_slope_sp <- cor.test(dat_slope$norm_pos, dat_slope$slope, method="spearman", exact=FALSE)
    spearman_results[[length(spearman_results) + 1]] <- tibble(
      dataset = nm, metric = "Growth-Dependent Slope",
      estimate = ct_slope_sp$estimate, p_value = ct_slope_sp$p.value
    )
  }
}

if (length(cor_results) > 0) {
  cor_res_df <- bind_rows(cor_results) %>%
    mutate(
      sig = p_to_stars_ns(p_value),
      species = get_species_from_ds(dataset),
      # Create clean factor levels for the x-axis alignment
      species = factor(species, levels = c("E. coli", "B. subtilis", "V. natriegens")),
      dataset_clean = factor(strip_labeller_dataset(dataset), levels = strip_labeller_dataset(c(ecoli_ds, bsub_ds, vnat_ds)))
    )
  
  write_csv(cor_res_df, file.path("tables", "correlation_continuous_position.csv"))
  
  # Plot
  p_cor <- ggplot(cor_res_df, aes(x = dataset_clean, y = estimate, color = metric)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.5) +
    
    # Pointrange draws the estimate dot and error bars representing the 95% CI
    geom_pointrange(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.6), size = 0.8) +
    
    # Position significance stars cleanly above or below the error bars
    geom_text(aes(label = sig, 
                  y = ifelse(estimate > 0, conf_high, conf_low),
                  vjust = ifelse(estimate > 0, -0.5, 1.5)), 
              position = position_dodge(width = 0.6), size = 5, show.legend = FALSE, family = FONT_FAMILY) +
    
    facet_grid(. ~ species, scales = "free_x", space = "free_x") +
    
    scale_color_manual(values = c("Mean Expression" = "#d7191c", "Growth-Dependent Slope" = "#2c7bb6")) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(
      title = "Correlation with Continuous Gene Position (norm_pos)", 
      subtitle = "Pearson correlation estimate \u00B1 95% Confidence Interval",
      y = "Correlation Coefficient (r)", 
      x = NULL, 
      color = "Predictor"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path("figs", "correlation_continuous_position_WIDE.pdf"), p_cor, width = 14, height = 7)
}

if (length(spearman_results) > 0) {
  spearman_res_df <- bind_rows(spearman_results) %>%
    mutate(
      sig = p_to_stars_ns(p_value),
      species = get_species_from_ds(dataset),
      species = factor(species, levels = c("E. coli", "B. subtilis", "V. natriegens")),
      dataset_clean = factor(strip_labeller_dataset(dataset), levels = strip_labeller_dataset(c(ecoli_ds, bsub_ds, vnat_ds))),
      metric = factor(metric, levels = c("Mean Expression", "Growth-Dependent Slope"))
    )
  
  p_spearman <- ggplot(spearman_res_df, aes(x = dataset_clean, y = estimate, fill = metric)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sig, y = estimate + sign(estimate)*0.02), 
              position = position_dodge(width = 0.7), size = 5, show.legend = FALSE) +
    facet_grid(. ~ species, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c("Mean Expression" = "#d7191c", "Growth-Dependent Slope" = "#2c7bb6")) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(
      title = "Spearman Correlation with Continuous Gene Position (norm_pos)", 
      y = "Spearman's rho", x = NULL, fill = "Predictor"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path("figs", "correlation_continuous_position_SPEARMAN.pdf"), p_spearman, width = 14, height = 7)
}

# ===============================================
# Env-pair Wilcoxon test
# Results in Fig.S1-S3 are produced by analyses where no filter of gene expression (gene_filter == "all")
# or filter of expression fold change (fc_mode == "unfiltered") is applied.
# ===============================================
ln_expr_threshold <- log(1e-5)
cutoff_grid <- tibble(rel = c(0.00, 0.20, 0.50), abs = c(0.00, 0.10, 0.25))
FC_MIN <- 0.1; FC_MAX <- 10.0

task12_results_accumulator <- list()

run_task12 <- function(assay_lab, ds_set) {
  message("\n[Task 1.2] Running for assay: ", assay_lab)
  
  for (k in seq_len(nrow(cutoff_grid))) {
    rel_cut <- cutoff_grid$rel[k]; abs_cut <- cutoff_grid$abs[k]
    cut_lab <- paste0("rel", rel_cut*100, "_abs", abs_cut)
    
    for (gene_filter in c("all", "high")) {
      for (fc_mode in c("unfiltered", "filtered")) {
        
        rows <- list()
        for (ds_name in ds_set) {
          env_meta <- env_meta_list[[ds_name]]
          ds_gene_stats <- ds_gene_stats_list_full[[ds_name]] 
          if (is.null(env_meta) || is.null(ds_gene_stats)) next
          
          df <- DATASETS_FULL[[ds_name]]
          
          base <- ds_gene_stats %>%
            left_join(df %>% select(gene, type) %>% distinct(), by = "gene") %>%
            filter(!is.na(type), is.finite(mean_ln), is.finite(slope))
          
          if (gene_filter == "high") base <- base %>% filter(mean_ln > ln_expr_threshold)
          if (nrow(base) == 0) next
          
          env_ord <- env_meta %>% arrange(mu_num)
          expr_mat <- df %>% select(gene, type, all_of(env_ord$env_col)) %>%
            semi_join(base %>% select(gene), by = "gene")
          
          for (i in 1:(nrow(env_ord) - 1)) {
            for (j in (i+1):nrow(env_ord)) {
              slow_col <- env_ord$env_col[i]; fast_col <- env_ord$env_col[j]
              mu_slow  <- env_ord$mu_num[i];  mu_fast  <- env_ord$mu_num[j]
              
              d_abs <- abs(mu_fast - mu_slow)
              d_rel <- if (is.finite(mu_slow) && mu_slow != 0) abs(mu_fast - mu_slow)/abs(mu_slow) else Inf
              if (!(d_abs > abs_cut && d_rel > rel_cut)) next
              
              pairs <- expr_mat %>%
                transmute(gene, type, slow = exp(.data[[slow_col]]), fast = exp(.data[[fast_col]]), ratio = fast/slow) %>%
                filter(is.finite(slow), is.finite(fast))%>%
                # Remove genes with 0 expression (pseudo-count 1e-8) in BOTH conditions ---
                filter(!(slow <= 1.01e-8 & fast <= 1.01e-8))
              
              if (fc_mode == "filtered") pairs <- pairs %>% filter(ratio >= FC_MIN, ratio <= FC_MAX)
              
              do_test <- function(sub) {
                if(nrow(sub) == 0) return(c(p=NA, n=0, med_log2fc=NA))
                
                # Calculate the effect size for the volcano plot x-axis
                med_fc <- median(log2(sub$ratio), na.rm=TRUE) 
                
                if(all(sub$fast - sub$slow == 0)) return(c(p=1, n=nrow(sub), med_log2fc=med_fc))
                
                # Changed to two-sided Wilcoxon test
                w <- try(wilcox.test(sub$fast, sub$slow, paired=TRUE, alternative="two.sided", exact=FALSE), silent=TRUE)
                c(p = if(inherits(w, "try-error")) NA else w$p.value, n = nrow(sub), med_log2fc = med_fc)
              }
              
              res_o <- do_test(pairs %>% filter(type == "oriC"))
              res_t <- do_test(pairs %>% filter(type == "ter"))
              
              common_info <- tibble(assay=assay_lab, dataset=ds_name, cutoff=cut_lab, 
                                    gene_filter=gene_filter, fc_mode=fc_mode, n_total=nrow(pairs))
              
              # Append the median log2FC to the results output
              rows[[length(rows)+1]] <- bind_cols(common_info, tibble(group="oriC", p=res_o['p'], n_genes=res_o['n'], med_log2fc=res_o['med_log2fc']))
              rows[[length(rows)+1]] <- bind_cols(common_info, tibble(group="ter",  p=res_t['p'], n_genes=res_t['n'], med_log2fc=res_t['med_log2fc']))
            }
          }
        }
        
        if (length(rows) > 0) {
          res <- bind_rows(rows)
          out_tag <- paste(assay_lab, gene_filter, fc_mode, cut_lab, sep = "_")
          write_csv(res, file.path("tables", paste0("task1_2_pairs_", out_tag, ".csv")))
          
          fractions <- res %>% filter(is.finite(p)) %>%
            group_by(assay, cutoff, gene_filter, fc_mode, dataset, group) %>%
            summarise(m_pairs=n(), frac_raw=mean(p<0.05), frac_bonf=mean(p<(0.05/n())), .groups="drop")
          
          write_csv(fractions, file.path("tables", paste0("task1_2_fractions_", out_tag, ".csv")))
          task12_results_accumulator[[out_tag]] <<- res
        }
      }
    }
  }
}

proteomics_ds <- names(DATASETS)[assay_map[names(DATASETS)] == "proteomics"]
run_task12("proteomics", proteomics_ds)

# Meta-Volcano
plot_task12_volcano <- function(res_all, setting_str, out_pdf) {
  
  # --- Set maximum limits for both axes ---
  MAX_LOGP <- 10 
  MAX_FC <- 1   
  
  df_plot <- res_all %>% 
    filter(is.finite(p), is.finite(med_log2fc)) %>%
    mutate(
      species = get_species_from_ds(as.character(dataset)),
      species = factor(species, levels = c("E. coli", "B. subtilis", "V. natriegens")),
      logp = -log10(p),
      
      # Cap the Y-axis (p-values)
      plot_logp = ifelse(logp > MAX_LOGP, MAX_LOGP, logp),
      
      # Cap the X-axis (fold-change) symmetrically ---
      plot_log2fc = case_when(
        med_log2fc > MAX_FC ~ MAX_FC,
        med_log2fc < -MAX_FC ~ -MAX_FC,
        TRUE ~ med_log2fc
      ),
      
      clean_dataset = strip_labeller_dataset(as.character(dataset))
    )
  
  # Calculate Bonferroni threshold PER SPECIES 
  bonf_lines <- df_plot %>%
    group_by(species) %>%
    summarise(
      m = n(),
      x_bonf = -log10(0.05 / m),
      .groups = "drop"
    )
  
  # Check if any points actually exceeded the limits in the raw data
  has_pos_cap <- any(res_all$med_log2fc > MAX_FC, na.rm = TRUE)
  has_neg_cap <- any(res_all$med_log2fc < -MAX_FC, na.rm = TRUE)
  
  # Update AES to use the new plot_log2fc for the X-axis
  p <- ggplot(df_plot, aes(x = plot_log2fc, y = plot_logp, color = group, shape = clean_dataset)) +
    # Reference lines
    geom_vline(xintercept = 0, linetype = "solid", color = "grey80", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +
    geom_hline(data = bonf_lines, aes(yintercept = x_bonf), linetype = "dotted", color = "black", linewidth = 0.8) +
    
    # Volcano points
    geom_point(size = 2.5, alpha = 0.8) +
    
    # Facet by species (free_x keeps V. natriegens zoomed in cleanly)
    facet_wrap(~ species, scales = "free_x", ncol = 3) +
    
    # Aesthetics mapping
    scale_color_manual(
      values = c("oriC" = "#2c7bb6", "ter" = "#d7191c"),
      name = "Region"
    ) +
    scale_shape_manual(
      values = c(1, 2, 0, 5, 6, 3, 4, 8), 
      name = "Dataset"
    ) +
    
    # Format Y-axis to handle the capped points
    scale_y_continuous(
      limits = c(0, MAX_LOGP * 1.05), 
      breaks = seq(0, MAX_LOGP, by = 5),
      labels = function(x) ifelse(x == MAX_LOGP, paste0("\u2265 ", MAX_LOGP), x) 
    ) +
    
    # Format X-axis to label the capped points dynamically ---
    scale_x_continuous(
      labels = function(x) {
        dplyr::case_when(
          is.na(x) ~ "",
          # Only add the >= sign if the break is at MAX_FC AND dots were actually capped
          x >= MAX_FC & has_pos_cap ~ paste0("\u2265 ", MAX_FC),
          x <= -MAX_FC & has_neg_cap ~ paste0("\u2264 -", MAX_FC),
          TRUE ~ as.character(x)
        )
      }
    ) +
    
    labs(
      title = paste("Meta-Volcano Plot", setting_str),
      x = expression(Median~log[2]("Fast" / "Slow"~Expression)),
      y = expression(-log[10](italic(p)*"-value"))
    ) +
    theme_pub() +
    theme(
      legend.position = "right",
      legend.box = "vertical"
    )
  
  ggsave(out_pdf, p, width = 14, height = 6, device = grDevices::cairo_pdf)
}

# Execute Volcano Plots
if (length(task12_results_accumulator) > 0) {
  for (tag in names(task12_results_accumulator)) {
    res_df <- task12_results_accumulator[[tag]]
    parts <- str_split(tag, "_")[[1]] 
    gf <- parts[2]; fm <- parts[3]; cl <- paste(parts[4:length(parts)], collapse = "_")
    out_pdf <- file.path("figs", paste0("task1_2_volcano_COMBINED_", gf, "_", fm, "_", cl, ".pdf"))
    plot_task12_volcano(res_df, tag, out_pdf)
  }
}

# ===============================================
#  Beta regression models
# ===============================================
task21_results_accumulator <- list()

run_task21 <- function(assay_lab, ds_set) {
  message("\n[Task 2.1] Running Model 3 for assay: ", assay_lab)
  all_rows <- list()
  
  for (ds_name in ds_set) {
    df <- DATASETS[[ds_name]]
    ds_gene_stats <- ds_gene_stats_list[[ds_name]]
    if (is.null(ds_gene_stats)) next
    
    base <- df %>% select(gene, norm_pos, type) %>%
      left_join(ds_gene_stats %>% select(gene, mean_ln, slope), by = "gene") %>%
      filter(is.finite(norm_pos), is.finite(mean_ln), is.finite(slope))
    
    for (gene_filter in c("all", "high")) {
      dat <- base
      if (gene_filter == "high") dat <- dat %>% filter(mean_ln > ln_expr_threshold)
      if (nrow(dat) < 10) next
      
      dat <- dat %>% mutate(z_slope = safe_z(slope), z_mean = safe_z(mean_ln))
      
      cor_p <- tryCatch(cor.test(dat$z_slope, dat$z_mean)$estimate, error=function(e) NA)
      vif   <- if(!is.na(cor_p)) 1/(1-cor_p^2) else NA
      
      fit_full <- try(betareg(norm_pos ~ z_slope + z_mean, data = dat), silent=TRUE)
      if (inherits(fit_full, "try-error")) next
      
      s <- summary(fit_full)
      cf <- s$coefficients$mean
      
      p_diff <- tryCatch(linearHypothesis(fit_full, "z_slope = z_mean")$`Pr(>Chisq)`[2], error=function(e) NA)
      
      efron_r2 <- function(m, y) { 1 - sum((y - fitted(m))^2) / sum((y - mean(y))^2) }
      r2_full <- efron_r2(fit_full, dat$norm_pos)
      
      all_rows[[length(all_rows)+1]] <- tibble(
        assay=assay_lab, dataset=ds_name, gene_filter=gene_filter, n_genes=nrow(dat),
        term=c("z_slope", "z_mean"),
        estimate=cf[c("z_slope","z_mean"), "Estimate"],
        se=cf[c("z_slope","z_mean"), "Std. Error"],
        p=cf[c("z_slope","z_mean"), "Pr(>|z|)"],
        p_wald_diff=p_diff, R2_efron_full=r2_full, vif=vif
      )
    }
  }
  
  if (length(all_rows) > 0) {
    out <- bind_rows(all_rows)
    write_csv(out, file.path("tables", paste0("task2_1_model3_", assay_lab, ".csv")))
    task21_results_accumulator[[assay_lab]] <<- out
  }
}

run_task21("proteomics", proteomics_ds)

plot_task21_combined <- function(out21, gene_filter_val, out_pdf) {
  
  # Group datasets: E. coli then B. subtilis. Remove V. natriegens for individual plots.
  plot_levels <- c(ecoli_ds, bsub_ds)
  
  df_plot <- out21 %>% 
    filter(gene_filter == gene_filter_val) %>%
    filter(dataset %in% plot_levels) %>%
    mutate(
      dataset = factor(dataset, levels = plot_levels),
      species = get_species_from_ds(as.character(dataset)),
      species = factor(species, levels = c("E. coli", "B. subtilis"))
    )
  
  # 1. Pre-calculate the y-position for the Wald bracket
  #    We need the max height of the bars+error or the stars to ensure the bracket clears them.
  wald_dat <- df_plot %>%
    group_by(species, dataset) %>%
    summarise(
      p_wald = unique(p_wald_diff), # Assuming p_wald is constant per group
      # Calculate the highest point of the graphic elements (bar + error + space for star)
      max_top = max(
        pmax(
          estimate + se,            # Top of positive error bar
          ifelse(estimate > 0, estimate + se + (0.05 * abs(estimate)), -Inf) # Est. space for top star
        ), 
        na.rm = TRUE
      ),
      # Set bracket y_pos slightly above that maximum
      y_pos = max_top * 1.2, 
      .groups = "drop"
    )
  
  # 2. Plotting Code
  p <- ggplot(df_plot, aes(x = term, y = estimate)) +
    geom_col(fill = "grey75", color = "grey25", width = 0.6) +
    geom_errorbar(aes(ymin = estimate - se, ymax = estimate + se), width = 0.2) +

    geom_text(
      aes(
        label = p_to_stars_ns(p), 
        # If positive, put at top of error bar; if negative, put at bottom
        y = ifelse(estimate > 0, estimate + se, estimate - se),
        # If positive, vjust -0.5 (up); if negative, vjust 1.5 (down)
        vjust = ifelse(estimate > 0, -0.5, 1.5)
      ), 
      size = 5, 
      family = FONT_FAMILY
    ) +
    
    # --- WALD BRACKET ---
    geom_segment(data = wald_dat, aes(x = 1, xend = 2, y = y_pos, yend = y_pos), inherit.aes = FALSE) +
    geom_text(data = wald_dat, aes(x = 1.5, y = y_pos, label = p_to_stars_ns(p_wald)), 
              vjust = -0.5, size = 5, family = FONT_FAMILY, inherit.aes = FALSE) +
    
    # Faceting and Scales
    ggh4x::facet_nested_wrap(vars(species, dataset), labeller = labeller(dataset = strip_labeller_dataset),
                             scales = "free_y", axes = "all", ncol = 3) +

    # Increased top expansion from 0.3 to 0.5 to add more vertical space
    scale_y_continuous(expand = expansion(mult = c(0.2, 0.5))) + 
    labs(title = paste("Task 2.1 Coefficients -", gene_filter_val), y = "Estimate") +
    theme_pub()
  
  ggsave(out_pdf, p, width = 11, height = 8.6, device = grDevices::cairo_pdf)
}

if (!is.null(task21_results_accumulator[["proteomics"]])) {
  out21 <- task21_results_accumulator[["proteomics"]]
  for (gf in c("all", "high")) {
    out_pdf <- file.path("figs", paste0("task2_1_model3_coeffs_COMBINED2_", gf, ".pdf"))
    plot_task21_combined(out21, gf, out_pdf)
  }
}

# ===============================================
# Meta-Analysis: E. coli, B. subtilis, V. nat
# ===============================================
message("\n[Meta-Analysis] Starting Consensus Analysis...")

get_r2 <- function(model, data) {
  # Retrieve the exact response vector used by the model (handles omitted NAs)
  yo <- model$y
  yp <- predict(model, type = "response")
  
  # Efron's Pseudo-R2
  1 - (sum((yo - yp)^2) / sum((yo - mean(yo))^2))
}

# 1. Helper to extract slope data from list of datasets
extract_meta_stats <- function(ds_names) {
  res_list <- lapply(ds_names, function(nm) {
    if (is.null(DATASETS[[nm]])) return(NULL) 
    
    df <- DATASETS[[nm]]
    env_meta <- env_meta_list[[nm]]
    if (is.null(env_meta) || nrow(env_meta) < 3) return(NULL)
    
    long_dat <- df %>%
      select(gene, type, norm_pos, all_of(env_meta$env_col)) %>%
      pivot_longer(cols = all_of(env_meta$env_col), names_to = "env_col", values_to = "ln_expr") %>%
      left_join(env_meta %>% select(env_col, mu_num), by = "env_col") %>%
      mutate(M = 0.619 * mu_num) %>%
      filter(is.finite(ln_expr), is.finite(M))
    
    long_dat %>%
      group_by(gene, type, norm_pos) %>%
      filter(n() >= 3, var(M) > 0) %>%
      summarise(
        n_obs = n(), mean_ln = mean(ln_expr),
        sd_ln = sd(ln_expr), # NEW
        mod = list(lm(ln_expr ~ M)), 
        .groups = "drop"
      ) %>%
      mutate(
        se_mean = sd_ln / sqrt(n_obs), # NEW
        slope = map_dbl(mod, ~coef(.x)[2]),
        se_slope = map_dbl(mod, ~summary(.x)$coefficients[2, "Std. Error"]),
        dataset = nm,
        gene = tolower(trimws(gene)) # Normalize gene names
      ) %>%
      select(-mod) %>%
      filter(is.finite(slope), is.finite(se_slope), is.finite(se_mean), se_mean > 0)%>%
      mutate(
       mean_ln_z = as.numeric(scale(mean_ln)),
        se_mean_z = se_mean / sd(mean_ln, na.rm = TRUE) # Scale SE by the same standard deviation
              )
  })
  bind_rows(res_list)
}


# 2. Helper to Calculate Consensus and Run Regression
run_consensus_pipeline <- function(ds_names, group_name, min_datasets = 1) {
  meta_df <- extract_meta_stats(ds_names)
  if (nrow(meta_df) == 0) return(NULL)
  
  # Consensus Calc
  cons_df <- meta_df %>%
    group_by(gene) %>%
    filter(n_distinct(dataset) >= min_datasets) %>%
    summarise(
      norm_pos = mean(norm_pos, na.rm = TRUE),
      weights_mean_sum = sum(1/se_mean^2, na.rm = TRUE),
      global_mean_ln = sum(mean_ln * (1/se_mean^2), na.rm = TRUE) / weights_mean_sum,
      weights_slope_sum = sum(1/se_slope^2, na.rm = TRUE),
      consensus_slope = sum(slope * (1/se_slope^2), na.rm = TRUE) / weights_slope_sum,
      .groups = "drop"
    )
  
  if (nrow(cons_df) < 10) return(NULL)
  
  # Standardize for Analysis
  analysis_df <- cons_df %>%
    filter(is.finite(consensus_slope), is.finite(global_mean_ln)) %>%
    mutate(
      norm_pos = norm_pos,
      scale_slope = as.numeric(scale(consensus_slope)), 
      scale_mean  = as.numeric(scale(global_mean_ln))
    )
  
  # --- MODELS ---
  # 1. Full Model
  mod_full <- try(betareg(norm_pos ~ scale_slope + scale_mean, data = analysis_df), silent=TRUE)
  if (inherits(mod_full, "try-error")) return(NULL)
  
  # --- SQUARED PEARSON CORRELATION (Explanatory Power) ---
  # 1. Get predicted values (response scale)
  obs_used  <- mod_full$y
  preds     <- fitted(mod_full)
  
  # Calculate Pearson correlation between Observed and Predicted
  pearson_r <- cor(obs_used, preds)
  
  # Square it to get the "Variance Explained" equivalent
  r2_val <- pearson_r^2
  
  # Stats
  p_diff <- tryCatch(linearHypothesis(mod_full, "scale_slope = scale_mean")$`Pr(>Chisq)`[2], error=function(e) NA)
  
  cor_p <- cor(analysis_df$scale_slope, analysis_df$scale_mean)

  vif_val <- as.numeric(1/(1-cor_p^2))
  
  cf <- summary(mod_full)$coefficients$mean
  
  # Calculate simple Pearson correlation (r) and p-value for each predictor vs Position
  cor_slope_mean <- cor.test(analysis_df$scale_slope, analysis_df$scale_mean, method = "pearson")
  
  # --- OUTPUT TIBBLE ---
  tibble(
    assay = "Consensus",
    dataset = group_name,
    n_genes = nrow(analysis_df),
    
    term = c("scale_slope", "scale_mean"),
    estimate = cf[c("scale_slope", "scale_mean"), "Estimate"],
    se = cf[c("scale_slope", "scale_mean"), "Std. Error"],
    p = cf[c("scale_slope", "scale_mean"), "Pr(>|z|)"],
    p_wald_diff = p_diff,
    
    R2_Pearson_Full = r2_val,
    vif = vif_val,
    pearson_r = cor_slope_mean$estimate,
    pearson_p = cor_slope_mean$p.value
  )
}

# 3. Run for 3 Groups
res_ecoli <- run_consensus_pipeline(ecoli_ds, "E. coli Consensus", min_datasets = 3)
res_bsub  <- run_consensus_pipeline(bsub_ds, "B. subtilis Consensus", min_datasets = 2) # Strict: need presence in both if using 2
res_vnat  <- run_consensus_pipeline(vnat_ds, "V. natriegens", min_datasets = 1)

# Combine
all_cons_res <- bind_rows(res_ecoli, res_bsub, res_vnat)

# Only proceed if we have results
if (nrow(all_cons_res) > 0) {
  write_csv(all_cons_res, file.path("tables", "consensus_model_results_with_var_part.csv"))
  print(head(all_cons_res))
  all_cons_res <- all_cons_res %>%
    mutate(dataset = factor(dataset, levels = c("E. coli Consensus", "B. subtilis Consensus", "V. natriegens")))
  
  # Wald Annotation
  wald_ann <- all_cons_res %>%
    group_by(dataset) %>%  # Changed from Group
    summarise(
      p_wald = unique(p_wald_diff),
      y_pos = max(estimate + 1.96 * se) + 0.1 * max(abs(estimate)),
      .groups = "drop"
    )
  
  # Plot
  p_cons <- ggplot(all_cons_res, aes(x = term, y = estimate)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_pointrange(aes(ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se), size = 0.8, color = "black") +
    
    # Wald Bracket
    geom_segment(data = wald_ann, aes(x = 1, xend = 2, y = y_pos, yend = y_pos), inherit.aes = FALSE) +
    geom_text(data = wald_ann, aes(x = 1.5, y = y_pos + 0.05 * max(abs(all_cons_res$estimate)), label = p_to_stars_ns(p_wald)), 
              vjust = -0.5, size = 6, family = if(exists("FONT_FAMILY")) FONT_FAMILY else "sans", inherit.aes = FALSE) +
    
    # Stars
    geom_text(aes(label = p_to_stars_ns(p), 
                  y = ifelse(estimate > 0, estimate + 1.96 * se, estimate - 1.96 * se),
                  vjust = ifelse(estimate > 0, -0.5, 1.5)), 
              size = 6, family = if(exists("FONT_FAMILY")) FONT_FAMILY else "sans") +
    
    facet_wrap(~dataset, scales = "free_y") +  # Changed from Group
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.3))) +
    labs(title = "Consensus Models: Drivers of Position", y = "Estimate (Standardized)", x = "Predictor") +
    theme_pub()
  
  ggsave("figs/consensus_beta_regression_ALL.pdf", p_cons, width = 12, height = 5, device = grDevices::cairo_pdf)
} else {
  message("No consensus models converged. Skipping figure.")
}

# ==============================================================================
# CONDITION DOWNSAMPLING ANALYSIS
# ==============================================================================
if (ENABLE_DOWNSAMPLING) {
message("\n[Downsampling] Starting Condition Downsampling Analysis...")

# 1. Helper: Extract stats using specific selected columns per dataset
extract_meta_stats_subset <- function(selected_cols_map) {
  # selected_cols_map is a list: names = dataset, value = vector of column names
  
  res_list <- lapply(names(selected_cols_map), function(nm) {
    if (is.null(DATASETS[[nm]])) return(NULL) 
    
    target_cols <- selected_cols_map[[nm]]
    # We strictly need >= 3 points to calculate a meaningful slope/error
    if (length(target_cols) < 3) return(NULL) 
    
    df <- DATASETS[[nm]]
    env_meta <- env_meta_list[[nm]] %>% filter(env_col %in% target_cols)
    
    if (nrow(env_meta) < 3) return(NULL)
    
    long_dat <- df %>%
      select(gene, type, norm_pos, all_of(env_meta$env_col)) %>%
      pivot_longer(cols = all_of(env_meta$env_col), names_to = "env_col", values_to = "ln_expr") %>%
      left_join(env_meta %>% select(env_col, mu_num), by = "env_col") %>%
      mutate(M = 0.619 * mu_num) %>%
      filter(is.finite(ln_expr), is.finite(M))
    
    long_dat %>%
      group_by(gene, type, norm_pos) %>%
      filter(n() >= 3, var(M) > 0) %>% # Ensure variance in growth rate exists
      summarise(
        n_obs = n(),
        mean_ln = mean(ln_expr),
        sd_ln = sd(ln_expr), 
        mod = list(lm(ln_expr ~ M)), 
        .groups = "drop"
      ) %>%
      mutate(
        se_mean = sd_ln / sqrt(n_obs),
        slope = map_dbl(mod, ~coef(.x)[2]),
        se_slope = map_dbl(mod, ~summary(.x)$coefficients[2, "Std. Error"]),
        dataset = nm,
        gene = tolower(trimws(gene)) # Normalize gene names for matching
      ) %>%
      select(-mod) %>%
      filter(is.finite(slope), is.finite(se_slope), is.finite(se_mean), se_mean > 0)
  })
  
  bind_rows(res_list)
}

# 2. Main Downsampling Simulation Function (Pre-Filtered)
run_downsampling_sim <- function(ds_list, group_label, step_size = 1, min_k = 5, n_reps = 100, min_datasets_req = 1) {
  
  # --- STEP A: DEFINE THE CORE UNIVERSE ---
  # Identify genes that meet the Consensus requirement (e.g., present in >= 3 datasets)
  # This uses the FULL datasets, ignoring the downsampling for a moment, to establish the list.
  all_genes_list <- lapply(ds_list, function(nm) {
    if(!is.null(DATASETS[[nm]])) return(tolower(trimws(DATASETS[[nm]]$gene)))
    return(NULL)
  })
  
  gene_counts <- table(unlist(all_genes_list))
  valid_genes <- names(gene_counts)[gene_counts >= min_datasets_req]
  
  message(sprintf("  %s: Restricting simulation to %d 'Core' genes (present in >= %d datasets)", 
                  group_label, length(valid_genes), min_datasets_req))
  
  # --- STEP B: BUILD GLOBAL CONDITION POOL ---
  condition_pool <- c()
  for (nm in ds_list) {
    if (!is.null(env_meta_list[[nm]])) {
      cols <- env_meta_list[[nm]]$env_col
      condition_pool <- c(condition_pool, paste0(nm, ":::", cols))
    }
  }
  
  max_conds <- length(condition_pool)
  if (max_conds < min_k) return(NULL)
  
  # Sequence of condition counts to test
  k_seq <- seq(min_k, max_conds, by = step_size)
  if (tail(k_seq, 1) != max_conds) k_seq <- c(k_seq, max_conds)
  k_seq <- unique(k_seq)
  
  results_list <- list()
  
  # --- STEP C: SIMULATION LOOP ---
  for (k in k_seq) {
    message(sprintf("    Processing k = %d / %d", k, max_conds))
    
    for (rep_i in 1:n_reps) {
      
      # 1. Randomly sample k conditions from the GLOBAL pool
      sampled_tokens <- sample(condition_pool, k)
      
      # 2. Reconstruct the map: Dataset -> [Selected Columns]
      curr_cols_map <- list()
      for (token in sampled_tokens) {
        parts <- str_split(token, ":::")[[1]]
        d_name <- parts[1]; c_name <- parts[2]
        if (is.null(curr_cols_map[[d_name]])) curr_cols_map[[d_name]] <- c(c_name)
        else curr_cols_map[[d_name]] <- c(curr_cols_map[[d_name]], c_name)
      }
      
      # 3. Extract Stats (Dataset-level slopes)
      meta_df <- extract_meta_stats_subset(curr_cols_map)
      
      if (is.null(meta_df) || nrow(meta_df) == 0) next
      
      # 4. FILTER: Keep ONLY the Core Genes identified in Step A
      meta_df <- meta_df %>% filter(gene %in% valid_genes)
      
      if (nrow(meta_df) == 0) next
      
      # 5. Consensus Calculation
      cons_df <- meta_df %>%
        group_by(gene) %>%
        summarise(
          norm_pos = mean(norm_pos, na.rm = TRUE),
          weights_mean_sum = sum(1/se_mean^2, na.rm = TRUE),
          global_mean_ln = sum(mean_ln * (1/se_mean^2), na.rm = TRUE) / weights_mean_sum,
          weights_slope_sum = sum(1/se_slope^2, na.rm = TRUE),
          consensus_slope = sum(slope * (1/se_slope^2), na.rm = TRUE) / weights_slope_sum,
          .groups = "drop"
        )
      
      # 6. Fit Beta Regression
      analysis_df <- cons_df %>%
        filter(is.finite(consensus_slope), is.finite(global_mean_ln), 
               norm_pos > 0, norm_pos < 1) %>%
        mutate(
          scale_slope = as.numeric(scale(consensus_slope)), 
          scale_mean  = as.numeric(scale(global_mean_ln))
        )
      
      fit <- try(betareg(norm_pos ~ scale_slope + scale_mean, data = analysis_df), silent=TRUE)
      
      if (!inherits(fit, "try-error")) {
        cf <- coef(fit)
        results_list[[length(results_list) + 1]] <- tibble(
          dataset = group_label,
          n_conditions = k,
          replicate = rep_i,
          n_genes_used = nrow(analysis_df),
          beta_slope = cf["scale_slope"],
          beta_mean = cf["scale_mean"],
          beta_diff = cf["scale_slope"] - cf["scale_mean"]
        )
      }
    }
  }
  
  bind_rows(results_list)
}

# 3. Run Simulation with MATCHED FILTERS
set.seed(46) 

# E. coli: 3 dataset minimum (Matches 'run_consensus_pipeline' call)
ds_res_ecoli <- run_downsampling_sim(ecoli_ds, "E. coli Consensus", step_size = 5, min_datasets_req = 3)

# B. subtilis: 2 dataset minimum (Matches 'run_consensus_pipeline' call)
ds_res_bsub  <- run_downsampling_sim(bsub_ds, "B. subtilis Consensus", step_size = 1, min_datasets_req = 2)

# V. natriegens: 1 dataset minimum (Matches 'run_consensus_pipeline' call)
ds_res_vnat  <- run_downsampling_sim(vnat_ds, "V. natriegens", step_size = 1, min_datasets_req = 1)

all_sim_res <- bind_rows(ds_res_ecoli, ds_res_bsub, ds_res_vnat)

# 4. Summarize & Plot
if (!is.null(all_sim_res) && nrow(all_sim_res) > 0) {
  
  # A. Calculate Mean and Standard Error (SE) across replicates
  # We pivot longer FIRST to calculate stats for all metrics (Slope, Mean, Diff) simultaneously
  sim_summary <- all_sim_res %>%
    pivot_longer(
      cols = c(beta_slope, beta_mean, beta_diff), 
      names_to = "metric_raw", 
      values_to = "val"
    ) %>%
    group_by(dataset, n_conditions, metric_raw) %>%
    summarise(
      mean_val = mean(val, na.rm=TRUE),
      se_val   = sd(val, na.rm=TRUE) / sqrt(n()), # Standard Error calculation
      .groups = "drop"
    ) %>%
    mutate(
      metric_label = case_when(
        metric_raw == "beta_slope" ~ "Slope (Growth)",
        metric_raw == "beta_mean"  ~ "Mean (Expression)",
        metric_raw == "beta_diff"  ~ "Difference (Slope - Mean)"
      ),
      dataset = factor(dataset, levels = c("E. coli Consensus", "B. subtilis Consensus", "V. natriegens"))
    )
  
  # Save Data
  write_csv(all_sim_res, file.path("tables", "downsampling_simulation_raw.csv"))
  write_csv(sim_summary, file.path("tables", "downsampling_simulation_summary.csv"))
  
  # Helper for Integer X-Axis Labels
  integer_breaks <- function(x) {
    unique(floor(pretty(seq(min(x), max(x)), n = 5)))
  }
  
  # B. Create Plot with Error Bars
  p_sim <- ggplot(sim_summary, aes(x = n_conditions, y = mean_val, color = metric_label)) +
    # Error Bars (Mean +/- SE)
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), 
                  width = 0.2, alpha = 0.6) +
    
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    
    facet_wrap(~dataset, scales = "free_x") +
    
    # Force Integer X-Axis
    scale_x_continuous(breaks = integer_breaks) +
    
    scale_color_manual(
      values = c("Slope (Growth)" = "#2c7bb6", 
                 "Mean (Expression)" = "#d7191c", 
                 "Difference (Slope - Mean)" = "black"),
      name = "Coefficient"
    ) +
    
    labs(
      title = "Effect of Increasing Condition Count on Model Coefficients",
      subtitle = "Mean \u00B1 Standard Error (SE) of 100 replicates. Gene universe fixed to Core Consensus.",
      x = "Number of Conditions Sampled",
      y = "Standardized Beta Coefficient"
    ) +
    theme_pub() +
    theme(legend.position = "bottom")
  
  ggsave("figs/downsampling_effect_size_evolution.pdf", p_sim, width = 12, height = 6)
  print(p_sim)
  
} else {
  message("Downsampling simulation produced no results.")
}
} else {
  message("\nDownsampling analysis is DISABLED.")
}
# ===============================================
# 4. Angle Variability Analysis
# ===============================================
message("\n[Angle Analysis] Running slope stability analysis...")

if (!requireNamespace("ggbreak", quietly = TRUE)) {
  message("Warning: 'ggbreak' package is required. Plotting standard histogram.")
  USE_GGBREAK <- FALSE
} else {
  library(ggbreak)
  USE_GGBREAK <- TRUE
}

# --- 1. Prepare Data ---
slopes_ecoli <- extract_meta_stats(ecoli_ds) %>% mutate(assay = "E. coli")
slopes_bsub  <- extract_meta_stats(bsub_ds)  %>% mutate(assay = "B. subtilis")
slopes <- bind_rows(slopes_ecoli, slopes_bsub)

stats_df <- slopes %>%
  filter(is.finite(slope)) %>%
  mutate(angle_deg = atan(slope) * 180 / pi) %>%
  group_by(assay, type, gene) %>%
  summarize(
    n = n(),
    mean_angle = mean(angle_deg),
    sd_angle = sd(angle_deg),
    sign_consistent = all(slope > 0) | all(slope < 0),
    .groups = "drop"
  ) %>%
  filter(n >= 2, is.finite(sd_angle))

# Save Summary
write_csv(
  stats_df %>% group_by(assay) %>% summarise(mean_sd=mean(sd_angle), n_consistent=sum(sign_consistent)),
  file.path("tables", "angle_variability_stats_summary.csv")
)

# --- 2. Data Binning ---
cutoff <- 60
bw <- 2
overflow_pos <- cutoff + 4

plot_data_full <- stats_df %>%
  mutate(
    bin_floor = floor(sd_angle / bw) * bw,
    x_pos = ifelse(sd_angle >= cutoff, overflow_pos, bin_floor + bw/2),
    is_overflow = sd_angle >= cutoff
  ) %>%
  group_by(assay, x_pos, is_overflow, sign_consistent) %>%
  summarise(real_count = n(), .groups = "drop")

# --- 3. Plotting Function ---
create_species_hist <- function(sp_name) {
  
  # Filter Data
  sp_data <- plot_data_full %>% filter(assay == sp_name)
  
  # Calculate Heights
  counts_summary <- sp_data %>%
    group_by(is_overflow) %>%
    summarise(max_val = max(tapply(real_count, x_pos, sum)), .groups="drop")
  
  max_normal <- max(counts_summary$max_val[counts_summary$is_overflow == FALSE], na.rm = TRUE)
  outlier_val <- max(counts_summary$max_val[counts_summary$is_overflow == TRUE], na.rm = TRUE)
  
  if(!is.finite(max_normal)) max_normal <- 0
  if(!is.finite(outlier_val)) outlier_val <- 0
  
  # Define Axis Breaks: 
  # We want standard ticks for the bottom, and ONLY the outlier value for the top.
  # This trick removes all other labels above the break.
  pretty_ticks <- pretty(c(0, max_normal), n = 4)
  custom_breaks <- sort(unique(c(pretty_ticks, outlier_val)))
  
  # Base Plot
  p <- ggplot(sp_data, aes(x = x_pos, y = real_count, fill = sign_consistent)) +
    geom_col(width = bw, position = "stack") +
    
    # Add Dashed Line for the Outlier Bar
    geom_segment(data = subset(sp_data, is_overflow),
                 aes(x = x_pos, xend = 0, y = outlier_val, yend = outlier_val),
                 linetype = "dashed", color = "black", inherit.aes = FALSE) +
    
    scale_x_continuous(
      breaks = c(seq(0, cutoff, by = 10), overflow_pos), 
      labels = c(seq(0, cutoff, by = 10), paste0(">", cutoff))
    ) +
    
    # Use Custom Breaks to force the label at outlier_val
    scale_y_continuous(breaks = custom_breaks, expand = expansion(mult = c(0, 0.05))) +
    
    scale_fill_manual(
      values = c("TRUE" = "#2c7bb6", "FALSE" = "#d7191c"),
      name = "Sign Consistency",
      labels = c("TRUE" = "Stable (Same Sign)", "FALSE" = "Unstable (Sign Switched)")
    ) +
    
    labs(
      title = sp_name,
      y = "Count",
      x = "SD of Angle (Degrees)" # Always show X title
    ) +
    
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(color = "black"),
      # Explicitly remove right axis elements
      axis.line.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.text.y.right = element_blank()
    )
  
  # Apply Break
  if (USE_GGBREAK && (outlier_val > 1.5 * max_normal)) {
    break_start <- max_normal * 1.1
    break_end   <- outlier_val * 0.9
    
    p <- p + 
      scale_y_break(c(break_start, break_end), scales = 1, space = 0.5, ticklabels = NULL) +
      theme(
        # Ensure right axis is gone in the broken plot too
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank()
      )
  }
  
  return(p)
}

# --- 4. Generate & Combine ---
p_ecoli <- create_species_hist("E. coli")
p_bsub  <- create_species_hist("B. subtilis")

# Combine using patchwork
# We stack them, but this time we might want B.subtilis on top or bottom.
# Typically E. coli is top, but user order in previous plots was B. sub top. 
# We'll stick to B. sub top.
final_plot <- p_bsub / p_ecoli + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Stability of Growth-Dependent Expression",
    subtitle = paste0("SD of Angle (arctan(slope)). Axis break applied for >", cutoff, " bin.")
  )

ggsave(
  filename = file.path("figs", "angle_sd_stability_break.pdf"),
  plot = final_plot,
  device = grDevices::cairo_pdf,
  width = 8, height = 10
)

message("Angle stability plot saved to figs/angle_sd_stability_break.pdf")




# ==============================================================================
# PROMOTER STRENGTH ANALYSIS (E. COLI) 
# ==============================================================================
message("\n[Promoter Analysis] Processing Promoter Strength vs Position (Shared/Unique Logic)...")

if (file.exists("LB_promoter.rds") && file.exists("M9_promoter.rds")) {
  
  # -------------------------------------------------------
  # 1. Build Consensus Dataset (Matched to Deleted Section Logic)
  # -------------------------------------------------------
  process_species_consensus_internal <- function(ds_list) {
    extract_meta_stats(ds_list) %>%
      group_by(gene) %>%
      summarise(
        norm_pos = mean(norm_pos, na.rm = TRUE),
        weights_mean_sum = sum(1/se_mean^2, na.rm = TRUE),
        global_mean_ln   = sum(mean_ln * (1/se_mean^2), na.rm = TRUE) / weights_mean_sum,
        #global_mean_ln = mean(mean_ln_z, na.rm = TRUE),
        weights_slope_sum = sum(1/se_slope^2, na.rm = TRUE),
        consensus_slope   = sum(slope * (1/se_slope^2), na.rm = TRUE) / weights_slope_sum,
        .groups = "drop"
      ) %>%
      filter(is.finite(consensus_slope), is.finite(global_mean_ln), is.finite(norm_pos))
  }
  
  ecoli_consensus <- process_species_consensus_internal(ecoli_ds) %>%
    mutate(gene = tolower(trimws(gene)))
  
  # -------------------------------------------------------
  # 2. Process Promoters: Identify Shared vs Unique
  # -------------------------------------------------------
  lb_raw <- readRDS("LB_promoter.rds") %>%
    mutate(gene = tolower(trimws(mapped_gene)), medium = "LB") %>%
    select(gene, start = left_boundary_coordinate, end = right_boundary_coordinate, 
           peak = peak_coordinate, activity = peak_activity, id = name)
  
  m9_raw <- readRDS("M9_promoter.rds") %>%
    mutate(gene = tolower(trimws(mapped_gene)), medium = "M9") %>%
    select(gene, start = left_boundary_coordinate, end = right_boundary_coordinate, 
           peak = peak_coordinate, activity = peak_activity, id = name)
  
  # A. Find Overlapping (Shared) Promoters
  # Join by gene, then filter for coordinate overlap
  shared_raw <- inner_join(lb_raw, m9_raw, by = "gene", suffix = c("_LB", "_M9"), relationship = "many-to-many") %>%
    filter(start_LB <= end_M9 & end_LB >= start_M9)
  
  shared_promoters <- shared_raw %>%
    mutate(
      LB_str = log(activity_LB+1),
      M9_str = log(activity_M9+1),
      Mean_Str = (LB_str + M9_str) / 2,
      # "Use mean peak coordinates for shared promoters"
      Final_Peak = (peak_LB + peak_M9) / 2, 
      type = "Shared"
    ) %>%
    select(gene, LB_str, M9_str, Mean_Str, Final_Peak, type, id_LB, id_M9)
  
  # B. Find LB-Only Promoters (No overlap in M9)
  lb_unique <- lb_raw %>%
    filter(!id %in% shared_raw$id_LB) %>%
    mutate(
      LB_str = log(activity+1),
      M9_str = log(0 + 1), # "assume strength is 0 in the other medium"
      Mean_Str = (LB_str + M9_str) / 2,
      Final_Peak = peak,
      type = "LB_Only"
    ) %>%
    select(gene, LB_str, M9_str, Mean_Str, Final_Peak, type)
  
  # C. Find M9-Only Promoters (No overlap in LB)
  m9_unique <- m9_raw %>%
    filter(!id %in% shared_raw$id_M9) %>%
    mutate(
      LB_str = log(0 +1),
      M9_str = log(activity + 1),
      Mean_Str = (LB_str + M9_str) / 2,
      Final_Peak = peak,
      type = "M9_Only"
    ) %>%
    select(gene, LB_str, M9_str, Mean_Str, Final_Peak, type)
  
  # -------------------------------------------------------
  # 3. Select Best Promoter Per Gene
  # -------------------------------------------------------
  all_candidates <- bind_rows(
    select(shared_promoters, -id_LB, -id_M9),
    lb_unique,
    m9_unique
  )
  
  # "Use the promoter with the highest mean peak strength"
  best_promoters <- all_candidates %>%
    group_by(gene) %>%
    arrange(desc(Mean_Str)) %>%
    slice(1) %>% # Takes top 1 per gene
    ungroup()

  # -------------------------------------------------------
  # 4. Merge with Consensus & Calculate Moving Averages
  # -------------------------------------------------------
  merged_data <- ecoli_consensus %>%
    inner_join(best_promoters, by = "gene") %>%
    mutate(
      Diff_Promoter = LB_str - M9_str
    ) %>%
    arrange(norm_pos)
  
  message("Calculating correlations between promoter strength and position...")
  
  cor_stats_prom <- merged_data %>%
    pivot_longer(
      cols = c(LB_str, M9_str, Mean_Str, Diff_Promoter),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    group_by(Metric) %>%
    summarise(
      Spearman_rho = cor(Value, norm_pos, method = "spearman", use = "complete.obs"),
      p_value = cor.test(Value, norm_pos, method = "spearman", use = "complete.obs")$p.value,
      n_genes = n(),
      .groups = "drop"
    ) 
  
  cor_stats_global <- ecoli_consensus %>%
     pivot_longer(cols = c(global_mean_ln, consensus_slope),names_to = "Metric", values_to = "Value") %>%
        group_by(Metric) %>%
         summarise(
          Spearman_rho = cor(Value, norm_pos, method = "spearman", use ="complete.obs"),
            p_value = cor.test(Value, norm_pos, method = "spearman", use ="complete.obs")$p.value,
           n_genes = n(),
            .groups = "drop"
         )
  
  # Combine and sort
  cor_stats <- bind_rows(cor_stats_prom, cor_stats_global) %>% arrange(desc(abs(Spearman_rho)))
  
  print(as.data.frame(cor_stats))
  write_csv(cor_stats, file.path("tables", "expression_promoter_strength_correlations.csv"))
  
  window_size <- 100

  # -------------------------------------------------------
  # 5. Plotting (Scatter + Dual Axis + Regression)
  # -------------------------------------------------------
  
  # --- PANEL A ---
  min_expr <- min(ecoli_consensus$global_mean_ln)
  max_expr <- max(ecoli_consensus$global_mean_ln)
  
  min_prom <- min(merged_data$Mean_Str)
  max_prom <- max(merged_data$Mean_Str)
  
  q_expr <- quantile(ecoli_consensus$global_mean_ln, probs = c(0.01, 0.99), na.rm = TRUE)
  min_expr <- q_expr[1]
  max_expr <- q_expr[2]
  
  q_prom <- quantile(merged_data$Mean_Str, probs = c(0.01, 0.99), na.rm = TRUE)
  min_prom <- q_prom[1]
  max_prom <- q_prom[2]
  
  scale_factor_1 <- (max_expr - min_expr) / (max_prom - min_prom)
  shift_1 <- min_expr
  
  m_expr <- summary(lm(global_mean_ln ~ norm_pos, data = ecoli_consensus))
  m_prom <- summary(lm(Mean_Str ~ norm_pos, data = merged_data))
  eq_expr <- sprintf("Expr: y = %.2fx %s %.2f (p = %.1e)", m_expr$coefficients[2,1], ifelse(m_expr$coefficients[1,1]<0, "-", "+"), abs(m_expr$coefficients[1,1]), m_expr$coefficients[2,4])
  eq_prom <- sprintf("Prom: y = %.2fx %s %.2f (p = %.1e)", m_prom$coefficients[2,1], ifelse(m_prom$coefficients[1,1]<0, "-", "+"), abs(m_prom$coefficients[1,1]), m_prom$coefficients[2,4])
  
  p1 <- ggplot() +
           geom_point(data = ecoli_consensus, aes(x = norm_pos, y = global_mean_ln, 
                              color = "Mean Expression"), alpha = 0.05, stroke = 0) +
          geom_smooth(data = ecoli_consensus, aes(x = norm_pos, y = global_mean_ln,
                              color = "Mean Expression"), method = "lm", linewidth =1.2, se = FALSE) +
 
    geom_point(data = merged_data, aes(x = norm_pos, y = (Mean_Str - min_prom) * scale_factor_1 
                                       + shift_1, color = "Mean Promoter"), alpha =0.05, stroke = 0) +
        geom_smooth(data = merged_data, aes(x = norm_pos, y = (Mean_Str - min_prom) * scale_factor_1 + shift_1, 
                                    color = "Mean Promoter"), method ="lm", linewidth = 1.2, se = FALSE) +
           scale_y_continuous(name = "Mean Expression (ln)",sec.axis = sec_axis(~ ((. - shift_1) / scale_factor_1) 
                       + min_prom,name = "Mean Promoter Strength")) +
           scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("oriC", "mid","ter")) +
         scale_color_manual(values = c("Mean Expression" = "black", "Mean Promoter" = "#0072B2")) +
     coord_cartesian(ylim = c(min_expr, max_expr)) +    
     annotate("text", x = 0.5, y = max_expr, label = eq_expr, vjust = 1, color = "black", fontface = "bold") +
     annotate("text", x = 0.5, y = max_expr - (max_expr - min_expr)*0.06, label = eq_prom, vjust = 1, color = "#0072B2", fontface = "bold") +
     labs(title = "Panel A: Expression Level & Promoter Strength", x =NULL) +
     theme_pub() +
     theme(legend.position = "top", legend.title = element_blank())
  
  # --- PANEL B ---
  min_slope <- min(ecoli_consensus$consensus_slope)
  max_slope <- max(ecoli_consensus$consensus_slope)
  
  min_diff <- min(merged_data$Diff_Promoter)
  max_diff <- max(merged_data$Diff_Promoter)
  
  q_slope <- quantile(ecoli_consensus$consensus_slope, probs = c(0.01,0.99), na.rm = TRUE)
  min_slope <- q_slope[1]
  max_slope <- q_slope[2]
  
  q_diff <- quantile(merged_data$Diff_Promoter, probs = c(0.01, 0.99),na.rm = TRUE)
  min_diff <- q_diff[1]
  max_diff <- q_diff[2]
 
  scale_factor_diff <- (max_slope - min_slope) / (max_diff - min_diff)
  shift_2 <- min_slope
   
  m_slope <- summary(lm(consensus_slope ~ norm_pos, data = ecoli_consensus))
  m_diff <- summary(lm(Diff_Promoter ~ norm_pos, data = merged_data))
  eq_slope <- sprintf("Slope: y = %.2fx %s %.2f (p = %.1e)", m_slope$coefficients[2,1], ifelse(m_slope$coefficients[1,1]<0, "-", "+"), abs(m_slope$coefficients[1,1]), m_slope$coefficients[2,4])
  eq_diff <- sprintf("Diff: y = %.2fx %s %.2f (p = %.1e)", m_diff$coefficients[2,1], ifelse(m_diff$coefficients[1,1]<0, "-", "+"), abs(m_diff$coefficients[1,1]), m_diff$coefficients[2,4])
  
  p2 <- ggplot() +geom_point(data = ecoli_consensus, aes(x = norm_pos, y = consensus_slope,
                                                         color = "Growth Slope"), alpha = 0.05, stroke = 0) +
       geom_smooth(data = ecoli_consensus, aes(x = norm_pos, y = consensus_slope, color = "Growth Slope"), 
                   method = "lm", linewidth = 1.2, se = FALSE) +

       geom_point(data = merged_data, aes(x = norm_pos, y = (Diff_Promoter - min_diff) * scale_factor_diff + 
                                            shift_2, color = "Diff (LB-M9)"), alpha = 0.05, stroke = 0) +
        geom_smooth(data = merged_data, aes(x = norm_pos, y = (Diff_Promoter - min_diff) * scale_factor_diff + 
                                shift_2, color = "Diff (LB-M9)"), method = "lm", linewidth = 1.2, se = FALSE) +
        scale_y_continuous(name = "Growth Slope (k)",
            sec.axis = sec_axis(~ ((. - shift_2) / scale_factor_diff) + min_diff, name = "Promoter Strength Diff (LB - M9)")
        ) +
        scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("oriC", "mid",  "ter")) +
        scale_color_manual(values = c("Growth Slope" = "black", "Diff (LB-M9)"= "#D55E00")) +    
      coord_cartesian(ylim = c(min_slope, max_slope)) +
         annotate("text", x = 0.5, y = max_slope, label = eq_slope, vjust = 1, color = "black", fontface = "bold") +
         annotate("text", x = 0.5, y = max_slope - (max_slope - min_slope)*0.06, label = eq_diff, vjust = 1, color = "#D55E00", fontface = "bold") +
         labs(title = "Panel B: Growth Slope & Promoter Diff", x = "Normalized Position") +
       theme_pub() +
        theme(legend.position = "top", legend.title = element_blank())
  
  prom_reg_stats <- tibble(
    Panel = c("A", "A", "B", "B"),
    Response = c("Expression", "Promoter Strength", "Growth Slope", "Promoter Diff"),
    Intercept = c(m_expr$coefficients[1,1], m_prom$coefficients[1,1], m_slope$coefficients[1,1], m_diff$coefficients[1,1]),
    Slope = c(m_expr$coefficients[2,1], m_prom$coefficients[2,1], m_slope$coefficients[2,1], m_diff$coefficients[2,1]),
    P_Value = c(m_expr$coefficients[2,4], m_prom$coefficients[2,4], m_slope$coefficients[2,4], m_diff$coefficients[2,4]),
    Equation = c(eq_expr, eq_prom, eq_slope, eq_diff)
  )
  write_csv(prom_reg_stats, file.path("tables", "promoter_strength_regression_stats.csv"))
  
  final_prom_plot <- p1 / p2
  ggsave("figs/promoter_strength_dual_axis_scatter.pdf", final_prom_plot, width = 8, height = 10)
  
  message("Promoter analysis complete. Table saved to tables/promoter_strength_correlations.csv. Plot saved to figs/promoter_strength_dual_axis_shared.pdf")
} else {
  message("Promoter RDS files missing.")
}