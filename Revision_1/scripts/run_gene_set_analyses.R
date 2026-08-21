# Five-gene-set summaries, TT attenuation, and core/accessory interaction models.
# Sourced once by final_gene_position_gene_sets.R during the all-gene run.

gene_set_labels <- c(
  all = "All eligible genes",
  core = "PanKB core",
  accessory = "PanKB accessory (including rare)",
  transcription_translation = "RNA polymerase + ribosomal protein",
  other = "Other genes"
)
dataset_order <- names(DATASETS)

gene_set_membership <- function(data, gene_set) {
  if (gene_set == "all") return(rep(TRUE, nrow(data)))
  if (gene_set == "transcription_translation") return(data$is_transcription_translation)
  if (gene_set == "other") return(data$is_other_gene)
  data$binary_core_accessory_main == gene_set
}

analysis_base <- purrr::imap(
  DATASETS,
  ~ .x |>
    dplyr::select(
      gene, norm_pos, type, dataset, species, source_strain,
      reference_identifier, binary_core_accessory_main,
      is_transcription_translation, is_other_gene
    ) |>
    dplyr::left_join(
      ds_gene_stats_list[[.y]] |> dplyr::select(gene, mean_ln, slope),
      by = "gene",
      relationship = "many-to-one"
    )
)

# QC counts for the five gene-set universes, calculated without fitting models.
gene_set_counts <- tidyr::crossing(
  dataset_name = dataset_order,
  gene_set = names(gene_set_labels)
) |>
  purrr::pmap_dfr(function(dataset_name, gene_set) {
    data <- analysis_base[[dataset_name]]
    keep <- gene_set_membership(data, gene_set) %in% TRUE &
      is.finite(data$norm_pos) & is.finite(data$mean_ln) & is.finite(data$slope)
    retained <- data[keep, , drop = FALSE]

    tibble::tibble(
      dataset = dataset_name,
      species = dplyr::first(data$species),
      gene_set = gene_set_labels[[gene_set]],
      mapping_confidence_version = "main",
      row_count = nrow(retained),
      unique_gene_count = dplyr::n_distinct(retained$gene)
    )
  }) |>
  dplyr::arrange(mapping_confidence_version, gene_set, match(dataset, dataset_order))

# The only retained supp_* outputs are two multipage PDFs, one page per gene set.
figure_data <- purrr::imap_dfr(
  analysis_base,
  function(data, dataset_name) {
    purrr::imap_dfr(
      gene_set_labels,
      function(label, gene_set) {
        data |>
          dplyr::filter(gene_set_membership(data, gene_set) %in% TRUE) |>
          dplyr::transmute(
            dataset_label = unname(name_map[dataset_name]),
            gene_set = label,
            norm_pos,
            lnE = mean_ln,
            k = slope
          )
      }
    )
  }
) |>
  dplyr::mutate(
    dataset_label = factor(dataset_label, levels = unname(name_map[dataset_order])),
    gene_set = factor(gene_set, levels = unname(gene_set_labels))
  )

supp_position_plot <- function(response, gene_set_label) {
  plot_data <- figure_data |>
    dplyr::filter(
      gene_set == gene_set_label,
      is.finite(norm_pos),
      is.finite(.data[[response]])
    )
  panel_stats <- plot_data |>
    dplyr::group_by(dataset_label) |>
    dplyr::group_modify(~ {
      fit <- stats::lm(
        stats::reformulate("norm_pos", response = response),
        data = .x
      )
      coefficients <- summary(fit)$coefficients
      slope <- coefficients["norm_pos", "Estimate"]
      intercept <- coefficients["(Intercept)", "Estimate"]
      y_limits <- stats::quantile(
        .x[[response]],
        probs = c(0.01, 0.99),
        names = FALSE
      )
      tibble::tibble(
        slope,
        intercept,
        y_min = y_limits[[1]],
        y_max = y_limits[[2]],
        formula_annotation = sprintf(
          "paste(italic(r), ' = ', %s, ', ', italic(P), ' = ', %s)",
          formatC(stats::cor(.x$norm_pos, .x[[response]], method = "pearson"), format = "g", digits = 3),
          formatC(coefficients["norm_pos", "Pr(>|t|)"], format = "g", digits = 3)
        ),
        n_annotation = sprintf("paste(italic(n), ' = ', %d, '.')", nrow(.x))
      )
    }) |>
    dplyr::ungroup()
  display_data <- plot_data |>
    dplyr::left_join(
      panel_stats |>
        dplyr::select(dataset_label, y_min, y_max),
      by = "dataset_label"
    ) |>
    dplyr::filter(
      .data[[response]] >= y_min,
      .data[[response]] <= y_max
    )
  panel_ranges <- panel_stats |>
    dplyr::select(dataset_label, y_min, y_max) |>
    tidyr::pivot_longer(
      cols = c(y_min, y_max),
      names_to = NULL,
      values_to = "y_limit"
    )

  ggplot2::ggplot(display_data, ggplot2::aes(norm_pos, .data[[response]])) +
    ggplot2::geom_point(alpha = 0.05, stroke = 0) +
    ggplot2::geom_blank(
      data = panel_ranges,
      ggplot2::aes(x = 0.5, y = y_limit),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_abline(
      data = panel_stats,
      ggplot2::aes(slope = slope, intercept = intercept),
      inherit.aes = FALSE, linewidth = 0.8, colour = "black"
    ) +
    ggplot2::geom_text(
      data = panel_stats,
      ggplot2::aes(x = -Inf, y = -Inf, label = formula_annotation),
      inherit.aes = FALSE, parse = TRUE,
      hjust = -0.05, vjust = -1.8, size = 4.5
    ) +
    ggplot2::geom_text(
      data = panel_stats,
      ggplot2::aes(x = -Inf, y = -Inf, label = n_annotation),
      inherit.aes = FALSE, parse = TRUE,
      hjust = -0.05, vjust = -0.3, size = 4.5
    ) +
    ggplot2::facet_wrap(~dataset_label, ncol = 3, scales = "free_y", drop = FALSE) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("oriC", "mid", "ter")
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title = gene_set_label,
      subtitle = "Points are genes; black lines are descriptive ordinary least-squares fits",
      x = "Normalized gene position",
      y = response
    ) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
      strip.text = ggplot2::element_text(face = "bold", size = 14),
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      plot.title = ggplot2::element_text(size = 20, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 15),
      panel.spacing = grid::unit(1, "lines")
    )
}

for (response in c("lnE", "k")) {
  grDevices::pdf(
    file.path("figs", paste0("supp_", response, "_by_position_all_gene_sets.pdf")),
    width = 14, height = 8, onefile = TRUE
  )
  for (gene_set_label in unname(gene_set_labels)) {
    print(supp_position_plot(response, gene_set_label))
  }
  grDevices::dev.off()
}

# Inverse-variance consensus estimates used by TT and interaction analyses.
mode_value <- function(x) {
  if (length(x)) names(sort(table(x), decreasing = TRUE))[1] else NA_character_
}

make_consensus_data <- function(
  ds_names,
  dataset_label,
  species_label,
  min_datasets = 1L
) {
  annotations <- purrr::map_dfr(
    DATASETS[ds_names],
    ~ .x |>
      dplyr::transmute(
        gene = tolower(trimws(gene)),
        binary_core_accessory_main,
        is_transcription_translation,
        is_other_gene
      )
  ) |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      binary_core_accessory_main = mode_value(binary_core_accessory_main[!is.na(binary_core_accessory_main)]),
      is_transcription_translation = any(is_transcription_translation %in% TRUE),
      is_other_gene = any(is_other_gene %in% TRUE),
      .groups = "drop"
    )

  extract_meta_stats(ds_names) |>
    dplyr::group_by(gene) |>
    dplyr::filter(dplyr::n_distinct(dataset) >= min_datasets) |>
    dplyr::summarise(
      norm_pos = mean(norm_pos),
      mean_ln = weighted.mean(mean_ln, 1 / se_mean^2),
      slope = weighted.mean(slope, 1 / se_slope^2),
      .groups = "drop"
    ) |>
    dplyr::left_join(annotations, by = "gene") |>
    dplyr::mutate(
      reference_identifier = gene,
      dataset = dataset_label,
      species = species_label
    )
}

consensus_analysis_base <- list(
  ecoli_consensus = make_consensus_data(
    ecoli_ds, "E. coli consensus", "E. coli", min_datasets = 3L
  ),
  bsub_consensus = make_consensus_data(
    bsub_ds, "B. subtilis consensus", "B. subtilis", min_datasets = 2L
  ),
  vnat_consensus = make_consensus_data(
    vnat_ds, "Zhu (V. natriegens)", "V. natriegens"
  )
)

# Interaction models intentionally retain the union of genes observed in any
# dataset, so they are the sole consensus-analysis exception to the thresholds.
interaction_consensus_base <- list(
  ecoli_consensus = make_consensus_data(
    ecoli_ds, "E. coli consensus", "E. coli"
  ),
  bsub_consensus = make_consensus_data(
    bsub_ds, "B. subtilis consensus", "B. subtilis"
  )
)

# Raw TT attenuation and equal-size-removal permutation null.
run_tt_attenuation <- function(
  data_list,
  B = 10000,
  predictor_columns = c(z_slope = "slope", z_mean = "mean_ln")
) {
  summaries <- list()
  nulls <- list()

  for (dataset_name in names(data_list)) {
    dat <- data_list[[dataset_name]] |>
      dplyr::filter(
        is.finite(norm_pos),
        dplyr::if_all(
          dplyr::all_of(unname(predictor_columns)),
          is.finite
        )
      ) |>
      dplyr::mutate(
        is_tt = is_transcription_translation %in% TRUE
      )
    for (predictor in names(predictor_columns)) {
      dat[[predictor]] <- safe_z(dat[[predictor_columns[[predictor]]]])
    }
    n_total <- nrow(dat)
    n_tt <- sum(dat$is_tt)
    if (n_tt < 5 || n_total - n_tt < 10) next

    for (predictor in names(predictor_columns)) {
      x <- dat[[predictor]]
      y <- dat$norm_pos
      beta_all <- cov(x, y) / var(x)
      beta_nonTT <- cov(x[!dat$is_tt], y[!dat$is_tt]) / var(x[!dat$is_tt])
      delta_TT <- beta_all - beta_nonTT

      set.seed(947)
      delta_random <- replicate(B, {
        excluded <- sample.int(n_total, n_tt)
        beta_all - cov(x[-excluded], y[-excluded]) / var(x[-excluded])
      })
      Z_TT <- (delta_TT - mean(delta_random)) / stats::sd(delta_random)
      p_lower <- (sum(delta_random <= delta_TT) + 1) / (B + 1)
      p_upper <- (sum(delta_random >= delta_TT) + 1) / (B + 1)

      summaries[[length(summaries) + 1]] <- tibble::tibble(
        dataset = dataset_name,
        predictor,
        n_total,
        n_tt,
        beta_all,
        beta_nonTT,
        delta_TT,
        Z_TT,
        p_emp_lower = p_lower,
        p_emp_upper = p_upper,
        #p_emp = ifelse(delta_TT < 0, p_lower, p_upper)
        p_emp = p_lower
      )
      nulls[[length(nulls) + 1]] <- tibble::tibble(
        dataset = dataset_name, predictor, delta_random
      )
    }
  }
  list(summary = dplyr::bind_rows(summaries), null = dplyr::bind_rows(nulls))
}

tt_predictor_labels <- c(
  z_slope = "Growth rate slope (k)",
  z_mean = "Mean expression (lnE)",
  z_mean_promoter_strength = "Mean promoter strength"
)

plot_tt_attenuation <- function(
  result,
  title_suffix,
  filename,
  display_order,
  bottom_result = NULL,
  bottom_display_order = NULL
) {
  label_map <- c(
    name_map,
    ecoli_consensus = "E. coli Consensus",
    bsub_consensus = "B. subtilis Consensus",
    vnat_consensus = "V. natriegens"
  )

  build_panel <- function(panel_result, panel_display_order) {
    panel_display_order <- sub(" dataset$", "", panel_display_order)
    summary_data <- panel_result$summary |>
      dplyr::mutate(
        dataset_label = factor(
          sub(" dataset$", "", unname(label_map[dataset])),
          levels = panel_display_order
        ),
        predictor_label = factor(
          unname(tt_predictor_labels[predictor]),
          levels = unname(tt_predictor_labels)
        ),
        annotation = sprintf(
          "delta_TT = %.4f\nZ_TT = %.2f\nP = %.4f",
          delta_TT, Z_TT, p_emp
        )
      )
    null_data <- panel_result$null |>
      dplyr::mutate(
        dataset_label = factor(
          sub(" dataset$", "", unname(label_map[dataset])),
          levels = panel_display_order
        ),
        predictor_label = factor(
          unname(tt_predictor_labels[predictor]),
          levels = unname(tt_predictor_labels)
        )
      )

    ggplot2::ggplot(null_data, ggplot2::aes(delta_random)) +
      ggplot2::geom_histogram(bins = 50, fill = "#77BBDD", color = "grey30", alpha = 0.6) +
      ggplot2::geom_vline(
        data = summary_data, ggplot2::aes(xintercept = delta_TT),
        color = "#FF8899", linetype = "dashed", linewidth = 0.9
      ) +
      ggplot2::geom_text(
        data = summary_data,
        ggplot2::aes(x = Inf, y = Inf, label = annotation),
        vjust = 1.3, hjust = 1.1, size = 5, color = "black", fontface = "bold"
      ) +
      ggplot2::facet_grid(dataset_label ~ predictor_label, scales = "free") +
      ggplot2::theme_classic(base_size = 16) +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
        strip.text = ggplot2::element_text(face = "bold", size = 14),
        axis.text = ggplot2::element_text(size = 14),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        plot.title = ggplot2::element_text(size = 20, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 15),
        panel.border = ggplot2::element_rect(color = "grey70", fill = NA)
      ) +
      ggplot2::labs(
        x = "Change in positional regression coefficient after gene removal",
        y = "Permutation count"
      )
  }

  main_plot <- build_panel(result, display_order)
  if (is.null(bottom_result)) {
    plot <- main_plot +
      ggplot2::labs(
        title = paste0("Transcription/Translation (TT) Attenuation Permutation Test (", title_suffix, ")"),
        subtitle = "Observed raw delta_TT vs. delta_random from 10,000 equal-size gene removals"
      )
    height <- ifelse(length(display_order) > 5, 14, 9)
  } else {
    bottom_plot <- build_panel(bottom_result, bottom_display_order)
    bottom_row <- bottom_plot +
      patchwork::plot_spacer() +
      patchwork::plot_layout(widths = c(1, 1))
    main_plot <- main_plot +
      ggplot2::labs(
        title = paste0("Transcription/Translation (TT) Attenuation Permutation Test (", title_suffix, ")"),
        subtitle = "Observed raw delta_TT vs. delta_random from 10,000 equal-size gene removals"
      )
    plot <- main_plot / bottom_row +
      patchwork::plot_layout(
        heights = c(length(display_order), length(bottom_display_order))
      )
    height <- ifelse(length(display_order) > 5, 17, 12)
  }

  ggplot2::ggsave(file.path("figs", paste0(filename, ".pdf")), plot, width = 14, height = height)
  ggplot2::ggsave(file.path("figs", paste0(filename, ".png")), plot, width = 14, height = height, dpi = 300)
}

tt_individual <- run_tt_attenuation(analysis_base)
tt_consensus <- run_tt_attenuation(consensus_analysis_base)
plot_tt_attenuation(
  tt_individual, "Individual Datasets", "tt_attenuation_individual",
  unname(name_map[dataset_order])
)

# E. coli promoter summaries for the consensus interaction models.
read_promoters <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::filter(!is.na(reference_identifier)) |>
    dplyr::transmute(
      gene = tolower(trimws(reference_identifier)),
      start = left_boundary_coordinate,
      end = right_boundary_coordinate,
      peak = peak_coordinate,
      activity = peak_activity,
      id = name
    )
}

lb <- read_promoters(file.path(REVISION_ROOT, "data", "promoter", "mapped", "LB_promoter_identifiers.csv"))
m9 <- read_promoters(file.path(REVISION_ROOT, "data", "promoter", "mapped", "M9_promoter_identifiers.csv"))
shared <- dplyr::inner_join(lb, m9, by = "gene", suffix = c("_LB", "_M9"), relationship = "many-to-many") |>
  dplyr::filter(start_LB <= end_M9, end_LB >= start_M9)

promoter_candidates <- dplyr::bind_rows(
  shared |>
    dplyr::transmute(
      gene,
      LB_str = log(activity_LB + 1),
      M9_str = log(activity_M9 + 1)
    ),
  lb |>
    dplyr::filter(!id %in% shared$id_LB) |>
    dplyr::transmute(gene, LB_str = log(activity + 1), M9_str = 0),
  m9 |>
    dplyr::filter(!id %in% shared$id_M9) |>
    dplyr::transmute(gene, LB_str = 0, M9_str = log(activity + 1))
) |>
  dplyr::mutate(
    Mean_Str = (LB_str + M9_str) / 2,
    # Existing definition of growth dependence in promoter strength:
    # positive values indicate stronger promoter activity in LB than in M9.
    Diff_Promoter = LB_str - M9_str
  ) |>
  dplyr::slice_max(Mean_Str, n = 1, by = gene, with_ties = FALSE)

# Reference position and gene-class annotations for promoter interaction
# models. These are loaded independently of the expression datasets.
ecoli_pankb_annotations <- readr::read_csv(
  file.path(REVISION_ROOT, "data", "raw_data", "PanKB", "Ecoli_focal_gene_info.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(!is.na(original_locus_tag), original_locus_tag != "") |>
  dplyr::inner_join(
    readr::read_csv(
      file.path(REVISION_ROOT, "data", "raw_data", "PanKB", "Ecoli_gene_annotations.csv"),
      show_col_types = FALSE
    ),
    by = "gene"
  ) |>
  dplyr::arrange(dplyr::desc(original_exact_match)) |>
  dplyr::transmute(
    gene = tolower(trimws(original_locus_tag)),
    binary_core_accessory_main = dplyr::if_else(
      tolower(pangenomic_class) == "core", "core", "accessory"
    )
  ) |>
  dplyr::distinct(gene, .keep_all = TRUE)

ecoli_tt_genes <- readr::read_csv(
  file.path(REVISION_ROOT, "data", "reference_transcription_translation_genes.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(species == "Escherichia coli") |>
  dplyr::transmute(
    gene = tolower(trimws(reference_identifier)),
    is_transcription_translation = TRUE
  )

ecoli_interaction_annotations <- readr::read_csv(
  file.path(REVISION_ROOT, "data", "reference_gene_positions.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(species == "Escherichia coli", is.finite(norm_pos)) |>
  dplyr::transmute(
    gene = tolower(trimws(reference_identifier)),
    reference_identifier = gene,
    norm_pos
  ) |>
  dplyr::inner_join(ecoli_pankb_annotations, by = "gene") |>
  dplyr::left_join(ecoli_tt_genes, by = "gene") |>
  dplyr::mutate(
    is_transcription_translation = dplyr::coalesce(
      is_transcription_translation, FALSE
    )
  )

tt_promoter_consensus <- run_tt_attenuation(
  list(
    ecoli_consensus = dplyr::inner_join(
      consensus_analysis_base$ecoli_consensus,
      promoter_candidates,
      by = "gene"
    )
  ),
  predictor_columns = c(z_mean_promoter_strength = "Mean_Str")
)

plot_tt_attenuation(
  tt_consensus,
  "Consensus Datasets",
  "tt_attenuation_consensus",
  c("E. coli Consensus", "B. subtilis Consensus", "V. natriegens")
)

attach_interaction_annotations <- function(dataset_name) {
  extract_meta_stats(dataset_name) |>
    dplyr::select(gene, norm_pos, mean_ln, slope) |>
    dplyr::left_join(
      DATASETS[[dataset_name]] |>
        dplyr::transmute(
          gene = tolower(trimws(gene)),
          reference_identifier,
          binary_core_accessory_main,
          is_transcription_translation
        ),
      by = "gene",
      relationship = "many-to-one"
    )
}

interaction_dataset_labels <- c(
  schmidt_2016 = "Schmidt (E. coli)",
  peebo_2015 = "Peebo (E. coli)",
  li_2014_all = "Li (E. coli)",
  valgepea_2013_prot = "Valgepea (E. coli)",
  zhu_E_coli = "Zhu (E. coli)",
  Goelzer_2015 = "Goelzer (B. subtilis)",
  zhu_B_subtilis = "Zhu (B. subtilis)",
  zhu_V_natriegens = "Zhu (V. natriegens)"
)

# Expression outcomes use the unfiltered expression-consensus union. Promoter
# outcomes use every promoter gene with the annotations needed by the model.
ecoli_interaction_data <- dplyr::full_join(
  interaction_consensus_base$ecoli_consensus |>
    dplyr::select(gene, mean_ln, slope),
  promoter_candidates,
  by = "gene"
) |>
  dplyr::inner_join(ecoli_interaction_annotations, by = "gene")

interaction_datasets <- list(
  `E. coli consensus` = list(
    data = ecoli_interaction_data,
    species = "E. coli", type = "Consensus",
    outcomes = c("Mean_Str", "Diff_Promoter", "mean_ln", "slope")
  ),
  `B. subtilis consensus` = list(
    data = interaction_consensus_base$bsub_consensus,
    species = "B. subtilis", type = "Consensus",
    outcomes = c("mean_ln", "slope")
  )
)

for (dataset_name in names(interaction_dataset_labels)) {
  label <- interaction_dataset_labels[[dataset_name]]
  interaction_datasets[[label]] <- list(
    data = attach_interaction_annotations(dataset_name),
    species = get_species_from_ds(dataset_name),
    type = ifelse(dataset_name == "zhu_V_natriegens", "Consensus_Representative", "Individual"),
    outcomes = c("mean_ln", "slope")
  )
}

outcome_labels <- c(
  Mean_Str = "Mean promoter strength",
  Diff_Promoter = "Growth dependence in promoter strength",
  mean_ln = "Mean expression",
  slope = "Growth slope"
)

# OLS interaction: outcome ~ norm_pos * gene class, with core as reference.
fit_interaction_model <- function(data, outcome, dataset, species, dataset_type, exclude_tt) {
  dat <- data |>
    dplyr::filter(
      binary_core_accessory_main %in% c("core", "accessory"),
      is.finite(norm_pos),
      is.finite(.data[[outcome]])
    )
  if (exclude_tt) dat <- dat |> dplyr::filter(!is_transcription_translation)
  dat <- dat |>
    dplyr::mutate(
      binary_core_accessory_main =
        relevel(factor(binary_core_accessory_main), ref = "core")
    )

  fit <- stats::lm(
    stats::reformulate("norm_pos * binary_core_accessory_main", response = outcome),
    data = dat
  )
  beta <- stats::coef(fit)
  variance <- stats::vcov(fit)
  residual_df <- stats::df.residual(fit)
  t_critical <- stats::qt(0.975, residual_df)
  interaction_term <- "norm_pos:binary_core_accessory_mainaccessory"

  core_estimate <- beta[["norm_pos"]]
  core_se <- sqrt(variance["norm_pos", "norm_pos"])
  accessory_estimate <- core_estimate + beta[[interaction_term]]
  accessory_se <- sqrt(
    variance["norm_pos", "norm_pos"] + variance[interaction_term, interaction_term] +
      2 * variance["norm_pos", interaction_term]
  )
  interaction_estimate <- beta[[interaction_term]]
  interaction_se <- sqrt(variance[interaction_term, interaction_term])

  slope_table <- tibble::tibble(
    dataset, species, dataset_type,
    outcome_raw = outcome,
    outcome_label = outcome_labels[[outcome]],
    gene_class = c("Core", "Accessory"),
    filter_mode = ifelse(exclude_tt, "Excl. Trans/Transl (Sensitivity)", "Full (All Genes)"),
    estimate = c(core_estimate, accessory_estimate),
    std_error = c(core_se, accessory_se),
    residual_df,
    t_statistic = estimate / std_error,
    p_value = 2 * stats::pt(-abs(t_statistic), residual_df),
    confidence_interval_low = estimate - t_critical * std_error,
    confidence_interval_high = estimate + t_critical * std_error,
    n_obs = nrow(dat)
  )
  interaction_table <- tibble::tibble(
    dataset, species, dataset_type,
    outcome_raw = outcome,
    outcome_label = outcome_labels[[outcome]],
    filter_mode = ifelse(exclude_tt, "Excl. Trans/Transl (Sensitivity)", "Full (All Genes)"),
    term = interaction_term,
    estimate = interaction_estimate,
    std_error = interaction_se,
    residual_df,
    t_statistic = estimate / std_error,
    p_value = 2 * stats::pt(-abs(t_statistic), residual_df),
    confidence_interval_low = estimate - t_critical * std_error,
    confidence_interval_high = estimate + t_critical * std_error,
    n_obs = nrow(dat)
  )
  list(slopes = slope_table, interaction = interaction_table)
}

interaction_fits <- list()
for (dataset in names(interaction_datasets)) {
  info <- interaction_datasets[[dataset]]
  for (outcome in info$outcomes) {
    for (exclude_tt in c(FALSE, TRUE)) {
      interaction_fits[[length(interaction_fits) + 1]] <- fit_interaction_model(
        info$data, outcome, dataset, info$species, info$type, exclude_tt
      )
    }
  }
}

interaction_slopes <- purrr::map_dfr(interaction_fits, "slopes") |>
  dplyr::mutate(p_value_fdr = p.adjust(p_value, method = "fdr"))
interaction_effects <- purrr::map_dfr(interaction_fits, "interaction") |>
  dplyr::mutate(p_value_fdr = p.adjust(p_value, method = "fdr"))

interaction_theme <- function() {
  ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black", size = 14),
      axis.title = ggplot2::element_text(face = "bold", size = 16),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      strip.text = ggplot2::element_text(face = "bold", size = 14),
      plot.title = ggplot2::element_text(face = "bold", size = 20),
      plot.caption = ggplot2::element_text(size = 13, hjust = 0),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 14),
      legend.position = "bottom",
      panel.spacing = grid::unit(1.2, "lines")
    )
}

interaction_plot_data <- interaction_effects |>
  dplyr::mutate(
    significance = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    outcome_label = factor(outcome_label, levels = unname(outcome_labels)),
    model_mode = factor(
      ifelse(grepl("Excl", filter_mode), "Excl. TT Genes", "All Genes"),
      levels = c("All Genes", "Excl. TT Genes")
    )
  )

interaction_panel <- function(data, title, caption = NULL, facet_ncol = 2) {
  ggplot2::ggplot(data, ggplot2::aes(estimate, dataset, color = model_mode)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = confidence_interval_low, xmax = confidence_interval_high),
      orientation = "y", width = 0.35,
      position = ggplot2::position_dodge(width = 0.65)
    ) +
    ggplot2::geom_point(size = 2.8, position = ggplot2::position_dodge(width = 0.65)) +
    ggplot2::geom_text(
      ggplot2::aes(label = significance),
      position = ggplot2::position_dodge(width = 0.65),
      vjust = -0.6, size = 5.5, show.legend = FALSE, fontface = "bold"
    ) +
    ggplot2::facet_wrap(~outcome_label, scales = "free_x", ncol = facet_ncol) +
    ggplot2::scale_color_manual(values = c("All Genes" = "#77BBDD", "Excl. TT Genes" = "#FF5522"), name = "Model") +
    ggplot2::labs(
      title = title,
      x = "Position x gene-class interaction coefficient (95% finite-sample t CI)",
      y = NULL,
      caption = caption
    ) +
    interaction_theme()
}

caption <- paste0(
  "Interaction coefficient = accessory positional slope - core positional slope.",
  "\n",
  "Error bars and stars use finite-sample t inference with the model residual degrees of freedom:",
  "\n",
  "* p < 0.05, ** p < 0.01, *** p < 0.001."
)
consensus_levels <- c("E. coli consensus", "B. subtilis consensus", "Zhu (V. natriegens)")
individual_levels <- unname(interaction_dataset_labels[names(interaction_dataset_labels) != "zhu_V_natriegens"])
consensus_plot_outcomes <- c("mean_ln", "slope", "Mean_Str", "Diff_Promoter")
consensus_plot_outcome_labels <- unname(outcome_labels[consensus_plot_outcomes])

consensus_plot_data <- interaction_plot_data |>
  dplyr::filter(dataset %in% consensus_levels) |>
  dplyr::mutate(dataset = factor(dataset, levels = rev(consensus_levels)))
consensus_plot <- interaction_panel(
  consensus_plot_data |>
    dplyr::filter(outcome_raw %in% consensus_plot_outcomes) |>
    dplyr::mutate(
      outcome_label = factor(
        outcome_label,
        levels = consensus_plot_outcome_labels
      )
    ),
  "Consensus Datasets: Position x Gene-Class Interaction Effects",
  caption,
  facet_ncol = 4
)

individual_plot_data <- interaction_plot_data |>
  dplyr::filter(dataset %in% individual_levels) |>
  dplyr::mutate(dataset = factor(dataset, levels = rev(individual_levels)))
individual_plot <- interaction_panel(
  individual_plot_data,
  "Individual Datasets: Position x Gene-Class Interaction Effects",
  caption
)

ggplot2::ggsave("figs/fig_interaction_coefficients_forest_consensus.pdf", consensus_plot, width = 20, height = 7)
ggplot2::ggsave("figs/fig_interaction_coefficients_forest_consensus.png", consensus_plot, width = 20, height = 7, dpi = 300)
ggplot2::ggsave("figs/fig_interaction_coefficients_forest_individual.pdf", individual_plot, width = 13, height = 10)
ggplot2::ggsave("figs/fig_interaction_coefficients_forest_individual.png", individual_plot, width = 13, height = 10, dpi = 300)

readr::write_csv(gene_set_counts, "tables/gene_set_dataset_counts.csv")
readr::write_csv(tt_individual$summary, "tables/tt_attenuation_individual.csv")
readr::write_csv(
  dplyr::bind_rows(tt_consensus$summary, tt_promoter_consensus$summary),
  "tables/tt_attenuation_consensus.csv"
)
readr::write_csv(interaction_slopes, "tables/interaction_models_group_slopes_summary.csv")
readr::write_csv(interaction_effects, "tables/interaction_models_interaction_term_summary.csv")

message("Completed five-gene-set summaries, TT attenuation, and interaction models.")
