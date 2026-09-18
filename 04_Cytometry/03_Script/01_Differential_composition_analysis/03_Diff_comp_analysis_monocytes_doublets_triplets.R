############################# Monocyte singlets vs multiplets ##########################

# Loading the csv file
path_csv = file.path("/home/quatreville","03_BIOINFO_PROJECT", "01_MADBIO_Cameroon", "01_scRNASeq", "86_cytometry_Toulouse", "00_InputData","Mono_singlets_doublets_triplets.csv")
PATH_ANALYSIS_OUTPUT = file.path("/home/quatreville","03_BIOINFO_PROJECT", "01_MADBIO_Cameroon", "01_scRNASeq", "86_cytometry_Toulouse", "05_Output")
cyto_df <- read.csv(path_csv,
  sep = ";",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
                                        
# To add a condition column
cyto_df <- cyto_df %>%
  mutate(
    condition = case_when(
      grepl("^C", Individual) ~ "C",
      grepl("^S", Individual) ~ "S",
      TRUE ~ NA_character_
    ),
    condition = factor(condition, levels = c("C", "S"))
  )
  
# Check
table(cyto_df$Individual, cyto_df$condition)

# ---------------------------------------------------------
# Checks
# ---------------------------------------------------------



# Check number of cell types per individual
table(cyto_df$Individual)


# Check that every individual has the 3 categories
cyto_df %>%
    group_by(Individual) %>%
    summarise(
        n_cell_types = n_distinct(Cell_type),
        cell_types = paste(
            sort(unique(Cell_type)),
            collapse = ", "
        ),
        .groups = "drop"
    )


# ---------------------------------------------------------
# Put Individual and Cell_type in factors
# ---------------------------------------------------------

cyto_sccomp <- cyto_df %>%
    mutate(
        Individual = factor(Individual),
        Cell_type = factor(
            Cell_type,
            levels = c(
                "Mono_singlets",
                "Mono_doublets",
                "Mono_triplet"
            )
        )
    )


# ---------------------------------------------------------
# Check total counts per individual
# ---------------------------------------------------------

cyto_sccomp %>%
    group_by(Individual) %>%
    summarise(
        total_monocytes = sum(Count),
        .groups = "drop"
    )


# ---------------------------------------------------------
# Calculate proportions
# ---------------------------------------------------------

cyto_prop <- cyto_sccomp %>%
    group_by(Individual) %>%
    mutate(
        proportion = Count / sum(Count)
    ) %>%
    ungroup()


# ---------------------------------------------------------
# Summary of proportions
# ---------------------------------------------------------

cyto_prop %>%
    group_by(Cell_type, condition) %>%
    summarise(
        mean_proportion = mean(proportion),
        median_proportion = median(proportion),
        sd_proportion = sd(proportion),
        .groups = "drop"
    )


# ---------------------------------------------------------
# sccomp
# ---------------------------------------------------------

cyto_sccomp_fit <- sccomp_estimate(
    .data = cyto_sccomp,
    formula_composition = ~ condition,
    sample = "Individual",
    cell_group = "Cell_type",
    abundance = "Count"
)


# ---------------------------------------------------------
# Extract effects for condition
# ---------------------------------------------------------

cyto_effects <- cyto_sccomp_fit %>%
    dplyr::filter(
        parameter == "conditionS"
    ) %>%
    dplyr::select(
        Cell_type,
        factor,
        c_lower,
        c_effect,
        c_upper,
        c_rhat,
        c_ess_bulk,
        c_ess_tail
    )

cyto_effects


# ---------------------------------------------------------
# Test significance
# ---------------------------------------------------------

cyto_sccomp_test <- sccomp_test(
    cyto_sccomp_fit
)


cyto_sccomp_test %>%
    dplyr::filter(
        parameter == "conditionS"
    ) %>%
    print(n = 10)


# ---------------------------------------------------------
# Extract results for doublets and triplets
# ---------------------------------------------------------

cyto_sccomp_results <- cyto_sccomp_test %>%
    dplyr::filter(
        parameter == "conditionS",
        Cell_type %in% c(
            "Mono_doublets",
            "Mono_triplet"
        )
    ) %>%
    dplyr::select(
        Cell_type,
        factor,
        c_lower,
        c_effect,
        c_upper,
        c_pH0,
        c_FDR,
        c_rhat,
        c_ess_bulk,
        c_ess_tail
    )


cyto_sccomp_results


# ---------------------------------------------------------
# Save results
# ---------------------------------------------------------

write.csv(
    cyto_sccomp_results,
    file = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_mono_doublets_triplets_results.csv"
    ),
    row.names = FALSE
)


# =========================================================
# Visualization
# =========================================================


# ---------------------------------------------------------
# General visualization
# ---------------------------------------------------------

ggplot(
    cyto_prop,
    aes(
        x = condition,
        y = proportion,
        group = Individual
    )
) +
    geom_boxplot(
        outlier.shape = NA
    ) +
    geom_jitter(
        width = 0.1,
        height = 0
    ) +
    facet_wrap(
        ~ Cell_type,
        scales = "free_y"
    ) +
    theme_bw()


# =========================================================
# DOUBLETS
# =========================================================

doublets_plot <- cyto_prop %>%
    filter(
        Cell_type == "Mono_doublets"
    )


# Get sccomp FDR for doublets
fdr_doublets <- cyto_sccomp_test %>%
    filter(
        Cell_type == "Mono_doublets",
        parameter == "conditionS"
    ) %>%
    pull(c_FDR)


# Colors
cols <- c(
    "C" = "orange",
    "S" = "red"
)


# Dotplot
sccomp_doublets_dotplot <- ggplot(
    doublets_plot,
    aes(
        x = condition,
        y = proportion,
        color = condition
    )
) +

    # One dot = one biological individual
    geom_jitter(
        width = 0.3,
        height = 0,
        size = 5,
        alpha = 0.9
    ) +

    # Group mean
    stat_summary(
        fun = mean,
        geom = "crossbar",
        width = 0.4,
        linewidth = 0.8,
        color = "black",
        fatten = 0,
        alpha = 0.8
    ) +

    # sccomp FDR
    geom_signif(
        comparisons = list(
            c("C", "S")
        ),
        annotations = paste0(
            "FDR = ",
            formatC(
                fdr_doublets,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            doublets_plot$proportion,
            na.rm = TRUE
        ) * 1.08,
        tip_length = 0.01,
        textsize = 7,
        color = "black"
    ) +

    scale_color_manual(
        values = cols
    ) +

    theme_classic(
        base_family = "Arial",
        base_size = 16
    ) +

    labs(
        title = "Mono doublets",
        y = "Proportion",
        x = ""
    ) +

    theme(
        legend.position = "none",

        axis.text.y = element_text(
            size = 18,
            family = "Arial"
        ),

        axis.text.x = element_text(
            size = 18,
            family = "Arial"
        ),

        axis.title.y = element_text(
            size = 20,
            family = "Arial"
        ),

        axis.title.x = element_text(
            size = 20,
            family = "Arial"
        ),

        plot.title = element_text(
            size = 20,
            family = "Arial",
            hjust = 0.5
        )
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(0.05, 0.15)
        )
    )


print(sccomp_doublets_dotplot)


# Save doublets plot
ggsave(
    plot = sccomp_doublets_dotplot,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_mono_doublets_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    device = cairo_pdf,
    dpi = 600
)


# =========================================================
# TRIPLETS
# =========================================================

triplets_plot <- cyto_prop %>%
    filter(
        Cell_type == "Mono_triplet"
    )


# Get sccomp FDR for triplets
fdr_triplets <- cyto_sccomp_test %>%
    filter(
        Cell_type == "Mono_triplet",
        parameter == "conditionS"
    ) %>%
    pull(c_FDR)


# Dotplot
sccomp_triplets_dotplot <- ggplot(
    triplets_plot,
    aes(
        x = condition,
        y = proportion,
        color = condition
    )
) +

    # One dot = one biological individual
    geom_jitter(
        width = 0.3,
        height = 0,
        size = 5,
        alpha = 0.9
    ) +

    # Group mean
    stat_summary(
        fun = mean,
        geom = "crossbar",
        width = 0.4,
        linewidth = 0.8,
        color = "black",
        fatten = 0,
        alpha = 0.8
    ) +

    # sccomp FDR
    geom_signif(
        comparisons = list(
            c("C", "S")
        ),
        annotations = paste0(
            "FDR = ",
            formatC(
                fdr_triplets,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            triplets_plot$proportion,
            na.rm = TRUE
        ) * 1.08,
        tip_length = 0.01,
        textsize = 7,
        color = "black"
    ) +

    scale_color_manual(
        values = cols
    ) +

    theme_classic(
        base_family = "Arial",
        base_size = 16
    ) +

    labs(
        title = "Mono triplets",
        y = "Proportion",
        x = ""
    ) +

    theme(
        legend.position = "none",

        axis.text.y = element_text(
            size = 18,
            family = "Arial"
        ),

        axis.text.x = element_text(
            size = 18,
            family = "Arial"
        ),

        axis.title.y = element_text(
            size = 20,
            family = "Arial"
        ),

        axis.title.x = element_text(
            size = 20,
            family = "Arial"
        ),

        plot.title = element_text(
            size = 20,
            family = "Arial",
            hjust = 0.5
        )
    ) +

    scale_y_continuous(
        expand = expansion(
            mult = c(0.05, 0.15)
        )
    )


print(sccomp_triplets_dotplot)


# Save triplets plot
ggsave(
    plot = sccomp_triplets_dotplot,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_mono_triplets_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    device = cairo_pdf,
    dpi = 600
)









