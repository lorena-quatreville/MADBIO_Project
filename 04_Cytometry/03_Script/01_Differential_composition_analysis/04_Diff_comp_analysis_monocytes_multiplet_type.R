############################# Monocyte singlets vs multiplets ##########################

# Loading the csv file
path_csv = file.path("/home/quatreville","03_BIOINFO_PROJECT", "01_MADBIO_Cameroon", "01_scRNASeq", "86_cytometry_Toulouse", "00_InputData","Mono_multiplet_type.csv")
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


# Check that every individual has the 9 categories
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


# Check all category names
sort(unique(cyto_df$Cell_type))


# ---------------------------------------------------------
# Put Individual and Cell_type in factors
# ---------------------------------------------------------

cyto_sccomp <- cyto_df %>%
    mutate(
        Individual = factor(Individual),

        Cell_type = factor(
            Cell_type,
            levels = c(
                "NCmono_Bcell",
                "NCmono_Tcell",
                "NCmono_Tcell_Bcell",
                "Imono_Bcell",
                "Imono_Tcell",
                "Imono_Tcell_Bcell",
                "Cmono_Bcell",
                "Cmono_Tcell",
                "Cmono_Tcell_Bcell"
            )
        )
    )


# ---------------------------------------------------------
# Check Cell_type factor
# ---------------------------------------------------------

table(
    cyto_sccomp$Cell_type,
    useNA = "ifany"
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


# =========================================================
# sccomp
# =========================================================

cyto_sccomp_fit <- sccomp_estimate(
    .data = cyto_sccomp,
    formula_composition = ~ condition,
    sample = "Individual",
    cell_group = "Cell_type",
    abundance = "Count"
)


# ---------------------------------------------------------
# Extract effects for condition S vs C
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


# =========================================================
# Test significance
# =========================================================

cyto_sccomp_test <- sccomp_test(
    cyto_sccomp_fit
)


cyto_sccomp_test %>%
    dplyr::filter(
        parameter == "conditionS"
    ) %>%
    print(n = 20)


# ---------------------------------------------------------
# Extract result specifically for NCmono_Tcell
# ---------------------------------------------------------

cyto_sccomp_multiplet_type <- cyto_sccomp_test %>%
    dplyr::filter(
        parameter == "conditionS"
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




# ---------------------------------------------------------
# Save result
# ---------------------------------------------------------

write.csv(
    cyto_sccomp_multiplet_type,
    file = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_multiplet_type_results.csv"
    ),
    row.names = FALSE
)


# =========================================================
# Visualization
# =========================================================




# =========================================================
# NCmono_Tcell
# =========================================================

NCmono_Tcell_plot <- cyto_prop %>%
    filter(
        Cell_type == "NCmono_Tcell"
    )


# ---------------------------------------------------------
# Get sccomp FDR
# ---------------------------------------------------------

fdr_NCmono_Tcell <- cyto_sccomp_test %>%
    filter(
        Cell_type == "NCmono_Tcell",
        parameter == "conditionS"
    ) %>%
    pull(c_FDR)


# Check FDR
fdr_NCmono_Tcell


# ---------------------------------------------------------
# Colors
# ---------------------------------------------------------

cols <- c(
    "C" = "orange",
    "S" = "red"
)


# =========================================================
# Dotplot
# =========================================================

sccomp_NCmono_Tcell_dotplot <- ggplot(
    NCmono_Tcell_plot,
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
            "p = ",
            formatC(
                fdr_NCmono_Tcell,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            NCmono_Tcell_plot$proportion,
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
        title = "NCmono T-cell multiplets",
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


print(sccomp_NCmono_Tcell_dotplot)


# ---------------------------------------------------------
# Save plot
# ---------------------------------------------------------

ggsave(
    plot = sccomp_NCmono_Tcell_dotplot,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_NCmono_Tcell_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    device = cairo_pdf,
    dpi = 600
)






# =========================================================
# Cmono_Bcell
# =========================================================

Cmono_Bcell_plot <- cyto_prop %>%
    filter(
        Cell_type == "Cmono_Bcell"
    )


# ---------------------------------------------------------
# Get sccomp FDR
# ---------------------------------------------------------

fdr_Cmono_Bcell <- cyto_sccomp_test %>%
    filter(
        Cell_type == "Cmono_Bcell",
        parameter == "conditionS"
    ) %>%
    pull(c_FDR)


# Check FDR
fdr_Cmono_Bcell


# ---------------------------------------------------------
# Colors
# ---------------------------------------------------------

cols <- c(
    "C" = "orange",
    "S" = "red"
)


# =========================================================
# Dotplot
# =========================================================

sccomp_Cmono_Bcell_dotplot <- ggplot(
   Cmono_Bcell_plot,
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
            "p = ",
            formatC(
                fdr_Cmono_Bcell,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            Cmono_Bcell_plot$proportion,
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
        title = "Cmono B-cell multiplets",
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


print(sccomp_Cmono_Bcell_dotplot)


# ---------------------------------------------------------
# Save plot
# ---------------------------------------------------------

ggsave(
    plot = sccomp_Cmono_Bcell_dotplot,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_Cmono_Bcell_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    device = cairo_pdf,
    dpi = 600
)










# =========================================================
# NCmono_Tcell_Bcell
# =========================================================

NCmono_Tcell_Bcell_plot <- cyto_prop %>%
    filter(
        Cell_type == "NCmono_Tcell_Bcell"
    )


# ---------------------------------------------------------
# Get sccomp FDR
# ---------------------------------------------------------

fdr_NCmono_Tcell_Bcell <- cyto_sccomp_test %>%
    filter(
        Cell_type == "NCmono_Tcell_Bcell",
        parameter == "conditionS"
    ) %>%
    pull(c_FDR)


# Check FDR
fdr_NCmono_Tcell_Bcell


# ---------------------------------------------------------
# Colors
# ---------------------------------------------------------

cols <- c(
    "C" = "orange",
    "S" = "red"
)


# =========================================================
# Dotplot
# =========================================================

sccomp_NCmono_Tcell_Bcell_dotplot <- ggplot(
    NCmono_Tcell_Bcell_plot,
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
            "p = ",
            formatC(
                fdr_NCmono_Tcell_Bcell,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            NCmono_Tcell_Bcell_plot$proportion,
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
        title = "NCmono_Tcell_Bcell multiplets",
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


print(sccomp_NCmono_Tcell_Bcell_dotplot)


# ---------------------------------------------------------
# Save plot
# ---------------------------------------------------------

ggsave(
    plot = sccomp_NCmono_Tcell_Bcell_dotplot,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_NCmono_Tcell_Bcell_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    device = cairo_pdf,
    dpi = 600
)




