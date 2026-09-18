############################# Monocyte singlets vs multiplets ##########################

# Loading the csv file
path_csv = file.path("/home/quatreville","03_BIOINFO_PROJECT", "01_MADBIO_Cameroon", "01_scRNASeq", "86_cytometry_Toulouse", "00_InputData","Mono_singlet_multiplet.csv")
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

  
  # Chacking if all individuals get the same cell types
  table(cyto_df$Individual)
  
  # To check that total counts correspond to the FowJo totals
  cyto_df %>%
  group_by(Individual) %>%
  summarise(
    calculated_total = sum(Count),
    .groups = "drop"
  )

# Put Individual and Cell_type in factors
 cyto_sccomp <- cyto_df %>%
  mutate(
    Individual = factor(Individual),
    Cell_type = factor(Cell_type)
  ) 


# ---------------------------------------------------------
# Check total number of monocytes per individual
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
# Specifically calculate multiplet proportion
# ---------------------------------------------------------

multiplet_prop <- cyto_prop %>%
    filter(Cell_type == "Mono_multiplets")


multiplet_prop %>%
    select(
        Individual,
        condition,
        Count,
        proportion
    )


# ---------------------------------------------------------
# sccomp
# ---------------------------------------------------------

cyto_sccomp_fit <- sccomp_estimate(
    .data = cyto_sccomp,
    formula_composition = ~ condition,
    .sample = Individual,
    .cell_group = Cell_type,
    .abundance = Count
)


# ---------------------------------------------------------
# Extract condition effect
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
# Save sccomp results
# ---------------------------------------------------------

cyto_sccomp_results <- cyto_sccomp_test %>%
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


write.csv(
    cyto_sccomp_results,
    file = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_mono_singlet_multiplet_results.csv"
    ),
    row.names = FALSE
)


# ---------------------------------------------------------
# Mono multiplets - dotplot
# ---------------------------------------------------------

multiplets_plot <- cyto_prop %>%
    filter(
        Cell_type == "Mono_multiplets"
    )


# Extract sccomp FDR for multiplets
fdr_val <- cyto_sccomp_test %>%
    filter(
        Cell_type == "Mono_multiplets",
        parameter == "conditionS"
    ) %>%
    pull(c_FDR)


# Colors
cols <- c(
    "C" = "orange",
    "S" = "red"
)


# ---------------------------------------------------------
# Dotplot
# ---------------------------------------------------------

sccomp_multiplets_dotplot <- ggplot(
    multiplets_plot,
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
                fdr_val,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            multiplets_plot$proportion,
            na.rm = TRUE
        ) * 1.08,
        tip_length = 0.01,
        textsize = 7,
        color = "black"
    ) +

    # Colors
    scale_color_manual(
        values = cols
    ) +

    theme_classic(
        base_family = "Arial",
        base_size = 16
    ) +

    labs(
        title = "Monocyte multiplets",
        y = "Proportion of multiplets",
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


print(sccomp_multiplets_dotplot)


# ---------------------------------------------------------
# Save plot
# ---------------------------------------------------------

ggsave(
    plot = sccomp_multiplets_dotplot,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_mono_multiplets_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    device = cairo_pdf,
    dpi = 600
)












