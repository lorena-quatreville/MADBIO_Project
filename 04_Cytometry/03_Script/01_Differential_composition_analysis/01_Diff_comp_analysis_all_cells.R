

# Loading the csv file
path_csv = file.path("/home/quatreville","03_BIOINFO_PROJECT", "01_MADBIO_Cameroon", "01_scRNASeq", "86_cytometry_Toulouse", "00_InputData","Cytometry_all_celltype_without_negative_long_format.csv")
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

# Removing the "Total" rows
cyto_sccomp <- cyto_df %>%
  filter(Cell_type != "Total")
  
  # Chacking if all individuals get the same cell types
  table(cyto_sccomp$Individual)
  
  # To check that total counts correspond to the FowJo totals
  cyto_sccomp %>%
  group_by(Individual) %>%
  summarise(
    calculated_total = sum(Count),
    .groups = "drop"
  )

# Put Individual and Cell_type in factors
 cyto_sccomp <- cyto_sccomp %>%
  mutate(
    Individual = factor(Individual),
    Cell_type = factor(Cell_type)
  ) 


# Calculate proportions just to control visually the data
cyto_prop <- cyto_sccomp %>%
  group_by(Individual) %>%
  mutate(
    proportion = Count / sum(Count)
  ) %>%
  ungroup()

cyto_prop %>%
  group_by(Cell_type, condition) %>%
  summarise(
    mean_proportion = mean(proportion),
    median_proportion = median(proportion),
    .groups = "drop"
  )


# scComp
cyto_sccomp_fit <- sccomp_estimate(
    .data = cyto_sccomp,
    formula_composition = ~ condition,
    .sample = Individual,
    .cell_group = Cell_type,
    .abundance = Count
)

cyto_effects <- cyto_sccomp_fit %>%
    dplyr::filter(parameter == "conditionS") %>%
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

# Tester la significativité des effets
cyto_sccomp_test <- sccomp_test(cyto_sccomp_fit)
cyto_sccomp_test %>%
    dplyr::filter(parameter == "conditionS") %>%
    print(n = 18)

# Save
cyto_sccomp_results <- cyto_sccomp_test %>%
    dplyr::filter(parameter == "conditionS") %>%
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
        paste0("sccomp_condition_results.csv")
    ),
    row.names = FALSE
)


# To see if results are coherent (visualization)
library(ggplot2)

ggplot(
    cyto_prop,
    aes(x = condition, y = proportion, group = Individual)
) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0) +
    facet_wrap(~ Cell_type, scales = "free_y") +
    theme_bw()
    
    
    
    
# ---------------------------------------------------------
# Classical monocytes - sccomp dotplot
# ---------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggsignif)

# Extract proportions for non-classical monocytes
classicalmono_plot <- cyto_prop %>%
    filter(Cell_type == "Classical monocytes")

# Get sccomp FDR
fdr_val <- cyto_sccomp_test %>%
    filter(
        Cell_type == "Classical monocytes",
        parameter == "conditionS"
    ) %>%
    pull(c_FDR)

# Colors
cols <- c(
    "C" = "orange",
    "S" = "red"
)

# Dotplot
sccomp_classicalmono_dotplot <- ggplot(
    classicalmono_plot,
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
        comparisons = list(c("C", "S")),
        annotations = paste0(
            "p = ",
            formatC(
                fdr_val,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            classicalmono_plot$proportion,
            na.rm = TRUE
        ) * 1.08,
        tip_length = 0.01,
        textsize = 7,
        color = "black"
    ) +
    
    # Colors
    scale_color_manual(values = cols) +
    
    theme_classic(
        base_family = "Arial",
        base_size = 16
    ) +
    
    labs(
        title = "Classical monocytes",
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

print(sccomp_classicalmono_dotplot)

ggsave(
    plot = sccomp_classicalmono_dotplot,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(
            "sccomp_classical_monocytes_dotplot.pdf"
        )
    ),
    width = 8,
    height = 6,
    device = cairo_pdf,
    dpi = 600
)



################################################

# Monocytes groupés

###############################################






# ---------------------------------------------------------
# 1. Regroup Classical + Non-classical monocytes
# ---------------------------------------------------------

cyto_sccomp_mono <- cyto_sccomp %>%
    mutate(
        Cell_type = case_when(
            Cell_type %in% c(
                "Classical monocytes",
                "Non-classical monocytes"
            ) ~ "Monocytes",
            TRUE ~ as.character(Cell_type)
        )
    ) %>%
    group_by(Individual, condition, Cell_type) %>%
    summarise(
        Count = sum(Count),
        .groups = "drop"
    ) %>%
    mutate(
        Individual = factor(Individual),
        Cell_type = factor(Cell_type)
    )


# Check that monocytes have been correctly regrouped
cyto_sccomp_mono %>%
    filter(Cell_type == "Monocytes") %>%
    arrange(Individual)


# ---------------------------------------------------------
# 2. Calculate proportions
# ---------------------------------------------------------

cyto_prop_mono <- cyto_sccomp_mono %>%
    group_by(Individual) %>%
    mutate(
        proportion = Count / sum(Count)
    ) %>%
    ungroup()


# Optional check
cyto_prop_mono %>%
    filter(Cell_type == "Monocytes") %>%
    select(Individual, condition, Cell_type, Count, proportion) %>%
    arrange(Individual)


# ---------------------------------------------------------
# 3. scComp
# ---------------------------------------------------------

cyto_sccomp_fit_mono <- sccomp_estimate(
    .data = cyto_sccomp_mono,
    formula_composition = ~ condition,
    .sample = Individual,
    .cell_group = Cell_type,
    .abundance = Count
)


# ---------------------------------------------------------
# 4. Extract effects
# ---------------------------------------------------------

cyto_effects_mono <- cyto_sccomp_fit_mono %>%
    dplyr::filter(parameter == "conditionS") %>%
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

cyto_effects_mono


# ---------------------------------------------------------
# 5. Test significance
# ---------------------------------------------------------

cyto_sccomp_test_mono <- sccomp_test(cyto_sccomp_fit_mono)

cyto_sccomp_test_mono %>%
    dplyr::filter(parameter == "conditionS") %>%
    print(n = 30)


# ---------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------

cyto_sccomp_results_mono <- cyto_sccomp_test_mono %>%
    dplyr::filter(parameter == "conditionS") %>%
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
    cyto_sccomp_results_mono,
    file = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_condition_results_monocytes.csv"
    ),
    row.names = FALSE
)


# ---------------------------------------------------------
# 7. Visualization - Monocytes
# ---------------------------------------------------------

# Extract proportions for Monocytes
monocytes_plot <- cyto_prop_mono %>%
    filter(Cell_type == "Monocytes")


# Get sccomp FDR
fdr_val_mono <- cyto_sccomp_test_mono %>%
    filter(
        Cell_type == "Monocytes",
        parameter == "conditionS"
    ) %>%
    pull(c_FDR)


# Colors
cols <- c(
    "C" = "orange",
    "S" = "red"
)


# Dotplot
sccomp_monocytes_dotplot <- ggplot(
    monocytes_plot,
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
        comparisons = list(c("C", "S")),
        annotations = paste0(
            "p < ",
            formatC(
                fdr_val_mono + 0.001,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            monocytes_plot$proportion,
            na.rm = TRUE
        ) * 1.08,
        tip_length = 0.01,
        textsize = 7,
        color = "black"
    ) +

    # Colors
    scale_color_manual(values = cols) +

    theme_classic(
        base_family = "Arial",
        base_size = 16
    ) +

    labs(
        title = "Monocytes",
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


print(sccomp_monocytes_dotplot)


# ---------------------------------------------------------
# 8. Save plot
# ---------------------------------------------------------

ggsave(
    plot = sccomp_monocytes_dotplot,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_monocytes_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    device = cairo_pdf,
    dpi = 600
)















