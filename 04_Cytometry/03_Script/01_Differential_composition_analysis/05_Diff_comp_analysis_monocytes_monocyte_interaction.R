############################# Monocyte singlets vs multiplets ##########################

# Loading the csv file
path_csv = file.path("/home/quatreville","03_BIOINFO_PROJECT", "01_MADBIO_Cameroon", "01_scRNASeq", "86_cytometry_Toulouse", "00_InputData","Cmono_NCmono_Imono_interaction.csv")
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

# ---------------------------------------------------------
# Check condition
# ---------------------------------------------------------

table(
    cyto_df$Individual,
    cyto_df$condition
)


# ---------------------------------------------------------
# Checks
# ---------------------------------------------------------

# Check number of categories per individual
table(cyto_df$Individual)


# Check that every individual has the 6 categories
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
                "Cmono_singlets",
                "Cmono_multiplets",
                "NCmono_singlet",
                "NCmono_multiplet",
                "Imono_singlets",
                "Imono_multiplet"
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

set.seed(1234)

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
    ) %>%
    print(n = 20)


# ---------------------------------------------------------
# Extract results for all 6 categories
# ---------------------------------------------------------

cyto_sccomp_singlets_multiplets <- cyto_sccomp_test %>%
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
# Display results
# ---------------------------------------------------------

cyto_sccomp_singlets_multiplets


# =========================================================
# Save results
# =========================================================

write.csv(
    cyto_sccomp_singlets_multiplets,
    file = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_singlets_vs_multiplets_results.csv"
    ),
    row.names = FALSE
)




# =========================================================
# Visualization of singlet vs multiplet proportions
# Cmono and NCmono only
# =========================================================


# =========================================================
# Cmono
# =========================================================

# Build Cmono dataframe
df_Cmono <- cyto_prop %>%
    filter(
        Cell_type %in% c(
            "Cmono_singlets",
            "Cmono_multiplets"
        )
    ) %>%
    mutate(
        Interaction = case_when(
            Cell_type == "Cmono_singlets" ~ "Singlet",
            Cell_type == "Cmono_multiplets" ~ "Multiplet"
        ),
        Group = paste0(condition, "_", Interaction)
    )


# Retrieve sccomp FDR
p_singlet_C <- cyto_sccomp_singlets_multiplets$c_FDR[
    cyto_sccomp_singlets_multiplets$Cell_type == "Cmono_singlets"
]

p_multiplet_C <- cyto_sccomp_singlets_multiplets$c_FDR[
    cyto_sccomp_singlets_multiplets$Cell_type == "Cmono_multiplets"
]


# Colors
cols_Cmono <- c(
    "C_Singlet" = "#88CCEE",
    "C_Multiplet" = "#332288",
    "S_Singlet" = "#88CCEE",
    "S_Multiplet" = "#332288"
)


# Plot Cmono
dotplot_Cmono <- ggplot(
    df_Cmono,
    aes(
        x = Group,
        y = proportion,
        color = Group
    )
) +
   geom_jitter( width = 0.3, 
   height = 0, 
   size = 5, 
   alpha = 0.9 
   ) +
    
    stat_summary(
        fun = mean,
        geom = "crossbar",
        color = "black",
        width = 0.5,
        linewidth = 0.8
    ) +
    
    scale_color_manual(
        values = cols_Cmono
    ) +
    
    # Singlet comparison: C vs S
    geom_signif(
        comparisons = list(
            c("C_Singlet", "S_Singlet")
        ),
        annotations = paste0(
            "p < ",
            formatC(
                p_singlet_C + 0.001,
                format = "e" ,
                digits = 1
            )
        ),
        y_position = max(
            df_Cmono$proportion,
            na.rm = TRUE
        ) * 1.05,
        tip_length = 0.01,
        textsize = 10,
        color = "black"
    ) +
    
    # Multiplet comparison: C vs S
    geom_signif(
        comparisons = list(
            c("C_Multiplet", "S_Multiplet")
        ),
        annotations = paste0(
            "p = ",
            formatC(
                p_multiplet_C,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            df_Cmono$proportion,
            na.rm = TRUE
        ) * 1.15,
        tip_length = 0.01,
        textsize = 10,
        color = "black"
    ) +
    
    theme_classic(
        base_family = "Arial",
        base_size = 16
    ) +
    
    labs(
        title = NULL,
        x = "Condition",
        y = "Proportion of cells",
        color = "Group"
    ) +
    
    theme(
        legend.position = "right",
        
        legend.text = element_text(
            size = 20,
            family = "Arial"
        ),
        
        legend.title = element_blank(),
        
        axis.title.y = element_text(
            size = 25,
            family = "Arial"
        ),
        
        axis.title.x = element_text(
            size = 25,
            family = "Arial"
        ),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_text(
            size = 25,
            family = "Arial"
        ),
        
        plot.title = element_blank()
    ) +
    
    scale_y_continuous(
        expand = expansion(
            mult = c(0.05, 0.15)
        )
    )


print(dotplot_Cmono)


# Save Cmono PDF
ggsave(
    plot = dotplot_Cmono,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_Cmono_singlet_multiplet_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    dpi = 600,
    device = cairo_pdf
)





# =========================================================
# NCmono
# =========================================================

# Build NCmono dataframe
df_NCmono <- cyto_prop %>%
    filter(
        Cell_type %in% c(
            "NCmono_singlet",
            "NCmono_multiplet"
        )
    ) %>%
    mutate(
        Interaction = case_when(
            Cell_type == "NCmono_singlet" ~ "Singlet",
            Cell_type == "NCmono_multiplet" ~ "Multiplet"
        ),
        Group = paste0(condition, "_", Interaction)
    )


# Retrieve sccomp FDR
p_singlet_NC <- cyto_sccomp_singlets_multiplets$c_FDR[
    cyto_sccomp_singlets_multiplets$Cell_type == "NCmono_singlet"
]

p_multiplet_NC <- cyto_sccomp_singlets_multiplets$c_FDR[
    cyto_sccomp_singlets_multiplets$Cell_type == "NCmono_multiplet"
]


# Colors
cols_NCmono <- c(
    "C_Singlet" = "#E69F00",
    "C_Multiplet" = "red",
    "S_Singlet" = "#E69F00",
    "S_Multiplet" = "red"
)


# Plot NCmono
dotplot_NCmono <- ggplot(
    df_NCmono,
    aes(
        x = Group,
        y = proportion,
        color = Group
    )
) +
   geom_jitter( width = 0.3, 
   height = 0, 
   size = 5, 
   alpha = 0.9 
   ) +
    
    stat_summary(
        fun = mean,
        geom = "crossbar",
        color = "black",
        width = 0.5,
        linewidth = 0.8
    ) +
    
    scale_color_manual(
        values = cols_NCmono
    ) +
    
    # Singlet comparison: C vs S
    geom_signif(
        comparisons = list(
            c("C_Singlet", "S_Singlet")
        ),
        annotations = paste0(
            "p = ",
            formatC(
                p_singlet_NC,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            df_NCmono$proportion,
            na.rm = TRUE
        ) * 1.05,
        tip_length = 0.01,
        textsize = 10,
        color = "black"
    ) +
    
    # Multiplet comparison: C vs S
    geom_signif(
        comparisons = list(
            c("C_Multiplet", "S_Multiplet")
        ),
        annotations = paste0(
            "p = ",
            formatC(
                p_multiplet_NC,
                format = "e",
                digits = 1
            )
        ),
        y_position = max(
            df_NCmono$proportion,
            na.rm = TRUE
        ) * 1.15,
        tip_length = 0.01,
        textsize = 10,
        color = "black"
    ) +
    
    theme_classic(
        base_family = "Arial",
        base_size = 16
    ) +
    
    labs(
        title = NULL,
        x = "Condition",
        y = "Proportion of cells",
        color = "Group"
    ) +
    
    theme(
        legend.position = "right",
        
        legend.text = element_text(
            size = 20,
            family = "Arial"
        ),
        
        legend.title = element_blank(),
        
        axis.title.y = element_text(
            size = 25,
            family = "Arial"
        ),
        
        axis.title.x = element_text(
            size = 25,
            family = "Arial"
        ),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_text(
            size = 25,
            family = "Arial"
        ),
        
        plot.title = element_blank()
    ) +
    
    scale_y_continuous(
        expand = expansion(
            mult = c(0.05, 0.15)
        )
    )


print(dotplot_NCmono)


# Save NCmono PDF
ggsave(
    plot = dotplot_NCmono,
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        "sccomp_NCmono_singlet_multiplet_dotplot.pdf"
    ),
    width = 8,
    height = 6,
    dpi = 600,
    device = cairo_pdf
)




