#
#
# Christian Ayala
# Functions for 3_Differential_analysis.Rmd
#
#
# -------------------------------------------------------------------------
get_samples <- function(metadata.df, Treatment, value){
  # Get value to filter samples
  selector <- syms({{Treatment}})
  
  samples <- metadata.df %>% 
    filter((!!! selector) == value)
  
  # Get only sample names
  samples <- samples$Sample_code
  
  return(samples)
}

# -------------------------------------------------------------------------
get_vectors <- function(df, filter_by, value, get_col){
  # Column where the value will be filtered
  filter_col <- syms({{filter_by}})
  
  # Column that will be retrieved
  get <- syms({{get_col}})
  
  vector <- df %>% 
    filter((!!! filter_col) == value) %>% 
    pull((!!! get))
  
  return(vector)
}

# -------------------------------------------------------------------------
get_diff_table <- function(auc_matrix, control.sample_list, treatment.sample_list, log2_transformed = FALSE){
  
  # Get the AUC values per each sample and calculate the means per feature
  temp.df_control <- auc_matrix %>% 
    select(all_of(control.sample_list))
  control_means <- rowMeans(temp.df_control, na.rm = TRUE)
  
  temp.df_treatment <- auc_matrix %>% 
    select(all_of(treatment.sample_list))
  treatment_means <- rowMeans(temp.df_treatment, na.rm = TRUE)
  
  diff_table <- as.data.frame(cbind(control_means, treatment_means))
  
  if(log2_transformed == TRUE){
    diff_table <- diff_table %>% 
      mutate(log2FC = treatment_means - control_means)
  } else {
    diff_table <- diff_table %>% 
      mutate(ratio = treatment_means/control_means) %>% # get the control/treatment ratio
      mutate(log2FC = log2(ratio)) # calculate log2FC
  }
  
  rownames(diff_table) <- rownames(auc_matrix)
  
  # Initialize pvalues matrix
  pvalues <- data.frame(row.names = rownames(auc_matrix), pval = rep(0, length(rownames(auc_matrix))))
  
  #Calculate pvalue per each of the features
  for(i in 1:nrow(pvalues)){
    t.test <- t.test(as.numeric(temp.df_control[i,]), as.numeric(temp.df_treatment[i,]), paired = FALSE)
    pvalues$pval[i] <- t.test$p.value
  }
  pvalues <- pvalues %>% 
    rownames_to_column(var = 'FeatureID')
  
  diff_table <- diff_table %>% 
    rownames_to_column(var = 'FeatureID') %>% 
    left_join(pvalues, by = 'FeatureID') %>% 
    filter(!is.nan(pval))
  
  diff_table$pval.adj <- p.adjust(diff_table$pval, method = 'fdr')
  
  return(diff_table)
  
}

# -------------------------------------------------------------------------
plot_volcano <- function(df, 
                         log2FC, 
                         pval, 
                         log2FC.threshold, 
                         pval.threshold, 
                         search_label = FALSE,
                         label_ids_control = NA,
                         label_ids_treatment = NA){
  
  #Generate label for the plot
  
  significant_points <- df %>% 
    select(FeatureID, {{log2FC}}, {{pval}}) %>% 
    filter(abs({{log2FC}}) > log2FC.threshold,
           -log10({{pval}}) > -log10(pval.threshold)) %>% 
    pull(FeatureID)
  
  if(search_label){
    df <- df %>% 
      mutate(label = case_when(
        (FeatureID %in% label_ids_control) & 
          (FeatureID %in% label_ids_treatment) ~ 'Label in both',
        FeatureID %in% label_ids_control ~ 'Label in control',
        FeatureID %in% label_ids_treatment ~ 'Label in treatment',
        TRUE ~ 'Unlabeled'))
    
    shapes <- c('Label in both' = 16,
                'Label in control' = 15,
                'Label in treatment' = 17,
                'Unlabeled' = 1)
    
    plot <- df %>%
      mutate(color4plot = ifelse(FeatureID %in% significant_points, 'significant', 'non-significant')) %>% 
      ggplot(aes(x = {{log2FC}},
                 y = -log10({{pval}}),
                 shape = label)) +
      geom_point(aes(color = color4plot)) +
      scale_color_manual(values = c("significant" = 'red', 'non-significant' = 'black')) +
      geom_vline(xintercept = c(-{{log2FC.threshold}}, {{log2FC.threshold}}),
                 linetype = 'dotted',
                 size = 1,
                 color = 'blue') +
      geom_hline(yintercept = -log10({{pval.threshold}}),
                 linetype = 'dotted',
                 size = 1,
                 color = 'blue') +
      scale_shape_manual(values = shapes) +
      guides(color = 'none') +
      theme_bw() +
      labs(title = 'Volcano plot',
           x = expression("Log"[2]*" Fold Change"),
           y = expression("-Log"[10]*" pvalue")) +
      theme(plot.title = element_text(hjust = 0.5,
                                      face = 'bold'),
            plot.subtitle = element_text(hjust = 0.5,
                                         face = 'bold'))
  } else {
    plot <- df %>%
      mutate(color4plot = ifelse(FeatureID %in% significant_points, 'significant', 'non-significant')) %>% 
      ggplot(aes(x = {{log2FC}},
                 y = -log10({{pval}}))) +
      geom_point(aes(color = color4plot)) +
      scale_color_manual(values = c("significant" = 'red', 'non-significant' = 'black')) +
      geom_vline(xintercept = c(-{{log2FC.threshold}}, {{log2FC.threshold}}),
                 linetype = 'dotted',
                 size = 1,
                 color = 'blue') +
      geom_hline(yintercept = -log10({{pval.threshold}}),
                 linetype = 'dotted',
                 size = 1,
                 color = 'blue') +
      theme_bw() +
      labs(title = 'Volcano plot',
           x = expression("Log"[2]*" Fold Change"),
           y = expression("-Log"[10]*" pvalue")) +
      theme(plot.title = element_text(hjust = 0.5,
                                      face = 'bold'),
            plot.subtitle = element_text(hjust = 0.5,
                                         face = 'bold'),
            legend.position = 'none')
  }
  
  
  
  return(plot)
}

# -------------------------------------------------------------------------
plot_venn <- function(my_list, my_colors){
  venn(my_list,
       zcolor = {{my_colors}},
       ilcs = 1,
       sncs = 1)
}

# -------------------------------------------------------------------------
plot_vank <- function(df, color_by, facet_by = NULL, facet_by2 = NULL){
  ggplot(df,
         aes(x = O_to_C,
             y = H_to_C,
             color = {{color_by}})) +
    geom_point(size = 2) +
    scale_color_igv() +
    theme_bw() +
    labs(title = 'Van Krevelen Diagram',
         x = 'O/C',
         y = 'H/C') +
    theme(plot.title = element_text(face = 'bold',
                                    hjust = 0.5)) +
    facet_grid(rows = vars({{facet_by}}),
               cols = vars({{facet_by2}}))
  
}

# -------------------------------------------------------------------------
fix_pval_distribution <- function(diff_table){
  fdrtool_res <- fdrtool(diff_table$pval, statistic = 'pvalue') 
  diff_table$fdrtool_pvalue <- fdrtool_res$pval
  diff_table$fdrtool_qval <- fdrtool_res$qval
  
  return(diff_table)
}
