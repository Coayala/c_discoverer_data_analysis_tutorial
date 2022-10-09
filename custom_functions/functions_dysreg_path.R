#
#
# Christian Ayala
# Functions for 4_Dysregulated_Pathways.Rmd
#
#
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
plot_col <- function(df, my_x, my_y, color_by1, color_by2, my_colors, dodge = FALSE){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}})) +
    geom_col(aes(fill = {{color_by1}},
                 color = {{color_by2}}),
             size = 2,
             width = 0.75,
             position = ifelse(dodge == TRUE, 'dodge', 'stack')) +
    scale_fill_manual(values = my_colors) +
    scale_color_jama() +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
}