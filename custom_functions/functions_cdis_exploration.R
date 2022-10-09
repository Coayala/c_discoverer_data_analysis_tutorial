#
#
# Christian Ayala
# Functions for Labeled_analysis.Rmd
#
#
# -------------------------------------------------------------------------
find_duplicates <- function(df, ...){ #Based on get_dupes from the janitor package (https://github.com/sfirke/janitor/blob/master/R/get_dupes.R)
  ##Find elements that are duplicated in a given dataframe
  
  ## Get which columns will be used to group and find duplicates
  expr <- rlang::expr(c(...))
  pos <- tidyselect::eval_select(expr, data = df)
  
  ## Check if using specific columns or the whole dataframe to find duplication
  if(rlang::dots_n(...) == 0){ # whole dataframe
    col_names <- names(df)
    col_names <- rlang::syms(col_names)
  } else { # only specific columns
    col_names <- names(pos)
    col_names <- rlang::syms(col_names)
  }
  
  # Count and filter duplicated columns
  dup.df <- df %>% 
    add_count(!!! col_names, name = "Times_repeated") %>% 
    filter(Times_repeated > 1) %>% 
    select(Times_repeated, everything())
  
  
  return(dup.df)
}

solve_duplicates <- function(df, ...){
  # Function with the set the rules to solve duplicate compounds in the data matrix
  
  ## Obtain duplicated columns
  dup.df <- find_duplicates(df, ...)
  
  # Rules to solve duplicated columns based on type of information
  
  
}
# -------------------------------------------------------------------------
get_elements <- function(df){
  ## Obtain element list in the formula
  element_list <- str_c(df$Formula[!is.na(df$Formula)], collapse = ' ') %>% 
    str_remove_all('[:digit:]') %>% 
    str_split(' ')
  
  element_list <- unique(unlist(element_list))
  
  return(element_list)
}


# -------------------------------------------------------------------------
separate_formula <- function(df){
  # This function will split the Formula column
  # into columns with the number of each element
  
  ## Get the elements that make each of the compounds
  element_list <- get_elements(df)
  
  ## Get formula in a df were results will be stored
  result_df <- select(df, FeatureID, Formula) 
  
  el_cols <- map(element_list, function(x){
    temp_df <- result_df %>% 
      mutate(el = str_extract(Formula, paste0(x, '\\d{0,}'))) %>% 
      filter(!is.na(el)) %>% 
      mutate(el = ifelse(el == x, 1, str_remove(el, '[:alpha:]')),
             el = as.numeric(el)) %>% 
      rename({{x}} := el)
  })
  
  replace_list <- map(element_list, function(x) 0)
  names(replace_list) <- element_list
  
  final_df <- reduce(el_cols, full_join) %>% 
    replace_na(replace_list) %>% 
    left_join(df)
  
  
  return(final_df)
  
}

# -------------------------------------------------------------------------
calc_ratios_n_idxs <- function(df){
  # This function will calculate H/c and O/C ratios
  # as well as other thermodynamics index
  
  ## Get ratios
  df <- df %>% 
    mutate(H_to_C = H / C) %>% 
    mutate(O_to_C = O / C)
  
  ## Calculate thermodynamic indices
  df <- df %>% 
    mutate(NOSC = -((4*C + H - 3*N - 2* O + 5*P - 2*S) / C) + 4) %>% 
    mutate(GFE = 60.3 - 28.5*NOSC) %>% 
    mutate(DBE = 1 + 0.5 * (2*C - H + N + P)) %>% 
    mutate(DBE_O = DBE - O) %>% 
    mutate(AI = (1 + C - O - S - ((H + P + N) * 0.5)) / (C - O - S - N - P)) %>% 
    mutate(AI_mod = (1 + C - (O * 0.5) - S - ((H + P + N) * 0.5)) / (C - (O * 0.5) - S - N - P)) %>% 
    mutate(DBE_AI = 1 + C -O -S - (0.5 * (H + N + P)))
}

# -------------------------------------------------------------------------
calc_classes <- function(df){
  df <- df %>% 
    mutate(Class = case_when(
      between(O_to_C, 0, 0.3) & between(H_to_C, 1.5, 2.5) ~ 'Lipid',
      between(O_to_C, 0, 0.125) & between(H_to_C, 0.8, 1.5) ~ 'Unsaturated HC',
      between(O_to_C, 0, 0.95) & between(H_to_C, 0.2, 0.8) ~ 'Condensed HC',
      between(O_to_C, 0.3, 0.55) & between(H_to_C, 1.5, 2.3) ~ 'Protein',
      between(O_to_C, 0.55, 0.7) & between(H_to_C, 1.5, 2.2) ~ 'Amino Sugar',
      between(O_to_C, 0.7, 1.5) & between(H_to_C, 1.5, 2.5) ~ 'Carbohydrate',
      between(O_to_C, 0.125, 0.65) & between(H_to_C, 0.8, 1.5) ~ 'Lignin',
      between(O_to_C, 0.65, 1.1) & between(H_to_C, 0.8, 1.5) ~ 'Tannin', 
      TRUE ~ 'Other'))
             
             

  return(df)
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
plot_col <- function(df, my_x, my_y, color_by1, dodge = FALSE){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}})) +
    geom_col(aes(fill = {{color_by1}}),
             size = 2,
             width = 0.75,
             position = ifelse(dodge == TRUE, 'dodge', 'stack')) +
    scale_fill_jco() +
    scale_color_jama() +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
}

# -------------------------------------------------------------------------
plot_boxplot <- function(df, my_x, my_y, color_by, my_comparisons = NULL){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}},
             fill = {{color_by}})) +
    geom_boxplot() +
    scale_fill_jama() +
    stat_compare_means(comparisons = {{my_comparisons}},
                       method = 't.test',
                       label = 'p.signif') +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
}

# -------------------------------------------------------------------------
plot_venn <- function(my_list, my_colors){
  venn(my_list,
       zcolor = {{my_colors}},
       ilcs = 1,
       sncs = 1)
}

# -------------------------------------------------------------------------
plot_density <- function(df, my_x, color_by, facet_by = NULL, facet_by2 = NULL){
  ggplot(df,
         aes(x = {{my_x}},
             fill = {{color_by}})) +
    geom_density(alpha = 0.6) +
    scale_fill_jama() +
    theme_bw() +
    facet_grid(rows = vars({{facet_by}}),
               cols = vars({{facet_by2}})) +
    theme(axis.title.y = element_blank()) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
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