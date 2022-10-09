#
# Christian Ayala
# Functions for Data normalization and statistics
#
#
#
## The following normalization functions were adapted from the NormalyzerDe
## package at https://github.com/ComputationalProteomics/NormalyzerDE/blob/master/R/normMethods.R
# -------------------------------------------------------------------------
global.norm <- function(matrix, transform_data = TRUE){
  # This function will perform normalization based in the global AUC of each sample
  # and the median of such intensities across samples
  
  colsum <- colSums(matrix, na.rm = TRUE)
  colsum.median <- median(colsum)
  norm.matrix <- data.frame(matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colsum[col]) * colsum.median
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  
  return(norm.matrix)
}


# -------------------------------------------------------------------------
median.norm <- function(matrix, transform_data = TRUE){
  # This function will perform data normalization based in the  median AUC of each sample
  
  colmedian <- apply(matrix, 2, FUN = median, na.rm = TRUE)
  colmedian.mean <- mean(colmedian)
  norm.matrix <- data.frame(matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colmedian[col]) * colmedian.mean
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  return(norm.matrix)
} 

# -------------------------------------------------------------------------
mean.norm <- function(matrix, transform_data = TRUE){
  # This function will perform data normalization based in the  mean AUC of each sample
  
  colmean <- colMeans(matrix, na.rm = TRUE)
  colmean.mean <- mean(colmean)
  norm.matrix <- data.frame(matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colmean[col]) * colmean.mean
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  return(norm.matrix)
} 

# -------------------------------------------------------------------------
vsn.norm <- function(matrix){
  # This functions tries to adjust the data to the vsn normalization
  
  norm.matrix <- suppressMessages(vsn::justvsn(as.matrix(matrix)))
  norm.matrix <- as.data.frame(norm.matrix)
  return(norm.matrix)
}
# -------------------------------------------------------------------------
cycloess.norm <- function(matrix, missing_data = FALSE){
    # This functions tries to adjust the data to the vsn normalization
    
    if(missing_data){
      norm.matrix <- log2(matrix + 1)
    } else {
      norm.matrix <- log2(matrix)
    }
    norm.matrix <- limma::normalizeCyclicLoess(norm.matrix, method = 'fast')
    norm.matrix <- as.data.frame(norm.matrix)
    rownames(norm.matrix) <- rownames(matrix)
    return(norm.matrix)
}

# -------------------------------------------------------------------------

max.norm <- function(matrix, transform_data = TRUE){
  # This function will perform data normalization based in the  max AUC of each sample
  
  colmax <- apply(matrix, 2, FUN = max, na.rm = TRUE)
  norm.matrix <- data.frame(matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix)))
  for(col in 1:ncol(matrix)){
    norm.matrix[,col] <- (matrix[,col] / colmax[col])
  }
  colnames(norm.matrix) <- colnames(matrix)
  rownames(norm.matrix) <- rownames(matrix)
  if(transform_data == TRUE){
    norm.matrix <- log2(norm.matrix + 1)
  }
  return(norm.matrix)
} 

# -------------------------------------------------------------------------
plot_boxplot <- function(df, my_x, my_y){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}})) +
    geom_boxplot() +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
}

# -------------------------------------------------------------------------
normalize_by_all <- function(df){

  no_norm.df <- gather(df, key = 'SampleID', value = 'AUC')
  no_norm.plot <- plot_boxplot(no_norm.df, SampleID, AUC) +
    labs(title = 'No normalization') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())

  gi_norm.df <- global.norm(df)
  gi_norm.df <- gather(gi_norm.df, key = 'SampleID', value = 'AUC')
  gi_norm.plot <- plot_boxplot(gi_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by global AUC',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  

  mean_norm.df <- mean.norm(df)
  mean_norm.df <- gather(mean_norm.df, key = 'SampleID', value = 'AUC')
  mean_norm.plot <- plot_boxplot(mean_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by mean',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  

  median_norm.df <- median.norm(df)
  median_norm.df <- gather(median_norm.df, key = 'SampleID', value = 'AUC')
  median_norm.plot <- plot_boxplot(median_norm.df, SampleID, AUC) +
    labs(title = 'Normalization by median',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  

  vsn_norm.df <- vsn.norm(df)
  vsn_norm.df <- gather(vsn_norm.df, key = 'SampleID', value = 'AUC')
  vsn_norm.plot <- plot_boxplot(vsn_norm.df, SampleID, AUC) +
    labs(title = 'VSN',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  

  cycloess_norm.df <- cycloess.norm(df)
  cycloess_norm.df <- gather(cycloess_norm.df, key = 'SampleID', value = 'AUC')
  cycloess_norm.plot <- plot_boxplot(cycloess_norm.df, SampleID, AUC) +
    labs(title = 'LOESS normalization',
         y =  'Normalized AUC') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  

  all_norm.plot <- ggarrange(no_norm.plot, gi_norm.plot, mean_norm.plot,
                                median_norm.plot, vsn_norm.plot, cycloess_norm.plot,
                                nrow = 2,
                                ncol = 3)
  
  return(all_norm.plot)

}

#---------------------------------------------------------------------------
plot_cumvar <- function(eigen){
  dimensions <- c(1:dim(eigen)[1])
  cumvar <- ggplot(data=eigen, aes(x=as.factor(dimensions), y=cumulative.variance.percent/100, group=1)) +
    geom_line(size=1.2, color='black') +
    geom_point(size=3, color='black') +
    xlab('PC axes') + ylab('Amount of explained variance') +
    theme_bw() +
    theme(plot.title = element_text(face="bold", hjust = 0.5)) +
    labs(title = 'Cumulative variance plot')
  return(cumvar)
}
