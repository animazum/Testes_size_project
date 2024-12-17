

### Data analysis of the bootstrap PGLS results

library(ggplot2)
library(tidyr)
library(dplyr)

analyze_pgls_bootstrap <- function(file_path) {
  # Read the bootstrap summary
  results <- readRDS(file_path)
  
  # Get column names (statistics)
  stat_names <- colnames(results$mean)
  
  # Function to summarize statistics across genes
  summarize_stats <- function(stat_idx) {
    data.frame(
      Statistic = stat_names[stat_idx],
      Mean = mean(results$mean[, stat_idx], na.rm = TRUE),
      SD = mean(results$sd[, stat_idx], na.rm = TRUE),
      Lower_CI = mean(results$quantiles[1, , stat_idx], na.rm = TRUE),
      Upper_CI = mean(results$quantiles[2, , stat_idx], na.rm = TRUE),
      Genes_Significant = sum(results$quantiles[2, , stat_idx] < 0.05, na.rm = TRUE)
    )
  }
  
  # Create summary for each statistic
  summary_list <- lapply(1:ncol(results$mean), summarize_stats)
  summary_df <- do.call(rbind, summary_list)
  
  # Print basic summary
  cat("PGLS Bootstrap Analysis Summary:\n")
  cat("-------------------------------\n")
  cat(sprintf("Total number of genes analyzed: %d\n", nrow(results$mean)))
  cat(sprintf("Number of statistics per gene: %d\n\n", ncol(results$mean)))
  
  # Print detailed statistics
  cat("Averaged Statistics Across All Genes:\n")
  print(summary_df, digits = 3)
  
  # Create coefficient plot for averaged data
  coef_data <- summary_df[grep("Coefficient", summary_df$Statistic, ignore.case = TRUE),]
  
  p1 <- ggplot(coef_data, aes(x = Statistic, y = Mean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Average PGLS Coefficients with 95% Confidence Intervals",
         x = "",
         y = "Estimate") +
    theme(axis.text.y = element_text(size = 10))
  
  print(p1)
  
  # Create p-value plot for averaged data
  pval_data <- summary_df[grep("p_value|pval", summary_df$Statistic, ignore.case = TRUE),]
  
  p2 <- ggplot(pval_data, aes(x = Statistic, y = Mean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Average PGLS P-values with 95% Confidence Intervals",
         x = "",
         y = "P-value") +
    theme(axis.text.y = element_text(size = 10))
  
  print(p2)
  
  # Create gene-wise significance summary
  pval_cols <- grep("p_value|pval", stat_names, ignore.case = TRUE)
  sig_summary <- data.frame(
    Statistic = stat_names[pval_cols],
    Total_Genes = nrow(results$mean),
    Significant_Genes = sapply(pval_cols, function(i) 
      sum(results$quantiles[2, , i] < 0.05, na.rm = TRUE)),
    Percent_Significant = sapply(pval_cols, function(i) 
      100 * sum(results$quantiles[2, , i] < 0.05, na.rm = TRUE) / nrow(results$mean))
  )
  
  cat("\nGene-wise Significance Summary:\n")
  cat("-----------------------------\n")
  print(sig_summary, digits = 3)
  
  # Create histogram of p-values for each statistic
  for(i in pval_cols) {
    p3 <- ggplot(data.frame(pvalue = results$mean[, i]), aes(x = pvalue)) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
      theme_minimal() +
      labs(title = paste("Distribution of P-values for", stat_names[i]),
           x = "P-value",
           y = "Count")
    
    print(p3)
  }
  
  # Return both summary statistics and gene-wise significance
  return(list(
    summary_stats = summary_df,
    significance_summary = sig_summary
  ))
}


# Analyze the results
results_analysis <- analyze_pgls_bootstrap("PGLS_bootstrap_summary.rds")

# Access the summary statistics
summary_stats <- results_analysis$summary_stats

# Access the gene-wise significance summary
sig_summary <- results_analysis$significance_summary
