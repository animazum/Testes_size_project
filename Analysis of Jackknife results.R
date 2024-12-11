# Load or install required packages
required_packages <- c("rmarkdown", "ggplot2", "reshape2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(ggplot2)
library(reshape2)

generate_jackknife_pdf <- function(jackknife_object, alpha = 0.05, max_iterations = 1000, output_pdf = "jackknife_analysis.pdf") {
  # Validate input
  required_fields <- c("jackknife_results", "jackknife_means", "jackknife_vars", "failed_iterations")
  if (!is.list(jackknife_object) || !all(required_fields %in% names(jackknife_object))) {
    stop("Invalid jackknife_object structure.")
  }
  
  # Extract components
  results <- jackknife_object$jackknife_results
  failed <- jackknife_object$failed_iterations
  if (length(results) > max_iterations) {
    warning(sprintf("Truncating results to max_iterations (%d).", max_iterations))
    results <- results[1:max_iterations]
  }
  
  # Convert results to a dataframe
  results_df <- do.call(rbind, results)
  if (is.null(results_df) || nrow(results_df) == 0) stop("No valid results to analyze.")
  
  # Calculate summary statistics
  means <- colMeans(results_df, na.rm = TRUE)
  se <- apply(results_df, 2, sd, na.rm = TRUE) / sqrt(nrow(results_df))
  ci_lower <- apply(results_df, 2, quantile, alpha / 2, na.rm = TRUE)
  ci_upper <- apply(results_df, 2, quantile, 1 - alpha / 2, na.rm = TRUE)
  cv <- apply(results_df, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)
  z_scores <- scale(results_df, center = TRUE, scale = TRUE)
  influential_species <- which(apply(abs(z_scores), 1, max, na.rm = TRUE) > 2)
  
  # Create coefficient plot
  coef_data <- data.frame(Parameter = names(means), Estimate = means, SE = se, 
                          CI_lower = ci_lower, CI_upper = ci_upper)
  coef_plot <- ggplot(coef_data, aes(x = Parameter, y = Estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Coefficient Estimates with Confidence Intervals", y = "Estimate", x = "Parameter")
  
  # Create stability plot
  stability_data <- data.frame(Parameter = names(cv), CV = cv)
  stability_plot <- ggplot(stability_data, aes(x = Parameter, y = CV)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Parameter Stability Analysis", y = "Coefficient of Variation (%)", x = "Parameter")
  
  # Create density plot
  melted_results <- melt(results_df)
  dist_plot <- ggplot(melted_results, aes(x = value)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    facet_wrap(~variable, scales = "free") +
    theme_minimal() +
    labs(title = "Distribution of Parameter Estimates", x = "Value", y = "Density")
  
  # Save to PDF
  pdf(file = output_pdf, width = 8.5, height = 11)
  print(coef_plot)
  print(stability_plot)
  print(dist_plot)
  dev.off()
  
  cat(sprintf("PDF report generated: %s\n", output_pdf))
}

# Example usage
# Replace `jackknife_object` with your actual object.
generate_jackknife_pdf(jackknife_results, alpha = 0.05, output_pdf = "jackknife_report.pdf")


ggplot()


###############

bootstrap_analysis <- function(data, n_bootstrap = 1000, alpha = 0.05) {
  # Store bootstrap results
  bootstrap_results <- replicate(n_bootstrap, {
    sampled_data <- data[sample(1:nrow(data), replace = TRUE), ]
    # Calculate your statistic (e.g., correlation)
    cor(sampled_data$gene_families, sampled_data$trait, use = "complete.obs")
  })
  
  # Summary statistics
  mean_estimate <- mean(bootstrap_results)
  ci <- quantile(bootstrap_results, c(alpha / 2, 1 - alpha / 2))
  
  list(
    bootstrap_distribution = bootstrap_results,
    mean = mean_estimate,
    ci_lower = ci[1],
    ci_upper = ci[2]
  )
}

# Example usage
bootstrap_results <- bootstrap_analysis(results_df)
print(bootstrap_results$mean)
print(bootstrap_results$ci_lower)
print(bootstrap_results$ci_upper)

# Plot bootstrap distribution
library(ggplot2)
ggplot(data.frame(Value = bootstrap_results$bootstrap_distribution), aes(x = Value)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(title = "Bootstrap Distribution of Correlation", x = "Correlation", y = "Density")


