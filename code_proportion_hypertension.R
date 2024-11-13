library(dplyr)

# Read the cleaned data
data <- read.csv("stroke_data_clean.csv")
data <- data[data$gender != "Other", ]
################################
###Sample Size Calculation######
################################

# Set parameters
confidence_level <- 0.95
delta <- 0.015 #from paper
variance_guess <- 0.25 #the conservative guess 0.5*0.5
N <- nrow(data)
z_alpha_2 <- qnorm(1 - (1 - confidence_level)/2)

# Calculate initial sample size
n <- (z_alpha_2^2 * variance_guess) / delta^2

# Apply FPC
n_star <- n/(1 + n/N)
n_final <- ceiling(n_star)

# Calculate FPC
FPC <- 1 - n_final/N

cat("Required sample size with FPC:", n_final, "\n")

################################
###SRS##########################
################################

# Take SRS
set.seed(123)
sample_indices <- sample(1:N, n_final)
sample_data <- data[sample_indices, ]

# Calculate population means for auxiliary variables
pop_mean_glucose <- mean(data$avg_glucose_level)
pop_prop_glucose_high <- mean(data$avg_glucose_level >= 100)

# 1. Vanilla Estimator
p_hat <- mean(sample_data$hypertension)
var_p_hat <- FPC * (p_hat * (1 - p_hat))/(n_final - 1)
se_p_hat <- sqrt(var_p_hat)
ci_vanilla <- c(p_hat - z_alpha_2 * se_p_hat,
                p_hat + z_alpha_2 * se_p_hat)

# 2. Ratio Estimator with Continuous Glucose
R_cont <- mean(sample_data$hypertension)/mean(sample_data$avg_glucose_level)
p_ratio_cont <- R_cont * pop_mean_glucose

# Calculate variance for continuous ratio estimator
residuals_cont <- sample_data$hypertension -
  R_cont * sample_data$avg_glucose_level
var_residuals_cont <- 1 / (n_final - 1) * sum(residuals_cont^2)
var_ratio_cont <- FPC * var_residuals_cont/n_final
se_ratio_cont <- sqrt(var_ratio_cont)
ci_ratio_cont <- c(p_ratio_cont - z_alpha_2 * se_ratio_cont,
                   p_ratio_cont + z_alpha_2 * se_ratio_cont)

# 3. Ratio Estimator with Binary Glucose
sample_binary <- ifelse(sample_data$avg_glucose_level >= 100, 1, 0)
R_bin <- mean(sample_data$hypertension)/mean(sample_binary)
p_ratio_bin <- R_bin * pop_prop_glucose_high

# Calculate variance for binary ratio estimator
residuals_bin <- sample_data$hypertension - R_bin * sample_binary
var_residuals_bin <- 1/(n_final - 1) * sum(residuals_bin^2)
var_ratio_bin <- FPC * var_residuals_bin/n_final
se_ratio_bin <- sqrt(var_ratio_bin)
ci_ratio_bin <- c(p_ratio_bin - z_alpha_2 * se_ratio_bin,
                  p_ratio_bin + z_alpha_2 * se_ratio_bin)

# Print SRS results
cat("\nSimple Random Sampling Results:\n")
cat("\n1. Vanilla Estimator:\n")
cat("Estimated Proportion:", round(p_hat, 4), "\n")
cat("Standard Error:", round(se_p_hat, 4), "\n")
cat("95% CI: [", round(ci_vanilla[1], 4), ",", round(ci_vanilla[2], 4), "]\n")

cat("\n2. Ratio Estimator (Continuous Glucose):\n")
cat("Estimated Proportion:", round(p_ratio_cont, 4), "\n")
cat("Standard Error:", round(se_ratio_cont, 4), "\n")
cat("95% CI: [", round(ci_ratio_cont[1], 4), ",", round(ci_ratio_cont[2], 4), "]\n")

cat("\n3. Ratio Estimator (Binary Glucose):\n")
cat("Estimated Proportion:", round(p_ratio_bin, 4), "\n")
cat("Standard Error:", round(se_ratio_bin, 4), "\n")
cat("95% CI: [", round(ci_ratio_bin[1], 4), ",", round(ci_ratio_bin[2], 4), "]\n")

cat("True Population Proportion:", round(mean(data$hypertension), 4), "\n")

################################
###Stratified Sampling##########
################################

# Calculate stratum information by gender
strata_info <- data %>%
  group_by(smoking_status) %>%
  summarise(
    N_h = n(),
    prop = n()/N,
    p_h = mean(hypertension),
    var_h = p_h * (1 - p_h)
  )

# Calculate optimal allocation (assume equal variance for each strata)
strata_info$n_h <- round(n_final*strata_info$N_h/N)


# Stratified sampling function
stratified_sampling <- function(data, strata_info, seed = 123) {
  set.seed(seed)
  stratum_samples <- list()
  stratum_estimates <- list()

  for(i in 1:nrow(strata_info)) {
    # Get stratum data
    stratum_data <- data %>%
      filter(smoking_status == strata_info$smoking_status[i])

    # Sample size for this stratum
    n_h <- strata_info$n_h[i]
    N_h <- strata_info$N_h[i]

    # Draw sample
    stratum_samples[[i]] <- stratum_data %>%
      slice_sample(n = n_h)

    # Get population means for auxiliary variables in this stratum
    pop_mean_glucose_h <- mean(stratum_data$avg_glucose_level)
    pop_prop_glucose_high_h <- mean(stratum_data$avg_glucose_level >= 100)

    # 1. Vanilla Estimator
    p_hat_h <- mean(stratum_samples[[i]]$hypertension)
    var_p_hat_h <- (1 - n_h/N_h) * (p_hat_h * (1 - p_hat_h))/(n_h - 1)

    # 2. Ratio Estimator (Continuous)
    R_cont_h <- mean(stratum_samples[[i]]$hypertension) /
      mean(stratum_samples[[i]]$avg_glucose_level)
    ratio_est_cont_h <- R_cont_h * pop_mean_glucose_h

    residuals_cont_h <- stratum_samples[[i]]$hypertension -
      R_cont_h * stratum_samples[[i]]$avg_glucose_level
    var_residuals_cont_h <- 1 / (n_h - 1) * sum(residuals_cont_h^2)
    var_ratio_cont_h <- (1 - n_h/N_h) * var_residuals_cont_h/n_h

    # 3. Ratio Estimator (Binary)
    sample_binary_h <- ifelse(stratum_samples[[i]]$avg_glucose_level >= 100, 1, 0)
    R_bin_h <- mean(stratum_samples[[i]]$hypertension)/mean(sample_binary_h)
    ratio_est_bin_h <- R_bin_h * pop_prop_glucose_high_h

    residuals_bin_h <- stratum_samples[[i]]$hypertension -
      R_bin_h * sample_binary_h
    var_residuals_bin_h <- 1 / (n_h - 1) * sum(residuals_bin_h^2)
    var_ratio_bin_h <- (1 - n_h/N_h) * var_residuals_bin_h/n_h

    # Store estimates
    stratum_estimates[[i]] <- data.frame(
      stratum = strata_info$smoking_status[i],
      N_h = N_h,
      n_h = n_h,
      vanilla = p_hat_h,
      ratio_cont = ratio_est_cont_h,
      ratio_bin = ratio_est_bin_h,
      var_vanilla = var_p_hat_h,
      var_ratio_cont = var_ratio_cont_h,
      var_ratio_bin = var_ratio_bin_h
    )
  }

  # Combine stratum estimates
  stratum_estimates <- bind_rows(stratum_estimates)

  # Calculate overall estimates
  overall_est <- data.frame(
    method = c("Vanilla", "Ratio (Continuous)", "Ratio (Binary)"),
    estimate = c(
      sum(stratum_estimates$vanilla * strata_info$N_h) / N,
      sum(stratum_estimates$ratio_cont * strata_info$N_h) / N,
      sum(stratum_estimates$ratio_bin * strata_info$N_h) / N
    ),
    variance = c(
      sum((strata_info$N_h/N)^2 * stratum_estimates$var_vanilla),
      sum((strata_info$N_h/N)^2 * stratum_estimates$var_ratio_cont),
      sum((strata_info$N_h/N)^2 * stratum_estimates$var_ratio_bin)
    )
  )

  # Add standard errors and confidence intervals
  overall_est$se <- sqrt(overall_est$variance)
  overall_est$ci_lower <- overall_est$estimate - z_alpha_2 * overall_est$se
  overall_est$ci_upper <- overall_est$estimate + z_alpha_2 * overall_est$se

  list(
    overall_estimates = overall_est,
    stratum_estimates = stratum_estimates
  )
}

# Run stratified sampling analysis
results <- stratified_sampling(data, strata_info)

# Print stratified sampling results
cat("\nStratified Sampling Results:\n")
cat("\nOverall Estimates:\n")
print(data.frame(
  Method = results$overall_estimates$method,
  Estimate = round(results$overall_estimates$estimate, 4),
  SE = round(results$overall_estimates$se, 4),
  CI_Lower = round(results$overall_estimates$ci_lower, 4),
  CI_Upper = round(results$overall_estimates$ci_upper, 4)
))

cat("\nTrue Population Proportion:", round(mean(data$hypertension), 4), "\n")

cat("\nResults by Stratum:\n")
stratum_results <- results$stratum_estimates %>%
  select(stratum, N_h, n_h, vanilla, ratio_cont, ratio_bin) %>%
  mutate(across(where(is.numeric), ~round(., 4)))
print(stratum_results)
