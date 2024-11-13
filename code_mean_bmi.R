library("dplyr")
stroke_data <- read.csv("healthcare-dataset-stroke-data.csv")
head(stroke_data)

stroke_data[stroke_data == "N/A"] <- NA
stroke_data_nona <- na.omit(stroke_data)

stroke_data_clean <- stroke_data_nona[stroke_data_nona$smoking_status != "Unknown", ]
stroke_data_clean$bmi <- as.numeric(stroke_data_clean$bmi)
stroke_data_clean$glucose_high <- ifelse(stroke_data_clean$avg_glucose_level >= 100, 1, 0)


# Check the dimensions before and after cleaning
print("Original dimensions:")
dim(stroke_data)
print("Dimensions after removing NA values and unknown values in smoking_status:")
dim(stroke_data_clean)

# Set seed for reproducibility
set.seed(123)

# Calculate initial parameters for sample size calculation
N <- nrow(stroke_data_clean)
S2_guess <- var(stroke_data_clean$bmi)  # Using population variance as our guess
z_alpha_2 <- qnorm(0.975)  # 95% CI, so alpha/2 = 0.025
delta <- 0.5  # Desired half-width of CI

# Calculate initial sample size (ignoring FPC)
n <- (z_alpha_2^2 * S2_guess) / delta^2

# Adjust for FPC
n_star <- n / (1 + n/N)

# Round up to get final sample size
n_final <- ceiling(n_star)

# FPC calculation
FPC <- 1 - n_final/N

# Print results
cat("Initial sample size (n):", n, "\n")
cat("Final sample size with FPC (n*):", n_final, "\n")

################################
###SRS##########################
################################
sample_indices <- sample(1:N, n_final)
sample <- stroke_data_clean[sample_indices, ]

# Calculate population means and proportions
pop_mean_glucose <- mean(stroke_data_clean$avg_glucose_level)
pop_proportion_glucose_high <- mean(stroke_data_clean$glucose_high)

# 1. Regression Estimator
reg_model <- lm(bmi ~ avg_glucose_level, data = sample)
b <- coef(reg_model)[2]

# Calculate regression estimate
y_reg <- sample_mean_bmi + b * (pop_mean_glucose - sample_mean_glucose)

# Calculate standard error for regression estimate
residuals <- residuals(reg_model)
var_residuals <- 1/(n_final - 1) * sum(residuals^2)
var_reg <- (1 - n_final/N) * var_residuals/n_final
se_reg <- sqrt(var_reg)

# Calculate confidence interval for regression estimator
ci_lower <- y_reg - z_alpha_2 * se_reg
ci_upper <- y_reg + z_alpha_2 * se_reg

# 2. Ratio Estimator with Continuous Glucose
# Ratio estimate
y_ratio_continuous <- mean(sample$bmi) / mean(sample$avg_glucose_level) * pop_mean_glucose

# Calculate variance of ratio estimator
residuals_continuous <- sample$bmi - mean(sample$bmi) / mean(sample$avg_glucose_level) * sample$avg_glucose_level
var_residuals_continuous <- 1/(n_final - 1) * sum(residuals_continuous^2)
var_ratio_continuous <- FPC * var_residuals_continuous/n_final
se_ratio_continuous <- sqrt(var_ratio_continuous)

# CI for ratio estimator (continuous)
ci_ratio_continuous <- c(y_ratio_continuous - z_alpha_2 * se_ratio_continuous,
                         y_ratio_continuous + z_alpha_2 * se_ratio_continuous)

# 3. Ratio Estimator with Binary Glucose
# Ratio estimate
y_ratio_binary <- mean(sample$bmi)/mean(sample$glucose_high) * pop_proportion_glucose_high

# Calculate variance of ratio estimator
residuals_binary <- sample$bmi - mean(sample$bmi)/mean(sample$glucose_high) * sample$glucose_high
var_residuals_binary <- 1/(n_final - 1) * sum(residuals_binary^2)
var_ratio_binary <- FPC * var_residuals_binary/n_final
se_ratio_binary <- sqrt(var_ratio_binary)

# CI for ratio estimator (binary)
ci_ratio_binary <- c(y_ratio_binary - z_alpha_2 * se_ratio_binary,
                     y_ratio_binary + z_alpha_2 * se_ratio_binary)

# Print results
cat("\nResults for Mean BMI Estimation:\n")
cat("\n1. Regression Estimation Results:\n")
cat("Estimated Mean BMI:", round(y_reg, 3), "\n")
cat("Standard Error:", round(se_reg, 3), "\n")
cat("95% CI: [", round(ci_lower, 3), ",", round(ci_upper, 3), "]\n")

cat("\n2. Ratio Estimator with Continuous Glucose Level:\n")
cat("Estimate:", round(y_ratio_continuous, 3), "\n")
cat("Standard Error:", round(se_ratio_continuous, 3), "\n")
cat("95% CI: [", round(ci_ratio_continuous[1], 3), ",",
    round(ci_ratio_continuous[2], 3), "]\n")

cat("\n3. Ratio Estimator with Binary Glucose Level:\n")
cat("Estimate:", round(y_ratio_binary, 3), "\n")
cat("Standard Error:", round(se_ratio_binary, 3), "\n")
cat("95% CI: [", round(ci_ratio_binary[1], 3), ",",
    round(ci_ratio_binary[2], 3), "]\n")

cat("\nTrue Population Mean BMI:", round(mean(stroke_data_clean$bmi), 3), "\n")

################################
###Stratified Sampling##########
################################

data <- data %>%
  mutate(age_strata = case_when(
    age < 20 ~ "0-19",
    age >= 20 & age < 30 ~ "20-29",
    age >= 30 & age < 40 ~ "30-39",
    age >= 40 & age < 50 ~ "40-49",
    age >= 50 & age < 60 ~ "50-59",
    age >= 60 & age < 70 ~ "60-69",
    age >= 70 & age < 80 ~ "70-79",
    age >= 80 ~ "80+"
  ))

# Convert to factor with ordered levels
data$age_strata <- factor(data$age_strata,
                          levels = c("0-19", "20-29", "30-39", "40-49",
                                     "50-59", "60-69", "70-79", "80+"))

# Calculate new strata info
strata_info <- data %>%
  group_by(age_strata) %>%
  summarise(
    N_h = n(),
    prop = n()/nrow(data)
  )

# Calculate optimal allocation
strata_info$n_h <- round(n_final * strata_info$N_h/nrow(data))

# Print new strata information
print(strata_info)

# Stratified sampling with three estimators
stratified_sampling <- function(data, strata_info, seed = 123) {
  set.seed(seed)
  stratum_samples <- list()
  stratum_estimates <- list()

  for(i in 1:nrow(strata_info)) {
    # Get stratum data
    stratum_data <- data %>%
      filter(age_strata == strata_info$age_strata[i])

    # Sample size for this stratum
    n_h <- strata_info$n_h[i]
    N_h <- strata_info$N_h[i]

    # Draw sample
    stratum_samples[[i]] <- stratum_data %>%
    slice_sample(n = n_h)

    # Get population means for auxiliary variables in this stratum
    pop_mean_glucose_h <- mean(stratum_data$avg_glucose_level)
    pop_mean_glucose_binary_h <- mean(stratum_data$avg_glucose_level >= 100)

    # Sample means
    sample_mean_bmi_h <- mean(stratum_samples[[i]]$bmi)
    sample_mean_glucose_h <- mean(stratum_samples[[i]]$avg_glucose_level)

    # 1. Regression Estimator
    reg_model_h <- lm(bmi ~ avg_glucose_level, data = stratum_samples[[i]])
    b_h <- coef(reg_model_h)[2]
    reg_est_h <- sample_mean_bmi_h + b_h * (pop_mean_glucose_h - sample_mean_glucose_h)

    # Variance of regression estimator
    residuals_reg_h <- residuals(reg_model_h)
    var_residuals_reg_h <- 1/(n_h - 1)*sum(residuals_reg_h^2)
    var_reg_h <- (1 - n_h/N_h) * var_residuals_reg_h / n_h

    # 2. Ratio Estimator (Continuous)
    R_cont_h <- mean(stratum_samples[[i]]$bmi) / mean(stratum_samples[[i]]$avg_glucose_level)
    ratio_est_cont_h <- R_cont_h * pop_mean_glucose_h

    # Variance of ratio estimator (continuous)
    residuals_ratio_cont_h <- stratum_samples[[i]]$bmi -
      R_cont_h * stratum_samples[[i]]$avg_glucose_level
    var_residuals_ratio_cont_h <- 1/(n_h -1) * sum(residuals_ratio_cont_h^2)
    var_ratio_cont_h <- (1 - n_h/N_h) * var_residuals_ratio_cont_h / n_h

    # 3. Ratio Estimator (Binary)
    sample_binary_h <- ifelse(stratum_samples[[i]]$avg_glucose_level >= 100, 1, 0)
    R_bin_h <- mean(stratum_samples[[i]]$bmi) / mean(sample_binary_h)
    ratio_est_bin_h <- R_bin_h * pop_mean_glucose_binary_h

    # Variance of ratio estimator (binary)
    residuals_ratio_bin_h <- stratum_samples[[i]]$bmi -
      R_bin_h * sample_binary_h
    var_residuals_ratio_bin_h <- 1/(n_h -1) * sum(residuals_ratio_bin_h^2)
    var_ratio_bin_h <- (1 - n_h/N_h) * var_residuals_ratio_bin_h / n_h

    # Store estimates
    stratum_estimates[[i]] <- data.frame(
      stratum = strata_info$age_strata[i],
      N_h = N_h,
      n_h = n_h,
      regression = reg_est_h,
      ratio_cont = ratio_est_cont_h,
      ratio_bin = ratio_est_bin_h,
      var_reg = var_reg_h,
      var_ratio_cont = var_ratio_cont_h,
      var_ratio_bin = var_ratio_bin_h
    )
  }

  # Combine stratum estimates
  stratum_estimates <- bind_rows(stratum_estimates)

  # Calculate overall estimates
  overall_est <- data.frame(
    method = c("Regression", "Ratio (Continuous)", "Ratio (Binary)"),
    estimate = c(
      sum(stratum_estimates$regression * strata_info$N_h) / sum(strata_info$N_h),
      sum(stratum_estimates$ratio_cont * strata_info$N_h) / sum(strata_info$N_h),
      sum(stratum_estimates$ratio_bin * strata_info$N_h) / sum(strata_info$N_h)
    ),
    variance = c(
      sum((strata_info$N_h/sum(strata_info$N_h))^2 * stratum_estimates$var_reg),
      sum((strata_info$N_h/sum(strata_info$N_h))^2 * stratum_estimates$var_ratio_cont),
      sum((strata_info$N_h/sum(strata_info$N_h))^2 * stratum_estimates$var_ratio_bin)
    )
  )

  # Add standard errors and confidence intervals
  overall_est$se <- sqrt(overall_est$variance)
  overall_est$ci_lower <- overall_est$estimate - 1.96 * overall_est$se
  overall_est$ci_upper <- overall_est$estimate + 1.96 * overall_est$se

  list(
    overall_estimates = overall_est,
    stratum_estimates = stratum_estimates
  )
}

# Run the analysis
results <- stratified_sampling(data, strata_info)

# Print results
cat("\nStratified Sampling Results:\n")
cat("\nOverall Estimates:\n")
print(data.frame(
  Method = results$overall_estimates$method,
  Estimate = round(results$overall_estimates$estimate, 3),
  SE = round(results$overall_estimates$se, 3),
  CI_Lower = round(results$overall_estimates$ci_lower, 3),
  CI_Upper = round(results$overall_estimates$ci_upper, 3)
))

cat("\nTrue Population Mean BMI:", round(mean(data$bmi), 3), "\n")

cat("\nResults by Stratum:\n")
stratum_results <- results$stratum_estimates %>%
  select(stratum, N_h, n_h, regression, ratio_cont, ratio_bin) %>%
  mutate(across(where(is.numeric), ~round(., 3)))
print(stratum_results)

################################################################################################################################
################################################################################################################################
################################################################################################################################

