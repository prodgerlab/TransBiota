## TransBiota Project - Transfeminine microbiota analysis
## Determine associations between microbial taxa and other parameters like CS, BC, and Symptoms
## Elaborated by: Jorge Rojas-Vargas


## Load libraries

# Explore microbiome
library(MASS)  # For Box-Cox transformation
library(lmerTest) # For the ANOVA Type III test
library(performance) # Statistics R²
library(gamlss) # Zero-inflated beta distribution


# Load the data
load("combined_data_cluster_w_TGC.RData")



## Trying 4 approaches for Linear Mixed Models
# lmer model (no transformation)
# lmer model (logit transformation)
# lmer model (Arcsine square root transformation)
# lmer model (Box-Cox transformation)
# Alternative: Zero-adjusted beta regression with gamlss



### Defining a function to check lmer model diagnostics
check_model_fit <- function(model, title) {
  # Shapiro-Wilk test
  res <- residuals(model)
  sw_test <- shapiro.test(res)
  
  # QQ plot
  par(mfrow=c(2,2))
  qqnorm(res, main=paste(title, "- Normal Q-Q Plot"))
  qqline(res)
  
  # Residuals vs fitted
  plot(fitted(model), res,
       main=paste(title, "- Residuals vs Fitted"),
       xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red", lty=2)
  
  # Print Shapiro-Wilk results
  print(sw_test)
  
  # Add R-squared
  print(r2(model))
}

### Defining a function to check the Zero-adjusted beta regression model diagnostics
check_beta_model_fit <- function(model, title) {
  # Get residuals (normalized quantile residuals for gamlss)
  res <- residuals(model)
  
  # Shapiro-Wilk test
  sw_test <- shapiro.test(res)
  
  # Create diagnostic plots
  par(mfrow=c(2,2))
  
  # QQ plot
  qqnorm(res, main=paste(title, "- Normal Q-Q Plot"))
  qqline(res)
  
  # Residuals vs fitted
  plot(fitted(model), res,
       main=paste(title, "- Residuals vs Fitted"),
       xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red", lty=2)
  
  # Density plot of residuals
  #plot(density(res), main=paste(title, "- Density of Residuals"))
  #curve(dnorm(x, mean=mean(res), sd=sd(res)), add=TRUE, col="red")
  
  # ACF plot
  acf(res, main=paste(title, "- ACF of Residuals"))
  
  # Print Shapiro-Wilk results
  print(sw_test)
  
  # Print model summary statistics
  summary(model_beta)
}



### Performing all the transformations and adding to the combined_data_cluster table


# List of TC variables to transform
tc_vars <- c("TC1", "TC2", "TC3", "TC4")

# Loop through each TC variable
for (tc in tc_vars) {
  # Adding a small constant slightly smaller than the smallest non-zero value
  epsilon <- min(combined_data_cluster_w_TGC[[tc]][combined_data_cluster_w_TGC[[tc]] > 0], na.rm = TRUE) / 2
  combined_data_cluster_w_TGC[[paste0(tc, "_logit")]] <- log((combined_data_cluster_w_TGC[[tc]] + epsilon) / 
                                                               (1 - combined_data_cluster_w_TGC[[tc]] + epsilon))
  
  # Arcsine square root transformation
  combined_data_cluster_w_TGC[[paste0(tc, "_asin")]] <- asin(sqrt(combined_data_cluster_w_TGC[[tc]]))
  
  # Box-Cox transformation
  # First shift away from zero
  shift <- min(combined_data_cluster_w_TGC[[tc]][combined_data_cluster_w_TGC[[tc]] > 0], na.rm = TRUE) / 2
  tc_shifted <- combined_data_cluster_w_TGC[[tc]] + shift
  
  # Find optimal lambda
  bc <- boxcox(tc_shifted ~ 1)
  lambda <- bc$x[which.max(bc$y)]
  
  # Apply transformation
  if (abs(lambda) > 1e-5) {  # If lambda ≠ 0
    combined_data_cluster_w_TGC[[paste0(tc, "_box")]] <- (tc_shifted^lambda - 1) / lambda
  } else {  # If lambda ≈ 0
    combined_data_cluster_w_TGC[[paste0(tc, "_box")]] <- log(tc_shifted)
  }
}



### Performing the different calculations


## --------------------------------------------------------------------------------
## Basic lmer model
## --------------------------------------------------------------------------------

# Perform the modeling
model_basic <- lmer(TC1 ~ BC + CS + (1|participant), 
                    data = combined_data_cluster_w_TGC)

check_model_fit(model_basic, "model_basic")

## --------------------------------------------------------------------------------
## lmer model with logit transformation:
## --------------------------------------------------------------------------------

## Trying the model with logit transformed data
model_logit <- lmer(TC1_logit ~ BC + CS + (1|participant),
                    data = combined_data_cluster_w_TGC)

check_model_fit(model_logit, "model_logit")


## --------------------------------------------------------------------------------
## Zero-adjusted beta regression approach:
## This might be more appropriate since we're dealing with proportions.
## --------------------------------------------------------------------------------

## BEZI is a zero-inflated beta distribution
model_beta <- gamlss(TC1 ~ BC + CS + random(factor(participant)),
                     family = BEZI,
                     data = na.omit(combined_data_cluster_w_TGC))

check_beta_model_fit(model_beta, "model_beta")


## --------------------------------------------------------------------------------
##
## Arcsine square root transformation:
## This is traditionally used for proportions:
##
## --------------------------------------------------------------------------------

## Arcsine square root transformation
model_asin <- lmer(TC1_asin ~ BC + CS + (1|participant),
                   data = combined_data_cluster_w_TGC)

check_model_fit(model_asin, "model_asin")


## --------------------------------------------------------------------------------
##
## Box-Cox transformation:
## We can find the optimal λ parameter:
##
## --------------------------------------------------------------------------------

## Box-Cox transformation:
model_box <- lmer(TC1_box ~ BC + CS + (1|participant),
                  data = combined_data_cluster_w_TGC)

check_model_fit(model_box, "model_box")

summary(model_box)

anova(model_box, type = 3)

## --------------------------------------------------------------------------------
##
## Analysis of the results
##
## --------------------------------------------------------------------------------

### FOR TC1

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.5882 (good normality)
## - Logit transformation: p = 0.00003445 (violates normality)
## - Beta regression: p = 0.2199 (good normality)
## - Arcsine transformation: p = 0.6991 (good normality)
## - Box-Cox transformation: p = 0.6376 (good normality)

## R² values (where available):
## - No transformation: Conditional = 0.472, Marginal = 0.144
## - Logit: Conditional = 0.741, Marginal = 0.144
## - Arcsine: Conditional = 0.544, Marginal = 0.143
## - Box-Cox: Conditional = 0.543, Marginal = 0.141

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Good residual normality (p = 0.2199)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale
## - Shows significant effects that other models missed

## 3. **Key Findings from the Beta Regression**:

## Significant effects:
## - BC_4: β = -0.866, p < 0.001 (strong negative effect)
## - CS_1: β = -0.415, p < 0.001 (strong negative effect)

## Non-significant:
## - BC_2: β = -0.141, p = 0.394
## - BC_3: β = -0.391, p = 0.384

## These coefficients are on the logit scale and show that:
## - BC_4 is strongly associated with lower TC1 values
## - CS_1 is strongly associated with lower TC1 values




### FOR TC2

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.007701 (violates normality)
## - Logit transformation: p = 0.0009265 (violates normality)
## - Beta regression: p = 0.9795 (excellent normality)
## - Arcsine transformation: p = 0.7035 (good normality)
## - Box-Cox transformation: p = 0.5431 (good normality)

## R² values (where available):
## - No transformation: Conditional = 0.701, Marginal = 0.081
## - Logit: Conditional = 0.721, Marginal = 0.072
## - Arcsine: Conditional = 0.718, Marginal = 0.080
## - Box-Cox: Conditional = 0.715, Marginal = 0.070

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Best residual normality (p = 0.9795)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale
## - Shows significant effects that other models missed

## 3. **Key Findings from the Beta Regression**:

## Significant effects:
## - BC_2: β = 0.391, p = 0.014 (positive effect)
## - BC_3: β = 0.491, p = 0.003 (positive effect)
## - BC_4: β = 0.683, p < 0.001 (strong positive effect)
## - CS_1: β = -0.492, p < 0.001 (strong negative effect)
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## These coefficients are on the logit scale and show that:
## - BC_2 is associated with higher TC2 values
## - BC_3 is associated with higher TC2 values
## - BC_4 is strongly associated with higher TC2 values
## - CS_1 is strongly associated with lower TC2 values



### FOR TC3

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 8.772e-13 (violates normality)
## - Logit transformation: p = 0.0002133 (violates normality)
## - Beta regression: p = 0.002048 (violates normality)
## - Arcsine transformation: p = 0.000000001533 (violates normality)
## - Box-Cox transformation: p = 0.03343 (violates normality)

## R² values (where available):
## - No transformation: Conditional = 0.504, Marginal = 0.076
## - Logit: Conditional = 0.588, Marginal = 0.090
## - Arcsine: Conditional = 0.565, Marginal = 0.088
## - Box-Cox: Conditional = 0.645, Marginal = 0.097

## 2. **Best Model Selection**:
## - Keep the Beta regression

## 3. **Key Findings from the Beta Regression**:

## Significant effects:
## - BC_4: β = 0.880, p < 0.001 (strong positive effect)
## - CS_1: β = 0.389, p < 0.001 (strong positive effect)

## **Key Findings from the Box-Cox transformation**:

## Non-significant:
## - BC_2: β = 0.289, p = 0.209
## - BC_3: β = 0.009, p = 0.970

## These coefficients are on the logit scale and show that:
## - BC_4 is associated with higher TC1 values
## - CS_1 is associated with higher TC1 values


### FOR TC4

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.005082 (violates normality)
## - Logit transformation: p = 5.924e-11 (violates normality)
## - Beta regression: p = 0.4291 (excellent normality)
## - Arcsine transformation: p = 0.00005459 (violates normality)
## - Box-Cox transformation: p = 0.000437 (violates normality)

## R² values (where available):
## - No transformation: Conditional = 0.668, Marginal = 0.142
## - Logit: Conditional = 0.472, Marginal = 0.144
## - Arcsine: Conditional = 0.620, Marginal = 0.149
## - Box-Cox: Conditional = 0.630, Marginal = 0.146

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Best residual normality (p = 0.4291)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale
## - Shows significant effects that other models missed

## 3. **Key Findings from the Beta Regression**:

## Significant effects:
## - BC_3: β = -0.306, p = 0.050 (marginal effect)
## - BC_4: β = -0.634, p < 0.001 (strong negative effect)
## - CS_1: β = 0.447, p < 0.001 (strong positive effect)

## Non-significant:
## - BC_2: β = -0.265, p = 0.072

## These coefficients are on the logit scale and show that:
## - BC_2 is marginal associated with lower TC4 values
## - BC_3 is marginal associated with lower TC4 values
## - BC_4 is strongly associated with lower TC4 values
## - CS_1 is strongly associated with higher TC4 values





##
## Cross-validation for beta regression across TC1, TC2, TC3, TC4
##

# Set up cross-validation
set.seed(123)  # for reproducibility
n <- nrow(combined_data_cluster_w_TGC)
k <- 5  # number of folds

# Create indices for folds
fold_indices <- cut(seq(1, n), breaks = k, labels = FALSE)
fold_indices <- sample(fold_indices)  # randomize

# Initialize a list to store results for each TC
tc_cv_results <- list()

# Loop through TC1, TC2, TC3, and TC4
for (tc in c("TC1", "TC2", "TC3", "TC4")) {
  cat("\nCross-validation for", tc, ":\n")
  
  # Initialize results storage for this TC
  cv_results <- matrix(NA, nrow = k, ncol = 2)
  colnames(cv_results) <- c("RMSE", "MAE")
  
  # Perform k-fold cross-validation
  for (i in 1:k) {
    # Create training and test sets
    test_indices <- which(fold_indices == i)
    train_indices <- which(fold_indices != i)
    
    train_df <- combined_data_cluster_w_TGC[train_indices, ]
    test_df <- combined_data_cluster_w_TGC[test_indices, ]
    
    # Fit model
    tryCatch({
      model_train <- gamlss(as.formula(paste(tc, "~ BC + CS + random(factor(participant))")),
                            family = BEZI,
                            data = na.omit(train_df))
      
      # Make predictions
      pred <- predict(model_train, newdata = test_df, type = "response")
      
      # Calculate error metrics
      rmse <- sqrt(mean((test_df[[tc]] - pred)^2, na.rm = TRUE))
      mae <- mean(abs(test_df[[tc]] - pred), na.rm = TRUE)
      
      # Store results
      cv_results[i, ] <- c(rmse, mae)
    }, error = function(e) {
      cat("Error in fold", i, "for", tc, ":", conditionMessage(e), "\n")
    })
  }
  
  # Calculate summary statistics for this TC
  cv_summary <- data.frame(
    Metric = c("RMSE", "MAE"),
    Mean = c(mean(cv_results[, 1], na.rm = TRUE),
             mean(cv_results[, 2], na.rm = TRUE)),
    SD = c(sd(cv_results[, 1], na.rm = TRUE),
           sd(cv_results[, 2], na.rm = TRUE))
  )
  
  # Round to 4 decimal places
  cv_summary$Mean <- round(cv_summary$Mean, 4)
  cv_summary$SD <- round(cv_summary$SD, 4)
  
  # Store results for this TC
  tc_cv_results[[tc]] <- list(Summary = cv_summary, ByFold = cv_results)
  
  # Print results
  cat("\nCross-validation results for", tc, ":\n")
  print(cv_summary)
  
  cat("\nResults by fold for", tc, ":\n")
  rownames(cv_results) <- paste("Fold", 1:k)
  print(round(cv_results, 4))
}

# Access results for a specific TC, e.g., TC2:
tc_cv_results[["TC1"]]
tc_cv_results[["TC2"]]
tc_cv_results[["TC3"]]
tc_cv_results[["TC4"]]



##
## Odds analysis for beta regression across TC1, TC2, TC3, TC4
##

## Converting logit coefficients to odds ratios
# Extract fixed-effect coefficients for the mu parameter
fixed_coeffs <- coef(model_beta, what = "mu")  # Extract the mu coefficients
# Extract the variance-covariance matrix for the mu parameter
vcov_mu <- vcov(model_beta, type = "vcov", what = "mu")
# Compute standard errors by taking the square root of the diagonal of the variance-covariance matrix
mu_se <- sqrt(diag(vcov_mu))
# Extract only the fixed effects coefficients and SEs
fixed_coeffs <- model_beta$mu.coefficients[1:5]  # Excluding the random effect
fixed_se <- mu_se[1:5]  # Using only the first 5 standard errors

# Calculate odds ratios and CIs
odds_ratios <- exp(fixed_coeffs)
ci_lower <- exp(fixed_coeffs - 1.96 * fixed_se)
ci_upper <- exp(fixed_coeffs + 1.96 * fixed_se)
p_values <- 2 * (1 - pnorm(abs(fixed_coeffs/fixed_se)))

# Create odds ratio table
odds_table <- data.frame(
  Term = names(fixed_coeffs),
  Coefficient = fixed_coeffs,
  OR = odds_ratios,
  CI_lower = ci_lower,
  CI_upper = ci_upper,
  P_value = p_values,
  row.names = NULL
)

# Round numeric columns to 3 decimal places
odds_table[, 2:6] <- round(odds_table[, 2:6], 3)

print(odds_table)








## --------------------------------------------------------------------------------
##
## FOR SYMPTOMS
##
## --------------------------------------------------------------------------------


# Prepare data for comparing Symptoms individually (repeated samples)

# Define the variables of interest and the unique symptom levels
symptom_levels <- c("Bleeding", "Discharge", "Itching_burning", "Malodour", "Pain", "None")

# Create an empty data frame to store the large subsample
large_subsample <- data.frame()

# Outer loop to iterate over each symptom level
for (level in symptom_levels) {
  # Create a subset of the table with rows that have 1 for the current level
  subsample <- combined_data_cluster_w_TGC[combined_data_cluster_w_TGC[[level]] == 1, ]
  
  # Add a column to indicate the symptom group
  subsample$Symptom_Group <- level
  
  # Append the subsample to the large data frame
  large_subsample <- rbind(large_subsample, subsample)
}

# Ensure that 'Symptom' is a factor and set "None" as the reference
large_subsample$Symptom_Group <- factor(large_subsample$Symptom_Group, 
                                        levels = c("None", "Bleeding", "Discharge", 
                                                   "Itching_burning", "Malodour", "Pain"))


### Performing all the transformations and adding to the large_subsample table

# List of TC variables to transform
tc_vars <- c("TC1", "TC2", "TC3", "TC4")


# Loop through each TC variable
for (tc in tc_vars) {
  # Adding a small constant slightly smaller than the smallest non-zero value
  epsilon <- min(subsample_grande[[tc]][subsample_grande[[tc]] > 0], na.rm = TRUE) / 2
  subsample_grande[[paste0(tc, "_logit")]] <- log((subsample_grande[[tc]] + epsilon) / 
                                                    (1 - subsample_grande[[tc]] + epsilon))
  
  # Arcsine square root transformation
  subsample_grande[[paste0(tc, "_asin")]] <- asin(sqrt(subsample_grande[[tc]]))
  
  # Box-Cox transformation
  # First shift away from zero
  shift <- min(subsample_grande[[tc]][subsample_grande[[tc]] > 0], na.rm = TRUE) / 2
  tc_shifted <- subsample_grande[[tc]] + shift
  
  # Find optimal lambda
  bc <- boxcox(tc_shifted ~ 1)
  lambda <- bc$x[which.max(bc$y)]
  
  # Apply transformation
  if (abs(lambda) > 1e-5) {  # If lambda ≠ 0
    subsample_grande[[paste0(tc, "_box")]] <- (tc_shifted^lambda - 1) / lambda
  } else {  # If lambda ≈ 0
    subsample_grande[[paste0(tc, "_box")]] <- log(tc_shifted)
  }
}




### Performing the different calculations

## --------------------------------------------------------------------------------
## Basic lmer model
## --------------------------------------------------------------------------------

# Perform the modelling
model_basic <- lmer(TC1 ~ Symptom_Group + (1|participant), 
                    data = subsample_grande)

check_model_fit(model_basic, "model_basic")

## --------------------------------------------------------------------------------
## lmer model with logit transformation:
## --------------------------------------------------------------------------------

## Trying the model with logit transformed data
model_logit <- lmer(TC1_logit ~ Symptom_Group + (1|participant),
                    data = subsample_grande)

check_model_fit(model_logit, "model_logit")


## --------------------------------------------------------------------------------
## Zero-adjusted beta regression approach:
## --------------------------------------------------------------------------------

## BEZI is zero-inflated beta distribution
library(gamlss)
model_beta <- gamlss(TC1 ~ Symptom_Group + random(factor(participant)),
                     family = BEZI,
                     data = na.omit(subsample_grande))

check_beta_model_fit(model_beta, "model_beta")


## --------------------------------------------------------------------------------
## Arcsine square root transformation:
## --------------------------------------------------------------------------------

## Arcsine square root transformation
model_asin <- lmer(TC1_asin ~ Symptom_Group + (1|participant),
                   data = subsample_grande)

check_model_fit(model_asin, "model_asin")


## --------------------------------------------------------------------------------
## Box-Cox transformation:
## --------------------------------------------------------------------------------

## Box-Cox transformation:
model_box <- lmer(TC1_box ~ Symptom_Group + (1|participant),
                  data = subsample_grande)

check_model_fit(model_box, "model_box")

summary(model_box)

anova(model_box, type = 3)


## --------------------------------------------------------------------------------
## Analysis of the results
## --------------------------------------------------------------------------------

### FOR TC1

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.2893 (good normality)
## - Logit transformation: p = 0.000001931 (violates normality)
## - Beta regression: p = 0.4494 (good normality)
## - Arcsine transformation: p = 0.1275 (good normality)
## - Box-Cox transformation: p = 0.08825 (violates normality)

## R² values (where available):
## - No transformation: Conditional = 0.458, Marginal = 0.006
## - Logit: Conditional = 0.725, Marginal = 0.002
## - Arcsine: Conditional = 0.524, Marginal = 0.003
## - Box-Cox: Conditional = 0.523, Marginal = 0.003

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Best residual normality (p = 0.4494)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale

## 3. **Key Findings from the Beta Regression**:

## Significant effects:

## Non-significant:
## - Bleeding: β = -0.075, p = 0.704
## - Discharge: β = -0.150, p = 0.567
## - Itching: β = -0.012, p = 0.959
## - Malodour: β = -0.040, p = 0.790
## - Pain: β = -0.313, p = 0.561

## These coefficients are on the logit scale and show that:




### FOR TC2

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.01882 (violates normality)
## - Logit transformation: p = 0.001078 (violates normality)
## - Beta regression: p = 0.9642 (excellent normality)
## - Arcsine transformation: p = 0.8053 (good normality)
## - Box-Cox transformation: p = 0.5369 (violates normality)

## R² values (where available):
## - No transformation: Conditional = 0.674, Marginal = 0.023
## - Logit: Conditional = 0.707, Marginal = 0.026
## - Arcsine: Conditional = 0.693, Marginal = 0.025
## - Box-Cox: Conditional = 0.696, Marginal = 0.026

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Best residual normality (p = 0.9642)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale

## 3. **Key Findings from the Beta Regression**:

## Significant effects:
## - Malodour: β = -0.434, p = 0.00485 (strong negative effect)

## Non-significant:
## - Bleeding: β = 0.254, p = 0.0759 (marginal positive effect)
## - Discharge: β = -0.016, p = 0.930
## - Itching: β = -0.056, p = 0.785
## - Pain: β = -0.403, p = 0.296

#> p_values_fdr
#[1] 0.0242500 0.1897500 0.9300000 0.9300000 0.4933333

## These coefficients are on the logit scale and show that:
## - Malodour is strongly associated with lower TC2 values
## - Bleeding is marginally associated with higher TC2 values, NO after fdr adjustment



### FOR TC3

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 2.599e-13 (violates normality)
## - Logit transformation: p = 0.000002948 (violates normality)
## - Beta regression: p = 0.002107 (violates normality)
## - Arcsine transformation: p = 0.0000000001182 (violates normality)
## - Box-Cox transformation: p = 0.001149 (violates normality)

## R² values (where available):
## - No transformation: Conditional = 0.487, Marginal = 0.004
## - Logit: Conditional = 0.571, Marginal = 0.009
## - Arcsine: Conditional = 0.543, Marginal = 0.002
## - Box-Cox: Conditional = 0.637, Marginal = 0.011

## 2. **Best Model Selection**:
## Let's try BENZI and Box-Cox

## 3. **Key Findings from the Beta Regression**:

## No significant effects

## **Key Findings from the Box-Cox transformation**:

## No significant effects

## These coefficients are on the logit scale and show that:
## Symptoms are not associated with TC3 values


### FOR TC4

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.2004 (violates normality)
## - Logit transformation: p = 3.556e-11 (violates normality)
## - Beta regression: p = 0.6635 (excellent normality)
## - Arcsine transformation: p = 0.001729 (violates normality)
## - Box-Cox transformation: p = 0.02395 (violates normality)

## R² values (where available):
## - No transformation: Conditional = 0.660, Marginal = 0.042
## - Logit: Conditional = 0.497, Marginal = 0.064
## - Arcsine: Conditional = 0.624, Marginal = 0.051
## - Box-Cox: Conditional = 0.633, Marginal = 0.050

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Best residual normality (p = 0.6635)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale
## - Shows significant effects that other models missed

## 3. **Key Findings from the Beta Regression**:

## Significant effects:
## - Malodour: β = 0.46366, p = 0.000204 (strong positive effect)
## - Pain: β = 1.15195, p = 0.000811 (strong positive effect)

## Non-significant:
## - Bleeding: β = -0.041, p = 0.759
## - Discharge: β = 0.0765, p = 0.705 
## - Itching: β = 0.046, p = 0.824


## These coefficients are on the logit scale and show that:
## - Malodour is strongly associated with higher TC4 values
## - Pain is strongly associated with higher TC4 values






## --------------------------------------------------------------------------------
##
## FOR AGE, TSVy and PH
##
## --------------------------------------------------------------------------------



## --------------------------------------------------------------------------------
## Basic lmer model
## --------------------------------------------------------------------------------

# Perform the modelling
model_basic <- lmer(TC1 ~ age + TSVy + pH + (1|participant), 
                    data = combined_data_cluster_w_TGC)

check_model_fit(model_basic, "model_basic")

## --------------------------------------------------------------------------------
## lmer model with logit transformation:
## --------------------------------------------------------------------------------

## Trying the model with logit transformed data
model_logit <- lmer(TC1_logit ~ age + TSVy + pH + (1|participant),
                    data = combined_data_cluster_w_TGC)

check_model_fit(model_logit, "model_logit")


## --------------------------------------------------------------------------------
## Zero-adjusted beta regression approach:
## This might be more appropriate since we're dealing with proportions.
## --------------------------------------------------------------------------------

## BEZI is zero-inflated beta distribution
model_beta <- gamlss(TC1 ~ age + TSVy + pH + random(factor(participant)),
                     family = BEZI,
                     data = na.omit(combined_data_cluster_w_TGC))

check_beta_model_fit(model_beta, "model_beta")



## --------------------------------------------------------------------------------
##
## Arcsine square root transformation:
## This is traditionally used for proportions:
##
## --------------------------------------------------------------------------------

## Arcsine square root transformation
model_asin <- lmer(TC1_asin ~ age + TSVy + pH + (1|participant),
                   data = combined_data_cluster_w_TGC)

check_model_fit(model_asin, "model_asin")


## --------------------------------------------------------------------------------
##
## Box-Cox transformation:
## We can find the optimal λ parameter:
##
## --------------------------------------------------------------------------------

## Box-Cox transformation:
model_box <- lmer(TC1_box ~ age + TSVy + pH + (1|participant),
                  data = combined_data_cluster_w_TGC)

check_model_fit(model_box, "model_box")

summary(model_box)

anova(model_box, type = 3)

## --------------------------------------------------------------------------------
##
## Analysis of the results
##
## --------------------------------------------------------------------------------

### FOR TC1

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.643 (good normality)
## - Logit transformation: p = 0.0000004452 (violates normality)
## - Beta regression: p = 0.8557 (good normality)
## - Arcsine transformation: p = 0.7437 (good normality)
## - Box-Cox transformation: p = 0.7007 (good normality)

## R² values (where available):
## - No transformation: Conditional = 0.487, Marginal = 0.007
## - Logit: Conditional = 0.717, Marginal = 0.009
## - Arcsine: Conditional = 0.544, Marginal = 0.008
## - Box-Cox: Conditional = 0.542, Marginal = 0.008

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Good residual normality (p = 0.8557)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale
## - Shows significant effects that other models missed

## 3. **Key Findings from the Beta Regression**:

## Non-significant:
## - age: β = 0.0037, p = 0.4706
## - TSVy: β = -0.0135, p = 0.3330
## - pH: β = 0.1324, p = 0.0629 (marginally associated with higher TC1 values)



### FOR TC2

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.006567 (violates normality)
## - Logit transformation: p = 0.04432 (violates normality)
## - Beta regression: p = 0.8426 (excellent normality)
## - Arcsine transformation: p = 0.5363 (good normality)
## - Box-Cox transformation: p = 0.4943 (good normality)

## R² values (where available):
## - No transformation: Conditional = 0.691, Marginal = 0.002
## - Logit: Conditional = 0.753, Marginal = 0.005
## - Arcsine: Conditional = 0.713, Marginal = 0.002
## - Box-Cox: Conditional = 0.722, Marginal = 0.003

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Best residual normality (p = 0.8426)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale
## - Shows significant effects that other models missed

## 3. **Key Findings from the Beta Regression**:

## Non-significant:
## - age: β = -0.0047, p = 0.2857
## - TSVy: β = 0.0105, p = 0.3398
## - pH: β = -0.0263, p = 0.6608


### FOR TC3

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 1.511e-12 (violates normality)
## - Logit transformation: p = 0.0001147 (violates normality)
## - Beta regression: p = 0.00375 (excellent normality)
## - Arcsine transformation: p = 0.000000001994 (good normality)
## - Box-Cox transformation: p = 0.02174 (good normality)

## R² values (where available):
## - No transformation: Conditional = 0.465, Marginal = 0.022
## - Logit: Conditional = 0.579, Marginal = 0.004
## - Arcsine: Conditional = 0.713, Marginal = 0.002
## - Box-Cox: Conditional = 0.643, Marginal = 0.004

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Best residual normality
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale
## - Shows significant effects that other models missed

## 3. **Key Findings from the Beta Regression**:

## Non-significant:
## - age: β = -0.0018, p = 0.7899
## - TSVy: β = 0.0092, p = 0.5849
## - pH: β = -0.1427, p = 0.9198




### FOR TC4

## 1. **Model Diagnostics Comparison**:

## Normality (Shapiro-Wilk test p-values):
## - No transformation: p = 0.02434 (violates normality)
## - Logit transformation: p = 5.638e-11 (violates normality)
## - Beta regression: p = 0.3688 (excellent normality)
## - Arcsine transformation: p = 0.0001176 (violates normality)
## - Box-Cox transformation: p = 0.0009294 (violates normality)

## R² values (where available):
## - No transformation: Conditional = 0.668, Marginal = 0.013
## - Logit: Conditional = 0.441, Marginal = 0.021
## - Arcsine: Conditional = 0.613, Marginal = 0.015
## - Box-Cox: Conditional = 0.623, Marginal = 0.013

## 2. **Best Model Selection**:
## The zero-inflated beta regression (BEZI) appears to be the most appropriate model because:
## - Best residual normality (p = 0.3688)
## - Handles zero values appropriately
## - Directly models proportions
## - Provides meaningful coefficients on the logit scale
## - Shows significant effects that other models missed

## 3. **Key Findings from the Beta Regression**:

## Significant effects:
## - pH: β = 0.134, p = 0.0263 (positive effect)

## Non-significant:
## - age: β = -0.002464, p = 0.5812
## - TSVy: β = -0.01208, p = 0.3072



