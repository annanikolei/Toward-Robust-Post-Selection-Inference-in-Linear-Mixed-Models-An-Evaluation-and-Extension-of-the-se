################################################################################
# Figure: Omitted Variable Bias with interaction 
# first- and second-order incorrectness
################################################################################

# ---------------------------------------------------------------------------
# Data generation with omitted variable bias (OVB)
# ---------------------------------------------------------------------------

# Generate continuous predictor X
X <- rnorm(n = 100000, mean = 0, sd = 10)

# Generate binary treatment variable D
D <- rbinom(n = length(X), 1, 0.5)

# True data-generating process (includes interaction term X*D)
# y = -3.5*X + 7*D*X  (no main effect of D)
y.interaction.without.errors <- -3.5*X + 0*D + 7*D*X

# Theoretical explanation of OVB:
# Cov(X, X*D) = E(X)*Cov(X,D) + E(D)*Var(X) = 0.5 * 100
# Bias in linear model slope:
# beta1 = -beta3 * Cov(X, X*D) / Var(X) = -7 * 50 / 100 = -3.5

# Add noise to the outcome variable
y.noisy <- y.interaction.without.errors + rnorm(length(D), 0, 10)

# Create datasets (true and noisy)
data <- data.frame(x = X,
                   y = y.interaction.without.errors)

noisy_data <- data.frame(x = X, y = y.noisy)

# Fit misspecified linear model (omitting interaction term)
approx_model <- lm(y.interaction.without.errors ~ X)

# ---------------------------------------------------------------------------
# Plot: First-order incorrectness (mean structure misspecification)
# ---------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

# Visualization of true model vs linear approximation
first_order_assumption_plot <- ggplot(data, aes(x = x, y = y)) +
  
  # True data (invisible / background)
  geom_point(color = "white", size = 0.1) +   
  
  # Noisy observed data
  geom_point(data = noisy_data, aes(x = x, y = y), alpha = 0.3, size = 0.5, color = "gray") +  
  
  # True regression line for D = 0
  geom_abline(intercept = 0,
              slope = -3.5, color = "mediumpurple", size = 1) +
  
  # True regression line for D = 1 (slope changes due to interaction)
  geom_abline(intercept = 0,
              slope = +7 + -3.5, color = "mediumpurple", size = 1) +
  
  # Labels for the two groups
  annotate("text", x = 19, y = 45, label = "d = 1", color = "mediumpurple", hjust = 0) +
  annotate("text", x = 19, y = -95, label = "d = 0", color = "mediumpurple", hjust = 0) +
  
  # Restrict x-axis for better visualization
  coord_cartesian(xlim = c(-15, 15)) +
  
  # Linear approximation (misspecified model)
  geom_abline(intercept = coefficients(approx_model)[1],
              slope = coefficients(approx_model)[2], color = "red", size = 1) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 13)
  ) +
  
  labs(title = "",
       x = "x",
       y = "y")

# ---------------------------------------------------------------------------
# Residual plots: Second-order incorrectness (variance structure)
# ---------------------------------------------------------------------------

# Residuals from true interaction model vs noisy data
residuals.interaction <- y.interaction.without.errors - y.noisy
data.residuals.interaction <- data.frame(x = X, residuals = residuals.interaction)

# Predicted values from misspecified linear model
y.linear <- coefficients(approx_model)[1] + coefficients(approx_model)[2] * X

# Residuals from misspecified model
residuals.linear <- y.linear - y.noisy

# Create dataframe
data.residuals.linear <- data.frame(x = X, residuals = residuals.linear)

library(ggplot2)
library(patchwork)

# Common y-axis limits for comparability
ymin <- min(c(data.residuals.interaction$residuals,
              data.residuals.linear$residuals))
ymax <- max(c(data.residuals.interaction$residuals,
              data.residuals.linear$residuals))

# Residuals for correctly specified model
interaction_residuals_plot <-
  ggplot(data.residuals.interaction, aes(x = x, y = residuals)) +
  geom_point(alpha = 0.3, color = "mediumpurple") +
  geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(ymin, ymax), xlim = c(-15, 15)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = "x", y = "Residuals")

# Residuals for misspecified linear model
workingModel_residuals_plot <-
  ggplot(data.residuals.linear, aes(x = x, y = residuals)) +
  geom_point(alpha = 0.3, color = "red") +
  geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(ymin, ymax), xlim = c(-15, 15)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(x = "x", y = "Residuals")

# Combine residual plots
combined_plot <- interaction_residuals_plot + workingModel_residuals_plot
combined_plot

# ---------------------------------------------------------------------------
# Simulation study: Sampling distributions
# ---------------------------------------------------------------------------

# Data-generating process for simulation
dgp_SKV <- function(){
  
  X <- rnorm(n = 70, mean = 0, sd = 10)
  D <- rbinom(n = length(X), 1, 0.5)
  y.interaction.without.errors <- -3.5*X + 0*D + 7*D*X
  y <- y.interaction.without.errors + rnorm(length(D), 0, 10)
  
  data <- data.frame(x = X,
                     xz = X - mean(X),  # centered X (used later for covariance calculations)
                     d = D,
                     y = y)
  
  return(data)
}

# One simulation run
one_simulation_SKV <- function(){
  
  data = dgp_SKV()
  
  ## True model (correct specification)
  lm_true <- lm(y ~ x + d + x * d, data = data)
  b_x_true <- coefficients(lm_true)["x"]
  stderr_x_true <- coef(summary(lm_true))["x", "Std. Error"]
  res_stderror_true <- summary(lm_true)$sigma
  
  ## Working model (misspecified: omitted interaction)
  lm_wm <- lm(y ~ x, data = data)
  b_x_wm <- coefficients(lm_wm)["x"]
  stderr_x_wm <- coef(summary(lm_wm))["x", "Std. Error"]
  res_stderror_wm <- summary(lm_wm)$sigma
  
  # Return key statistics
  return(c(Interaction = b_x_true, Linear = b_x_wm, 
           stderr_x_true = stderr_x_true, stderr_x_wm = stderr_x_wm,
           res_stderror_true = res_stderror_true, res_stderror_wm = res_stderror_wm,
           sd.residuals.linear = sd(data$y),
           cov.residuals.linear.squared = cov(data$y^2, data$xz^2)))
}

# Run simulation study
simulation_study_SKV <- function(n_iter = 10000){
  
  results <- lapply(1:n_iter, function(i_iter){
    one_simulation_SKV()
  })
  
  df_results <- do.call(rbind, results)
  colnames(df_results)[c(1,2)] <- c("Interaction","Linear")
  
  return(df_results)
}

# Execute large simulation
SimulationStudy_SKV <- simulation_study_SKV(100000)

# Summary statistics (used for theoretical comparisons)
mean(SimulationStudy_SKV[,"res_stderror_true"])
mean(SimulationStudy_SKV[,"stderr_x_true"])
sd(SimulationStudy_SKV[,"Interaction"])

mean(SimulationStudy_SKV[,"sd.residuals.linear"]) 
mean(SimulationStudy_SKV[,"res_stderror_wm"])

mean_cov_e_xz <- mean(SimulationStudy_SKV[,"cov.residuals.linear.squared"])

sd(SimulationStudy_SKV[,"Linear"])

# Theoretical variance approximation (adjusted)
sqrt(1/162 *((mean(SimulationStudy_SKV[,"res_stderror_wm"])^2 * 
                var(rep(seq(-20,20,0.5), each = 2)) + mean_cov_e_xz) /
               var(rep(seq(-20,20,0.5), each = 2))^2))

# Standard formula under (incorrect) assumptions
mean(SimulationStudy_SKV[,"stderr_x_wm"])
sqrt(1/162 *((mean(SimulationStudy_SKV[,"res_stderror_wm"])^2 / 
                var(rep(seq(-20,20,0.5), each = 2)))))

# ---------------------------------------------------------------------------
# Plot sampling distributions
# ---------------------------------------------------------------------------

library(ggplot2)
library(tidyr)

# Convert results to dataframe
df <- as.data.frame(SimulationStudy_SKV)[,c(1,2)]

# Reshape to long format for ggplot
df_long <- pivot_longer(df, cols = everything(),
                        names_to = "Model",
                        values_to = "Estimated_Regression_Weight_for_x")

# Histogram of sampling distributions
SKV_plot <- ggplot(df_long, aes(x = Estimated_Regression_Weight_for_x, fill = Model)) +
  geom_histogram(alpha = 0.5, 
                 bins = 70, 
                 position = "identity") +
  theme_minimal(base_size = 14) +  
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 13)
  ) +
  labs(
    x = "Estimated Regression Weight for x",
    y = "Frequency") +
  scale_fill_manual(
    values = c("Interaction" = "mediumpurple", "Linear" = "red")
  )

# ---------------------------------------------------------------------------
# Combine all plots
# ---------------------------------------------------------------------------

correctness_assumptions_plots <- 
  first_order_assumption_plot / 
  SKV_plot /
  (interaction_residuals_plot + workingModel_residuals_plot)

# Save final figure
ggsave("correctness_assumptions_plots.png", 
       plot = correctness_assumptions_plots, 
       dpi = 900, width = 8, height = 8, units = "in")