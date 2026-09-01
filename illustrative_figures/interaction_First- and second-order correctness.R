################################################################################
# Figure: Omitted Variable Bias with interaction 
# first- and second-order incorrectness
################################################################################

library(lmtest)
# ---------------------------------------------------------------------------
# Data generation with omitted variable bias (OVB)
# ---------------------------------------------------------------------------

# Generate fixed predictor X
X <- rep(seq(-20,20,1), each = 1000)

# Generatefixed binary treatment variable D
D <- rep(c(0,1),length.out=length(X))

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
  
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 20)
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
  theme_minimal(base_size = 20) +
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
  theme_minimal(base_size = 20) +
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
  
  # Generate fixed predictor X
  X <- rep(seq(-20,20,1), each = 2)
  D <- rep(c(0,1),length.out=length(X))
  y.interaction.without.errors <- -3.5*X + 0*D + 7*D*X
  y <- y.interaction.without.errors + rnorm(length(D), 0, 10)
  
  data <- data.frame(x = X,
                     xz = X - mean(X),  # centered X (used later for covariance calculations)
                     d = D,
                     y = y)
  
  return(data)
}

# helper function to calculate standard errors under known residual cov matrix
kappasquare_lm <- function(
    COV_estimator,
    which_effects,
    X_A_work,
    R_list){
  
  kappa_square_lm_list <- 
    setNames(vector("list",length(COV_estimator)),
             COV_estimator) 
  
  # vt:= e_j'*(X_A'*X_A)_inv' * X_A' 
  X_A_workt <- t(X_A_work)
  X_At.X_A_inv <- solve(X_A_workt %*% X_A_work)
  
  for(iestimator in COV_estimator) {
    X_At.R_list.X_A <- X_A_workt %*% R_list[[iestimator]] %*% X_A_work
    kappasquare_lm_all_coeffs <- X_At.X_A_inv %*% X_At.R_list.X_A %*% X_At.X_A_inv
    
    kappa_square_lm_list[[iestimator]] <- kappasquare_lm_all_coeffs
  }
  
  
  return(kappa_square_lm_list)
  
}



# One simulation run
one_simulation_SKV <- function(){
  
  data = dgp_SKV()
  
  ## True model (correct specification)
  vcov_lm_true <- 
    kappasquare_lm(
      COV_estimator = "supplied",
      which_effects = which(colnames(as.matrix(model.matrix(lm(y ~ x + d + x * d, data = data))))=="x"),
      X_A_work = as.matrix(model.matrix(lm(y ~ x + d + x * d, data = data))),
      R_list = list(supplied = diag(x=100,ncol=nrow(data),nrow=nrow(data))  ))
  
  interaction_output <- coeftest(lm(y ~ x + d + x * d, data = data), 
                         vcov=function(x) as.matrix(vcov_lm_true[["supplied"]]))
  
  b_x_true <- interaction_output["x","Estimate"]
  stderr_x_true <- interaction_output["x", "Std. Error"]
  
  ## Working model (misspecified: omitted interaction)
  vcov_lm_linear <- 
    kappasquare_lm(
      COV_estimator = "supplied",
      which_effects = which(colnames(as.matrix(model.matrix(lm(y ~ x , data = data))))=="x"),
      X_A_work = as.matrix(model.matrix(lm(y ~ x , data = data))),
      R_list = list(supplied = diag(x=100,ncol=nrow(data),nrow=nrow(data))  ))
  
  linear_output <- coeftest(lm(y ~ x , data = data), 
                                 vcov=function(x) as.matrix(vcov_lm_linear[["supplied"]]))
  
  b_x_wm <- linear_output["x","Estimate"]
  stderr_x_wm <- linear_output["x", "Std. Error"]
  
  
  ## True model (estimated residual cov matrix)
  lm_true <- lm(y ~ x + d + x * d, data = data)
  b_x_true_est <- coefficients(lm_true)["x"]
  stderr_x_true_est <- coef(summary(lm_true))["x", "Std. Error"]
 
  ## Working model (estimated residual cov matrix)
  lm_wm <- lm(y ~ x, data = data)
  b_x_wm_est <- coefficients(lm_wm)["x"]
  stderr_x_wm_est <- coef(summary(lm_wm))["x", "Std. Error"]
  
  # Return key statistics
  return(c(Interaction = b_x_true, Linear = b_x_wm, 
           stderr_x_true = stderr_x_true, stderr_x_wm = stderr_x_wm,
           b_x_true_est = b_x_true_est,  b_x_wm_est = b_x_wm_est,
           stderr_x_true_est=stderr_x_true_est,
           stderr_x_wm_est=stderr_x_wm_est))
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
mean(SimulationStudy_SKV[,"stderr_x_true"])
sd(SimulationStudy_SKV[,"Interaction"])

mean(SimulationStudy_SKV[,"stderr_x_wm"])
sd(SimulationStudy_SKV[,"Linear"])
# correct estimation of standard error

# with estimated residual covariance matrix
mean(SimulationStudy_SKV[,"stderr_x_true_est"])
sd(SimulationStudy_SKV[,"b_x_true_est.x"])

mean(SimulationStudy_SKV[,"stderr_x_wm_est"]) 
sd(SimulationStudy_SKV[,"b_x_wm_est.x"])
# overestimation of empirical standard error


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
  theme_minimal(base_size = 20) +  
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 20)
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
  #SKV_plot /
  (interaction_residuals_plot + workingModel_residuals_plot)

# Save final figure
ggsave("correctness_assumptions_plots.png", 
       plot = correctness_assumptions_plots, 
       dpi = 300, width = 10, height = 5, units = "in")

