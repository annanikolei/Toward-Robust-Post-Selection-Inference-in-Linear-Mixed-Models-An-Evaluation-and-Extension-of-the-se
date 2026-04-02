################################################################################
# Figure: data-generating model and working model in Simulation Study 2
# quadratic relationship vs linear relationship
################################################################################

library(lme4)
library(whitening)

# Total number of observations (balanced design: m groups with equal size)
n_obs = 40^2 
m = 40  # number of level-2 units (groups)

## Create level-2 ID variable (grouping variable)
l2.id <- rep(1:m, each = n_obs/m)

## Number of fixed effects (number of predictors)
p <- 4 # X1, X2, X3, X4

## Generation of X -----------------------------------------------------------

# Generate half of the data from multivariate normal distribution
X.l1_per_l2_half <- mvtnorm::rmvnorm(n = n_obs/m/2, 
                                     mean = rep(0, p), 
                                     sigma = diag(p),
                                     checkSymmetry = F)

# Create symmetric data around the mean to ensure balanced structure
# This helps enforce symmetry and reduces skewness
X.l1_per_l2 <- rbind(X.l1_per_l2_half, 
                     -X.l1_per_l2_half + mean(X.l1_per_l2_half)) 

# Whitening transformation:
# removes correlations and standardizes variance across predictors
X.l1_per_l2.white <- whiten(X.l1_per_l2, center = TRUE) 

# Replicate level-1 data for all level-2 units
X.l1 <- do.call(rbind, replicate(m, X.l1_per_l2.white, simplify = FALSE))

# Standardize variables to variance = 1 and shift mean to 2
# Result: E(Xj) = 2 and Var(Xj) = 1 for all predictors
X.l1 <- X.l1 / apply(X.l1, 2, sd)
X.l1 <- X.l1 + 2 - colMeans(X.l1)

# Add intercept column (column of ones)
X <- cbind(1, X.l1)
colnames(X) <- c("intercept","X1","X2","X3","X4")

# Create quadratic terms (squared predictors)
X.quadratic <- cbind(X, X.l1^2)

# Assign meaningful column names
colnames(X.quadratic) <- c("intercept","X1","X2","X3","X4",
                           "X1_squared","X2_squared","X3_squared","X4_squared")

# Reorder columns: each predictor followed by its squared term
X.quadratic <- X.quadratic[,c("intercept",
                              "X1","X1_squared",
                              "X2","X2_squared",
                              "X3","X3_squared",
                              "X4","X4_squared")]

library(ggplot2)

# True parameter vector (data-generating process)
beta  <- c(1,0,-0.5,0,0.5,-4,1,4,-1)

# ---------------------------------------------------------------------------
# 1. NOISE
# ---------------------------------------------------------------------------

# Select predictor (quadratic term of X3, depending on column order)
X_noise <- X.quadratic[, 6]

# True relationship (quadratic function using corresponding beta values)
y_true_noise <- X.quadratic[, c(6, 7)] %*% beta[6:7]

# Create dataframe for plotting
df_noise <- data.frame(
  x = X_noise,
  y = y_true_noise
)

# Fit models:
# (1) correctly specified quadratic model
# (2) misspecified linear model
quad_fit_noise <- lm(y ~ x + I(x^2), data = df_noise)
lin_fit_noise <- lm(y ~ x, data = df_noise)

# Create grid for smooth prediction lines
x_vals <- seq(min(X_noise), max(X_noise), length.out = 200)
pred_df_noise <- data.frame(x = x_vals)

# Generate predictions from both models
pred_df_noise$y_quad <- predict(quad_fit_noise, newdata = pred_df_noise)
pred_df_noise$y_lin <- predict(lin_fit_noise, newdata = pred_df_noise)

# Plot for NOISE case
Noise_plot <- ggplot(df_noise, aes(x = x, y = y)) +
  geom_line(data = pred_df_noise, aes(x = x, y = y_quad), color = "mediumblue", size = 1.2) +
  geom_line(data = pred_df_noise, aes(x = x, y = y_lin), color = "deeppink", size = 1.2) +
  geom_point(size = 2, color = "cornflowerblue") +
  labs(title = "Working Noise", x = "x", y = "y") +
  
  # Annotate true functional form
  annotate("text", x = 1.5, y = -2, 
           label = "x^2 - 4*x + 1", 
           parse = TRUE, 
           color = "blue", 
           size = 5) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 13)
  )

# ---------------------------------------------------------------------------
# 2. SIGNAL
# ---------------------------------------------------------------------------

# Select predictor (X1 and its squared term)
X_signal <- X.quadratic[, 2]

# True signal-generating relationship
y_true_signal <- X.quadratic[, c(2, 3)] %*% beta[2:3]

# Create dataframe
df_signal <- data.frame(
  x = X_signal,
  y = y_true_signal
)

# Fit models (quadratic vs linear)
quad_fit_signal <- lm(y ~ x + I(x^2), data = df_signal)
lin_fit_signal <- lm(y ~ x, data = df_signal)

# Prediction grid
x_vals_signal <- seq(min(X_signal), max(X_signal), length.out = 200)
pred_df_signal <- data.frame(x = x_vals_signal)

# Predictions
pred_df_signal$y_quad <- predict(quad_fit_signal, newdata = pred_df_signal)
pred_df_signal$y_lin <- predict(lin_fit_signal, newdata = pred_df_signal)

# Plot for SIGNAL case
Signal_plot <- ggplot(df_signal, aes(x = x, y = y)) +
  geom_line(data = pred_df_signal, aes(x = x, y = y_quad), color = "mediumblue", size = 1.2) +
  geom_line(data = pred_df_signal, aes(x = x, y = y_lin), color = "deeppink", size = 1.2) +
  geom_point(size = 2, color = "cornflowerblue") +
  labs(title = "Working Signal", x = "x", y = "y") +
  
  # Annotate true model
  annotate("text", x = 0.5, y = -2, 
           label = "-0.5*x^2 + 1", 
           parse = TRUE, 
           color = "blue", 
           size = 5) +
  
  theme_minimal(base_size = 14) +  
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 13)
  )

library(patchwork)

# Combine both plots side-by-side
DataGenerating_WM_plots <- Noise_plot + Signal_plot

# Save combined plot as high-resolution image
ggsave("DataGenerating_WM_plots.png", 
       plot = DataGenerating_WM_plots, 
       dpi = 900, width = 8, height = 4, units = "in")