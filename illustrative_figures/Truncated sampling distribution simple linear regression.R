################################################################################
# Figure: Truncated sampling distribution / distribution of test statistic
################################################################################

library(ggplot2)

# Define the range for the x-axis (grid for plotting the distribution)
x_vals <- seq(-4, 4, length.out = 1000)

# Create data frame for standard normal density
df <- data.frame(
  x = x_vals,
  y = dnorm(x_vals)
)

# Degrees of freedom setup (based on sample size N and model dimension)
N = 500

# Truncation thresholds (selection region)
# Observations between these values will later be removed (set to zero density)
trunc_lower <- qt(0.1, df = N-2)
trunc_upper <- qt(0.9, df = N-2)

# Standard (correct) critical values without selection
lower <- qt(0.05, df = N-2)
upper <- qt(0.95, df = N-2)

# ---------------------------------------------------------------------------
# Plot 1: Distribution WITHOUT selection
# ---------------------------------------------------------------------------

# Data frame for vertical lines (used to display correct critical values)
vlines <- data.frame(
  xintercept = c(lower, upper),
  label = "Correct Critical Values"
)

# Plot of the full (untruncated) distribution
normal_dis <- ggplot(df, aes(x, y)) +
  
  # Filled density area
  geom_area(fill = "lightgray", alpha = 0.5) +
  
  # Density curve
  geom_line(size = 0.5) +
  
  # Add vertical lines for critical values
  geom_vline(data = vlines, aes(xintercept = xintercept, color = label),
             show.legend = TRUE) +
  
  # Customize legend and colors
  scale_color_manual(name = "", values = c("Correct Critical Values" = "forestgreen")) +
  
  # Label x-axis (null hypothesis point)
  scale_x_continuous(breaks = c(0), labels = c(expression(rho == 0))) + 
  
  # Fix y-axis limits
  scale_y_continuous(limits = c(0, max(df$y)))+
  
  labs(
    title = "Distribution of Test Statistic without Selection",
    x = "Test Statistic Values",
    y = "Density"
  ) +
  
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
# Apply truncation (selection effect)
# ---------------------------------------------------------------------------

# Create new density where values inside the selection interval are removed
# (density set to zero in [trunc_lower, trunc_upper])
df$y_trunc <- ifelse(df$x >= trunc_lower & df$x <= trunc_upper, 0, df$y)

# Adjusted ("correct") critical values after truncation:
# These account for the reduced probability mass due to selection
lower_truncated <- qt(0.05 * 2 * pnorm(trunc_lower), df = N-2)
upper_truncated <- -qt(0.05 * 2 * pnorm(trunc_lower), df = N-2)

# Combine original (incorrect) and adjusted (correct) critical values
vlines_combined <- data.frame(
  xintercept = c(lower, upper, lower_truncated, upper_truncated),
  label = c("Incorrect Critical Values", "Incorrect Critical Values",
            "Correct Critical Values", "Correct Critical Values"),
  linetype = c("solid", "solid", "dashed", "dashed")
)

# ---------------------------------------------------------------------------
# Plot 2: Distribution AFTER selection (truncated distribution)
# ---------------------------------------------------------------------------

truncated_dis <- ggplot(df, aes(x, y = y_trunc)) +
  
  # Filled truncated density
  geom_area(fill = "lightgray", alpha = 0.5) +
  
  # Density curve
  geom_line(size = 0.5) +
  
  # Add both incorrect (solid) and correct (dashed) critical values
  geom_vline(
    data = vlines_combined,
    aes(xintercept = xintercept, color = label, linetype = linetype),
    show.legend = TRUE
  ) +
  
  # Color mapping for correct vs incorrect critical values
  scale_color_manual(
    name = "",
    values = c("Correct Critical Values" = "forestgreen",
               "Incorrect Critical Values" = "purple")
  ) +
  
  # Use linetype as specified in the data frame
  scale_linetype_identity() +  
  
  # Axis formatting
  scale_x_continuous(breaks = c(0), labels = c(expression(rho == 0))) + 
  scale_y_continuous(limits = c(0, max(df$y)))+
  
  labs(
    title = "Distribution of Test Statistic after Selection",
    x = "Test Statistic Values",
    y = "Density"
  ) +
  
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

# Combine both plots vertically (before vs after selection)
truncation_plots <- normal_dis / truncated_dis

# Save the figure as high-resolution image
ggsave("truncation_plots.png", 
       plot = truncation_plots, 
       dpi = 900, width = 8, height = 5, units = "in")