library(ggplot2)
library(dplyr)

# Parameters for the normal distribution
mean <- 0
sd <- 1
x_vals <- seq(-4, 4, length.out = 1000)

# Compute density values of the normal distribution
density <- dnorm(x_vals, mean, sd)

# Compute quantiles (thresholds for regions)
q10 <- qnorm(0.1, mean, sd)
q90 <- qnorm(0.90, mean, sd)
q97 <- qnorm(0.98, mean, sd)

# Create base dataframe with x values and corresponding density
df <- data.frame(
  x = x_vals,
  y = density
)

# Define regions based on quantiles:
# low1: lower 10%, low2: between 90% and 98%, high: above 98%, mid: everything else
df$region2 <- with(df, ifelse(x < q10, "low1",
                              ifelse(x > q90 & x < q97, "low2",
                                     ifelse(x > q97, "high", "mid"))))

# Define colors for regions (used later in plots)
fill_colors2 <- c("low1" = "darkblue", "low2" = "darkblue", "high" = "cornflowerblue")

# Remove the "mid" region (only keep tail regions for histogram-like visualization)
df_nomid <- subset(df, region2 != "mid")

# Define bins for grouping (small width for fine resolution)
bin_width <- 0.1
bins <- seq(min(df_nomid$x), max(df_nomid$x), by = bin_width)

# Aggregate data into bins and compute mean x and density per bin
df_hist <- df_nomid %>%
  mutate(bin = cut(x, breaks = bins, include.lowest = TRUE)) %>%
  group_by(bin, region2) %>%
  summarise(x = mean(x), y = mean(y), .groups = "drop")

# Maximum density value (used for scaling points)
max_density <- max(df$y)

max_points <- 30  # Maximum number of stacked points per bin (controls resolution)

# Convert density into number of points per bin (scaled)
df_hist <- df_hist %>%
  mutate(points = ceiling(y / max_density * max_points))

set.seed(42)

# Determine jitter width depending on distance from center
# (more spread in the tails, less near the center)
x_min <- min(df_hist$x)
x_max <- max(df_hist$x)
max_jitter <- 0.00015

get_jitter_width <- function(x) {
  dist_to_center <- abs(x - (x_min + (x_max - x_min)/2))
  max_dist <- (x_max - x_min) / 2
  jitter_width <- max_jitter * (dist_to_center / max_dist)
  return(jitter_width)
}

# Generate stacked points with controlled jitter
df_points <- df_hist %>%
  rowwise() %>%
  do({
    n <- .$points           # number of points in this bin
    x_val <- .$x            # bin center
    region <- .$region2     # region label
    y_bin <- .$y            # density value
    
    jitter_width <- get_jitter_width(x_val)
    
    # Apply uniform random horizontal jitter within [-jitter_width, jitter_width]
    x_jittered <- x_val + runif(n, -jitter_width, jitter_width)
    
    # Create vertically stacked points from 0 up to density value
    y_points <- seq(0, y_bin, length.out = n)
    
    data.frame(x = x_jittered, region2 = region, y_point = y_points)
  }) %>%
  ungroup()

# Histogram-like point plot
histplot <- ggplot(df_points, aes(x = x, y = y_point, color = region2)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = fill_colors2, guide = "none") +
  scale_x_continuous(breaks = c(0, q97), labels = c(expression(rho), expression(t[obs]))) +
  labs(x = "Test Statistic Values", y = "Frequency") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_cartesian(ylim = c(-0.03, max_density * 1.05))+
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 13)
  )

# --- Recompute distribution setup (same as above, reused for second plot) ---

mean <- 0
sd <- 1
x_vals <- seq(-4, 4, length.out = 1000)

density <- dnorm(x_vals, mean, sd)

q10 <- qnorm(0.1, mean, sd)
q90 <- qnorm(0.9, mean, sd)
q97 <- qnorm(0.98, mean, sd)

df <- data.frame(
  x = x_vals,
  y = density
)

# Define regions again (same logic as before)
df$region2 <- with(df, ifelse(x < q10, "low1",
                              ifelse(x > q90 & x < q97, "low2",
                                     ifelse(x > q97, "high", "mid"))))

# Include "mid" region color for area plot
fill_colors2 <- c("low1" = "darkblue", "low2" = "darkblue", "mid" = "red", "high" = "cornflowerblue")

# Draw random samples:
# combination of histogram points and additional normally distributed points
set.seed(123)
samples <- c(df_points$x, rnorm(500, mean=mean(c(q10,q90)),sd=0.8))

# Assign each sample to a region (same logic as above)
samples_region <- ifelse(samples < q10, "low1",
                         ifelse(samples > q90 & samples < q97, "low2",
                                ifelse(samples > q97, "high", "mid")))

# Create dataframe for sample points (plotted slightly below curve)
df_samples <- data.frame(x = samples, y = -0.01, region2 = samples_region)

# Plot normal distribution with colored regions and sample points
normalplot <- ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "black") +
  geom_area(data = subset(df, region2 == "mid"), aes(fill = region2), alpha = 0.5) +
  geom_area(data = subset(df, region2 == "low1"), aes(fill = region2), alpha = 0.5) +
  geom_area(data = subset(df, region2 == "low2"), aes(fill = region2), alpha = 0.5) +
  geom_area(data = subset(df, region2 == "high"), aes(fill = region2), alpha = 0.5) +
  geom_point(data = df_samples, aes(x = x, y = y, color = region2), size = 2) +
  scale_fill_manual(values = fill_colors2, guide = "none") +
  scale_color_manual(values = fill_colors2, guide = "none") +
  scale_x_continuous(breaks = c(0, q97),
                     labels = c(expression(rho), expression(t[obs]))) +
  labs(
    title = "",
    x = "",
    y = "Density"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 13)
  )+
  coord_cartesian(ylim = c(-0.03, max(df$y)))

library(dplyr)
library(ggplot2)

# Colors for fraction plot
fill_colors2 <- c("low1" = "darkblue", "low2" = "darkblue", "high" = "cornflowerblue")

# Numerator: only "high" region points
numerator_df <- df_samples %>%
  filter(region2 == "high")

# Denominator: all tail regions (low1, low2, high)
denominator_df <- df_samples %>%
  filter(region2 %in% c("low1", "low2", "high")) %>%
  arrange(factor(region2, levels = c("low1", "low2", "high")))

# Helper function to arrange points in a horizontal line (like a count visualization)
create_inline_points <- function(data, y_val, point_spacing = 0.6) {
  n <- nrow(data)
  total_width <- point_spacing * (n - 1)
  x_start <- -total_width / 2  # center around x = 0
  data$x <- x_start + seq(0, by = point_spacing, length.out = n)
  data$y <- y_val
  data
}

# Create numerator and denominator lines
numerator_line <- create_inline_points(numerator_df, y_val = 1.1,  point_spacing = 0.4)
denominator_line <- create_inline_points(denominator_df, y_val = 0.9, point_spacing = 0.4)

# Combine both into one dataframe
bruch_df <- bind_rows(numerator_line, denominator_line)

# Positions for "#" labels (counts)
x_num_hash <- min(numerator_line$x) - 0.4
x_den_hash <- min(denominator_line$x) - 0.4

# Fraction-style plot (visualizing ratio)
zetaplot <- ggplot(bruch_df, aes(x = x, y = y, color = region2)) +
  geom_point(size = 2) +
  
  # Fraction line
  geom_segment(aes(x = min(x) - 0.2, xend = max(x) + 0.2, y = 1, yend = 1),
               inherit.aes = FALSE, size = 0.8) +
  
  # Left-side formula label
  annotate("text", x = min(bruch_df$x) - 1, y = 1, 
           label = expression(zeta(rho) == ""), parse = TRUE,
           size = 6, hjust = 1) +
  
  # "#" symbols indicating counts
  annotate("text", x = x_num_hash, y = 1.1, label = "#", size = 6, fontface = "bold") +
  annotate("text", x = x_den_hash, y = 0.9, label = "#", size = 6, fontface = "bold") +
  
  scale_color_manual(values = fill_colors2, guide = "none") +
  coord_fixed(ratio = 10, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 40))

library(patchwork)

# Combine all plots vertically
normalplot / histplot / zetaplot +
  plot_layout(
    widths = c(1),
    heights = c(2,2,1)
  )

# Combine only main plots
MCplots <- normalplot / histplot +
  plot_layout(
    widths = c(1),
    heights = c(2,2,1)
  )

MCplots

# Save plot to file (high resolution)
ggsave("MC_plots.png", 
       plot = MCplots, 
       dpi = 900, width = 9, height = 7, units = "in")