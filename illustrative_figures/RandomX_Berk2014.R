library(ggplot2)

set.seed(123)

f <- function(x) {
  0.5 * x + 0.8 * x^2
}

x_pop <- seq(-2, 2, length.out = 200)
y_pop <- f(x_pop)
pop_data <- data.frame(x = x_pop, y = y_pop)

x_blue <- runif(12, -2, 0.2)
y_blue <- f(x_blue) 
blue_data <- data.frame(x = x_blue, y = y_blue, group = "Blue Sample")

x_red <- runif(12, 0, 2)
y_red <- f(x_red) 
red_data <- data.frame(x = x_red, y = y_red, group = "Brown Sample")

lm_pop  <- lm(y ~ x, data = pop_data)

pop_data$y_lin <- predict(lm_pop, newdata = pop_data)
pop_data$group_curve <- "True Nonlinear Mean"
pop_data$group_lin   <- "Population Linear Approximation"

RandomX <- ggplot() +
  # True non-linear
  geom_line(data = pop_data,
            aes(x, y, color = "True Nonlinear Mean"),
            linetype = "dashed", size = 1) +
  
  # Linear Population-Approximation 
  geom_line(data = pop_data,
            aes(x, y_lin, color = "Population Linear Approximation"),
            size = 1.2) +
  
  # Blue sample + Regression
  geom_point(data = blue_data,
             aes(x, y, color = "Blue Sample"), size = 2) +
  geom_smooth(data = blue_data,
              aes(x, y, color = "Blue Sample"),
              method = "lm", se = FALSE, size = 1) +
  
  # Brown sample + Regression
  geom_point(data = red_data,
             aes(x, y, color = "Brown Sample"), size = 2) +
  geom_smooth(data = red_data,
              aes(x, y, color = "Brown Sample"),
              method = "lm", se = FALSE, size = 1) +
  
  scale_color_manual(
    name = NULL,  
    values = c(
      "Blue Sample" = "blue",
      "Brown Sample" = "brown",
      "True Nonlinear Mean" = "black",
      "Population Linear Approximation" = "gray80"
    )
  ) + labs(
    title = "",
    x = "x",
    y = "y"
  ) +theme_minimal(base_size =20)+
  theme(legend.position = "right") +
  guides(color = guide_legend(ncol = 1))

ggsave(
  "plot_RandomX_Appendix.png",
  plot = RandomX,
  dpi = 300,
  width = 9,
  height = 3.5,
  units = "in"
)

