## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

source("SimStudy2 - Discussion - BayesianConditional_helper.R")
source("SimStudy2 - Discussion - BayesianConditional_dgp.R")
source("SimStudy2 - Discussion - BayesianConditional_Selection_Process.R")
source("SimStudy2 - Discussion - BayesianConditional_one_simulation.R")
source("SimStudy2 - Discussion - BayesianConditional_simulation_study.R")

setwd("SimStudy2_Results")
SimulationStudy_SanityCheck <-
  SimSanity_simulation_study(  n_iter = 1000, #number of replications per condition
                               n_obs = 40^2, # number of observations
                               m = 40, # number of level-2 units
                               R_squared = c(0.6), # pseudo coefficient of determination
                               D = matrix(c(2,0.2,
                                            0.2,1.5), nrow=2, byrow=T),
                               perspective = c("marginal","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                               robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                               B_MCsamples = 500, # number of MC samples in selective inference process
                               remove_singular_fits = F,
                               interesting_vars = c("X2"), # for which variables do we want to do post-selection inference
                               interesting_working_model_vars = c("X1","X2"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                               checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  )

saveRDS(SimulationStudy_SanityCheck, 
file=paste0("SimStudy2_Discussion_Bayesian_conditional",".RDS"))


## QQ-plots --------------------------------------------------------------------
setwd("SimStudy2_Results")
SimulationStudy_SanityCheck <- readRDS("SimStudy2_Discussion_Bayesian_conditional.RDS")
SimSanity_Results <- SimulationStudy_SanityCheck$post_selection_pvalues

SimSanity_Results$coefficient <- ifelse(SimSanity_Results$coefficient == "X1", "L1-noise",'"L1-noise (random slope)"~x[2]')

library(ggplot2)
library(dplyr)

results_all_conditions_df <- SimSanity_Results  %>%
  mutate(
    label_facet = paste0(
      "Selected Fixed Effects = ", working_model
    )
  )

results_all_conditions_df <- results_all_conditions_df[results_all_conditions_df$perspective %in% c("conditional","marginal"),]

estimator_colors <- c(
  "supplied"            = "#009E73",  
  "ModelEstimate"       =  "#8C510A"  
  
)

results_all_conditions_df$perspective <- 
  factor(results_all_conditions_df$perspective)

results_all_conditions_df$perspective <- 
  relevel(results_all_conditions_df$perspective, "marginal")

levels(results_all_conditions_df$perspective) <- 
  c("Marginal (eff)", "Conditional")


plot_conditional_no_misspecification <- ggplot(results_all_conditions_df, 
                         aes(
                           sample = pval_two,
                           color  = estimator)) +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  theme_minimal(base_size = 18) +  # Basis-Schriftgröße
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA),  
    strip.text = element_text(color = "black"),
    legend.position = "bottom",       # <--- Hier Legende unten setzen
    legend.box = "vertical"
  ) +
  geom_qq(  
    distribution = stats::qunif, size = 0.8, geom = "line") +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )+
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator"
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  scale_color_manual(
    values = estimator_colors,
    labels = function(x) parse(text = c(
      supplied                 = "supplied",
      ModelEstimate               = "ModelEstimate"
    )[x])
  )
plot_conditional_no_misspecification

ggsave("plot_conditional_no_misspecification.png", 
       plot = plot_conditional_no_misspecification, 
       dpi = 900, width = 5, height = 6, units = "in")

