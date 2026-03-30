setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /SimSanityCheck_Problems_with_bayesian_conditional")
source("SimSanity_helper.R")
source("SimSanity_dgp.R")
source("SimSanity_SelectionProcess.R")
source("SimSanity_one_simulation.R")
source("SimSanity_simulation_study.R")

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

SimulationStudy_SanityCheck <- readRDS("SimStudy2_Discussion_Bayesian_conditional.RDS")
SimSanity_Results <- SimulationStudy_SanityCheck$post_selection_pvalues

## QQ-plots --------------------------------------------------------------------

SimSanity_Results$coefficient <- ifelse(SimSanity_Results$coefficient == "X1", "L1-noise",'"L1-noise (random slope)"~x[2]')

library(ggplot2)
library(dplyr)

# Hilfsspalte für lesbare Beschriftung
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

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Figures/Figures project 1")
ggsave("plot_conditional_no_misspecification.png", 
       plot = plot_conditional_no_misspecification, 
       dpi = 900, width = 5, height = 6, units = "in")

################################################################################
# Omitted random slope 
################################################################################

source("SimSanity_OmittedRandomSlope_dgp.R")
source("SimSanity_OmittedRandomSlope_Selection_Process.R")

## Normal

start <- Sys.time()
SimulationStudy_SanityCheck_normalstderrors <-
  #suppressMessages(
    SimSanity_simulation_study(  n_iter = 500, #number of replications per condition
                                 n_obs = 40^2, # number of observations
                                 m = 40, # number of level-2 units
                                 R_squared = c(0.6), # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,
                                              0.2,1.5), nrow=2, byrow=T),
                                 perspective = c("marginal"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 B_MCsamples = 500, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X2","X3"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
   # )
  )
end <- Sys.time()
saveRDS(SimulationStudy_SanityCheck_normalstderrors, 
        file=paste0("OMITTEDSLOPE_SanityCheck_n=40_R_squared=0.6_marg_notRobust_MCconvergencecheck_1000iterations",".RDS"))


## Robust

start <- Sys.time()
SimulationStudy_SanityCheck_CR2 <-
  suppressMessages(
    SimSanity_simulation_study(  n_iter = 500, #number of replications per condition
                                 n_obs = 40^2, # number of observations
                                 m = 40, # number of level-2 units
                                 R_squared = c(0.6), # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,
                                              0.2,1.5), nrow=2, byrow=T),
                                 perspective = c("marginal"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 B_MCsamples = 500, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X2","X3"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
    )
  )
end <- Sys.time()
saveRDS(SimulationStudy_SanityCheck_CR2, 
        file=paste0("OMITTEDSLOPE_SanityCheck_n=40_R_squared=0.6_marg_CR2_MCconvergencecheck_500iterations",".RDS"))





## QQ-plots --------------------------------------------------------------------

SimSanity_Results <- SimulationStudy_SanityCheck_normalstderrors$post_selection_pvalues

SimSanity_Results <- SimulationStudy_SanityCheck_CR2$post_selection_pvalues


library(ggplot2)
library(dplyr)

# Hilfsspalte für lesbare Beschriftung
results_all_conditions_df <- SimSanity_Results  %>%
  mutate(
    label_facet = paste0(
      "R² = ", R_squared_pseudo,
      ", N = ", n_obs,
      ", M = ", m,
      ", coeff = ", coefficient,
      ", persp = ", perspective,
      ", workingModel = ", working_model
    )
  )

# QQ-Plot gegen Uniformverteilung
ggplot(results_all_conditions_df, aes(sample = pval_two, color = estimator)) +
  stat_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey") + 
  facet_wrap(~ label_facet, labeller = label_wrap_gen(width = 35)) +
  labs(
    title = "Simulation Study Std.errors",
    subtitle = "CR2 std.errors",# robust standard errors ?
    x = "Expected p-values",
    y = "Observed p-values",
    color = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9) 
  )



################################################################################
# Noise and signal with and without random slope 
################################################################################

source("SimSanity_signal_noise_with(out)_rs_dgp.R")
source("SimSanity_signal_noise_with(out)_rs_Selection_Process.R")

## Normal

start <- Sys.time()
SimulationStudy_SanityCheck_signal_noise_with_out_rs <-
  #suppressMessages(
  SimSanity_simulation_study(  n_iter = 150, #number of replications per condition
                               n_obs = 40^2, # number of observations
                               m = 40, # number of level-2 units
                               D = matrix(c(2,0.2,0.2,
                                            0.2,1.5,0.2,
                                            0.2,0.2,1.5), nrow=3, byrow=T),
                               R_squared = c(0.55), # pseudo coefficient of determination
                               perspective = c("marginal","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                               robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                               B_MCsamples = 500, # number of MC samples in selective inference process
                               remove_singular_fits = F,
                               interesting_vars = c("X1","X2","X3","X4"), # for which variables do we want to do post-selection inference
                               interesting_working_model_vars = c("X1","X2","X3","X4"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                               checkConvergence_for_check_Fun = T # convergence check for MC sampling?
                               # )
  )
end <- Sys.time()
saveRDS(SimulationStudy_SanityCheck_signal_noise_with_out_rs, 
        file=paste0("X1X2X3X4signal_noise_with(out)_rs_SanityCheck_n=40_R_squared=0.6_marg_notRobust_MCconvergencecheck_150iterations",".RDS"))


## QQ-plots ----------------------

SimSanity_Results <- SimulationStudy_SanityCheck_signal_noise_with_out_rs$post_selection_pvalues


library(ggplot2)
library(dplyr)

# Hilfsspalte für lesbare Beschriftung
results_all_conditions_df <- SimSanity_Results  %>%
  mutate(
    label_facet = paste0(
      "R² = ", R_squared_pseudo,
      ", N = ", n_obs,
      ", M = ", m,
      ", coeff = ", coefficient,
      ", persp = ", perspective,
      ", workingModel = ", working_model
    )
  )

# QQ-Plot gegen Uniformverteilung
ggplot(results_all_conditions_df, aes(sample = pval_two, color = estimator)) +
  stat_qq(distribution = stats::qunif) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey") + 
  facet_wrap(~ label_facet, labeller = label_wrap_gen(width = 35)) +
  labs(
    title = "Simulation Study Correctly specified model",
    subtitle = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9) 
  )


## Coverage rate ----------------------

SimSanity_Results$truebeta <- ifelse(SimSanity_Results$coefficient %in% c("X1","X2"),4,0)

SimSanity_Results$covered <- 
  ifelse(SimSanity_Results$cil < SimSanity_Results$truebeta & SimSanity_Results$ciu > SimSanity_Results$truebeta,1,0)



library(dplyr)

coverage <- SimSanity_Results %>%
  group_by(perspective, estimator, coefficient) %>%
  summarise(prop_covered = mean(covered, na.rm = TRUE), .groups = "drop")

coverage[coverage$perspective!="naive",]



################################################################################
# non-diagonal residual covariance matrix R, i.e., heteroscedasticity induced
# by specification of R
################################################################################

source("SimSanity_Rnondiagonal_dgp.R")
source("SimSanity_Rnondiagonal_Selection_Process.R")

## Normal

start <- Sys.time()
SimulationStudy_SanityCheck_nondiagonalR <-
  #suppressMessages(
  SimSanity_simulation_study(  n_iter = 5, #number of replications per condition
                               n_obs = 40^2, # number of observations
                               m = 40, # number of level-2 units
                               R_squared = c(0.6), # pseudo coefficient of determination
                               D = NULL,
                               perspective = c("marginal"), # marginal and/or conditional (for robust std. errors we only need marginal)
                               robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                               B_MCsamples = 500, # number of MC samples in selective inference process
                               remove_singular_fits = F,
                               interesting_vars = c("X2","X3"), # for which variables do we want to do post-selection inference
                               interesting_working_model_vars = c("X1","X2","X3"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                               checkConvergence_for_check_Fun = T # convergence check for MC sampling?
                               # )
  )
end <- Sys.time()
saveRDS(SimulationStudy_SanityCheck_nondiagonalR , 
        file=paste0("nondiagonalR_SanityCheck_n=40_R_squared=0.6_marg_notRobust_MCconvergencecheck_500iterations",".RDS"))

start <- Sys.time()
SimulationStudy_SanityCheck_nondiagonalR_Rsupplied <-
  #suppressMessages(
  SimSanity_simulation_study(  n_iter = 500, #number of replications per condition
                               n_obs = 40^2, # number of observations
                               m = 40, # number of level-2 units
                               R_squared = c(0.6), # pseudo coefficient of determination
                               D = NULL,
                               perspective = c("marginal"), # marginal and/or conditional (for robust std. errors we only need marginal)
                               robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                               B_MCsamples = 500, # number of MC samples in selective inference process
                               remove_singular_fits = F,
                               interesting_vars = c("X2","X3"), # for which variables do we want to do post-selection inference
                               interesting_working_model_vars = c("X1","X2","X3"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                               checkConvergence_for_check_Fun = T, # convergence check for MC sampling?
                               R_diag = F
                               # )
  )
end <- Sys.time()
saveRDS(SimulationStudy_SanityCheck_nondiagonalR_Rsupplied , 
        file=paste0("nondiagonalR_Rsupplied_SanityCheck_n=40_R_squared=0.6_marg_notRobust_MCconvergencecheck_1000iterations",".RDS"))


