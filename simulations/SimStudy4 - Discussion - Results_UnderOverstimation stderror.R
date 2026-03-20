#####################################################################################
# Results Discussion - Under-/Overestimation of Standard errors - Simulation Study 4
####################################################################################

## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

## Load Simulation Study 4 files
source("SimStudy4_helper.R")
source("SimStudy4_dgp.R")
source("SimStudy4_SelectionProcess.R")
source("SimStudy4_one_simulation.R")
source("SimStudy4_simulation_study.R")
source("SimStudy4 - Discussion - one_simulation_UnderOverstimation stderror.R")
source("SimStudy4 - Discussion - simulation_study_UnderOverstimation stderror.R")


R_squared = 0.3
n_obs = 40^2 ; m = 40

# load fixed predictor matriced
dgp_X_output <- readRDS("dgp_X_output_n_obs=1600_m=40.RDS")

# load list_supplied_D_work
load(paste0("~/SimulationResults_supplied_D/SIM4_list_supplied_D_work_fullrandom=(Intercept)X2X4_n=400_Rsquared=",R_squared,"_nrsamples=10000.RData"))



################################################################################
## Run the Simulation ##
################################################################################

setwd("~/SimStudy4_Results")

# Number of iterations
N <- 1000


save_file <- paste0(
  "Discussion_SimulationStudy4_underestimation_stderror.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy4_Discussion_Results_simulation_study(
        n_iter = 1, 
        n_obs = n_obs,
        m = m,
        R_squared = R_squared,
        l2.id = dgp_output$l2.id, # id variable for level-2 units
        X=dgp_X_output$X, # X-matrix 
        X.quadratic=dgp_X_output$X.quadratic, # X-matrix (including squared variables)
        Z=dgp_X_output$Z, # design matrix corresponding to X
        D =  matrix(c(2,0.2,0.1,
                      0.2,1.5,0.2,
                      0.1,0.2,1.5),
                    nrow = 3, byrow=T) , # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
        robust_std.errors = "CR2" , # (which) robust std.err in marginal perspective
        alpha = 0.000001,
        remove_singular_fits = F, #remove singular fit 
        interesting_working_model_ranef = c("(Intercept)","X2"),
        interesting_fixed_effects = c("X4")
      )
    },
    error = function(e) {
      warning("Fehler in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  results_list[[i]] <- res
  
  saveRDS(results_list, save_file)
}


save_file <- paste0(
  "Discussion_SimulationStudy4_overestimation_stderror.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy4_Discussion_Results_simulation_study(
        n_iter = 1, 
        n_obs = n_obs,
        m = m,
        R_squared = R_squared,
        l2.id = dgp_output$l2.id, # id variable for level-2 units
        X=dgp_X_output$X, # X-matrix 
        X.quadratic=dgp_X_output$X.quadratic, # X-matrix (including squared variables)
        Z=dgp_X_output$Z, # design matrix corresponding to X
        D =  matrix(c(2,0.2,0.1,
                      0.2,1.5,0.2,
                      0.1,0.2,1.5),
                    nrow = 3, byrow=T) , # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
        robust_std.errors = "CR2" , # (which) robust std.err in marginal perspective
        alpha = 10^-31,
        remove_singular_fits = F, #remove singular fit 
        interesting_working_model_ranef = c("(Intercept)","X2","X4"),
        interesting_fixed_effects = c("X4")
      )
    },
    error = function(e) {
      warning("Fehler in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  results_list[[i]] <- res
  
  saveRDS(results_list, save_file)
}





################################################################################
## QQ-Plots, hist, standard error ##
################################################################################
setwd("~/SimStudy4_Results")

SimulationStudy4_underestimation <- readRDS("Discussion_SimulationStudy4_underestimation_stderror.RDS")
SimulationStudy4_overestimation <- readRDS("Discussion_SimulationStudy4_overestimation_stderror.RDS")


sim_list <- list(
  underestimation = SimulationStudy4_underestimation,
  overestimation  = SimulationStudy4_overestimation
)


## QQ-plots --------------------------------------------------------------------

for(sim_name in names(sim_list)){
  
  SimStudy4_Discussion_Results <- do.call(
    rbind,
    lapply(sim_list[[sim_name]], function(x) x[[1]])
  )


SimStudy4_Discussion_Results$pval_type <- ifelse(SimStudy4_Discussion_Results$method%in% c("Naive"),"Naive","Standard (no selection)")

SimStudy4_Discussion_Results$estimator <- "ModelEstimate_CR2"

SimStudy4_Discussion_Results$coefficient <- "L1-noise~(random~slope)~x[4]"
SimStudy4_Discussion_Results$working_model <-  "(Intercept)X2"

library(ggplot2)
library(dplyr)


results_all_conditions_df <- SimStudy4_Discussion_Results  %>%
  mutate(
    label_facet = paste0(
      #"R² = ", R_squared_pseudo,
      "Selected/Working Random Effects = ", working_model
      
    )
  )


results_all_conditions_df$estimator <- factor(results_all_conditions_df$estimator)


results_all_conditions_df$pval_type <-relevel(as.factor(results_all_conditions_df$pval_type),"Standard (no selection)")

#-------------------------
# Manual color palette
# -----------------------------
estimator_colors <- c(
  "fullCandidate_supplied" ="black",
  "supplied"            = "#009E73",  # Grün
  "ModelEstimate"       =  "#8C510A",  # braun
  "ModelEstimate_CR2"       =  "#8C510A",  # braun
  "ICM"                 = "#D55E00",  # Orange
  "VarY"                =  "#0072B2",  # Blau
  "ModelEstimate_SIGMA" = "#56B4E9",  # Hellblau
  "ModelEstimate_R"     = "#E69F00",  # Gelb-Orange
  "ICM_SIGMA"           = "#CC79A7",  # Violett
  "ICM_R"               = "#6A3D9A",  # Dunkellila
  # deutlich dunkleres Türkis (petrol)
  "fullCandidate"           = "white"
  
)

# -----------------------------
# Plot
# -----------------------------


plot_SimStudy4_Discussion <- ggplot(results_all_conditions_df, 
                                                   aes(
                                                     sample = pvals,
                                                     color  = estimator,
                                                     linetype = pval_type)) +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  theme_minimal(base_size = 18) +  # Basis-Schriftgröße
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA),  
    strip.text = element_text(color = "black"),
    legend.position = "bottom",       # <--- Hier Legende unten setzen
    legend.box = "vertical"
  ) + 
   facet_grid(
    . ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )+
  geom_qq(  
    distribution = stats::qunif, size = 0.8, geom = "line") +
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  scale_color_manual(
    values = estimator_colors,
    labels = function(x) parse(text = c(
      fullCandidate_supplied = "fullCandidate (supplied)",
      supplied                 = "supplied",
      ModelEstimate               = "ModelEstimate",
      ModelEstimate_CR2            = "ModelEstimate(CR2)",
      ICM                         = "ICM",
      VarY                        = "VarY",
      ModelEstimate_SIGMA         = "ModelEstimate[Sigma]",
      ModelEstimate_R             = "ModelEstimate[R]",
      ICM_SIGMA                   = "ICM[Sigma]",
      ICM_R                       = "ICM[R]",
      fullCandidate               = "fullCandidate"
    )[x])
  )
plot_SimStudy4_Discussion

ggsave(
  paste0("plot_SimStudy4_Discussion", sim_name, ".png"),
  plot = plot,
  dpi = 900,
  width = 6,
  height = 6,
  units = "in"
)

}


## hist and standard errors-------------------------------------
sim_list <- list(
  underestimation = SimulationStudy4_underestimation,
  overestimation  = SimulationStudy4_overestimation
)

for(sim_name in names(sim_list)){
  
  SimStudy4_Discussion_Results <- do.call(
    rbind,
    lapply(sim_list[[sim_name]], function(x) x[[1]])
  )
  
  results_all_conditions_df <- SimStudy4_Discussion_Results %>%
    rename(
      pval_two = pvals,
      estimator = method
    ) %>%
    mutate(label_facet = coefficient)
  
  
  # -----------------------------
  # Histogram Estimates
  # -----------------------------
  
  plot_est <- ggplot(results_all_conditions_df,
                     aes(x = est)) +
    geom_histogram(
      bins = 30,
      color = "black",
      fill = "lightgrey"
    ) +
    facet_grid(coefficient ~ estimator) +
    labs(
      title = "Distribution of Estimates",
      x = "Estimate",
      y = "Frequency"
    ) +
    theme_bw()
  
  ggsave(
    paste0("hist_est_", sim_name, ".png"),
    plot_est,
    width = 8,
    height = 6,
    dpi = 600
  )
  
  
  # -----------------------------
  # Histogram Standard Errors
  # -----------------------------
  
  plot_se <- ggplot(results_all_conditions_df,
                    aes(x = se)) +
    geom_histogram(
      bins = 30,
      color = "black",
      fill = "lightgrey"
    ) +
    facet_grid(coefficient ~ estimator) +
    labs(
      title = "Distribution of Std. Errors",
      x = "Std. Error",
      y = "Frequency"
    ) +
    theme_bw()
  
  ggsave(
    paste0("hist_se_", sim_name, ".png"),
    plot_se,
    width = 8,
    height = 6,
    dpi = 600
  )
  
  
  # -----------------------------
  # Histogram p-values
  # -----------------------------
  
  plot_p <- ggplot(results_all_conditions_df,
                   aes(x = pval_two)) +
    geom_histogram(
      bins = 30,
      color = "black",
      fill = "lightgrey"
    ) +
    facet_grid(coefficient ~ estimator) +
    labs(
      title = "Distribution of p-values",
      x = "p-value",
      y = "Frequency"
    ) +
    theme_bw()
  
  ggsave(
    paste0("hist_pvals_", sim_name, ".png"),
    plot_p,
    width = 8,
    height = 6,
    dpi = 600
  )
  
  
  # -----------------------------
  # Summaries
  # -----------------------------
  
  summary_est <- results_all_conditions_df %>%
    group_by(coefficient, estimator) %>%
    summarise(
      mean_est = mean(est, na.rm = TRUE),
      sd_est   = sd(est, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_se <- results_all_conditions_df %>%
    group_by(coefficient, estimator) %>%
    summarise(
      mean_se = mean(se, na.rm = TRUE),
      sd_se   = sd(se, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(paste("Summary estimates:", sim_name))
  print(summary_est)
  
  print(paste("Summary SE:", sim_name))
  print(summary_se)
  
}


