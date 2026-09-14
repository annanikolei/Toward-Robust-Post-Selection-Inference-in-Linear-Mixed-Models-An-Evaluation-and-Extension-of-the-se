################################################################################
# Discussion - Naive Post-selection Inference with supplied - Simualtion Study 4
################################################################################
## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

## Load Simulation Study 4 files
source("SimStudy4_helper.R")
source("SimStudy4_dgp.R")
source("SimStudy4_SelectionProcess.R")
source("SimStudy4 - Discussion - one_simulation_Naive supplied.R")
source("SimStudy4 - Discussion - simulation_study_Naive supplied.R")


R_squared = 0.3
n_obs = 40^2 ; m = 40

# load fixed predictor matriced
dgp_X_output <- readRDS("dgp_X_output_Sim4_n_obs=1600_m=40.RDS")

################################################################################
## Run the Simulation ##
################################################################################

setwd("SimStudy4_Results")


# Number of planned iterations
N <- 1000

save_file <- paste0(
  "Discussion_SimulationStudy4_marginal(eff)_Naive_supplied_500MC.RDS"
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
      SimStudy4_Naive_supplied_simulation_study(n_iter = 1,
                                                l2.id = dgp_X_output$l2.id, 
                                                n_obs = dgp_X_output$n_obs, # number of observations
                                                m = dgp_X_output$m, # number of level-2 units
                                                X = dgp_X_output$X , # X-matrix 
                                                X.quadratic = dgp_X_output$X.quadratic, # X-matrix (including squared variables)
                                                Z = dgp_X_output$Z, # design matrix corresponding to X
                                                perspective = c("marginal"),
                                                R_squared = 0.3, # pseudo coefficient of determination
                                                D = matrix(c(2,0.2,0.1,
                                                             0.2,1.5,0.2,
                                                             0.1,0.2,1.5),
                                                           nrow = 3, byrow=T) , # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                                remove_singular_fits = F, #remove singular fit 
                                                alpha = 0.000001, # alpha for significance selection
                                                interesting_working_model_ranef = c("(Intercept)","X2"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                                                interesting_fixed_effects = c("X3","X4"))
    },
    error = function(e) {
      warning("error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  results_list[[i]] <- res
  
  saveRDS(results_list, save_file)
}


save_file <- paste0(
  "Discussion_SimulationStudy4_marginal(not_eff)_Naive_supplied_500MC.RDS"
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
      SimStudy4_Naive_supplied_simulation_study(n_iter = 1,
                 l2.id = dgp_X_output$l2.id, 
                 n_obs = dgp_X_output$n_obs, # number of observations
                 m = dgp_X_output$m, # number of level-2 units
                 X = dgp_X_output$X , # X-matrix 
                 X.quadratic = dgp_X_output$X.quadratic, # X-matrix (including squared variables)
                 Z = dgp_X_output$Z, # design matrix corresponding to X
                 perspective = c("marginal_not_efficient"),
                 R_squared = 0.3, # pseudo coefficient of determination
                 D = matrix(c(2,0.2,0.1,
                              0.2,1.5,0.2,
                              0.1,0.2,1.5),
                            nrow = 3, byrow=T) , # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                 remove_singular_fits = F, #remove singular fit 
                 alpha = 0.000001, # alpha for significance selection
                 interesting_working_model_ranef = c("(Intercept)","X2"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                 interesting_fixed_effects = c("X3","X4"))
    },
    error = function(e) {
      warning("error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  results_list[[i]] <- res
  
  saveRDS(results_list, save_file)
}


################################################################################
## QQ-Plots
################################################################################
setwd("SimStudy4_Results")

Discussion_SimulationStudy4_Naive_supplied <- c(readRDS("~/Downloads/selective-inference-lmm-simulation-main/simulations/SimStudy4_Results/Discussion_SimulationStudy4_marginal(not_eff)_Naive_supplied_500MC.RDS"),
                                                readRDS("~/Downloads/selective-inference-lmm-simulation-main/simulations/SimStudy4_Results/Discussion_SimulationStudy4_marginal(eff)_Naive_supplied_500MC.RDS"))
SimulationStudy4_Naive_supplied <- do.call(rbind, lapply(Discussion_SimulationStudy4_Naive_supplied, `[[`, 1))

SimulationStudy4_Naive_supplied$pval_type <- "Naive"
SimulationStudy4_Naive_supplied$estimator <- "supplied"
SimulationStudy4_Naive_supplied$working_model <- "(Intercept)X2"

SimulationStudy4_Naive_supplied$coefficient <-
  ifelse(SimulationStudy4_Naive_supplied$coefficient == "X3",
         "L1-noise~x[3]",
         "L1-noise~(random~slope)~x[4]")


library(ggplot2)
library(dplyr)


results_all_conditions_df <- SimulationStudy4_Naive_supplied  %>%
  mutate(
    label_facet = paste0(
      #"R² = ", R_squared_pseudo,
      "Selected Random Effects = ", working_model
      
    )
  )

results_all_conditions_df$perspective <-
  factor(results_all_conditions_df$perspective)


results_all_conditions_df$perspective <- factor(
  results_all_conditions_df$perspective,
  
  levels = c(
    "marginal",
    "marginal_not_efficient"
  ),
  
  labels = c(
    "Marginal\n(eff)",
    "Marginal\n(not eff)"
  )
)

results_all_conditions_df$coefficient <- factor(
  results_all_conditions_df$coefficient,
  levels = c("L1-noise~x[3]", "L1-noise~(random~slope)~x[4]")
)

# -----------------------------
# Manual color palette
# -----------------------------

estimator_colors <- c(
  "fullCandidate_supplied" ="black",
  "supplied"            = "#009E73",  
  "ModelEstimate"       =  "#8C510A",  
  "ICM"                 = "#D55E00",  
  "VarY"                =  "#0072B2", 
  "ModelEstimate_SIGMA" = "#56B4E9",  
  "ModelEstimate_R"     = "#E69F00", 
  "ICM_SIGMA"           = "#CC79A7",  
  "ICM_R"               = "#6A3D9A",  
  "fullCandidate"           = "white"
  
)

# -----------------------------
# Plot
# -----------------------------


plot_SimStudy4_Naive_supplied<- ggplot(results_all_conditions_df, 
                                aes(
                                  sample = pval_two,
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
      ICM                         = "ICM",
      VarY                        = "VarY",
      ModelEstimate_SIGMA         = "ModelEstimate[Sigma]",
      ModelEstimate_R             = "ModelEstimate[R]",
      ICM_SIGMA                   = "ICM[Sigma]",
      ICM_R                       = "ICM[R]"
    )[x])
  )


plot_SimStudy4_Naive_supplied

ggsave("plot_SimStudy4_Naive_supplied.png", 
       plot = plot_SimStudy4_Naive_supplied, 
       dpi = 300, width = 9, height = 6.5, units = "in")





