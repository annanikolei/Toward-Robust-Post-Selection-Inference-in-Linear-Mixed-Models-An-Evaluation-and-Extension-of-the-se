## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

## Load Simulation Study 1 files
source("SimStudy3_helper.R")
source("RandomX_SimStudy3_dgp.R")
source("SimStudy3_SelectionProcess.R")
source("RandomX_SimStudy3_one_simulation.R")
source("RandomX_SimStudy3_simulation_study.R")


R_squared = 0.3
R_squared = "great" # this is the condition with the small level-1 resiudal variance
n_obs = 40^2 ; m = 40


########################
# Run Simulations
########################

########################
# no_HC_correction
########################

## all perspectives no_HC_correction
#### X1234
setwd("RandomX_SimStudy3_Results")
N <- 1000


save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}



########################
# HC3 correction
########################

#### X1234
N <- 1000

save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag =F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}

SimulationStudy3_no_HC_X1234 <- readRDS(paste0("RandomX_SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_",R_squared,"_1000MC.RDS"))
SimulationStudy3_HC3_X1234 <- readRDS(paste0("RandomX_SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_",R_squared,"_1000MC.RDS"))
SimulationStudy3 <- c(SimulationStudy3_no_HC_X1234, 
                      SimulationStudy3_HC3_X1234)
saveRDS(SimulationStudy3,paste0("RandomX_SimStudy3_workingModel=X1234_R2=",R_squared,"_1000MC.RDS"))

################################################################################
## QQ-Plots, Convergence ##
################################################################################

setwd("RandomX_SimStudy3_Results")
SimulationStudy3 <- readRDS(paste0("RandomX_SimStudy3_workingModel=X1234_R2=",R_squared,"_1000MC.RDS"))
SimStudy3_Results <- do.call(rbind, lapply(SimulationStudy3 ,
                                           function(x) x$post_selection_pvalues)
)


SimStudy3_Results$pval_type <- ifelse(SimStudy3_Results$estimator%in% c("naive"),"Naive", ifelse(SimStudy3_Results$estimator%in% c("fullCandidate","fullCandidate_supplied"),"Standard (no selection)","Selfmade"))
SimStudy3_Results$estimator <- ifelse(SimStudy3_Results$estimator=="naive","ModelEstimate",SimStudy3_Results$estimator)
SimStudy3_Results$perspective <- ifelse(SimStudy3_Results$estimator %in% c("ModelEstimate_SIGMA_replaced","ModelEstimate_R_replaced",
                                                                           "ICM_SIGMA_replaced","ICM_R_replaced"),"marginal_not_efficient_HC3", SimStudy3_Results$perspective  )                                                                        



new_persp <- c(
  "marginal",
  "marginal_not_efficient",
  "marginal_not_efficient_HC3",
  "conditional",
  "fullCandidate"
)
naive_rows <- SimStudy3_Results[SimStudy3_Results$perspective == "naive", ]
naive_expanded <- naive_rows[rep(seq_len(nrow(naive_rows)), each = length(new_persp)), ]
naive_expanded$perspective <- rep(new_persp, times = nrow(naive_rows))
SimStudy3_Results <- SimStudy3_Results[SimStudy3_Results$perspective != "naive", ]
SimStudy3_Results <- rbind(SimStudy3_Results, naive_expanded)

# supplied is the same if HC3 estimator is used => duplicate this
supplied_rows <- SimStudy3_Results[SimStudy3_Results$estimator == "supplied" & SimStudy3_Results$perspective=="marginal_not_efficient", ]
supplied_expanded <- supplied_rows[rep(seq_len(nrow(supplied_rows)), each = length("marginal_not_efficient_HC3")), ]
supplied_expanded$perspective <- rep("marginal_not_efficient_HC3", times = nrow(supplied_rows))

SimStudy3_Results <- rbind(SimStudy3_Results, supplied_expanded )


SimStudy3_Results$coefficient <-
  ifelse(SimStudy3_Results$coefficient == "X3",
         "L1-noise~x[3]",
         "L1-noise~(random~slope)~x[4]")


library(ggplot2)
library(dplyr)


results_all_conditions_df <- SimStudy3_Results  %>%
  mutate(
    label_facet = paste0(
      #"R² = ", R_squared_pseudo,
      "Selected Fixed = ", working_model
      
    )
  )

results_all_conditions_df$perspective <- factor(results_all_conditions_df$perspective)
results_all_conditions_df$perspective <- factor(
  results_all_conditions_df$perspective,
  levels = c("conditional", "marginal", "marginal_not_efficient","marginal_not_efficient_HC3","fullCandidate"),
  labels = c("Conditional", "Marginal (eff)", "Marginal (not eff)","Marginal (not eff) HC3","fullCandidate")
)

results_all_conditions_df$estimator <- factor(results_all_conditions_df$estimator)
results_all_conditions_df$estimator <- factor(
  results_all_conditions_df$estimator,
  levels = c("fullCandidate_supplied","supplied","ModelEstimate", "ICM","VarY",
             "ModelEstimate_SIGMA_replaced","ModelEstimate_R_replaced",
             "ICM_SIGMA_replaced","ICM_R_replaced",
             "fullCandidate"),
  labels = c("fullCandidate_supplied","supplied","ModelEstimate", "ICM","VarY",
             "ModelEstimate_SIGMA","ModelEstimate_R",
             "ICM_SIGMA","ICM_R",
             "fullCandidate")
)

results_all_conditions_df$pval_type<-relevel(as.factor(results_all_conditions_df$pval_type),"Selfmade")
results_all_conditions_df <- results_all_conditions_df[results_all_conditions_df$pval_type=="Selfmade",] #remove naive
results_all_conditions_df$coefficient <- factor(
  results_all_conditions_df$coefficient,
  levels = c("L1-noise~x[3]", "L1-noise~(random~slope)~x[4]")
)

results_all_conditions_df <- 
  results_all_conditions_df[results_all_conditions_df$estimator %in% c("supplied","ModelEstimate", "ICM","VarY",
                                                                       "ModelEstimate_SIGMA","ModelEstimate_R"),]


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


results_all_conditions_df_SimStudy3 <- 
  results_all_conditions_df[results_all_conditions_df$perspective %in% c("Marginal (eff)","Marginal (not eff)","Marginal (not eff) HC3"),]



plot_SimStudy3 <- ggplot(results_all_conditions_df_SimStudy3, 
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

plot_SimStudy3
ggsave(paste0("plot_RandomX_SimStudy3_R2=",R_squared,".png"), 
       plot = plot_SimStudy3, 
       dpi = 900, width = 9, height = 10, units = "in")

