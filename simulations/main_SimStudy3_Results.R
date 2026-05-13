## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

## Load Simulation Study 1 files
source("SimStudy3_helper.R")
source("SimStudy3_dgp.R")
source("SimStudy3_SelectionProcess.R")
source("SimStudy3_one_simulation.R")
source("SimStudy3_simulation_study.R")


R_squared = 0.3
n_obs = 40^2 ; m = 40

# load fixed predictor matriced
dgp_X_output <- readRDS("dgp_X_output_n_obs=1600_m=40.RDS")

# load list_supplied_D_work
load(paste0("SimulationResults_supplied_D/SIM3_list_supplied_D_work_full=X1X2X3X4_interesting=X3X4_n=400_Rsquared=",R_squared,"_nrsamples=10000.RData"))


########################
# Run Simulations
########################

########################
# no_HC_correction
########################

## all perspectives no_HC_correction
#### X1234
setwd("SimStudy3_Results")
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_", 0.3, "_1000MC.RDS"
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
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
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


#### X123
N <- 1000

    
save_file <- paste0(
  "SimStudy3_workingModel=X123_no_HC_correction_allperspectives_R2_", 0.3, "_1000MC.RDS"
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
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3"),
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

     


#### X124
N <- 1000

    
save_file <- paste0(
  "SimStudy3_workingModel=X124_no_HC_correction_allperspectives_R2_", 0.3, "_1000MC.RDS"
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
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X4"),
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
  "SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_", 0.3, "_1000MC.RDS"
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
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 R_diag =F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
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

     


#### X123
N <- 1000

    
save_file <- paste0(
  "SimStudy3_workingModel=X123_HC3_marginal(not_eff)_R2_", 0.3, "_1000MC.RDS"
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
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 R_diag =F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3"),
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

     

#### X124
N <- 1000

    
save_file <- paste0(
  "SimStudy3_workingModel=X124_HC3_marginal(not_eff)_R2_", 0.3, "_1000MC.RDS"
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
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 R_diag =F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X4"),
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


################################################################################
## QQ-Plots, Convergence ##
################################################################################
setwd("SimStudy3_Results")

  SimulationStudy3_no_HC_X1234 <- readRDS("SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_0.3_1000MC.RDS")
  SimulationStudy3_HC3_X1234 <- readRDS("SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_0.3_1000MC.RDS")
    
  SimulationStudy3_no_HC_X123 <- readRDS("SimStudy3_workingModel=X123_no_HC_correction_allperspectives_R2_0.3_1000MC.RDS")
  SimulationStudy3_HC3_X123 <- readRDS("SimStudy3_workingModel=X123_HC3_marginal(not_eff)_R2_0.3_1000MC.RDS")
    
  SimulationStudy3_no_HC_X124 <- readRDS("SimStudy3_workingModel=X124_no_HC_correction_allperspectives_R2_0.3_1000MC.RDS")
  SimulationStudy3_HC3_X124 <- readRDS("SimStudy3_workingModel=X124_HC3_marginal(not_eff)_R2_0.3_1000MC.RDS")
  
SimulationStudy3 <- c(SimulationStudy3_no_HC_X1234, 
                        SimulationStudy3_HC3_X1234,
                        SimulationStudy3_no_HC_X123,
                        SimulationStudy3_HC3_X123, 
                        SimulationStudy3_no_HC_X124,
                        SimulationStudy3_HC3_X124)
                     
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

ggsave("plot_SimStudy3.png", 
       plot = plot_SimStudy3, 
       dpi = 900, width = 12.4, height = 10, units = "in")


results_all_conditions_df_SimStudy3_Appendix <- 
  results_all_conditions_df[results_all_conditions_df$perspective %in% c("Conditional"),]

plot_SimStudy3_Appendix <- ggplot(results_all_conditions_df_SimStudy3_Appendix, 
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
      ICM_R                       = "ICM[R]",
      fullCandidate               = "fullCandidate"
    )[x])
  )
plot_SimStudy3_Appendix

ggsave("plot_SimStudy3_Appendix.png", 
       plot = plot_SimStudy3_Appendix, 
       dpi = 900, width = 12.2, height = 6, units = "in")




### convergence issues
# initial model fitting

SimStudy3_convergence <- do.call(rbind, lapply(SimulationStudy3 ,
                                               function(x) x$convergence_singular_info)
)



colMeans(SimStudy3_convergence[sapply(SimStudy3_convergence, is.numeric)], na.rm = F)

0.0006666667/(0.0006666667+189.3731666667)  # less than 0.001%

library(dplyr)

View(SimStudy3_convergence %>%
  group_by(working_model) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))))

0.001/(0.001 + 199.3790) # less than 0.001%
0.0005/70.0790 # less than 0.001%

## monte carlo resampling process

plot_df <- results_all_conditions_df %>%
  filter(
    perspective %in% c("Marginal (eff)","Marginal (not eff) HC3","Marginal (not eff)"), 
    pval_type != "Naive"
  ) %>%
  mutate(
    estimator = factor(estimator),
    pval_type = factor(pval_type),
    perspective = factor(perspective),
    label_facet = factor(label_facet),
    coefficient = factor(coefficient)
  )

plot_MCResampling_SimStudy3 <-
  ggplot(
    plot_df,
    aes(
      x ="",
      y = total_resamplings_due_to_convergence,
      color = estimator,
      group = interaction(pval_type, estimator)
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.alpha = 0.25,
    fill = NA
  ) +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )+
  labs(
    x = "",
    y = "Number of MC Resamplings due to convergence",
    color = "Estimator"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    panel.spacing.y = unit(1.5, "lines"),  
    panel.spacing.x = unit(1, "lines")     
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
      ICM_R                       = "ICM[R]",
      fullCandidate               = "fullCandidate"
    )[x])
  )


plot_MCResampling_SimStudy3
ggsave("plot_MCResampling_SimStudy3.png", 
       plot = plot_MCResampling_SimStudy3, 
       dpi = 900, width = 9, height = 7.5, units = "in")





plot_df <- results_all_conditions_df %>%
  filter(
    perspective %in% c("Conditional"), 
    pval_type != "Naive"
  ) %>%
  mutate(
    estimator = factor(estimator),
    pval_type = factor(pval_type),
    perspective = factor(perspective),
    label_facet = factor(label_facet),
    coefficient = factor(coefficient)
  )

plot_MCResampling_SimStudy3_Appendix <-
  ggplot(
    plot_df,
    aes(
      x ="",
      y = total_resamplings_due_to_convergence,
      color = estimator,
      group = interaction(pval_type, estimator)
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.alpha = 0.25,
    fill = NA
  ) +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )+
  labs(
    x = "",
    y = "Number of MC Resamplings due to convergence",
    color = "Estimator"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    panel.spacing.y = unit(1.5, "lines"),  
    panel.spacing.x = unit(1, "lines")     
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
      ICM_R                       = "ICM[R]",
      fullCandidate               = "fullCandidate"
    )[x])
  )


plot_MCResampling_SimStudy3_Appendix
ggsave("plot_MCResampling_SimStudy3_Appendix.png", 
       plot = plot_MCResampling_SimStudy3_Appendix, 
       dpi = 900, width = 9, height = 7.5, units = "in")

     
