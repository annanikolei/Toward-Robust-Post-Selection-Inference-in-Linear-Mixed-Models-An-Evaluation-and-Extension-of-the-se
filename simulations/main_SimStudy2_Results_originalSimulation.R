## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

## Load Simulation Study 2 files
source("SimStudy2_helper.R")
source("SimStudy2_dgp.R")
source("SimStudy2_SelectionProcess.R")
source("SimStudy2_one_simulation.R")
source("SimStudy2_simulation_study.R")


R_squared = 0.3
n_obs = 40^2 ; m = 40

# load fixed predictor matriced
dgp_X_output <- readRDS("dgp_X_output_n_obs=1600_m=40.RDS")

# load list_supplied_D_work
load(paste0("SimulationResults_supplied_D/SIM2_list_supplied_D_work_full=X1X2X3X4_interesting=X3X4_n=400_Rsquared=",R_squared,"_nrsamples=10000.RData"))


########################
# Run Simulations
########################
setwd("SimStudy2_Results")

## All perspectives together 
N <- 1000

###########################
# no HC correction
###########################

# bayesian corrected 
save_file <- paste0(
  "SimStudy2_workingModel=X1234_bayesian_corrected_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian_corrected"),
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = T,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

save_file <- paste0(
  "SimStudy2_workingModel=X123_bayesian_corrected_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian_corrected"),
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = T,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

save_file <- paste0(
  "SimStudy2_workingModel=X124_bayesian_corrected_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian_corrected"),
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = T,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X4"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}


# bayesian conditional, marginal (eff + not eff)
save_file <- paste0(
  "SimStudy2_workingModel=X1234_HomoTrue_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian"),
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = T,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}


save_file <- paste0(
  "SimStudy2_workingModel=X123_HomoTrue_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian"),
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = T,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

save_file <- paste0(
  "SimStudy2_workingModel=X124_HomoTrue_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian"),
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = T,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X4"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}


################
# HC3 correction
################

save_file <- paste0(
  "SimStudy2_workingModel=X1234_HomoFalse_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 R_diag = F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = FALSE,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

save_file <- paste0(
  "SimStudy2_workingModel=X123_HomoFalse_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 R_diag = F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = FALSE,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

save_file <- paste0(
  "SimStudy2_workingModel=X124_HomoFalse_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy2_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 dgp_X_output = dgp_X_output,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 R_diag = F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 R_diag = FALSE,
                                 list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X4"), # c("X1","X2","X3","X4")
                                 checkConvergence_for_check_Fun = T
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

SimulationStudy2_no_HC_X1234 <- readRDS("SimStudy2_workingModel=X1234_HomoTrue_allperspectives_R2_0.3_1000MC.RDS")
SimulationStudy2_HC3_X1234 <- readRDS("SimStudy2_workingModel=X1234_HomoFalse_allperspectives_R2_0.3_1000MC.RDS")

SimulationStudy2_no_HC_X123 <- readRDS("SimStudy2_workingModel=X123_HomoTrue_allperspectives_R2_0.3_1000MC.RDS")
SimulationStudy2_HC3_X123 <- readRDS("SimStudy2_workingModel=X123_HomoFalse_allperspectives_R2_0.3_1000MC.RDS")

SimulationStudy2_no_HC_X124 <- readRDS("SimStudy2_workingModel=X124_HomoTrue_allperspectives_R2_0.3_1000MC.RDS")
SimulationStudy2_HC3_X124 <- readRDS("SimStudy2_workingModel=X124_HomoFalse_allperspectives_R2_0.3_1000MC.RDS")

SimulationStudy2_bayesian_corrected_X1234 <- readRDS("~/Downloads/selective-inference-lmm-simulation-main/simulations/SimStudy2_Results/SimStudy2_workingModel=X1234_bayesian_corrected_R2_0.3_1000MC.RDS")
SimulationStudy2_bayesian_corrected_X123 <- readRDS("~/Downloads/selective-inference-lmm-simulation-main/simulations/SimStudy2_Results/SimStudy2_workingModel=X123_bayesian_corrected_R2_0.3_1000MC.RDS")
SimulationStudy2_bayesian_corrected_X124 <- readRDS("~/Downloads/selective-inference-lmm-simulation-main/simulations/SimStudy2_Results/SimStudy2_workingModel=X124_bayesian_corrected_R2_0.3_1000MC.RDS")

SimulationStudy2 <- c(SimulationStudy2_no_HC_X1234, 
                      SimulationStudy2_HC3_X1234,
                      SimulationStudy2_no_HC_X123,
                      SimulationStudy2_HC3_X123, 
                      SimulationStudy2_no_HC_X124,
                      SimulationStudy2_HC3_X124)

SimStudy2_Results <- do.call(rbind, lapply(SimulationStudy2 ,
                                           function(x) x$post_selection_pvalues)
)

SimStudy2_convergence <- do.call(rbind, lapply(SimulationStudy2 ,
                                               function(x) x$convergence_singular_info)
)

SimulationStudy2_bayesian_corrected <- c(
  SimulationStudy2_bayesian_corrected_X1234,
  SimulationStudy2_bayesian_corrected_X123,
  SimulationStudy2_bayesian_corrected_X124)

SimStudy2_Results_bayesian_corrected <- do.call(rbind, lapply(SimulationStudy2_bayesian_corrected ,
                                           function(x) x$post_selection_pvalues)
)

SimStudy2_convergence_bayesian_corrected <- do.call(rbind, lapply(SimulationStudy2_bayesian_corrected ,
                                               function(x) x$convergence_singular_info)
)




saveRDS(SimStudy2_Results,"SimStudy2_pvalues_MCresamplings_allconditions_R2_0.3_1000MC.RDS")
saveRDS(SimStudy2_convergence, "SimStudy2_convergence_singular_info_allconditions_R2_0.3_1000MC.RDS")
saveRDS(SimStudy2_Results_bayesian_corrected,"SimStudy2_pvalues_bayesian_corrected_R2_0.3_1000MC.RDS")
saveRDS(SimStudy2_convergence_bayesian_corrected, "SimStudy2_convergence_singular_info_bayesian_corrected_R2_0.3_1000MC.RDS")


################################################################################
## QQ-Plots, Convergence ##
################################################################################
setwd("SimStudy2_Results")
SimStudy2_Results <- rbind(readRDS("SimStudy2_pvalues_bayesian_corrected_R2_0.3_1000MC.RDS"),
                           readRDS("SimStudy2_pvalues_MCresamplings_allconditions_R2_0.3_1000MC.RDS"))

SimStudy2_Results$pval_type <- ifelse(SimStudy2_Results$estimator%in% c("naive"),"Naive", ifelse(SimStudy2_Results$estimator%in% c("fullCandidate","fullCandidate_supplied"),"Standard (no selection)","Selfmade"))
SimStudy2_Results$estimator <- ifelse(SimStudy2_Results$estimator=="naive","ModelEstimate",SimStudy2_Results$estimator)
SimStudy2_Results$perspective <- ifelse(SimStudy2_Results$estimator %in% c("ModelEstimate_SIGMA_replaced","ModelEstimate_R_replaced",
                                                                           "ICM_SIGMA_replaced","ICM_R_replaced"),"marginal_not_efficient_HC3", SimStudy2_Results$perspective  )                                                                        

SimStudy2_Results$perspective <- ifelse(SimStudy2_Results$estimator %in% c("ModelEstimate","supplied","ICM","VarY") &  SimStudy2_Results$perspective=="conditional",
                                      paste0(SimStudy2_Results$perspective,"_uncorrected"), SimStudy2_Results$perspective )                                                                        

SimStudy2_Results$perspective <- ifelse(SimStudy2_Results$estimator %in% c("supplied_bayesian_corrected","ModelEstimate_bayesian_corrected","ICM_bayesian_corrected","VarY_bayesian_corrected") &  SimStudy2_Results$perspective=="conditional",
                                        paste0(SimStudy2_Results$perspective,"_corrected"), SimStudy2_Results$perspective )   

SimStudy2_Results$estimator <- sub("_bayesian_corrected$", "", SimStudy2_Results$estimator)


new_persp <- c(
  "marginal",
  "marginal_not_efficient",
  "marginal_not_efficient_HC3",
  "conditional_uncorrected",
  "conditional_corrected",
  "fullCandidate"
)
naive_rows <- SimStudy2_Results[SimStudy2_Results$perspective == "naive", ]
naive_expanded <- naive_rows[rep(seq_len(nrow(naive_rows)), each = length(new_persp)), ]
naive_expanded$perspective <- rep(new_persp, times = nrow(naive_rows))
SimStudy2_Results <- SimStudy2_Results[SimStudy2_Results$perspective != "naive", ]
SimStudy2_Results <- rbind(SimStudy2_Results, naive_expanded)

# supplied is the same if HC3 estimator is used => duplicate this
supplied_rows <- SimStudy2_Results[SimStudy2_Results$estimator == "supplied" & SimStudy2_Results$perspective=="marginal_not_efficient", ]
supplied_expanded <- supplied_rows[rep(seq_len(nrow(supplied_rows)), each = length("marginal_not_efficient_HC3")), ]
supplied_expanded$perspective <- rep("marginal_not_efficient_HC3", times = nrow(supplied_rows))

SimStudy2_Results <- rbind(SimStudy2_Results, supplied_expanded )


SimStudy2_Results$coefficient <-
  ifelse(SimStudy2_Results$coefficient == "X3",
         "L1-noise~x[3]",
         "L1-noise~(RS)~x[4]")


library(ggplot2)
library(dplyr)


results_all_conditions_df <- SimStudy2_Results  %>%
  mutate(
    label_facet = paste0(
      #"R² = ", R_squared_pseudo,
      "Selected Fixed = ", working_model
      
    )
  )

results_all_conditions_df$label_facet <- gsub(
  "Selected Fixed = X",
  "Sel. FE = X",
  results_all_conditions_df$label_facet
)

results_all_conditions_df$label_facet <- gsub(
  "X(?=[0-9])",
  "",
  results_all_conditions_df$label_facet,
  perl = TRUE
)

results_all_conditions_df$label_facet <- sub(
  "Sel\\. FE = ",
  "Sel. FE = X",
  results_all_conditions_df$label_facet
)




results_all_conditions_df$perspective <- factor(results_all_conditions_df$perspective)
results_all_conditions_df$perspective <- factor(
  results_all_conditions_df$perspective,
  levels = c(
    "conditional_uncorrected",
    "conditional_corrected",
    "marginal",
    "marginal_not_efficient",
    "marginal_not_efficient_HC3",
    "fullCandidate"
  ),
  labels = c(
    "Conditional\n(uncorrected)",
    "Conditional\n(corrected)",
    "Marginal\n(efficient)",
    "Marginal\n(not eff)",
    "Marginal\n(not eff) HC3",
    "fullCandidate"
  )
)

results_all_conditions_df$estimator <- factor(results_all_conditions_df$estimator)
results_all_conditions_df$estimator <- factor(
  results_all_conditions_df$estimator,
  levels = c("fullCandidate_supplied","supplied","ModelEstimate", "ICM", "VarY",
             "ModelEstimate_SIGMA_replaced","ModelEstimate_R_replaced",
             "ICM_SIGMA_replaced","ICM_R_replaced", 
             "fullCandidate"),
  labels = c("fullCandidate_supplied","supplied","ModelEstimate", "ICM","VarY",
             "ModelEstimate_SIGMA","ModelEstimate_R",
             "ICM_SIGMA","ICM_R",
             "fullCandidate"))

results_all_conditions_df$pval_type<-relevel(as.factor(results_all_conditions_df$pval_type),"Selfmade")
results_all_conditions_df$coefficient <- factor(
  results_all_conditions_df$coefficient,
  levels = c("L1-noise~x[3]", "L1-noise~(RS)~x[4]")
)

results_all_conditions_df <- 
  results_all_conditions_df[
    results_all_conditions_df$perspective %in% 
      c(
        "Conditional\n(uncorrected)",
        "Conditional\n(corrected)",
        "Marginal\n(efficient)",
        "Marginal\n(not eff) HC3",
        "Marginal\n(not eff)"
      ),
  ]
# -----------------------------
# Manual color palette
# -----------------------------
estimator_colors <- c(
  "fullCandidate_supplied"              = "black",
  "supplied"                            = "#009E73",  
  "ModelEstimate"                       = "#8C510A",  
  "ICM"                                 = "#D55E00",  
  "VarY"                                = "#0072B2", 
  "ModelEstimate_SIGMA"                 = "#56B4E9",  
  "ModelEstimate_R"                     = "#E69F00", 
  "ICM_SIGMA"                           = "#CC79A7",  
  "ICM_R"                               = "#6A3D9A",  
  "fullCandidate"                       = "white"
  
)

# -----------------------------
# Plot
# -----------------------------


results_all_conditions_df_SimStudy2 <- 
  results_all_conditions_df[results_all_conditions_df$perspective %in% c( "Conditional\n(uncorrected)",
                                                                          "Conditional\n(corrected)",
                                                                          "Marginal\n(efficient)",
                                                                          "Marginal\n(not eff)"),]


plot_SimStudy2 <- ggplot(results_all_conditions_df_SimStudy2, 
                         aes(
                           sample = pval_two,
                           color  = estimator,
                           linetype = pval_type)) +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  theme_minimal(base_size = 18) +
  theme(
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    strip.text = element_text(
      color = "black",
      size = 19,
      lineheight = 0.8
    ),
    legend.text = element_text(size = 19),
    legend.title = element_text(size = 19),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 19)
  )+
  geom_qq(  
    distribution = stats::qunif, size = 1, geom = "line") +
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
      supplied                         = "supplied",
      ModelEstimate                    = "ModelEstimate",
      ICM                               = "ICM",
      VarY                              ="VarY",
      ModelEstimate_SIGMA               = "ModelEstimate[Sigma]",
      ICM_SIGMA                         = "ICM[Sigma]",
      ModelEstimate_R                   = "ModelEstimate[R]",
      ICM_R                             = "ICM[R]"
    )[x])
  )

plot_SimStudy2

ggsave("plot_SimStudy2.png", 
       plot = plot_SimStudy2, 
       dpi = 900, width = 16, height = 9, units = "in")


# results_all_conditions_df_SimStudy2_Appendix <- 
#   results_all_conditions_df[results_all_conditions_df$perspective %in% c("Marginal\n(not eff) HC3"),]
# 
# plot_SimStudy2_Appendix <- ggplot(results_all_conditions_df_SimStudy2_Appendix, 
#                                   aes(
#                                     sample = pval_two,
#                                     color  = estimator,
#                                     linetype = pval_type)) +
#   geom_abline(slope = 1, intercept = 0, color = "grey") +
#   theme_minimal(base_size = 18) +  # Basis-Schriftgröße
#   theme(
#     strip.background = element_rect(fill = "lightgrey", color = NA),  
#     strip.text = element_text(color = "black"),
#     legend.position = "bottom",       # <--- Hier Legende unten setzen
#     legend.box = "vertical"
#   ) +
#   geom_qq(  
#     distribution = stats::qunif, size = 0.8, geom = "line") +
#   facet_grid(
#     perspective ~ label_facet + coefficient,
#     labeller = labeller(
#       coefficient = label_parsed,
#       label_facet = label_wrap_gen(width = 35)
#     )
#   )+
#   labs(
#     title = "",
#     x = "Expected p-values",
#     y = "Observed p-values",
#     color = "Estimator",
#     linetype = "p-value type"
#   ) +
#   guides(
#     color = guide_legend(order = 1),
#     linetype = guide_legend(order = 2)
#   ) +
#   scale_color_manual(
#     values = estimator_colors,
#     labels = function(x) parse(text = c(
#       fullCandidate_supplied = "fullCandidate (supplied)",
#       supplied                 = "supplied",
#       ModelEstimate               = "ModelEstimate",
#       ICM                         = "ICM",
#       VarY                        = "VarY",
#       ModelEstimate_SIGMA         = "ModelEstimate[Sigma]",
#       ModelEstimate_R             = "ModelEstimate[R]",
#       ICM_SIGMA                   = "ICM[Sigma]",
#       ICM_R                       = "ICM[R]",
#       fullCandidate               = "fullCandidate"
#     )[x])
#   )
# plot_SimStudy2_Appendix
# 
# ggsave("plot_SimStudy2_Appendix.png", 
#        plot = plot_SimStudy2_Appendix, 
#        dpi = 900, width = 13, height = 6, units = "in")
# 
# 


### convergence issues
# initial model fitting

SimStudy2_convergence <- rbind(readRDS("SimStudy2_convergence_singular_info_allconditions_R2_0.3_1000MC.RDS"),
                           readRDS("SimStudy2_convergence_singular_info_bayesian_corrected_R2_0.3_1000MC.RDS"))

colMeans(SimStudy2_convergence[sapply(SimStudy2_convergence, is.numeric)], na.rm = F)

13.1075833/(13.1075833+66.7419167 ) #  0.165

library(dplyr)

View(SimStudy2_convergence %>%
       group_by(working_model) %>%
       summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))))

36.15550/(36.15550+
            182.4002
) # [1] 0.165
2.143000/(2.143000+
          12.077500)# [1] 0.15
1.02425/(1.02425+
           5.7480)#[1] 0.15


0.2020000/
3.59775
# [1] 0.056
1.20675/
20.543 #[1]0.05874264
0.10050/
  1.71175#[1] 0.05871184



## monte carlo resampling process
plot_df <- results_all_conditions_df %>%
  filter(
    perspective %in% c("Conditional\n(uncorrected)",
                       "Conditional\n(corrected)",
                       "Marginal\n(efficient)",
                       "Marginal\n(not eff)"), 
    pval_type != "Naive"
  ) %>%
  mutate(
    estimator = factor(estimator),
    pval_type = factor(pval_type),
    perspective = factor(perspective),
    label_facet = factor(label_facet),
    coefficient = factor(coefficient)
  )

plot_MCResampling_SimStudy2 <-
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
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    strip.text = element_text(
      color = "black",
      size = 19,
      lineheight = 0.8
    ),
    legend.text = element_text(size = 19),
    legend.title = element_text(size = 19),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 19)
  )+
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


plot_MCResampling_SimStudy2
ggsave("plot_MCResampling_SimStudy2.png", 
       plot = plot_MCResampling_SimStudy2, 
       dpi = 900, width = 15, height = 9, units = "in")





# plot_df <- results_all_conditions_df %>%
#   filter(
#     perspective %in% c("Marginal\n(not eff) HC3"), 
#     pval_type != "Naive"
#   ) %>%
#   mutate(
#     estimator = factor(estimator),
#     pval_type = factor(pval_type),
#     perspective = factor(perspective),
#     label_facet = factor(label_facet),
#     coefficient = factor(coefficient)
#   )
# 
# plot_MCResampling_SimStudy2_Appendix <-
#   ggplot(
#     plot_df,
#     aes(
#       x ="",
#       y = total_resamplings_due_to_convergence,
#       color = estimator,
#       group = interaction(pval_type, estimator)
#     )
#   ) +
#   geom_boxplot(
#     position = position_dodge(width = 0.8),
#     width = 0.7,
#     outlier.alpha = 0.25,
#     fill = NA
#   ) +
#   facet_grid(
#     perspective ~ label_facet + coefficient,
#     labeller = labeller(
#       coefficient = label_parsed,
#       label_facet = label_wrap_gen(width = 35)
#     )
#   )+
#   labs(
#     x = "",
#     y = "Number of MC Resamplings due to convergence",
#     color = "Estimator"
#   ) +
#   theme_minimal(base_size = 18) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.background = element_rect(fill = "lightgrey", color = NA),
#     strip.text = element_text(color = "black", size = 14),
#     legend.position = "bottom",
#     panel.spacing.y = unit(1.5, "lines"),  
#     panel.spacing.x = unit(1, "lines")     
#   ) +
#   scale_color_manual(
#     values = estimator_colors,
#     labels = function(x) parse(text = c(
#       fullCandidate_supplied = "fullCandidate (supplied)",
#       supplied                 = "supplied",
#       ModelEstimate               = "ModelEstimate",
#       ICM                         = "ICM",
#       VarY                        = "VarY",
#       ModelEstimate_SIGMA         = "ModelEstimate[Sigma]",
#       ModelEstimate_R             = "ModelEstimate[R]",
#       ICM_SIGMA                   = "ICM[Sigma]",
#       ICM_R                       = "ICM[R]",
#       fullCandidate               = "fullCandidate"
#     )[x])
#   )
# 
# 
# plot_MCResampling_SimStudy2_Appendix
# ggsave("plot_MCResampling_SimStudy2_Appendix.png", 
#        plot = plot_MCResampling_SimStudy2_Appendix, 
#        dpi = 900, width = 13, height = 6, units = "in")
# 


