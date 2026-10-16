################################################################################
# Simulation Study 3: Results
################################################################################

message("Start this script from the repository root.")

## Install and load necessary packages
source("simulations/Project1_install_packages.R")

setwd("Selfmade")
source("selfmade**.R")

setwd("../simulations")

## Load Simulation Study 3 files
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


RUN_SIMULATIONS <- FALSE 
results_dir <- "SimStudy3_Results"

if (RUN_SIMULATIONS && !dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

if (!dir.exists(results_dir)) {
  stop("Results folder not found. Extract SimStudy3_Results.zip first.")
}

setwd(results_dir)

if (RUN_SIMULATIONS) {
########################
# Run Simulations
########################
#######################
# CR2 correction
#######################

# efficient marginal perspective CR2 X1234
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X1234_marginal(eff)CR2_R2_", 0.3, "_1000MC.RDS"
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
                                 perspective = c("marginal"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
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

## eff marginal perspectives CR2 X123
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X123_marginal(eff)_CR2_R2_", 0.3, "_1000MC.RDS"
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
                                 perspective = c("marginal"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag = TRUE,
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

## eff marginal perspectives CR2 X124
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X124_marginal(eff)_CR2_R2_", 0.3, "_1000MC.RDS"
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
                                 perspective = c("marginal"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag = TRUE,
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

## not-efficient marginal perspective CR2 X1234
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X1234_marginal(not_eff)CR2_R2_", 0.3, "_1000MC.RDS"
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
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag = TRUE,
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


## not-eff marginal perspectives CR2 X123
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X123_marginal(not_eff)_CR2_R2_", 0.3, "_1000MC.RDS"
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
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag = TRUE,
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

## not eff marginal perspectives CR2 X124
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X124_marginal(not_eff)_CR2_R2_", 0.3, "_1000MC.RDS"
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
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag = TRUE,
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
# no_HC_correction
########################

## conditional bayesian corrected
#### X1234
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X1234_bayesian_corrected_R2_", 0.3, "_1000MC.RDS"
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
                                 perspective = c("conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian_corrected"),
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

## conditional bayesian corrected
#### X123
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X123_bayesian_corrected_R2_", 0.3, "_1000MC.RDS"
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
                                 perspective = c("conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian_corrected"),
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

## conditional bayesian corrected
#### X124
N <- 1000

save_file <- paste0(
  "SimStudy3_workingModel=X124_bayesian_corrected_R2_", 0.3, "_1000MC.RDS"
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
                                 perspective = c("conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 bayesian= c("bayesian_corrected"),
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


## all perspectives no_HC_correction
#### X1234
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
                                 bayesian= c("bayesian_corrected"),
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
                                 bayesian= c("bayesian"),
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
                                 bayesian= c("bayesian"),
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

}
################################################################################
## QQ-Plots, Convergence ##
################################################################################


SimulationStudy3_no_HC_X1234 <- readRDS("SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_0.3_1000MC.RDS")
SimulationStudy3_HC3_X1234 <- readRDS("SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_0.3_1000MC.RDS")

SimulationStudy3_no_HC_X123 <- readRDS("SimStudy3_workingModel=X123_no_HC_correction_allperspectives_R2_0.3_1000MC.RDS")
SimulationStudy3_HC3_X123 <- readRDS("SimStudy3_workingModel=X123_HC3_marginal(not_eff)_R2_0.3_1000MC.RDS")

SimulationStudy3_no_HC_X124 <- readRDS("SimStudy3_workingModel=X124_no_HC_correction_allperspectives_R2_0.3_1000MC.RDS")
SimulationStudy3_HC3_X124 <- readRDS("SimStudy3_workingModel=X124_HC3_marginal(not_eff)_R2_0.3_1000MC.RDS")

SimulationStudy3_bayesian_corrected_X1234 <- readRDS("SimStudy3_workingModel=X1234_bayesian_corrected_R2_0.3_1000MC.RDS")
SimulationStudy3_bayesian_corrected_X123 <- readRDS("SimStudy3_workingModel=X123_bayesian_corrected_R2_0.3_1000MC.RDS")
SimulationStudy3_bayesian_corrected_X124 <- readRDS("SimStudy3_workingModel=X124_bayesian_corrected_R2_0.3_1000MC.RDS")


SimulationStudy3 <- c(SimulationStudy3_bayesian_corrected_X1234,
                      SimulationStudy3_bayesian_corrected_X123,
                      SimulationStudy3_bayesian_corrected_X124,
                      SimulationStudy3_no_HC_X1234, 
                      SimulationStudy3_HC3_X1234,
                      SimulationStudy3_no_HC_X123,
                      SimulationStudy3_HC3_X123, 
                      SimulationStudy3_no_HC_X124,
                      SimulationStudy3_HC3_X124)

SimStudy3_Results <- do.call(rbind, lapply(SimulationStudy3 ,
                                           function(x) x$post_selection_pvalues)
)


SimulationStudy3_marginaleff_CR2_X1234 <- readRDS("SimStudy3_workingModel=X1234_marginal(eff)CR2_R2_0.3_1000MC.RDS")
SimulationStudy3_marginaleff_CR2_X123 <- readRDS("SimStudy3_workingModel=X123_marginal(eff)_CR2_R2_0.3_1000MC.RDS")
SimulationStudy3_marginaleff_CR2_X124 <- readRDS("SimStudy3_workingModel=X124_marginal(eff)_CR2_R2_0.3_1000MC.RDS")
SimulationStudy3_marginal_noteff_CR2_X1234 <- readRDS("SimStudy3_workingModel=X1234_marginal(not_eff)CR2_R2_0.3_1000MC.RDS")
SimulationStudy3_marginal_noteff_CR2_X123 <- readRDS("SimStudy3_workingModel=X123_marginal(not_eff)_CR2_R2_0.3_1000MC.RDS")
SimulationStudy3_marginal_noteff_CR2_X124 <- readRDS("SimStudy3_workingModel=X124_marginal(not_eff)_CR2_R2_0.3_1000MC.RDS")

SimulationStudy3_marginalCR2 <- c(SimulationStudy3_marginaleff_CR2_X1234,
                                  SimulationStudy3_marginaleff_CR2_X123,
                                  SimulationStudy3_marginaleff_CR2_X124,
                                  SimulationStudy3_marginal_noteff_CR2_X1234,
                                  SimulationStudy3_marginal_noteff_CR2_X123,
                                  SimulationStudy3_marginal_noteff_CR2_X124)
  
SimStudy3_marginalCR2_Results <- do.call(rbind, lapply(SimulationStudy3_marginalCR2 ,
                                                       function(x) x$post_selection_pvalues)
)
SimStudy3_marginalCR2_Results$perspective <- ifelse(SimStudy3_marginalCR2_Results$perspective %in% c("marginal","marginal_not_efficient","naive"),
                                                    paste0(SimStudy3_marginalCR2_Results$perspective,"_CR2"), SimStudy3_marginalCR2_Results$perspective )    
SimStudy3_marginalCR2_Results <- SimStudy3_marginalCR2_Results[SimStudy3_marginalCR2_Results$estimator != "supplied",]


SimStudy3_Results <- rbind(SimStudy3_marginalCR2_Results,SimStudy3_Results)

SimStudy3_Results$pval_type <- ifelse(SimStudy3_Results$estimator%in% c("naive"),"Naive", ifelse(SimStudy3_Results$estimator%in% c("fullCandidate","fullCandidate_supplied"),"Standard (no selection)","Selfmade"))
SimStudy3_Results$estimator <- ifelse(SimStudy3_Results$estimator=="naive","ModelEstimate",SimStudy3_Results$estimator)
SimStudy3_Results$perspective <- ifelse(SimStudy3_Results$estimator %in% c("ModelEstimate_SIGMA_replaced","ModelEstimate_R_replaced",
                                                                           "ICM_SIGMA_replaced","ICM_R_replaced"),"marginal_not_efficient_HC3", SimStudy3_Results$perspective  )                                                                        

SimStudy3_Results$perspective <- ifelse(SimStudy3_Results$estimator %in% c("ModelEstimate","supplied","ICM","VarY") &  SimStudy3_Results$perspective=="conditional",
                                        paste0(SimStudy3_Results$perspective,"_uncorrected"), SimStudy3_Results$perspective )                                                                        

SimStudy3_Results$perspective <- ifelse(SimStudy3_Results$estimator %in% c("supplied_bayesian_corrected","ModelEstimate_bayesian_corrected","ICM_bayesian_corrected","VarY_bayesian_corrected") &  SimStudy3_Results$perspective=="conditional",
                                        paste0(SimStudy3_Results$perspective,"_corrected"), SimStudy3_Results$perspective )  

SimStudy3_Results$estimator <- sub("_bayesian_corrected$", "", SimStudy3_Results$estimator)

new_persp <- c(
  "marginal",
  "marginal_not_efficient",
  "marginal_CR2",
  "marginal_not_efficient_CR2",
  "marginal_not_efficient_HC3",
  "conditional_uncorrected",
  "conditional_corrected",
  "fullCandidate"
)
naive_rows <- SimStudy3_Results[SimStudy3_Results$perspective == "naive", ]
naive_expanded <- naive_rows[rep(seq_len(nrow(naive_rows)), each = length(new_persp)), ]
naive_expanded$perspective <- rep(new_persp, times = nrow(naive_rows))
SimStudy3_Results <- SimStudy3_Results[SimStudy3_Results$perspective != "naive", ]
SimStudy3_Results <- rbind(SimStudy3_Results, naive_expanded)

# supplied is the same if HC3/CR2 estimator is used => duplicate this
# supplied results from marginal_not_efficient
supplied_rows <- SimStudy3_Results[
  SimStudy3_Results$estimator == "supplied" &
    SimStudy3_Results$perspective == "marginal_not_efficient", 
]

# Copy as marginal_not_efficient_HC3
supplied_expanded <- supplied_rows
supplied_expanded$perspective <- "marginal_not_efficient_HC3"

# Copy as marginal_not_efficient_CR2
supplied_expanded_CR2_noteff <- supplied_rows
supplied_expanded_CR2_noteff$perspective <- "marginal_not_efficient_CR2"


# supplied results from marginal
supplied_rows_eff <- SimStudy3_Results[
  SimStudy3_Results$estimator == "supplied" &
    SimStudy3_Results$perspective == "marginal", 
]

# Copy as marginal_CR2
supplied_expanded_CR2_eff <- supplied_rows_eff
supplied_expanded_CR2_eff$perspective <- "marginal_CR2"


# Add all copies to dataset
SimStudy3_Results <- rbind(
  SimStudy3_Results,
  supplied_expanded,
  supplied_expanded_CR2_eff,
  supplied_expanded_CR2_noteff
)
SimStudy3_Results$coefficient <-
  ifelse(SimStudy3_Results$coefficient == "X3",
         "L1-noise~x[3]",
         "L1-noise~(RS)~x[4]")


library(ggplot2)
library(dplyr)


results_all_conditions_df <- SimStudy3_Results  %>%
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
    "marginal_CR2",
    "marginal_not_efficient_CR2",
    "marginal_not_efficient_HC3",
    "fullCandidate"
  ),
  labels = c(
    "Conditional\n(uncorrected)",
    "Conditional\n(corrected)",
    "Marginal\n(efficient)",
    "Marginal\n(not eff)",
    "Marginal\n(eff) CR2 ", 
    "Marginal\n(not eff) CR2 ",
    "Marginal\n(not eff) HC3",
    "fullCandidate"
  )
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
             "fullCandidate"))

results_all_conditions_df$pval_type<-relevel(as.factor(results_all_conditions_df$pval_type),"Selfmade")
results_all_conditions_df$coefficient <- factor(
  results_all_conditions_df$coefficient,
  levels = c( "L1-noise~x[3]",
              "L1-noise~(RS)~x[4]")
)


# -----------------------------
# Manual color palette
# -----------------------------

library(ggplot2)
library(cowplot)

# -----------------------------
# Estimator colors
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
# Data
# -----------------------------

results_all_conditions_df_SimStudy3 <- 
  results_all_conditions_df[
    results_all_conditions_df$perspective %in% 
      c(
        "Marginal\n(efficient)",
        "Marginal\n(not eff)",
        "Marginal\n(eff) CR2 ", 
        "Marginal\n(not eff) CR2 ",
        "Marginal\n(not eff) HC3",
        "Conditional\n(corrected)"
      ),
  ]


# -----------------------------
# Base plot
# -----------------------------

plot_SimStudy3 <- ggplot(
  results_all_conditions_df_SimStudy3,
  aes(
    sample = pval_two,
    color = estimator,
    linetype = pval_type
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "grey"
  ) +
  geom_qq(
    distribution = stats::qunif,
    size = 1.2,
    geom = "line"
  ) +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
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
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
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

plot_SimStudy3

ggsave("plot_SimStudy3.png", 
       plot = plot_SimStudy3, 
       dpi = 300, width = 16, height = 12, units = "in")



SimStudy3_convergence <- do.call(rbind, lapply(c(SimulationStudy3,SimulationStudy3_marginalCR2) ,
                                               function(x) x$convergence_singular_info)
)



colMeans(SimStudy3_convergence[sapply(SimStudy3_convergence, is.numeric)], na.rm = F)

5.331911e-04/( 5.331911e-04 +1.901220e+02  )  # less than 0.001%

library(dplyr)

View(SimStudy3_convergence %>%
       group_by(working_model) %>%
       summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))))

0.0004/(0.0004+ 201) 
0.0011995202/(0.0011995202+365.6619)	# less than 0.001%
0.0004/
70 # less than 0.001%

## monte carlo resampling process

plot_df <- results_all_conditions_df_SimStudy3 %>%
  filter(
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


plot_MCResampling_SimStudy3
ggsave("plot_MCResampling_SimStudy3.png", 
       plot = plot_MCResampling_SimStudy3, 
       dpi = 300, width = 16, height = 12, units = "in")



# bayesian conditional (uncorrected)
results_all_conditions_df_SimStudy3_Appendix <- 
  results_all_conditions_df[
    results_all_conditions_df$perspective %in% 
      c(
         "Conditional\n(uncorrected)"
      ),
  ]


# -----------------------------
# Base plot
# -----------------------------

plot_SimStudy3_Appendix <- ggplot(
  results_all_conditions_df_SimStudy3_Appendix,
  aes(
    sample = pval_two,
    color = estimator,
    linetype = pval_type
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "grey"
  ) +
  geom_qq(
    distribution = stats::qunif,
    size = 1.2,
    geom = "line"
  ) +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
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
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
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

plot_SimStudy3_Appendix

ggsave("plot_SimStudy3_Appendix.png", 
       plot = plot_SimStudy3_Appendix, 
       dpi = 300, width = 16, height = 4, units = "in")




plot_df <- results_all_conditions_df_SimStudy3_Appendix %>%
  filter(
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


plot_MCResampling_SimStudy3_Appendix
ggsave("plot_MCResampling_SimStudy3_Appendix.png", 
       plot = plot_MCResampling_SimStudy3_Appendix, 
       dpi = 300, width =16, height = 6.5, units = "in")











     
