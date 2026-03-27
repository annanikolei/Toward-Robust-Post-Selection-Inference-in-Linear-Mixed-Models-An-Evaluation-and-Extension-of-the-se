setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")
source("SimStudy3_helper.R")
source("SimStudy3_dgp.R")
source("SimStudy3_SelectionProcess.R")
source("SimStudy3_one_simulation.R")
source("SimStudy3_simulation_study.R")

# load list with supplied random effects covariance matrices of all working models
# that are possibly interesting for Simulation Study 2: list_supplied_D_work
# with D = matrix(c(2,0.2,0.1,
#                   0.2,1.5,0.2,
#                   0.1,0.2,1.5),nrow=3, byrow=T)
R_squared = 0.3
load(paste0("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimulationResults_supplied_D/SIM3_list_supplied_D_work_full=X1X2X3X4_interesting=X3X4_n=400_Rsquared=",R_squared,"_nrsamples=10000.RData"))

n_obs = 40^2 ; m = 40
dgp_X_output <- readRDS(paste0("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/dgp_X_output_n_obs=",n_obs,"_m=",m,".RDS"))

########################
# normal standard errors
########################

## all perspectives no_HC_correction
#### X1234
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_new")
N <- 1000

# Speicherort
save_file <- paste0(
  "SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

# Vorhandene Datei laden, wenn vorhanden (weiterrechnen statt überschreiben)
if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Weiter bei Iteration: ", start_iter)
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
      warning("Fehler in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

message("FERTIG!")



#### X123
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_new")
N <- 1000

# Speicherort
save_file <- paste0(
  "SimStudy3_workingModel=X123_no_HC_correction_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

# Vorhandene Datei laden, wenn vorhanden (weiterrechnen statt überschreiben)
if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Weiter bei Iteration: ", start_iter)
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
      warning("Fehler in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

message("FERTIG!")


#### X124
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_new")
N <- 1000

# Speicherort
save_file <- paste0(
  "SimStudy3_workingModel=X124_no_HC_correction_allperspectives_R2_", 0.3, "_1000MC.RDS"
)

# Vorhandene Datei laden, wenn vorhanden (weiterrechnen statt überschreiben)
if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Weiter bei Iteration: ", start_iter)
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
      warning("Fehler in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  # Ergebnis speichern, NULL wird auch abgelegt um nichts zu verlieren
  results_list[[i]] <- res
  
  # Jede Iteration sofort sichern
  saveRDS(results_list, save_file)
}

message("FERTIG!")










SimulationStudy3_X123 <-
  #suppressMessages(
    SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                               n_obs = n_obs, # number of observations
                               m = m, # number of level-2 units
                               R_squared = R_squared, # pseudo coefficient of determination
                               dgp_X_output = dgp_X_output,
                               D = matrix(c(2,0.2,0.1,
                                            0.2,1.5,0.2,
                                            0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                               perspective = c("marginal","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                               robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                               B_MCsamples = 1000, # number of MC samples in selective inference process
                               remove_singular_fits = F,
                               list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                               interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                               interesting_working_model_vars = c("X1","X2","X3"),
                               checkConvergence_for_check_Fun = T # convergence check for MC sampling?
    ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(SimulationStudy3_X123, 
       file=paste0("2_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")


SimulationStudy3_X124 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(SimulationStudy3_X124, 
        file=paste0("Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")


SimulationStudy3_X124_2 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(SimulationStudy3_X124_2, 
        file=paste0("2_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")


SimulationStudy3_X1234 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X3","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(SimulationStudy3_X1234, 
        file=paste0("2_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")

# hist(SimulationStudy3$estimates_df[SimulationStudy3$estimates_df$predictor=="X3",1])
# mean(SimulationStudy3$estimates_df[SimulationStudy3$estimates_df$predictor=="X3",1]) #-0.004400027


# marginal not efficient -------------------------------------------------------
not_efficient_2_RdiagT_SimulationStudy3_X123 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             R_diag =T,
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X3"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(not_efficient_2_RdiagT_SimulationStudy3_X123, 
        file=paste0("not_efficient_2_RdiagT_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")


not_efficient_2_RdiagF_SimulationStudy3_X123 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
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
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X3"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(not_efficient_2_RdiagF_SimulationStudy3_X123, 
        file=paste0("not_efficient_2_RdiagF_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")




not_efficient_2_RdiagT_SimulationStudy3_X124 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             R_diag = T,
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(not_efficient_2_RdiagT_SimulationStudy3_X124, 
        file=paste0("not_efficient_2_RdiagT_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")

not_efficient_2_RdiagF_SimulationStudy3_X124 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             R_diag = F,
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(not_efficient_2_RdiagF_SimulationStudy3_X124, 
        file=paste0("not_efficient_2_RdiagF_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")




not_efficient_2_RdiagT_SimulationStudy3_X1234 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             R_diag = T,
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X3","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(not_efficient_2_RdiagT_SimulationStudy3_X1234, 
        file=paste0("not_efficient_2_RdiagT_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")

not_efficient_2_RdiagF_SimulationStudy3_X1234 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             R_diag = F,
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X3","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(not_efficient_2_RdiagF_SimulationStudy3_X1234, 
        file=paste0("not_efficient_2_RdiagF_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")


## HC3 standard errors
HC3not_efficient_RdiagF_SimulationStudy3_X123 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 1000, #number of replications per condition
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
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(HC3not_efficient_RdiagF_SimulationStudy3_X123, 
        file=paste0("HC3not_efficient_RdiagF_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_1000",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")


HC3not_efficient_RdiagF_SimulationStudy3_X124 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 1000, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                             cadjust = T,
                             R_diag = F,
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(HC3not_efficient_RdiagF_SimulationStudy3_X124, 
        file=paste0("HC3not_efficient_RdiagF_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_1000",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")


HC3not_efficient_RdiagF_SimulationStudy3_X1234 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                             cadjust = T,
                             R_diag = F,
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X3","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(HC3not_efficient_RdiagF_SimulationStudy3_X1234, 
        file=paste0("HC3not_efficient_RdiagF_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")

HC3not_efficient_2_RdiagF_SimulationStudy3_X1234 <-
  #suppressMessages(
  SimStudy3_simulation_study(n_iter = 500, #number of replications per condition
                             n_obs = n_obs, # number of observations
                             m = m, # number of level-2 units
                             R_squared = R_squared, # pseudo coefficient of determination
                             dgp_X_output = dgp_X_output,
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                             robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                             HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                             cadjust = T,
                             R_diag = F,
                             B_MCsamples = 1000, # number of MC samples in selective inference process
                             remove_singular_fits = F,
                             list_supplied_D_work = list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
                             interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                             interesting_working_model_vars = c("X1","X2","X3","X4"),
                             checkConvergence_for_check_Fun = T # convergence check for MC sampling?
  ) 
#)

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X")
saveRDS(HC3not_efficient_2_RdiagF_SimulationStudy3_X1234, 
        file=paste0("HC3not_efficient_2_RdiagF_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=",R_squared,"_marg+cond_notRobust_MCconvergencecheck_500",".RDS"))
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified")


# gemeinsame Auswertung mit QQ-Plots -------------------------------------------

SimulationStudy3_notRobust_X1X2X4 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_notRobust_X1X2X3X4 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS") 
SimulationStudy3_notRobust_X1X2X3 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_notRobust_X1X2X4_2 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/2_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_notRobust_X1X2X3X4_2 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/2_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_notRobust_X1X2X3_2 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/2_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")

HC0SimulationStudy3_not_efficient_RdiagF_X1X2X4 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_RdiagF_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_not_efficient_RdiagT_X1X2X4 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_RdiagT_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
HC0SimulationStudy3_not_efficient_RdiagF_X1X2X3 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_RdiagF_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_not_efficient_RdiagT_X1X2X3 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_RdiagT_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
HC0SimulationStudy3_not_efficient_RdiagF_X1X2X3X4 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_RdiagF_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_not_efficient_RdiagT_X1X2X3X4 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_RdiagT_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")

HC0SimulationStudy3_not_efficient_RdiagF_X1X2X4_2  <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_2_RdiagF_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_not_efficient_RdiagT_X1X2X4_2  <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_2_RdiagT_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
HC0SimulationStudy3_not_efficient_RdiagF_X1X2X3_2  <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_2_RdiagF_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_not_efficient_RdiagT_X1X2X3_2  <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_2_RdiagT_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
HC0SimulationStudy3_not_efficient_RdiagF_X1X2X3X4_2  <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_2_RdiagF_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
SimulationStudy3_not_efficient_RdiagT_X1X2X3X4_2  <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/not_efficient_2_RdiagT_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")

HC3SimulationStudy3_not_efficient_RdiagF_X1X2X4 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/HC3not_efficient_RdiagF_Xfixed_X1X2X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_1000.RDS")
HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/HC3not_efficient_RdiagF_Xfixed_X1X2X3_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_1000.RDS")
HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3X4<- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/HC3not_efficient_RdiagF_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")
HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3X4_2 <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fixed_X/HC3not_efficient_2_RdiagF_Xfixed_X1X2X3X4_SimulationStudy3_n=40_R_squared=0.3_marg+cond_notRobust_MCconvergencecheck_500.RDS")

SimulationStudy3_fullCandidateTruth <- readRDS("~/Documents/Forschung/Post-selection inference/Projekt 1 /Simulation Study 3 - random effects structure misspecified/SimStudy3_Results_fullCandidatesupplied/SimulationStudy3_fullCandidate_RdiagTRUE_n=40_R2_0.3_1000MC.RDS")
SimulationStudy3_fullCandidate <- do.call(rbind, lapply(SimulationStudy3_fullCandidateTruth  ,
                      function(x) x$post_selection_pvalues)
)

SimulationStudy3_fullCandidate_info <- do.call(rbind, lapply(SimulationStudy3_fullCandidateTruth  ,
                                                        function(x) x$convergence_singular_info)
)

SimulationStudy3_fullCandidate$working_model <-  SimulationStudy3_fullCandidate_info$working_model<- "X1X2X3X4"
SimulationStudy3_fullCandidate$estimator <-  SimulationStudy3_fullCandidate_info$estimator <-"fullCandidate_supplied"


library(dplyr)


SimStudy3_Results <- bind_rows(HC3SimulationStudy3_not_efficient_RdiagF_X1X2X4$post_selection_pvalues,
                           HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3$post_selection_pvalues,
                           HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3X4$post_selection_pvalues,
                           # HC0SimulationStudy3_not_efficient_RdiagF_X1X2X4$post_selection_pvalues,
                           # HC0SimulationStudy3_not_efficient_RdiagF_X1X2X3$post_selection_pvalues,
                           # HC0SimulationStudy3_not_efficient_RdiagF_X1X2X3X4$post_selection_pvalues,
                           SimulationStudy3_not_efficient_RdiagT_X1X2X4$post_selection_pvalues,
                           SimulationStudy3_not_efficient_RdiagT_X1X2X3$post_selection_pvalues,
                           SimulationStudy3_not_efficient_RdiagT_X1X2X3X4$post_selection_pvalues,
                           
                           # HC0SimulationStudy3_not_efficient_RdiagF_X1X2X4_2$post_selection_pvalues,
                           # HC0SimulationStudy3_not_efficient_RdiagF_X1X2X3_2$post_selection_pvalues,
                           # HC0SimulationStudy3_not_efficient_RdiagF_X1X2X3X4_2$post_selection_pvalues,
                           HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3X4_2$post_selection_pvalues,
                           SimulationStudy3_not_efficient_RdiagT_X1X2X4_2$post_selection_pvalues,
                           SimulationStudy3_not_efficient_RdiagT_X1X2X3_2$post_selection_pvalues,
                           SimulationStudy3_not_efficient_RdiagT_X1X2X3X4_2$post_selection_pvalues,
                           
                           SimulationStudy3_notRobust_X1X2X4$post_selection_pvalues,
                           SimulationStudy3_notRobust_X1X2X3$post_selection_pvalues,
                           SimulationStudy3_notRobust_X1X2X3X4$post_selection_pvalues,
                           SimulationStudy3_notRobust_X1X2X4_2$post_selection_pvalues,
                           SimulationStudy3_notRobust_X1X2X3_2$post_selection_pvalues,
                           SimulationStudy3_notRobust_X1X2X3X4_2$post_selection_pvalues,
                           SimulationStudy3_fullCandidate)

# table(SimulationStudy3_notRobust$post_selection_pvalues$working_model)
table(SimStudy3_Results$working_model)

SimStudy3_Results$pval_type <- ifelse(SimStudy3_Results$estimator%in% c("naive"),"Naive", ifelse(SimStudy3_Results$estimator%in% c("fullCandidate","fullCandidate_supplied"),"Standard (no selection)","Post-Selective"))

SimStudy3_Results$estimator <- ifelse(SimStudy3_Results$estimator=="naive","ModelEstimate",SimStudy3_Results$estimator)
SimStudy3_Results$perspective <- ifelse(SimStudy3_Results$estimator %in% c("ModelEstimate_SIGMA_replaced","ModelEstimate_R_replaced",
                                                                           "ICM_SIGMA_replaced","ICM_R_replaced"),"marginal_not_efficient_HC3", SimStudy3_Results$perspective  )                                                                        

# Perspektiven, die naive ersetzen sollen
new_persp <- c(
  "marginal",
  "marginal_not_efficient",
  "marginal_not_efficient_HC3",
  "conditional",
  "fullCandidate"
)

# naive-Zeilen auswählen
naive_rows <- SimStudy3_Results[SimStudy3_Results$perspective == "naive", ]

# vervierfachen
naive_expanded <- naive_rows[rep(seq_len(nrow(naive_rows)), each = length(new_persp)), ]

# perspective ersetzen
naive_expanded$perspective <- rep(new_persp, times = nrow(naive_rows))

# ursprüngliche naive-Zeilen entfernen
SimStudy3_Results <- SimStudy3_Results[SimStudy3_Results$perspective != "naive", ]

# alles wieder zusammenführen
SimStudy3_Results <- rbind(SimStudy3_Results, naive_expanded)

supplied_rows <- SimStudy3_Results[SimStudy3_Results$estimator == "supplied" & SimStudy3_Results$perspective=="marginal_not_efficient", ]

# verdoppeln
supplied_expanded <- supplied_rows[rep(seq_len(nrow(supplied_rows)), each = length("marginal_not_efficient_HC3")), ]

# perspective ersetzen
supplied_expanded$perspective <- rep("marginal_not_efficient_HC3", times = nrow(supplied_rows))

SimStudy3_Results <- rbind(SimStudy3_Results, supplied_expanded )


SimStudy3_Results$coefficient<-
  ifelse(SimStudy3_Results$coefficient == "X3", "L1-noise","L1-noise (random slope)")


## QQ-plots --------------------------------------------------------------------

library(ggplot2)
library(dplyr)


# Hilfsspalte für lesbare Beschriftung
results_all_conditions_df <- SimStudy3_Results  %>%
  mutate(
    label_facet = paste0(
      #"R² = ", R_squared_pseudo,
      "Selected Model = ", working_model
      
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

results_all_conditions_df$pval_type<-relevel(as.factor(results_all_conditions_df$pval_type),"Post-Selective")


warning("oder für alle Estiamtors und alle selected models?")
# results_all_conditions_df <- 
#   results_all_conditions_df[!(results_all_conditions_df$estimator %in% c("ICM Covariance Matrix","VarY Covariance Matrix")),]
# results_all_conditions_df <- 
#   results_all_conditions_df[(results_all_conditions_df$working_model %in% c("X1X2X3X4","fullCandidate")),]


results_all_conditions_df_w.o.fullC <- 
  results_all_conditions_df[results_all_conditions_df$perspective!="fullCandidate",]

results_all_conditions_df_fullC <- 
  results_all_conditions_df[results_all_conditions_df$perspective=="fullCandidate",]

# QQ-Plot gegen Uniformverteilung
ggplot(
  results_all_conditions_df,
  aes(
    sample   = pval_two,
    color    = estimator,
    linetype    = pval_type
  )
) +
  stat_qq(distribution = stats::qunif, geom = "line") +   # Linien statt Punkte
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey") + 
  facet_grid(perspective ~ label_facet + coefficient, 
             labeller = label_wrap_gen(width = 35)) +
  labs(
    title = paste0("Simulation Study 3 - normal standard errors - Rsquared = ", R_squared),
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

# -----------------------------
# Manual color palette
# -----------------------------
estimator_colors <- c(
  "fullCandidate_supplied" ="black",
  "supplied"            = "#009E73",  # Grün
  "ModelEstimate"       =  "#8C510A",  # braun
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


results_all_conditions_df_SimStudy3 <- 
  results_all_conditions_df_w.o.fullC[results_all_conditions_df_w.o.fullC$perspective!="Conditional",]

# für Konferenz
# results_all_conditions_df_SimStudy3 <- 
#   results_all_conditions_df_w.o.fullC[results_all_conditions_df_w.o.fullC$perspective %in% c("Marginal (not eff) HC3","Marginal (not eff)") &
#                                       results_all_conditions_df_w.o.fullC$working_model=="X1X2X3X4"&
#                                       results_all_conditions_df_w.o.fullC$estimator %in% c("ModelEstimate","supplied","ModelEstimate_SIGMA","ModelEstimate_R"),]
# 

ggplot(
  results_all_conditions_df_SimStudy3,
  aes(
    sample   = pval_two,
    color    = estimator,
    linetype    = pval_type
  )
) +
  stat_qq(distribution = stats::qunif, geom = "line") +   # Linien statt Punkte
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey") + 
  facet_grid(perspective ~ label_facet + coefficient, 
             labeller = label_wrap_gen(width = 35)) +
  scale_color_manual(values = estimator_colors) +
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
  ) +
  #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14)
  )


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
    labeller = label_wrap_gen(width = 35)
  ) +
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
plot_SimStudy3

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Figures/Figures project 1")
ggsave("plot_SimStudy3.png", 
       plot = plot_SimStudy3, 
       dpi = 900, width = 9, height = 10, units = "in")


results_all_conditions_df_SimStudy3_Appendix <- 
  results_all_conditions_df_w.o.fullC[results_all_conditions_df_w.o.fullC$perspective %in% c("Conditional"),]

ggplot(
  results_all_conditions_df_SimStudy3_Appendix,
  aes(
    sample   = pval_two,
    color    = estimator,
    linetype    = pval_type
  )
) +
  stat_qq(distribution = stats::qunif, geom = "line") +   # Linien statt Punkte
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey") + 
  facet_grid(perspective ~ label_facet + coefficient, 
             labeller = label_wrap_gen(width = 35)) +
  scale_color_manual(values = estimator_colors) +
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
  ) +
  #guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9)
  )
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
    labeller = label_wrap_gen(width = 35)
  ) +
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

setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /Figures/Figures project 1")
ggsave("plot_SimStudy3_Appendix.png", 
       plot = plot_SimStudy3_Appendix, 
       dpi = 900, width = 12.2, height = 6, units = "in")




### convergence issues

# initial model fitting

convergenceinfoX123 <-
  bind_rows(SimulationStudy3_notRobust_X1X2X3_2$convergence_singular_info,
            SimulationStudy3_notRobust_X1X2X3$convergence_singular_info,
            SimulationStudy3_not_efficient_RdiagT_X1X2X3_2$convergence_singular_info,
            SimulationStudy3_not_efficient_RdiagT_X1X2X3$convergence_singular_info,
            HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3$convergence_singular_info,
            )

convergenceinfoX1234 <-
  bind_rows(SimulationStudy3_notRobust_X1X2X3X4_2$convergence_singular_info,
            SimulationStudy3_notRobust_X1X2X3X4$convergence_singular_info,
            SimulationStudy3_not_efficient_RdiagT_X1X2X3X4_2$convergence_singular_info,
            SimulationStudy3_not_efficient_RdiagT_X1X2X3X4$convergence_singular_info,
            HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3X4$convergence_singular_info,
            HC3SimulationStudy3_not_efficient_RdiagF_X1X2X3X4_2$convergence_singular_info
            
  )          

convergenceinfoX124 <-
  bind_rows(HC3SimulationStudy3_not_efficient_RdiagF_X1X2X4$convergence_singular_info,
            SimulationStudy3_not_efficient_RdiagT_X1X2X4$convergence_singular_info,
            SimulationStudy3_not_efficient_RdiagT_X1X2X4_2$convergence_singular_info,
            SimulationStudy3_notRobust_X1X2X4$convergence_singular_info,
            SimulationStudy3_notRobust_X1X2X4_2$convergence_singular_info,
  )

convergenceinfo_fullCandidate <- SimulationStudy3_fullCandidate_info

colMeans(convergenceinfoX123[sapply(convergenceinfoX123, is.numeric)], na.rm = F)
colMeans(convergenceinfoX1234[sapply(convergenceinfoX1234, is.numeric)], na.rm = F)
colMeans(convergenceinfoX124[sapply(convergenceinfoX124, is.numeric)], na.rm = F)
colMeans(convergenceinfo_fullCandidate[sapply(convergenceinfo_fullCandidate, is.numeric)], na.rm = F)

0.0003333333/(0.0003333333+199.5936666667) 
0.0006666667/ (0.0006666667+370.2566666667)
0.0003333333/( 0.0003333333+2.8810000000 )
0/(0+1)

0.0000000000/70.1176666667 
0.0000000000/130.7403333333
0.0000000000/1.0226666667 

## monte carlo resampling process

library(dplyr)
library(tidyr)
library(ggplot2)

summary_stats_df <- results_all_conditions_df_w.o.fullC %>%
  group_by(
    estimator,
    pval_type,
    perspective,
    label_facet,
    coefficient,
    R_squared_pseudo,
    n_obs,
    m
  ) %>%
  summarise(
    q75  = quantile(total_resamplings_due_to_convergence, 0.75, na.rm = TRUE),
    q95  = quantile(total_resamplings_due_to_convergence, 0.95, na.rm = TRUE),
    q99 = quantile(total_resamplings_due_to_convergence, 0.99, na.rm = TRUE),
    mean = mean(total_resamplings_due_to_convergence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c( q75, q95,q99,  mean),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  mutate(
    statistic = factor(
      statistic,
      levels = c( "q75", "q95", "q99", "mean"),
      labels = c( "q.75", "q.95", "q.99",  "Mean")
    )
  )

plot_df <- summary_stats_df %>%
  filter(perspective != "Conditional" & pval_type!="Naive")

plot_MCResampling_SimStudy3 <- ggplot(
  plot_df,
  aes(
    x     = statistic,
    y     = value,
    color = estimator,
    shape = pval_type,
    group = interaction(estimator, pval_type)
  )
) +
  # Linien: nur 25%, 50%, 75% innerhalb des gefilterten Datensatzes
  geom_line(
    data = plot_df %>% filter(statistic != "Mean"),
    linewidth = 0.8,
    alpha = 0.8
  ) +
  # Punkte: alle Punkte inklusive Mean
  geom_point(
    size = 3,
    position = position_dodge(width = 0.2)
  ) +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = label_wrap_gen(width = 35)
  ) +
  labs(
    x = "",
    y = "MC Resamplings due to convergence",
    color = "Estimator",
    shape = "p-value type"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(color = "black", size = 16),
    legend.position = "bottom",
    legend.box = "vertical"
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


plot_df <- summary_stats_df %>%
  filter(perspective == "Conditional" )

plot_MCResampling_SimStudy3_Appendix <- ggplot(
  plot_df,
  aes(
    x     = statistic,
    y     = value,
    color = estimator,
    group = interaction(estimator, pval_type)
  )
) +
  # Linien: nur 25%, 50%, 75% innerhalb des gefilterten Datensatzes
  geom_line(
    data = plot_df %>% filter(statistic != "Mean"),
    linewidth = 0.8,
    alpha = 0.8
  ) +
  # Punkte: alle Punkte inklusive Mean
  geom_point(
    size = 3,
    position = position_dodge(width = 0.2)
  ) +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = label_wrap_gen(width = 35)
  ) +
  labs(
    x = "",
    y = "MC Resamplings due to convergence",
    color = "Estimator",
    shape = "p-value type"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(color = "black", size = 16),
    legend.position = "bottom",
    legend.box = "vertical"
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


