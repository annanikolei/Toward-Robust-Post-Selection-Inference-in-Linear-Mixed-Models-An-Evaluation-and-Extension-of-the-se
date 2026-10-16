################################################################################
# Simulation Study 3: Simulation of supplied working covariance matrix 
#                     for conditional perspective 
################################################################################
# Start this script from the repository root.

setwd("simulations")

## Install and load all necessary packages
source("Project1_install_packages.R")

## packages and functions
source("SimStudy3_helper.R")
source("SimStudy3_dgp.R")
source("SimStudy3_get_D_work_supplied.R")

## helper functions ------------------------------------------------------------

generate_variable_combinations <- function(variables) {
  combinations <- list()
  
  for (i in 1:length(variables)) {
    combo_list <- combn(variables, i, simplify = FALSE)
    for (combo in combo_list) {
      name <- paste0(combo, collapse = "")
      combinations[[name]] <- combo
    }
  }
  
  return(combinations)
}


## Simulation function ---------------------------------------------------------

## full_variable_set: vector with names of variables included in full candidate 
##                     model.
#
## interesting_vars: vector of (fixed effects) variables of which at least one 
#                    variable needs to be included in the working model.
#
## R_squared: vector of R_squared values 
#
## number_samples: number of samples for simulation


Simulation_function_supplied_D <- function(full_variable_set = c("X1","X2","X3","X4"),
                                           interesting_vars = c("X3","X4"),
                                           random_effects_str_working_model = " + (1 | l2.id)",
                                           sample_size = 400,
                                           R_squared,
                                           number_samples,
                                           D = matrix(c(2,0.2,0.1,
                                                        0.2,1.5,0.2,
                                                        0.1,0.2,1.5),nrow=3, byrow=T)
){
  
  
  ## Create list of variables in working models we want to investigate
  
  all_combinations_of_variables_in_working_models <- 
    generate_variable_combinations(full_variable_set)
  
  combination_of_variables_in_working_models <-
    Filter(function(x) any(interesting_vars %in% x), 
           all_combinations_of_variables_in_working_models)
  
  
  ## Create fixed predictor matrix for given (huge) sample size
  
  dgp_X_output <- SimStudy3_dgp_X( n_obs = sample_size^2, m = sample_size )
  
  X <- dgp_X_output$X
  X.quadratic <- dgp_X_output$X.quadratic
  Z <- dgp_X_output$Z
  l2.id <- dgp_X_output$l2.id
  
  
  ## Simulate working covariance matrix D by averaging over the estimated 
  ## working random effects covariance of the samples
  
  list_mean_D_work_estimated <- suppressMessages(
    
    setNames(lapply(R_squared, function(iR_squared){
      
      setNames( lapply(combination_of_variables_in_working_models, function(iworking_vars){
        
        D_work_supplied_output <-
          get_D_work_supplied(
            variables_in_working_model = iworking_vars,
            random_effects_str_working_model = random_effects_str_working_model,
            n_obs = sample_size^2,
            l2.id = l2.id,
            m = sample_size,
            X = X,
            X.quadratic = X.quadratic,
            D = D,
            Z = Z,
            R_squared = iR_squared,
            number_samples = number_samples,
            remove_singular_fits = T
          )
        
      }), names(combination_of_variables_in_working_models)  )
      
    }), paste0("R_squared = ", R_squared) ) )
  
  
  return(list_mean_D_work_estimated)
  
  
}


## Simulation Results------------------------------------------------------------------

plan(multisession) 

Simulation_supplied_D <- 
  Simulation_function_supplied_D(full_variable_set = c("X1","X2","X3","X4"),
                                 interesting_vars = c("X3","X4"),
                                 random_effects_str_working_model = " + (1 | l2.id)",
                                 sample_size = 400,
                                 R_squared = 0.3,
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T),
                                 number_samples = 10000)



list_supplied_D_work <- lapply(Simulation_supplied_D, function(model_list) {
  lapply(model_list, function(x) x$mean_D_work_estimated)
})
list_convergence <- lapply(Simulation_supplied_D[["R_squared = 0.3"]], function(x) x$counter_convergence_issues)
list_singular<- lapply(Simulation_supplied_D[["R_squared = 0.3"]], function(x) x$counter_singular_issues)

results_dir <- "SimulationResults_supplied_D"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
setwd(results_dir)

save(list_supplied_D_work, 
     file = "SIM3_list_supplied_D_work_full=X1X2X3X4_interesting=X3X4_n=400_Rsquared=0.3_nrsamples=10000.RData")

save(list_convergence, 
     file = "SIM3_convergenceIssues_full=X1X2X3X4_interesting=X3X4_n=400_Rsquared=0.3_nrsamples=10000.RData")

save(list_singular, 
     file = "SIM3_singularIssues_full=X1X2X3X4_interesting=X3X4_n=400_Rsquared=0.3_nrsamples=10000.RData")


