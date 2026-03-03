################################################################################
# Simulation Study 1: For each level of design factors and their combination, 
#                     performs simulation runs
################################################################################

# install packages -------------------------------------------------------------
install.packages("lme4")
installed.packages("lmerTest")
install.packages("Matrix")
install.packages("mvtnorm")

# load functions and packages --------------------------------------------------
library(lme4)
library(lmerTest)
library(Matrix)
library(mvtnorm)




################### simulation study function ##################################

SimStudy1_simulation_study <- 
  function(n_iter, # number of replications per condition
           n_obs = c(990,4830), # number of observations
           ICC_pseudo = c(0.05, 0.25), # pseudo intraclass-correlation
           checkConvergence_for_check_Fun, # T or F (should convergence checks for MC sampling be conducted)
           perspective,
           ignore_singular_fits = TRUE, # FALSE if singular fit-samples should be removed
           seed = NULL
  ){
    
    ## simulation conditions -------------------------------------------------------
    design <- expand.grid(i_iter = 1:n_iter,
                          n_obs = n_obs,
                          ICC_pseudo = ICC_pseudo)
    
    ## start simulation ------------------------------------------------------------
    if (ignore_singular_fits) {
      
      results <- lapply(1:nrow(design), function(i_design) {
        
        convergence_okay <- FALSE
        n_nonconverged <- 0  # Counter for failed convergence attempts
        
        while (!convergence_okay) {
          result <- tryCatch(
            {
              SimStudy1_one_simulation(
                n_obs = design$n_obs[i_design],
                perspective = perspective,
                ICC_pseudo = design$ICC_pseudo[i_design],
                checkConvergence_for_check_Fun,
                seed = seed
              )
            },
            error = function(e) {
              message(sprintf(" Fehler bei Design-Zeile %s: %s", i_design, e$message))
              return(NULL)  
            }
          )
          
          
          if (is.null(result)) {
            return(list(
              pvals_and_MCresamplings = paste0("NA_seed=",seed),
              n_nonconverged = paste0("NA_seed=",seed)
            ))
          }
          
          if((i_design-1) %% 200 == 0){print(result);print(i_design)}
          
          convergence_okay <- result$convergence_okay
          if (!convergence_okay) {
            n_nonconverged <- n_nonconverged + 1
            seed = 4000 + seed + n_nonconverged
          }
        }
        
        return(list(
          pvals_and_MCresamplings = result$pvals_and_MCresamplings,
          n_nonconverged = n_nonconverged
        ))
      })
      
      retry_counts <- sapply(results, function(res) res$n_nonconverged)
      results_list <- lapply(results, function(res) res$pvals_and_MCresamplings)
      
    } else {
      
      results <- lapply(1:nrow(design), function(i_design) {
        
        convergence_okay <- FALSE
        singular_fit <- TRUE
        n_retries <- 0  
        
        while (!convergence_okay || singular_fit) {
          
          result <- tryCatch(
            {
              SimStudy1_one_simulation(
                n_obs = design$n_obs[i_design],
                perspective = perspective,
                ICC_pseudo = design$ICC_pseudo[i_design],
                checkConvergence_for_check_Fun
              )
            },
            error = function(e) {
              message(sprintf("Error in Design-Row %s: %s", i_design, e$message))
              return(NULL)
            }
          )
          
          if (is.null(result)) {
            return(list(
              pvals_and_MCresamplings = NA,
              n_retries = NA
            ))
          }
          
          convergence_okay <- result$convergence_okay
          singular_fit <- result$singular_fit
          
          if (!convergence_okay || singular_fit) {
            n_retries <- n_retries + 1
          }
        }
        
        list(
          pvals_and_MCresamplings = result$pvals_and_MCresamplings,
          n_retries = n_retries
        )
      })
      
      retry_counts <- sapply(results, function(res) res$n_retries)
      results_list <- lapply(results, function(res) res$pvals_and_MCresamplings)
    }
    
    ## output ----------------------------------------------------------------------
    return(list( design = design,
                 pvalues_and_MCresamplings = results_list,
                 retry_counts = retry_counts))
  }
