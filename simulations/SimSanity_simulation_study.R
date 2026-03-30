################################################################################
# Sanity Check: Simulation Study - SimStudy 2 - Discussion - Bayesian conditional
################################################################################
library(dplyr)
library(purrr)
library(stringr)


######################
## Simulation Study ##
######################

SimSanity_simulation_study <- 
  function(n_iter, #number of replications per condition
           n_obs = 40^2, # number of observations
           m = 40, # number of level-2 units
           R_squared = c(0.1,0.6), # pseudo coefficient of determination
           D,
           perspective = c("marginal","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
           robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
           B_MCsamples = 1000, # number of MC samples in selective inference process
           remove_singular_fits = F,
           interesting_vars = c("X2"), # for which variables do we want to do post-selection inference
           interesting_working_model_vars = c("X1","X2"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
           checkConvergence_for_check_Fun = T, # T or F (should convergence checks for MC sampling be conducted)
           R_diag = TRUE
  ){
    
    
    ## simulation conditions ---------------------------------------------------
    
    design <- expand.grid(i_iter = 1:n_iter,
                          R_squared = R_squared)
    
    
    ## Simulation --------------------------------------------------------------
    
        results <- sapply(1:nrow(design), function(i_design) {
          
          
          # Set the maximum number of retries in case of an error
          max_retries <- 10000
          attempt <- 1
          result <- NULL
          
          while (attempt <= max_retries && is.null(result)) {
            # Try to execute the function and catch any errors
            result <- tryCatch({
              
              
              SanCheck_one_simulation( n_obs = n_obs,
                                       m = m,
                                       R_squared = design$R_squared[i_design], 
                                       D = D,
                                       perspective = perspective, 
                                       robust_std.errors_marginal = robust_std.errors_marginal,
                                       B_MCsamples = B_MCsamples,
                                       remove_singular_fits = remove_singular_fits,
                                       interesting_vars = interesting_vars,
                                       interesting_working_model_vars = interesting_working_model_vars,
                                       checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
                                       R_diag = R_diag) 
              
            }, error = function(e) {
              # Error handling: print a message and increment the attempt counter
              message(paste("Error at i_design =", i_design, "on attempt", attempt, ":", e$message))
              attempt <- attempt + 1
              return(NULL)  # Return NULL to indicate failure
              
            })
            
            # If the result is still NULL after the retries, proceed with the next iteration
            if (is.null(result)) {
              message(paste("Max retries reached for i_design =", i_design))
            }
          }
          
          if((i_design-1) %% 50 == 0){print(result);print(i_design)}
          
          # If the function still failed after max retries, return NA
          if (is.null(result)) {
            result <- NA
          }
          
          return(result)
          
        }, simplify = FALSE)
        
        
        ## prepare output
        
        names(results) <- paste0("R_squared = ",design$R_squared)
        
        post_selection_pvalues_df <- imap_dfr(results, function(res, r2) {
          res$post_selection_pvalues %>%
            mutate(R_squared_pseudo = r2)
        })
        post_selection_pvalues_df$n_obs <- rep(paste0("n_obs=",n_obs),
                                               nrow(post_selection_pvalues_df))
        post_selection_pvalues_df$m <- rep(paste0("m=",m),
                                               nrow(post_selection_pvalues_df))
        
        convergence_singular_info_df <- imap_dfr(results, function(res, r2) {
          res$convergence_singular_info %>%
            mutate(R_squared_pseudo = r2)
        })
        convergence_singular_info_df$n_obs <- rep(paste0("n_obs=",n_obs),
                                               nrow(convergence_singular_info_df))
        convergence_singular_info_df$m <- rep(paste0("m=",m),
                                           nrow(convergence_singular_info_df))
        
        
        
    post_selection_pvalues_all_conditions_df <- 
                    post_selection_pvalues_df[, c("perspective",
                                                  "estimator",
                                                  "coefficient",
                                                  "pval_two",
                                                  "cil",
                                                  "ciu",
                                                  "cil.smalleralpha",
                                                  "ciu.smalleralpha",
                                                  "total_resamplings_due_to_convergence",
                                                  "R_squared_pseudo",
                                                  "n_obs",
                                                  "m",
                                                  "working_model")]
    
    convergence_singular_info_all_conditions_df <-
      convergence_singular_info_df[,c("convergence_failures_full_model",
                                                     "convergence_okay_full_model",
                                                     "convergence_failures_selected_model",
                                                     "tries_selected_model",
                                                     "full_model_singular",
                                                     "selected_model_singular",
                                                     "R_squared_pseudo",
                                                     "n_obs",
                                                     "m",
                                                     "working_model")]
    
    
    
    ## output
    return( 
      list( post_selection_pvalues = post_selection_pvalues_all_conditions_df,
      convergence_singular_info = convergence_singular_info_all_conditions_df)
      )
    
  }


