################################################################################
# Simulation Study 3: Simulation Study 
################################################################################
library(dplyr)
library(purrr)
library(stringr)


######################
## Simulation Study ##
######################

SimStudy3_simulation_study <- 
  function(n_iter, #number of replications per condition
           n_obs = 40^2, # number of observations
           m = 40, # number of level-2 units
           dgp_X_output, # SimStudy3_dgp_X( n_obs = n_obs, m = m)
           R_squared = 0.1, # pseudo coefficient of determination
           D = matrix(c(2,0.2,0.1,
                        0.2,1.5,0.2,
                        0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
           perspective = c("marginal","conditional"), # marginal and/or conditional (for robust std. errors we only need marginal)
           bayesian= c("bayesian"),
           robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
           R_diag =T,
           HC_type = "HC0", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
           cadjust = F, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
           B_MCsamples = 1000, # number of MC samples in selective inference process
           remove_singular_fits = F,
           list_supplied_D_work, # list that includes supplied working random effects covariance matrices for all interesting working models
           interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
           interesting_working_model_vars = NULL, # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
           checkConvergence_for_check_Fun # T or F (should convergence checks for MC sampling be conducted)
           
  ){
    
    
    ## simulation conditions ---------------------------------------------------
    
    # how many different (fixed) X have to be generated
    design_X <- expand.grid(n_obs = n_obs, m = m)
    
    # additional simulation conditions
    design <- expand.grid(i_iter = 1:n_iter,
                          R_squared = R_squared)
    
    
    ## Simulation --------------------------------------------------------------
    
    results_all_conditions <-
      sapply(1:nrow(design_X), function(idesignX){
        
        n_obs <- design_X$n_obs[idesignX]
        m <- design_X$m[idesignX]
        
        if((n_obs/m) %% 2 != 0) warning("n_obs/m/2 needs to be whole number")
        
        X <- dgp_X_output$X
        X.quadratic <- dgp_X_output$X.quadratic
        Z <- dgp_X_output$Z
        l2.id <- dgp_X_output$l2.id
        
        ## start simulation over all conditions with same m and n_obs -------
        
        results <- sapply(1:nrow(design), function(i_design) {
          
          # Set the maximum number of retries in case of an error
          max_retries <- 10000
          attempt <- 1
          result <- NULL
          
          while (attempt <= max_retries && is.null(result)) {
            # Try to execute the function and catch any errors
            result <- tryCatch({
              
              
              SimStudy3_one_simulation(l2.id = l2.id, 
                                       n_obs = n_obs, 
                                       m = m, 
                                       X = X, 
                                       X.quadratic = X.quadratic, 
                                       Z = Z, 
                                       R_squared = design$R_squared[i_design], 
                                       D = D, 
                                       list_supplied_D_work =  list_supplied_D_work,
                                       perspective = perspective, 
                                       bayesian= bayesian,
                                       R_diag = R_diag,
                                       HC_type = HC_type, # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                       cadjust = cadjust, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                       robust_std.errors_marginal = robust_std.errors_marginal,
                                       B_MCsamples = B_MCsamples,
                                       remove_singular_fits = remove_singular_fits,
                                       interesting_vars = interesting_vars,
                                       interesting_working_model_vars = interesting_working_model_vars,
                                       checkConvergence_for_check_Fun = checkConvergence_for_check_Fun) 
              
            }, error = function(e) {
              # Error handling: print a message and increment the attempt counter
              message(paste("Error at i_design =", i_design, "on attempt", attempt, ":", e$message))
              attempt <- attempt + 1
              return(NULL)  # Return NULL to indicate failure
              
            })
            
            if((i_design-1) %% 100 == 0){print(result);print(i_design)}
            
            # If the result is still NULL after the retries, proceed with the next iteration
            if (is.null(result)) {
              message(paste("Max retries reached for i_design =", i_design))
            }
          }
          
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
        
        convergence_singular_info_df <- imap_dfr(results, function(res, r2) {
          res$convergence_singular_info %>%
            mutate(R_squared_pseudo = r2)
        })
        
        estimates_df <- imap_dfr(results, function(res, r2) {
          res$estimates %>%
            mutate(R_squared_pseudo = r2)
        })
        
        result_list <- list(post_selection_pvalues = post_selection_pvalues_df,
                            convergence_singular_info = convergence_singular_info_df,
                            estimates_df = estimates_df) 
        
        return(result_list)
        
      }, simplify = FALSE)
    
    
    ## prepare output
    names(results_all_conditions) <- 
      paste0("n_obs=", design_X$n_obs, "_m=", design_X$m)
    
    post_selection_pvalues_all_conditions_df <- 
      imap_dfr(results_all_conditions, function(res, cond_name) {
        condition_parts <- str_split(cond_name, "_m=", simplify = TRUE)
        n_obs_value <- str_remove(condition_parts[1], "n_obs=")
        m_value <- condition_parts[2]
        
        res$post_selection_pvalues %>%
          mutate(
            n_obs = as.numeric(n_obs_value),
            m = as.numeric(m_value),
            condition = cond_name
          )
      })
    
    convergence_singular_info_all_conditions_df <- 
      imap_dfr(results_all_conditions, function(res, cond_name) {
        condition_parts <- str_split(cond_name, "_m=", simplify = TRUE)
        n_obs_value <- str_remove(condition_parts[1], "n_obs=")
        m_value <- condition_parts[2]
        
        res$convergence_singular_info %>%
          mutate(
            n_obs = as.numeric(n_obs_value),
            m = as.numeric(m_value),
            condition = cond_name
          )
      })
    
    estimates_all_conditions_df <- 
      imap_dfr(results_all_conditions, function(res, cond_name) {
        condition_parts <- str_split(cond_name, "_m=", simplify = TRUE)
        n_obs_value <- str_remove(condition_parts[1], "n_obs=")
        m_value <- condition_parts[2]
        
        res$estimates_df %>%
          mutate(
            n_obs = as.numeric(n_obs_value),
            m = as.numeric(m_value),
            condition = cond_name
          )
      })
    
    
    post_selection_pvalues_all_conditions_df <- 
      post_selection_pvalues_all_conditions_df[,c("perspective",
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
      convergence_singular_info_all_conditions_df[,c("convergence_failures_full_model",
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
    return( list(
      post_selection_pvalues = post_selection_pvalues_all_conditions_df,
      convergence_singular_info = convergence_singular_info_all_conditions_df,
      estimates_df = estimates_all_conditions_df
    ) )
    
  }



