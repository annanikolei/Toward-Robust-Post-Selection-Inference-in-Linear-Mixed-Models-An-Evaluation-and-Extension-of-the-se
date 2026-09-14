#############################################################################################
# simulation study Discussion - Naive Post-selection with supplied - Simulation Study 4
#############################################################################################
library(dplyr)
library(purrr)
library(stringr)


######################
## Simulation Study ##
######################


SimStudy4_Naive_supplied_simulation_study <- 
  function(n_iter, #number of replications per condition
           l2.id, # id variable for level-2 units
           n_obs = 40^2, # number of observations
           m = 40, # number of level-2 units
           X, # X-matrix 
           X.quadratic, # X-matrix (including squared variables)
           Z, # design matrix corresponding to X
           R_squared = 0.3, # pseudo coefficient of determination
           D = matrix(c(2,0.2,0.1,
                        0.2,1.5,0.2,
                        0.1,0.2,1.5),
                      nrow = 3, byrow=T) , # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
           perspective = c("marginal_not_efficient"),
           checkConvergence_for_check_Fun = TRUE,
           remove_singular_fits = F, #remove singular fit 
           alpha = 0.000001, # alpha for significance selection
           interesting_working_model_ranef = c("(Intercept)","X2"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
           interesting_fixed_effects = c("X3","X4")
  ){
    
    
    ## simulation conditions ---------------------------------------------------
    
    # additional simulation conditions
    design <- expand.grid(i_iter = 1:n_iter,
                          R_squared = R_squared)
    
    
    ## Simulation --------------------------------------------------------------
    
    results <- sapply(1:nrow(design), function(i_design) {
      
      # Set the maximum number of retries in case of an error
      max_retries <- 100000000
      attempt <- 1
      result <- NULL
      
      while (attempt <= max_retries && is.null(result)) {
        # Try to execute the function and catch any errors
        result <- tryCatch({
          
          
          Naive_supplied_one_simulation(l2.id = l2.id, 
                                       n_obs = n_obs, 
                                       m = m, 
                                       X = X, 
                                       X.quadratic = X.quadratic, 
                                       Z = Z, 
                                       perspective=perspective,
                                       R_squared = design$R_squared[i_design], 
                                       D = D, 
                                       alpha = alpha,
                                       checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
                                       remove_singular_fits = remove_singular_fits,
                                       interesting_working_model_ranef = interesting_working_model_ranef,
                                       interesting_fixed_effects = interesting_fixed_effects 
          ) 
          
        }, error = function(e) {
          # Error handling: print a message and increment the attempt counter
          message(paste("Error at i_design =", i_design, "on attempt", attempt, ":", e$message))
          attempt <- attempt + 1
          return(NULL)  # Return NULL to indicate failure
          
        })
        
        #if((i_design-1) %% 10 == 0){print(result);print(i_design)}
        
        # If the result is still NULL after the retries, proceed with the next iteration
        if (is.null(result)) {
          message(paste("Max retries reached for i_design =", i_design))
        }
        
      }
      
      return(result)
    }, simplify = FALSE)
    
    return(results)
    
  }
