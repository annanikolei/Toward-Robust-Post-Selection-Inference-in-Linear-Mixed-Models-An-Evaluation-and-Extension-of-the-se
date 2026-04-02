################################################################################
# Simulation Study 2 - Function to obtain D_work_supplied
################################################################################

library(lme4)
library(performance)
library(future.apply)

## The conditional perspective does not require the random effects covariance
## matrix G (or D) to be known. Instead, it assumes that the working random 
## effects covariance matrix is known—that is, the covariance matrix of the 
## random effects in the specific working model, as if the model were fitted 
## to the entire population.
## Since it is impossible to fit LMMs with extremely large sample sizes, we 
## decided to simulate many samples (with fixed X) from the population model
## and fit the working model to obtain the random effects covariance matrix.
## G_work_supplied then is the mean random effects covariance matrix

get_D_work_supplied <- function(variables_in_working_model,
                                random_effects_str_working_model = " + (1 + X2 + X4 | l2.id)",
                                l2.id,
                                n_obs,
                                m,
                                X,
                                X.quadratic,
                                D,
                                Z,
                                R_squared,
                                number_samples,
                                remove_singular_fits = T) {
  
  list_results <- future_sapply(1:number_samples, function(isample) {
    
    tries <- 0
    convergence_failures <- 0
    singular_failures <- 0
    
    result <- tryCatch({
      convergence_okay <- FALSE
      singular_okay <- FALSE
      
      while (!(convergence_okay && singular_okay)) {
        dgp_output <- SimStudy2_dgp_y(l2.id = l2.id, 
                                      n_obs = n_obs,
                                      m = m, 
                                      X = X, 
                                      X.quadratic = X.quadratic, 
                                      Z = Z, 
                                      R_squared = R_squared, 
                                      D = D)
        
        dat <- data.frame(y = dgp_output$y, X, l2.id = l2.id)
        
        formula_working_model <- as.formula(
          paste0("y ~ ",
                 paste(variables_in_working_model, collapse = " + "),
                 random_effects_str_working_model)
        )
        
        fit_working_model <- lmer(formula_working_model, REML = FALSE, data = dat)
        
        convergence_okay <- check_convergence(fit_working_model)
        singular_okay <- if (!remove_singular_fits) TRUE else !check_singularity(fit_working_model)
        
        tries <- tries + 1
        if (!convergence_okay) {
          convergence_failures <- convergence_failures + 1
        }
        if (!singular_okay) {
          singular_failures <- singular_failures + 1 
        }
      }
      
      D_work_est <- as.matrix(VarCorr(fit_working_model)[["l2.id"]])
      return(list(D_work_est = D_work_est,
                  convergence_failures = convergence_failures,
                  singular_failures = singular_failures,
                  total_tries = tries))
      
    }, error = function(e) {
      message(paste("error for sample", isample, ":", e$message))
      return(NULL)
    })
    
    return(result)
    
  }, future.seed = TRUE, simplify = FALSE)
  
  # Filter NULL results
  valid_results <- Filter(Negate(is.null), list_results)
  
  # Extract D matrices
  D_matrices <- lapply(valid_results, function(res) res$D_work_est)
  mean_D_work_estimated <- Reduce("+", D_matrices) / length(D_matrices)
  
  # Count convergence/singular fit issues
  total_convergence_failures <- sum(sapply(valid_results, function(res) res$convergence_failures))
  total_singular_failures <- sum(sapply(valid_results, function(res) res$singular_failures))
  total_tries <- sum(sapply(valid_results, function(res) res$total_tries))
  
  return(list(
    mean_D_work_estimated = mean_D_work_estimated,
    counter_convergence_issues = total_convergence_failures,
    counter_singular_issues = total_singular_failures,
    counter_tries = total_tries
  ))
}


