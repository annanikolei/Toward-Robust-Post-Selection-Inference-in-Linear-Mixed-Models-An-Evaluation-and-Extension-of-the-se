################################################################################
# Simulation Study 1: single simulation run that generates data, analyzes the 
#                     data, and obtains parameters of interest
################################################################################


################### one_simulation function ####################################

SimStudy1_one_simulation <- function(n_obs, # number of observations
                                     ICC_pseudo, # pseudo intraclass-correlation
                                     perspective, 
                                     seed = NULL,
                                     checkConvergence_for_check_Fun # T or F (should convergence checks for MC sampling be conducted)
){
  
  
  ## Generate data by calling the dgp function -----------------------------------
  
  set.seed(seed)
  dgp_output <- SimStudy1_dgp(n_obs = n_obs, ICC_pseudo = ICC_pseudo)
  dat <- data.frame(y = dgp_output$y,dgp_output$X, l2.id = dgp_output$l2.id)
  
  
  ## Fit/Select full model (without selection) -----------------------------------
  
  mod <- lmer(y ~ 1 + X1.l1 + X2.l1 + X.l2 + (1 + X1.l1 | l2.id),
              data = dat,
              REML = FALSE)
  naive_convergence_okay <- naive_convergence_okay <- tryCatch(
    {
      performance::check_convergence(mod)
    },
    error = function(e) {
      FALSE
    }
  )
  
  singular_fit <- performance::check_singularity(mod)
  
  ## naive selective inference ---------------------------------------------------
  naive_pvals <- summary(mod)$coefficients[2:4,"Pr(>|t|)"]
  
  
  ## valid selective inference via numerical approximation -----------------------
  
  # since no selection occured, selection event always matches initial selection 
  # event
  checkFun.no_sel <- function(yb){ TRUE }
  
  # checkFun_checkConv <- function(yb){ 
  #   
  #   newdat <- data.frame(y = yb, dgp_output$X, l2.id = dgp_output$l2.id)
  #   new_full_model <- lmer(y ~ 1 + X1.l1 + X2.l1 + X.l2 + (1 + X1.l1 | l2.id),
  #                          data = newdat,
  #                          REML = FALSE)
  #   full_convergence_okay = check_convergence(new_full_model)
  #   
  #   if(full_convergence_okay){
  #     in_inference_region <- TRUE
  #   }
  #   
  #   if(!full_convergence_okay){
  #     in_inference_region <- NA
  #   }
  #   
  #   
  #   return( list(in_inference_region = in_inference_region,
  #                convergence_okay = full_convergence_okay)  )
  #   
  # }
  
  checkFun_checkConv <- function(yb){ 
    
    newdat <- data.frame(y = yb, dgp_output$X, l2.id = dgp_output$l2.id)
    
    new_full_model <- tryCatch(
      lmer(y ~ 1 + X1.l1 + X2.l1 + X.l2 + (1 + X1.l1 | l2.id),
           data = newdat,
           REML = FALSE),
      error = function(e) NULL
    )
    
    # If model failed to fit → not converged
    if (is.null(new_full_model)) {
      return(list(
        in_inference_region = NA,
        convergence_okay = FALSE
      ))
    }
    
    # Now safely check convergence
    full_convergence_okay <- tryCatch(
      check_convergence(new_full_model),
      error = function(e) FALSE
    )
    
    in_inference_region <- if (isTRUE(full_convergence_okay)) TRUE else NA
    
    return(list(
      in_inference_region = in_inference_region,
      convergence_okay = isTRUE(full_convergence_okay)
    ))
  }
  
  
  
  ## --- Default: set everything to NA -----------------------------------------
  
  selfmade_marginal_pvals <- NA
  selfmade_marginal_not_efficient_pvals <- NA
  selfmade_conditional_bayes_pvals <- NA
  selfmade_conditional_not_bayes_pvals <- NA
  
  selfmade_marginal_total_resamplings_due_to_convergence <- NA
  selfmade_marginal_not_efficient_total_resamplings_due_to_convergence <- NA
  selfmade_conditional_bayes_total_resamplings_due_to_convergence <- NA
  selfmade_conditional_not_bayes_total_resamplings_due_to_convergence <- NA
  
  
  if(naive_convergence_okay){
    
    selfmade_output <-
      mocasin(mod = mod, 
              checkFun = checkFun.no_sel,
              checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
              checkFun_checkConv = checkFun_checkConv,
              this_y = dat$y, 
              B_MCsamples = 1000, # number of MC-samples
              MC_distribution = "proposal",
              bayesian = c("not_bayesian", "bayesian"),
              which_marginal = 2:4, # selinf for X1.l1, X2.l1, X.l2
              which_conditional = 2:4, # selinf for X1.l1, X2.l1, X.l2
              rho_H0 = 0, 
              perspective = perspective, 
              robust_std.errors_marginal = F,
              marginal_COV_estimator = c("supplied"), 
              marginal_not_efficient_COV_estimator = c("supplied"),
              conditional_COV_estimator = c("supplied"), 
              diagF_not_efficient_SIGMA_or_R = NULL,
              symmetric_selection = T, 
              alpha = 0.05, 
              allow_smaller_alpha = T, 
              multiple = 7, 
              R_diag = TRUE, 
              Sigma_supplied = NULL,
              Z_supplied = dgp_output$Z,
              G_supplied = dgp_output$G,
              Lambda_supplied = NULL,
              Lambdat_supplied = NULL,
              G_work_supplied = dgp_output$G, 
              sigma_square_epsilon_supplied = dgp_output$sigma_square_epsilon)  
    
      
      ## ---- marginal ----------------------------------------------------------
      if ("marginal" %in% perspective) {
        
        selfmade_marginal_pvals <- 
          sapply(selfmade_output$marginal[[1]],
                 function(x) x[[1]]$pval_two)
        
        selfmade_marginal_total_resamplings_due_to_convergence <- 
          sapply(selfmade_output$marginal[[1]],
                 function(x) x[[1]]$total_resamplings_due_to_convergence)
        
        names(selfmade_marginal_pvals) <- names(naive_pvals)
        names(selfmade_marginal_total_resamplings_due_to_convergence) <- names(naive_pvals)
      }
      
      
      ## ---- marginal_not_efficient -------------------------------------------
      if ("marginal_not_efficient" %in% perspective) {
        
        selfmade_marginal_not_efficient_pvals <- 
          sapply(selfmade_output$marginal_not_efficient[[1]],
                 function(x) x[[1]]$pval_two)
        
        selfmade_marginal_not_efficient_total_resamplings_due_to_convergence <- 
          sapply(selfmade_output$marginal_not_efficient[[1]],
                 function(x) x[[1]]$total_resamplings_due_to_convergence)
        
        names(selfmade_marginal_not_efficient_pvals) <- names(naive_pvals)
        names(selfmade_marginal_not_efficient_total_resamplings_due_to_convergence) <- names(naive_pvals)
      }
      
      
      ## ---- conditional -------------------------------------------------------
      if ("conditional" %in% perspective) {
        
        selfmade_conditional_bayes_pvals <- 
          sapply(selfmade_output$conditional[[1]][[2]],
                 function(x) x[[1]]$pval_two)
        
        selfmade_conditional_not_bayes_pvals <- 
          sapply(selfmade_output$conditional[[1]][[1]],
                 function(x) x[[1]]$pval_two)
        
        selfmade_conditional_bayes_total_resamplings_due_to_convergence <- 
          sapply(selfmade_output$conditional[[1]][[2]],
                 function(x) x[[1]]$total_resamplings_due_to_convergence)
        
        selfmade_conditional_not_bayes_total_resamplings_due_to_convergence <- 
          sapply(selfmade_output$conditional[[1]][[1]],
                 function(x) x[[1]]$total_resamplings_due_to_convergence)
        
        names(selfmade_conditional_bayes_pvals) <- names(naive_pvals)
        names(selfmade_conditional_not_bayes_pvals) <- names(naive_pvals)
        
        names(selfmade_conditional_bayes_total_resamplings_due_to_convergence) <- names(naive_pvals)
        names(selfmade_conditional_not_bayes_total_resamplings_due_to_convergence) <- names(naive_pvals)
      }
    
    
    }
    
    
    ## --- Output ----------------------------------------------------------------
    
    return(list(
      convergence_okay = naive_convergence_okay,
      singular_fit = singular_fit,
      
      pvals_and_MCresamplings = cbind(
        naive_pvals,
        
        selfmade_marginal_pvals,
        selfmade_marginal_not_efficient_pvals,
        selfmade_conditional_bayes_pvals,
        selfmade_conditional_not_bayes_pvals,
        
        selfmade_marginal_total_resamplings_due_to_convergence,
        selfmade_marginal_not_efficient_total_resamplings_due_to_convergence,
        selfmade_conditional_bayes_total_resamplings_due_to_convergence,
        selfmade_conditional_not_bayes_total_resamplings_due_to_convergence)
    ))
    
  
  
  
}

