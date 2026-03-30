################################################################################
# Sanity Check Simulation - one simulation
################################################################################
library(performance)
library(clubSandwich)
setwd("~/Documents/Forschung/Post-selection inference/Masterarbeit/selfmade package nachprogrammieren/GitHub_selective_inference")
source("Hauptfunktion (analog zur mocasin Funktion).R")
source("parallel_selective-inference-Funktionen.R")
setwd("~/Documents/Forschung/Post-selection inference/Projekt 1 /SimSanityCheck_Problems_with_bayesian_conditional")

SanCheck_one_simulation <- function(n_obs, # number of observations
                                     m, # number of level-2 units
                                     R_squared = 0.3, # pseudo coefficient of determination
                                     D, # 2x2 covariance matrix of random effects (u0,u2) per level2 unit
                                     perspective, # marginal and/or conditional (for robust std. errors we only need marginal)
                                     B_MCsamples = 1000, # number of MC-samples
                                     remove_singular_fits = F, #remove singular fit 
                                     interesting_vars = c("X2"), # for which variables do we want to make post-selective inference
                                     interesting_working_model_vars = c("X1","X2"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                                    robust_std.errors_marginal = FALSE,
                                    checkConvergence_for_check_Fun = T,
                                    R_diag = TRUE){
  
  convergence_failures_full_model <- 0
  convergence_okay_full_model <- 0 
  
  convergence_failures_selected_model <- 0
  tries_selected_model <- 0
  
  convergence_okay = F
  singular_okay = F
  selected_var_names = NULL
  working_equal_to_selected = F
  
  while(!(convergence_okay == T && 
          singular_okay == T &&
          any(interesting_vars %in% selected_var_names) &&
          working_equal_to_selected == T) ){
    
    ## Generate data (y and u) by calling the dgp function ---------------------
    
    dgp_output <- SanCheck_dgp(n_obs = n_obs,m = m,R_squared = R_squared,D = D)
    
    dat <- data.frame(y = dgp_output$y, dgp_output$X, l2.id = dgp_output$l2.id)
    
    ## Selection process -------------------------------------------------------
    
    # fit fullest candidate model
    full_model <- modFun(dat)
    
    # check convergence of lmer with performance-package
    convergence_okay = check_convergence(full_model)
    if (!convergence_okay) {
      convergence_failures_full_model <- convergence_failures_full_model + 1
    }
    
    
    if(convergence_okay){
      convergence_okay_full_model <- convergence_okay_full_model + 1
      
      # save first full model that has no convergence issues for inference of
      # full Candidate model. You can not use the full model of the selected model
      # since the selected model includes at least one of the noise variables. Hence
      # the coefficient for this noise variable is unexpectedly high in the selected model
      # but also in the corresponding full model
      if(convergence_okay_full_model == 1){fullCandidate_model <- full_model}
      
      # check if singular fits occured
      if (remove_singular_fits == F){singular_okay = T}
      else {singular_okay = !check_singularity(full_model)}
      
      
      if(singular_okay){
        # select model via stepwise backward elimination
        selected_model <- SelFun(full_model)
      
       
        
        # check convergence and singularity of selected model
        convergence_okay = check_convergence(selected_model)
        if (remove_singular_fits == F){singular_okay = T }
        else { singular_okay = !check_singularity(selected_model)}
        
        
        # names of selected variables
        selected_var_names <- extractSelFun_fixed(selected_model)
        
        if(is.null(interesting_working_model_vars)){
          working_equal_to_selected = TRUE
        }
        else if(!is.null(interesting_working_model_vars)){
          working_equal_to_selected = setequal(interesting_working_model_vars,selected_var_names)
        }
        
        if(any(interesting_vars %in% selected_var_names)) {
          
          tries_selected_model <- tries_selected_model + 1
          if (!convergence_okay) convergence_failures_selected_model <- convergence_failures_selected_model + 1
          
        }
      }
    }
  } # end while-loop
  
  
  ### Inference ------------------------------------------------------------------
  
  sel.noise <- match.arg(interesting_vars, selected_var_names, several.ok = T) 
  which_coeffs <- which(selected_var_names %in% sel.noise)  + 1 # +1 to account for Intercept                                                                               #?# bisher nur noise
  
  ## inference for full candidate model ----------------------------------------
  
  # We only perform inference for the full candidate model on samples 
  # that do not exhibit convergence issues. This is because these are 
  # the same samples used for post-selection inference. If convergence 
  # problems introduce any bias, it should be detectable—since in that 
  # case, the p-values of the working noise from the full candidate model 
  # would not follow a uniform distribution. 
  
  summ_fullCandidate_model <- summary(fullCandidate_model)
  fullCandidate_est <- summ_fullCandidate_model$coefficients[interesting_vars, "Estimate"]
  
  # do not use robust standard errors
  if (robust_std.errors_marginal == FALSE) {
    
    fullCandidate_pvals <- summ_fullCandidate_model$coefficients[interesting_vars, "Pr(>|t|)"]
    fullCandidate_se <- summ_fullCandidate_model$coefficients[interesting_vars, "Std. Error"]
    
    # Critical t-value (we used the 1.965 value instead)
    t_crit <- 1.965
    
    # Confidence interval bounds
    lower_fullCandidate <- fullCandidate_est - t_crit * fullCandidate_se
    upper_fullCandidate <- fullCandidate_est + t_crit * fullCandidate_se
  }
  
  # use robust standard errors
  if (robust_std.errors_marginal != FALSE) {
    
    vcovCR <- vcovCR(fullCandidate_model, cluster = dgp_output$l2.id, type = robust_std.errors_marginal)
    
    robust_test <- coef_test(fullCandidate_model,
                             vcov = vcovCR,
                             test = "Satterthwaite")
    
    fullCandidate_pvals <- robust_test[interesting_vars, "p_Satt"]
    
    est <- robust_test$beta[interesting_vars]
    se <- robust_test$SE[interesting_vars]
    df <- robust_test$df[interesting_vars]
    
    # Critical t-value
    t_crit <- qt(1 - 0.05 / 2, df = df)
    
    # Confidence interval bounds
    lower_fullCandidate <- est - t_crit * se
    upper_fullCandidate <- est + t_crit * se
  }
  
  # Assign names to the results
  names(fullCandidate_pvals) <- names(fullCandidate_est) <- 
    names(lower_fullCandidate) <- names(upper_fullCandidate) <- interesting_vars
  
  
  
  
  
  ## Post-selective inference --------------------------------------------------
  
  ## naive post-selection inference ## -----------------------------------------
  summ_selected_model <- summary(selected_model)
  naive_est <- summ_selected_model$coefficients[sel.noise, "Estimate"]
  
  # do not use robust standard errors
  if (robust_std.errors_marginal == FALSE) {
    
    naive_pvals <- summ_selected_model$coefficients[sel.noise, "Pr(>|t|)"]
    naive_se <- summ_selected_model$coefficients[sel.noise, "Std. Error"]
    
    # Critical t-value (we used the 1.965 value instead)
    t_crit <- 1.965
    
    # Confidence interval bounds
    lower <- naive_est - t_crit * naive_se
    upper <- naive_est + t_crit * naive_se
  }
  
  # use robust standard errors
  if (robust_std.errors_marginal != FALSE) {
    
    vcovCR <- vcovCR(selected_model, cluster = dgp_output$l2.id, type = robust_std.errors_marginal)
    
    robust_test <- coef_test(selected_model,
                             vcov = vcovCR,
                             test = "Satterthwaite")
    
    naive_pvals <- robust_test[sel.noise, "p_Satt"]
    
    est <- robust_test$beta[sel.noise]
    se <- robust_test$SE[sel.noise]
    df <- robust_test$df[sel.noise]
    
    # Critical t-value
    t_crit <- qt(1 - 0.05 / 2, df = df)
    
    # Confidence interval bounds
    lower <- est - t_crit * se
    upper <- est + t_crit * se
  }
  
  # Assign names to the results
  names(naive_pvals) <- names(naive_est) <- names(lower) <- names(upper) <- sel.noise
  
  
  
  
  ## selfmade ## ---------------------------------------------------------------
  
  ## checkFun: Function to check match with initial selection event ------------
  checkFun <- function(yb){ 
    
    newdat <- data.frame(y = yb, dgp_output$X, l2.id = dgp_output$l2.id)
    new_full_model <- modFun(newdat)
    new_selected_model <- SelFun(new_full_model)
    
    return(
      setequal( extractSelFun_fixed(new_selected_model), 
                selected_var_names) )
    
    
  }
  
  checkFun_checkConv <- function(yb){ 
    
    newdat <- data.frame(y = yb, dgp_output$X, l2.id = dgp_output$l2.id)
    new_full_model <- modFun(newdat)
    convergence_full_okay = check_convergence(new_full_model)
    
    if(convergence_full_okay){
      new_selected_model <- SelFun(new_full_model)
      in_inference_region <- setequal( extractSelFun_fixed(new_selected_model), 
                                       selected_var_names)
      convergence_okay <- check_convergence(new_selected_model)
      
    }
    
    if(!convergence_full_okay){
      in_inference_region <- NA
      convergence_okay = FALSE
    }
    
    
    return( list(in_inference_region = in_inference_region,
                 convergence_okay = convergence_okay)  )
    
  }
  
  ## extract mean estimated working covariance matrix G for selected model -----
  working_vars <- paste0(sort(selected_var_names), collapse = "")
  mean_G_work_estimated <- dgp_output$G
  
  ## Numerical Approximation ---------------------------------------------------
  if(R_diag == TRUE){
  marginal_COV_estimator = c("supplied","ModelEstimate")
  conditional_COV_estimator = c("supplied","ModelEstimate")
  }
  
  if(R_diag == FALSE){
    marginal_COV_estimator = c("supplied")
    conditional_COV_estimator = c("supplied")}
  
  selfmade_output <-
    mocasin(mod = selected_model, 
            checkFun = checkFun,
            checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
            checkFun_checkConv = checkFun_checkConv,
            this_y = dat$y, 
            B_MCsamples = B_MCsamples, # number of MC-samples
            MC_distribution = "proposal",
            bayesian = c("bayesian"),
            which_marginal = which_coeffs, 
            which_conditional = which_coeffs, 
            rho_H0 = 0, 
            perspective = perspective, 
            robust_std.errors_marginal = robust_std.errors_marginal,
            marginal_COV_estimator = marginal_COV_estimator, 
            conditional_COV_estimator = conditional_COV_estimator, 
            symmetric_selection = T, 
            alpha = 0.05, 
            allow_smaller_alpha = T, 
            multiple = 7, 
            R_diag = R_diag, 
            Sigma_supplied = NULL,
            Z_supplied = dgp_output$Z,
            G_supplied = dgp_output$G,
            Lambda_supplied = NULL,
            Lambdat_supplied = NULL,
            G_work_supplied = mean_G_work_estimated, 
            sigma_square_epsilon_supplied = 
              as.numeric(dgp_output$sigma_square_epsilon),
            R_supplied = dgp_output$R) 
  
  
  ## Create ouput --------------------------------------------------------------
  
  ## selfmade output
  selfmade_output_df <- extract_pval_two(selfmade_output)
  coefficient_num <- as.integer(sub("coefficient ", "",
                                    selfmade_output_df$coefficient))
  selfmade_output_df$coefficient <- selected_var_names[coefficient_num - 1] # -1 to account for Intercept
  
  ## naive inference output
  naive_df <- data.frame(
    perspective = "naive",
    estimator = "naive",
    coefficient = names(naive_pvals),
    pval_two = as.numeric(naive_pvals),
    cil = as.numeric(lower),
    ciu = as.numeric(upper),
    cil.smalleralpha = as.numeric(lower),
    ciu.smalleralpha = as.numeric(upper),
    total_resamplings_due_to_convergence = NA, # due to MC sampling
    stringsAsFactors = FALSE
  )
  
  ## naive inference output
  fullCandidate_df <- data.frame(
    perspective = "fullCandidate",
    estimator = "fullCandidate",
    coefficient = names(fullCandidate_pvals),
    pval_two = as.numeric(fullCandidate_pvals),
    cil = as.numeric(lower_fullCandidate),
    ciu = as.numeric(upper_fullCandidate),
    cil.smalleralpha = as.numeric(lower_fullCandidate),
    ciu.smalleralpha = as.numeric(upper_fullCandidate),
    total_resamplings_due_to_convergence = NA, # due to MC sampling
    stringsAsFactors = FALSE
  )
  
  ## combined output 
  post_selection_pvalues_withoutfullCandidate <- rbind(selfmade_output_df, naive_df)
  post_selection_pvalues <- rbind(post_selection_pvalues_withoutfullCandidate,
                                  fullCandidate_df)
  post_selection_pvalues$working_model <- 
    cbind( c(rep(working_vars,nrow(post_selection_pvalues_withoutfullCandidate)),
             rep("fullCandidate",nrow(fullCandidate_df))) )
  
  
  
  return(list(post_selection_pvalues = post_selection_pvalues,
              convergence_singular_info =
                data.frame(
                  convergence_failures_full_model = convergence_failures_full_model,
                  convergence_okay_full_model = convergence_okay_full_model,
                  convergence_failures_selected_model = convergence_failures_selected_model,
                  tries_selected_model = tries_selected_model,
                  full_model_singular = check_singularity(full_model),
                  selected_model_singular = check_singularity(selected_model),
                  working_model = working_vars)))
  
}



