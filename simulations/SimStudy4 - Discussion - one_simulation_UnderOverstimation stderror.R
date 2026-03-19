#####################################################################################################
# one simualtion function Discussion - Under-/Overestimation of Standard errors - Simulation Study 4
#####################################################################################################


SimStudy4_CR2_one_simulation <- function(l2.id, # id variable for level-2 units
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
                                         robust_std.errors = "CR2", # (which) robust std.err in marginal perspective
                                         remove_singular_fits = F, #remove singular fit 
                                         alpha = 0.000001, # alpha for significance selection
                                         interesting_working_model_ranef = c("(Intercept)","X2"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                                         interesting_fixed_effects = c("X3","X4")
){
  
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
          working_equal_to_selected == T) ){
    
    ## Generate data (y and u) by calling the dgp function ---------------------
    
    dgp_output <- SimStudy4_dgp_y(l2.id = l2.id, 
                                  n_obs = n_obs,
                                  m = m, 
                                  X = X, 
                                  X.quadratic = X.quadratic, 
                                  Z = Z, 
                                  R_squared = R_squared, 
                                  D = D )
    
    dat <- data.frame(y = dgp_output$y, X, l2.id = l2.id)
    
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
      # full Candidate model (you can not use the full model of the selected model)
      # since the selected model includes at least one of the noise variables. Hence
      # the coefficient for this noise variable is unexpectedly high in the selected model
      # but also in the corresponding full model
      
      if(convergence_okay_full_model == 1){no_selection_dat <- dat }
      
      # check if singular fits occured
      if (remove_singular_fits == F){singular_okay = T}
      else {singular_okay = !check_singularity(full_model)}
      
      
      if(singular_okay){
        # select model 
        selected_model <- SelFun(dat ,full_model,alpha = alpha)
        
        # check convergence and singularity of selected model
        convergence_okay = check_convergence(selected_model)
        
        if (remove_singular_fits == F){singular_okay = T }
        else { singular_okay = !check_singularity(selected_model)}
        
        
        # names of selected random effects
        selected_var_names <- extractSelFun_random(selected_model)
        #print(selected_var_names)
        ## check if selected model matches interesting working models
        if(is.null(interesting_working_model_ranef)){
          working_equal_to_selected = TRUE
        }
        else if(!is.null(interesting_working_model_ranef)){
          
          if(is.list(interesting_working_model_ranef)){ 
            working_equal_to_selected = any(sapply(interesting_working_model_ranef, function(x) setequal(x, selected_var_names)))
          }
          
          else if(!is.list(interesting_working_model_ranef)){ 
            working_equal_to_selected = setequal(interesting_working_model_ranef,selected_var_names)
          }
        }
        
        tries_selected_model <- tries_selected_model + 1
        if (!convergence_okay) convergence_failures_selected_model <- convergence_failures_selected_model + 1
        
      }
    }
  } # end while-loop
  
  
  
  
  # inference---------------------------------------------------------------------
  if(identical(interesting_fixed_effects, c("X3","X4"))){which_coeffs <- c(4,5)}
  if(identical(interesting_fixed_effects,c("X4"))){which_coeffs <- c(5)}
  
  ## naive post-selection inference ## -----------------------------------------
  summ_selected_model <- summary(selected_model)
  naive_est <- summ_selected_model$coefficients[ interesting_fixed_effects, "Estimate"]
  
  # do not use robust standard errors
  if (robust_std.errors == FALSE) {
    
    naive_pvals <- summ_selected_model$coefficients[ interesting_fixed_effects, "Pr(>|t|)"]
    naive_se <- summ_selected_model$coefficients[ interesting_fixed_effects, "Std. Error"]
    
    # Critical t-value (we used the 1.965 value instead)
    t_crit <- 1.965
    
    # Confidence interval bounds
    lower <- naive_est - t_crit * naive_se
    upper <- naive_est + t_crit * naive_se
  }
  
  # use robust standard errors
  if (robust_std.errors != FALSE) {
    
    vcovCR <- vcovCR(selected_model, cluster = dgp_output$l2.id, type = robust_std.errors)
    
    robust_test <- coef_test(selected_model,
                             vcov = vcovCR,
                             test = "Satterthwaite")
    
    naive_pvals <- robust_test[ interesting_fixed_effects, "p_Satt"]
    
    naive_est <- robust_test$beta[which_coeffs]
    naive_se <- robust_test$SE[which_coeffs]
    df <- robust_test$df_Satt[which_coeffs]
    
    # Critical t-value
    t_crit <- qt(1 - 0.05 / 2, df = df)
    
    # Confidence interval bounds
    lower <- naive_est - t_crit * naive_se
    upper <- naive_est + t_crit * naive_se
  }
  
  names(naive_pvals) <- names(naive_est) <- names(lower) <- names(upper) <-  interesting_fixed_effects
  
  df_naive <- data.frame(est=naive_est,pvals=naive_pvals,ci.l=lower,ci.u=upper, se= naive_se)
  df_naive$coefficient <-  interesting_fixed_effects
  df_naive$method = "Naive"
  
  
  ## no_selection inference ## -----------------------------------------
  
  if (identical(interesting_working_model_ranef,c("(Intercept)", "X2"))) {
    no_selection_model <- lmer( y ~ 1 + X1 + X2 + X3 + X4 + (1 + X2 | l2.id), 
                                REML = FALSE, data = no_selection_dat) }
  else if (identical(interesting_working_model_ranef,c("(Intercept)", "X2","X4"))) {
    no_selection_model <- lmer( y ~ 1 + X1 + X2 + X3 + X4 + (1 + X2 + X4| l2.id), 
                                REML = FALSE, data = no_selection_dat) }
  else{ no_selection_model <- lmer( y ~ 1 + X1 + X2 + X3 + X4 + (1 | l2.id), 
                                    REML = FALSE, data = no_selection_dat)}
  
  summ_no_selection_model <- summary(no_selection_model)
  no_selection_est <- summ_no_selection_model$coefficients[ interesting_fixed_effects, "Estimate"]
  
  # do not use robust standard errors
  if (robust_std.errors == FALSE) {
    
    no_selection_pvals <- summ_no_selection_model$coefficients[ interesting_fixed_effects, "Pr(>|t|)"]
    no_selection_se <- summ_no_selection_model$coefficients[ interesting_fixed_effects, "Std. Error"]
    
    # Critical t-value (we used the 1.965 value instead)
    t_crit <- 1.965
    
    # Confidence interval bounds
    lower <- no_selection_est - t_crit * no_selection_se
    upper <- no_selection_est + t_crit * no_selection_se
  }
  
  # use robust standard errors
  if (robust_std.errors != FALSE) {
    
    vcovCR <- vcovCR(no_selection_model, cluster = dgp_output$l2.id, type = robust_std.errors)
    
    robust_test <- coef_test(no_selection_model,
                             vcov = vcovCR,
                             test = "Satterthwaite")
    
    no_selection_pvals <- robust_test[ interesting_fixed_effects, "p_Satt"]
    
    no_selection_est <- robust_test$beta[which_coeffs]
    no_selection_se <- robust_test$SE[which_coeffs]
    df <- robust_test$df_Satt[which_coeffs]
    
    # Critical t-value
    t_crit <- qt(1 - 0.05 / 2, df = df)
    
    # Confidence interval bounds
    lower <- no_selection_est - t_crit * no_selection_se
    upper <- no_selection_est + t_crit * no_selection_se
  }
  
  names(no_selection_pvals) <- names(no_selection_est) <- names(lower) <- names(upper) <-  interesting_fixed_effects
  
  df_no_selection <- data.frame(est=no_selection_est,pvals=no_selection_pvals,ci.l=lower,ci.u=upper,se= no_selection_se)
  df_no_selection$coefficient <-  interesting_fixed_effects
  df_no_selection$method = "no_selection"
  
  
  return(rbind.data.frame(df_naive,df_no_selection))
  
}

