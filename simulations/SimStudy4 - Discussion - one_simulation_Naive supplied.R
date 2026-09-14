#####################################################################################################
# one simualtion function Discussion - Naive Post-Selection with supplied - Simulation Study 4
#####################################################################################################


Naive_supplied_one_simulation <- function(l2.id, # id variable for level-2 units
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
                                         checkConvergence_for_check_Fun = TRUE,
                                         remove_singular_fits = F, #remove singular fit 
                                         alpha = 0.000001, # alpha for significance selection
                                         perspective = c("marginal_not_efficient"),
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
  
  
  
  ## Naive Post-selective inference with supplied --------------------------------------------------
  
  
  ## selfmade ## ---------------------------------------------------------------
  
  checkFun <- function(yb){ 
    
    return(TRUE)
    
  }
  
  checkFun_checkConv <- function(yb){ 
    
    newdat <- data.frame(y = yb, X, l2.id)
    new_full_model <- modFun(newdat)
    convergence_full_okay = check_convergence(new_full_model)
    
    if(convergence_full_okay){
      new_selected_model <- SelFun( newdat, new_full_model)
      in_inference_region <- TRUE
      convergence_okay <- check_convergence(new_selected_model)
      
    }
    
    if(!convergence_full_okay){
      in_inference_region <- NA
      convergence_okay = FALSE
    }
    
    
    return( list(in_inference_region = in_inference_region,
                 convergence_okay = convergence_okay)  )
    
  }
  
  
  selfmade_output <-
    mocasin(mod = selected_model, 
            checkFun = checkFun,
            checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
            checkFun_checkConv = checkFun_checkConv,
            this_y = dat$y, 
            B_MCsamples = 500, # number of MC-samples
            MC_distribution = "centered",
            bayesian = c("bayesian"),
            which_marginal = which_coeffs, 
            which_conditional = which_coeffs, 
            rho_H0 = 0, 
            perspective = perspective, 
            robust_std.errors_marginal = F,
            marginal_COV_estimator = c("supplied"), 
            conditional_COV_estimator = c("ModelEstimate","VarY","supplied"), 
            marginal_not_efficient_COV_estimator = c("supplied"), 
            diagF_not_efficient_SIGMA_or_R = c("SIGMA_replaced", "R_replaced"),
            symmetric_selection = T, 
            alpha = 0.05, 
            allow_smaller_alpha = T, 
            multiple = 7, 
            R_diag = T, 
            HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
            cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
            Sigma_supplied = NULL,
            Z_supplied = Z,
            G_supplied = dgp_output$G,
            Lambda_supplied = NULL,
            Lambdat_supplied = NULL,
            G_work_supplied = NULL, 
            sigma_square_epsilon_supplied = 
              as.numeric(dgp_output$sigma_squared_epsilon),
            R_supplied = Matrix(diag(as.numeric(dgp_output$sigma_squared_epsilon),n_obs), sparse = T)) 
  
  
  ## Create ouput --------------------------------------------------------------
  
  ## selfmade output
  selfmade_output_df <- extract_pval_two(selfmade_output)
  coefficient_num <- as.integer(sub("coefficient ", "",
                                    selfmade_output_df$coefficient))
  selfmade_output_df$coefficient <- interesting_fixed_effects
  

  
  return(selfmade_output_df)
  
  
}
