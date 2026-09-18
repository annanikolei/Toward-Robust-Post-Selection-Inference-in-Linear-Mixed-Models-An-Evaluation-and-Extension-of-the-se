################################################################################
# Simulation Study 3-Discussion: Conditional Fixed RE (One Simulation function)
################################################################################

library(performance)
library(clubSandwich)  

####################### one_simulation function ################################

SimStudy3_one_simulation <- function(l2.id, # id variable for level-2 units
                                     n_obs, # number of observations
                                     m, # number of level-2 units
                                     X, # X-matrix 
                                     X.quadratic, # X-matrix (including squared variables)
                                     Z, # design matrix corresponding to X
                                     R_squared = 0.1, # pseudo coefficient of determination
                                     D = matrix(c(2,0.2,0.1,
                                                  0.2,1.5,0.2,
                                                  0.1,0.2,1.5),
                                                nrow = 3, byrow=T) , # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                     list_supplied_D_work,
                                     u_fixed = u_fixed,
                                     perspective = "conditional", 
                                     bayesian= c("bayesian_corrected"),
                                     B_MCsamples = 1000, # number of MC-samples
                                     remove_singular_fits = F, #remove singular fit 
                                     interesting_var, # "X3" or "X4"
                                     interesting_working_model_vars = NULL, # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                                     checkConvergence_for_check_Fun = T,
                                     target # rho 
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
          any(interesting_var %in% selected_var_names) &&
          working_equal_to_selected == T) ){
    
    ## Generate data (y and fixed u) by calling the dgp function ---------------------
   
    
    dgp_output <- SimStudy3_dgp_y(l2.id = l2.id, 
                                  n_obs = n_obs,
                                  m = m, 
                                  X = X, 
                                  X.quadratic = X.quadratic, 
                                  Z = Z, 
                                  R_squared = R_squared, 
                                  u = u_fixed$u,
                                  D = u_fixed$D )
    
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
       # print(selected_var_names)
        
        if(is.null(interesting_working_model_vars)){
          working_equal_to_selected = TRUE
        }
        else if(!is.null(interesting_working_model_vars)){
          working_equal_to_selected = setequal(interesting_working_model_vars,selected_var_names)
        }
        
        if(any(interesting_var %in% selected_var_names)) {
          
          tries_selected_model <- tries_selected_model + 1
          if (!convergence_okay) convergence_failures_selected_model <- convergence_failures_selected_model + 1
          
        }
      }
    }
  } # end while-loop
  
  ### Inference ------------------------------------------------------------------
  
  sel.noise <- match.arg(interesting_var,selected_var_names, several.ok = T) 
  which_coeffs <- which(selected_var_names %in% sel.noise)  + 1 # +1 to account for Intercept                                                                               #?# bisher nur noise
  
  
  ## Post-selective inference --------------------------------------------------

  
  
  ## selfmade ## ---------------------------------------------------------------
  
  checkFun <- function(yb){ 
    
    newdat <- data.frame(y = yb, X, l2.id)
    new_full_model <- modFun(newdat)
    new_selected_model <- SelFun(new_full_model)
    
    return(
      setequal( extractSelFun_fixed(new_selected_model), 
                selected_var_names) )
    
    
  }
  
  
  
  checkFun_checkConv <- function(yb){ 
    
    newdat <- data.frame(y = yb, X, l2.id)
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
  
  # extract mean estimated working covariance matrix G for selected model
  working_vars <- paste0(sort(selected_var_names), collapse = "")
  mean_D_work_estimated = 
    list_supplied_D_work[[paste0("R_squared = ",R_squared)]][[working_vars]]
  mean_G_work_estimated <- bdiag(replicate(m, mean_D_work_estimated, simplify = F))
  
  if(length(interesting_var)>1){stop("length(interesting_var) should be 1")}
  
  if(interesting_var=="X4"){
    if(target == "rho") { rhoH0 = u_fixed$mean_u[3]}
    else if(target =="rho_smile") { rhoH0 = u_fixed$mean_u[3]}
  }
  else if(interesting_var=="X3"){rhoH0=0}
 
  
  selfmade_output <-
    mocasin(mod = selected_model, 
            checkFun = checkFun,
            checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
            checkFun_checkConv = checkFun_checkConv,
            this_y = dat$y, 
            B_MCsamples = B_MCsamples, # number of MC-samples
            MC_distribution = "proposal",
            bayesian = bayesian,
            which_conditional = which_coeffs, 
            rho_H0 = rhoH0, 
            perspective = perspective, 
            conditional_COV_estimator = c("ModelEstimate","ICM","VarY","supplied"), 
            symmetric_selection = T, 
            alpha = 0.05, 
            allow_smaller_alpha = T, 
            multiple = 7, 
            R_diag = TRUE, 
            Sigma_supplied = NULL,
            Z_supplied = Z,
            G_supplied = dgp_output$G,
            Lambda_supplied = NULL,
            Lambdat_supplied = NULL,
            G_work_supplied = mean_G_work_estimated, 
            sigma_square_epsilon_supplied = 
              as.numeric(dgp_output$sigma_squared_epsilon),
            R_supplied = Matrix(diag(as.numeric(dgp_output$sigma_squared_epsilon),n_obs), sparse = T)) 
  
  
  
  ## Create ouput --------------------------------------------------------------
  
  ## selfmade output
  selfmade_output_df <- extract_pval_two(selfmade_output)
  coefficient_num <- as.integer(sub("coefficient ", "",
                                    selfmade_output_df$coefficient))
  selfmade_output_df$coefficient <- selected_var_names[coefficient_num - 1] # -1 to account for Intercept
  
  
  
  ## combined output 
  post_selection_pvalues_withoutfullCandidate <- selfmade_output_df
  post_selection_pvalues <- post_selection_pvalues_withoutfullCandidate
  post_selection_pvalues$working_model <- 
     rep(working_vars,nrow(post_selection_pvalues_withoutfullCandidate))
  
  estimate <- selfmade_output[[1]][[1]][[1]][[1]]$`ModelEstimate_bayesian_corrected_coefficient 4`[,"t_obs"]
  
  return(list(post_selection_pvalues = post_selection_pvalues,
              convergence_singular_info =
                data.frame(
                  convergence_failures_full_model = convergence_failures_full_model,
                  convergence_okay_full_model = convergence_okay_full_model,
                  convergence_failures_selected_model = convergence_failures_selected_model,
                  tries_selected_model = tries_selected_model,
                  full_model_singular = check_singularity(full_model),
                  selected_model_singular = check_singularity(selected_model),
                  working_model = working_vars),
              estimate=estimate)
  )
  
}







