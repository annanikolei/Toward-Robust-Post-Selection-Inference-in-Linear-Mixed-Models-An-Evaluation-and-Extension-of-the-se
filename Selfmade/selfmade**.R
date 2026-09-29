################################################################################
# Main function (analogous to the mocasin function)
################################################################################

source("helper_functions.R")
source("CR-Standard_errors.R")
source("CovarianceMatrix-Function.R")
source("parallel_selective-inference-Function.R")
source("Testvector-Functions.R")

mocasin <- function(mod, # fit of the selected model
                    balanced = TRUE, # equal number of Level-1 units per cluster
                    checkFun = NULL, # function for checking agreement
                                     # with the initial selection event
                    checkConvergence_for_check_Fun = F, # if TRUE, convergence is checked by
                    # checkFun and Tbs is resampled if convergence is not given
                    # then checkFun_checkConv must also be specified
                    checkFun_checkConv = NULL, # checkFun with an additional convergence check.
                    # Must return two arguments: one indicating agreement with the initial selection
                    # event and one indicating whether convergence is adequate
                    this_y, # vector of observed y-values
                    B_MCsamples = 1000, # number of MC samples
                    min_B_MCsamples = 0,
                    MC_distribution = NULL, # MC samples from
                                            # proposal distribution "proposal"
                                            # or from N(0,kappa_square)
                                            # i.e., from "centered"
                    bayesian = c("not_bayesian",
                                 "bayesian",
                                 "bayesian_corrected"), # Bayesian covariance matrix
                    which_marginal = NULL, # selective inference for which fixed effects
                    which_conditional = NULL, # selective inference for which fixed and random effects
                    rho_H0 = 0, # null hypothesis for all weights to be tested
                    perspective = c("marginal","conditional"), # conditional,
                                                               # marginal (efficient) or "marginal_not_efficient" or
                                                               # all together
                    robust_std.errors_marginal = F, # should robust standard errors be used
                                                    # for the marginal perspective?
                                                    # and if so, which ones?
                                                    # FALSE, "CR0", "CR2"
                    AdjMcCaffrey=TRUE, #boolean,if CR2 robust standard errors in not-efficient: robustness with adjument matrix by McCaffrey et al. (2001)
                    marginal_COV_estimator = c("ModelEstimate", 
                                               "ICM",
                                               "supplied"), # which estimators
                                                            # for the covariance matrix
                    marginal_not_efficient_COV_estimator = c("ModelEstimate", 
                                                             "ICM",
                                                             "supplied"), # which estimators
                                                                          # for the covariance matrix
                    diagF_not_efficient_SIGMA_or_R = c("SIGMA_replaced", "R_replaced"), # if Rdiag = F, the residual cross-product covariance matrix can be estimated using the ALM. Should it replace the entire SIGMA or only R?
                   
                    conditional_COV_estimator = c("ModelEstimate", 
                                                  "ICM",
                                                  "VarY",
                                                  "supplied"), # which estimators
                                                               # for the covariance matrix
                    HC_type = "HC0", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                    cadjust = F, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                    symmetric_selection = T, # symmetric selection event?
                    alpha = 0.05, # significance level
                    allow_smaller_alpha = T, # computation of larger confidence intervals?
                    multiple = 7, # regulates sequence width in pval_CI function
                                  # when generating possible rho values for CI
                                  # calculation
                    R_diag = TRUE, # R = sigma^2 * IN (homogeneous + independent) or not?                  
                    Sigma_supplied = NULL,
                    Z_supplied = NULL,
                    G_supplied = NULL,
                    Lambda_supplied = NULL,
                    Lambdat_supplied = NULL,
                    G_work_supplied = NULL, # conditional perspective supplied
                    sigma_square_epsilon_supplied = NULL,
                    R_supplied = NULL){
  

  #-----------------------------------------------------------------------------
  #### Extract important parameters from the fitted model mod
  #-----------------------------------------------------------------------------
  
  
  if(("supplied" %in% marginal_not_efficient_COV_estimator)& robust_std.errors_marginal!=FALSE){

    warning("If supplied, we assume that robust standard errors are not needed. Supplied is therefore excluded")

    marginal_not_efficient_COV_estimator <- 
      marginal_not_efficient_COV_estimator[marginal_not_efficient_COV_estimator != "supplied"]

  }
  
  if(("supplied" %in% marginal_COV_estimator)& robust_std.errors_marginal!=FALSE){
    
    warning("If supplied, we assume that robust standard errors are not needed. Supplied is therefore excluded")
    
    marginal_COV_estimator <- 
      marginal_COV_estimator[marginal_COV_estimator != "supplied"]
    
  }
  
  if(AdjMcCaffrey==FALSE & length(marginal_not_efficient_COV_estimator)>1){
    
    warning("Sigma does not play a relevant role here apart from the exact definition of dirV. Therefore, marginal_not_efficient_COV_estimator is set to ModelEstimate")
    
    marginal_not_efficient_COV_estimator <- c("ModelEstimate")
  }
  
  # Extract y-vector, if necessary
  # if(is.null(this_y)){this_y <- mod@resp}
  
  # Extract number of observations
  n <- length(this_y)
  number_of_cluster = ngrps(mod)[[1]]
  
  # Extract the working matrices X_A and Z_B from the model                          
  X_A_work <- getME(mod, "X")
  Z_B_work <- getME(mod, "Z") # in the marginal perspective, no random effects
                              # are selected anyway, therefore
                              # generally corresponding to Z_supplied, but
                              # we still do not provide this as Z_supplied,
                              # in case we want to investigate such a violation
                              # in a later simulation.
  
  # Extract the estimated working Cholesky factor Lambda from the model.
  # It is needed for all estimators in the conditional perspective and for ModelEstimate
  # in the marginal perspective, so we simply compute it in all cases                      
  # The following holds: G_est = sigma_square_epsilon * Lambda %*% t(Lambda).
  # https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf (Equation 4)
  Lambda_work_est <- getME(mod, "Lambda")
  Lambdat_work_est <- getME(mod, "Lambdat") # = t(Lambda_work_est)
  G_work_est <- sigma(mod)^2 * Lambda_work_est %*% Lambdat_work_est                 
  

  #-----------------------------------------------------------------------------
  #### Construct or estimate the (inverse) covariance matrices
  #-----------------------------------------------------------------------------
  # If one of the two perspectives, marginal or conditional, is not specified,
  # set the corresponding COV_estimator to NULL
  
  if(!("marginal" %in% perspective)){
    marginal_COV_estimator <- NULL
  }
  
  if(!("marginal_not_efficient" %in% perspective)){
    marginal_not_efficient_COV_estimator <- NULL
    
    if(R_diag==TRUE){ diagF_not_efficient_COV_estimator <- NULL  }
  }
  
  if(!("conditional" %in% perspective)){
    conditional_COV_estimator <- NULL
  }
  
  # Determine which estimators need to be computed
  if(!is.null(marginal_COV_estimator)){
   marginal_COV_estimator <- match.arg(marginal_COV_estimator,
                                      c("ModelEstimate", 
                                        "ICM",
                                        "supplied"),
                                         several.ok = T)}
  
  if(!is.null(marginal_not_efficient_COV_estimator)){
    marginal_not_efficient_COV_estimator <-
      match.arg(marginal_not_efficient_COV_estimator,
                                        c("ModelEstimate", 
                                          "ICM",
                                          "supplied"),
                                        several.ok = T)}
  
  if(!is.null( diagF_not_efficient_SIGMA_or_R )){
    diagF_not_efficient_SIGMA_or_R <-
      match.arg(diagF_not_efficient_SIGMA_or_R,
                c("SIGMA_replaced", "R_replaced"),
                several.ok = T)}
  
  if(!is.null(conditional_COV_estimator)){ 
   conditional_COV_estimator <- match.arg(conditional_COV_estimator,
                                          c("ModelEstimate",
                                            "ICM",
                                            "VarY",
                                            "supplied"),
                                             several.ok = T)}
  
  
  
  
  ## Preliminary: computation of R and R_invers for the different estimators--------------
  
  fit_ICM <- NULL
  # Sigma_ICM <- NULL
  
  if(R_diag == T){
    
    if(c("ICM") %in% marginal_COV_estimator |
       c("ICM") %in% marginal_not_efficient_COV_estimator |
       c("ICM") %in% conditional_COV_estimator){
      fit_ICM <- minMod(mod)}
    
  R_R_inv_list <- R_R_invers(mod = mod,
             n = n,
             fit_ICM = fit_ICM,
             this_y = this_y,
             marginal_COV_estimator = unique(c(marginal_COV_estimator, marginal_not_efficient_COV_estimator)),
             conditional_COV_estimator = conditional_COV_estimator,
             Sigma_supplied = Sigma_supplied,
             sigma_square_epsilon_supplied = sigma_square_epsilon_supplied)
  
  }
  
  else if(R_diag == F){
    
    
  ## marginal or conditional perspective
    # Working with a non-diagonal R is currently only possible for the supplied COV_estimator
    # for the marginal and conditional perspectives
    if(!any( c("ICM","ModelEstimate","VarY") %in% c(marginal_COV_estimator,conditional_COV_estimator) ) ){
      if(is.null(R_supplied) && any(c("marginal","conditional") %in% perspective)){ 
      stop("R_diag = FALSE indicates that a non-diagonal
       Level-1 residual covariance matrix should be used. Currently, however, this is
       only possible when a non-diagonal Level-1 residual covariance matrix is supplied
       for the marginal and conditional perspectives.")
      }
      
      else if(any(c("ICM","ModelEstimate","VarY") %in% c(marginal_COV_estimator,conditional_COV_estimator))){
        stop("R_diag = FALSE indicates that a non-diagonal
       Level-1 residual covariance matrix should be used. Currently, this is only
       possible for diagonal Level-1 residual covariance matrices or when a non-diagonal
       Level-1 residual covariance matrix is supplied. An error is also returned if,
       in addition to supplied, any other estimator is specified.")
      }
    
      
  ## marginal not efficient perspective 
    if("marginal_not_efficient" %in% perspective){
      
      if(robust_std.errors_marginal != FALSE){stop("robust estimator for marginal not efficient is either produced with R_diag=F or via robust standard errors")}
      
      # For the not-efficient approach, we can essentially use the ALM formulas
      # and, in contrast to the ALM, allow other matrices than sigmaepsilon^2*I_N
      # for Sigma by using the residuals from the ALM.
      # From my perspective, there are several possible ways in which the estimator
      # could be defined.
      # Either the entire Sigma is replaced by this residual cross-product matrix.
      # Or only R is replaced and the ZGZ' component is added.
      
      # For ModelEstimate, the residual cross-product matrix can either be used as SIGMA
      # or as R, in which case the working-model component ZGZ' is added.
      
      # If R were replaced for ICM, we would additionally need G from the minimal
      # model. However, the minimal model is not an LMM but an ALM y~1, so we would
      # either have to specify that only SIGMA can be replaced here or use G from the
      # general working model mod.
      
      # Therefore, if R_diag==F, we need to specify whether R or SIGMA
      # should be replaced by the residual cross-product matrix.
      
      
      # Elements other than "supplied"
      marginal_not_efficient_COV_estimator_no_supplied <- marginal_not_efficient_COV_estimator[!grepl("supplied", marginal_not_efficient_COV_estimator)]
      
      # Combinations only for non-supplied estimators
      combos <- as.vector(outer(marginal_not_efficient_COV_estimator_no_supplied, diagF_not_efficient_SIGMA_or_R, paste, sep = "_"))
      
      # Append supplied separately
      marginal_not_efficient_COV_estimator  <- c(combos, marginal_not_efficient_COV_estimator[grepl("supplied", marginal_not_efficient_COV_estimator)])
      
      
      if( (is.null(R_supplied)|is.null(Z_supplied)|is.null(G_supplied)) && 
           is.null(Sigma_supplied) && 
           "supplied" %in% marginal_not_efficient_COV_estimator){
        
        stop("diagF_not_efficient_SIGMA_or_R indicates that R,Z,G or Sigma was supplied.
             However, R,Z,G_supplied or Sigma_supplied is NULL.")
      }
         
    }
     
      R_R_inv_list <- NULL
      
      if(any(c("marginal","conditional") %in% perspective)){ 
        
      R_R_inv_list <- R_R_invers_nondiagonal(marginal_COV_estimator = marginal_COV_estimator,
                                             conditional_COV_estimator = conditional_COV_estimator,
                                             R_supplied = R_supplied) 
      
      }
      
      
      
    }
     
  
    
    
  }
  
  
  ##------------------------------------------------------------------------------ 
 
  ### Marginal perspective
  if (c("marginal") %in%  perspective) {
    SIGMA_invers <- 
      SIGMA_invers_function( marginal_COV_estimator =  marginal_COV_estimator,
                             Sigma_supplied = Sigma_supplied,
                             Z_supplied = Z_supplied,
                             G_supplied = G_supplied,
                             Lambda_supplied = Lambda_supplied,
                             Lambdat_supplied = Lambdat_supplied,
                             R_supplied = R_R_inv_list$R$supplied,
                             R_supplied_invers = R_R_inv_list$R_inv$supplied,
                             R_est = R_R_inv_list$R$ModelEstimate,
                             R_est_invers = R_R_inv_list$R_inv$ModelEstimate,
                             Z_B_work = Z_B_work,
                             G_work_est = G_work_est,
                             fit_ICM = fit_ICM,
                             R_ICM = R_R_inv_list$R$ICM,
                             R_ICM_invers = R_R_inv_list$R_inv$ICM,
                             n = n )
  }
  
  ### Marginal perspective (not efficient)
  if (c("marginal_not_efficient") %in%  perspective) {
    
  if(!balanced){ stop("The estimation of clustered errors assumes a balanced design.")}
    
  if(R_diag == T){
    SIGMA <- 
      SIGMA_function_diagT( marginal_not_efficient_COV_estimator =  marginal_not_efficient_COV_estimator,
                            Sigma_supplied = Sigma_supplied,
                            Z_supplied = Z_supplied,
                            G_supplied = G_supplied,
                            Lambda_supplied = Lambda_supplied,
                            Lambdat_supplied = Lambdat_supplied,
                            R_supplied = R_R_inv_list$R$supplied,
                            R_est = R_R_inv_list$R$ModelEstimate,
                            Z_B_work = Z_B_work,
                            G_work_est = G_work_est,
                            fit_ICM = fit_ICM,
                            R_ICM = R_R_inv_list$R$ICM,
                            n = n )
    }
    
    if(R_diag == F){
      SIGMA <- 
        SIGMA_function_diagF(  mod = mod,
                               marginal_not_efficient_COV_estimator,
                               G_work_est = G_work_est,
                               Z_B_work = Z_B_work,
                               R_supplied = R_supplied,
                               Z_supplied = Z_supplied,
                               G_supplied = G_supplied,
                               Sigma_supplied = Sigma_supplied,
                               balanced = balanced,
                               HC_type = HC_type, 
                               cadjust = cadjust, 
                               n = n,
                               number_of_cluster = number_of_cluster )
      }
    
  }
  
  ### Conditional perspective
  if (c("conditional") %in% perspective) {
    
    K_inv.Ct.R_inv_C_block <- K_inv.Ct.R_inv_list_function(
    conditional_COV_estimator = conditional_COV_estimator,
    X_A_work = X_A_work,
    Z_B_work = Z_B_work,
    G_work_est = G_work_est,
    R_supplied_invers = R_R_inv_list$R_inv$supplied,
    G_work_supplied = G_work_supplied,
    R_est_invers = R_R_inv_list$R_inv$ModelEstimate,
    R_ICM_invers = R_R_inv_list$R_inv$ICM,
    R_VarY_invers = R_R_inv_list$R_inv$VarY)

  }
  
  
  #-----------------------------------------------------------------------------
  #### Test vector v or vt, variance of the test statistic, and dirV
  #-----------------------------------------------------------------------------

  
  ## marginal efficient 
  if(c("marginal") %in% perspective) {
    # If which_marginal is not specified, test all coefficients                    
    if(is.null(which_marginal)){which_marginal <- 1:ncol(X_A_work)}
    
    testvec_marginal_list <- testvec_marginal(
    marginal_COV_estimator = marginal_COV_estimator,
    which_marginal = which_marginal,
    X_A_work = X_A_work,
    SIGMA_invers = SIGMA_invers,
    robust_std.errors_marginal = robust_std.errors_marginal,
    mod = mod,
    this_y = this_y)
    
  }
  
  ## marginal not efficient 
  if(c("marginal_not_efficient") %in% perspective) {
    # If which_marginal is not specified, test all coefficients                     
    if(is.null(which_marginal)){which_marginal <- 1:ncol(X_A_work)}
    
    testvec_marginal_not_efficient_list <- testvec_marginal_not_efficient(
      AdjMcCaffrey=AdjMcCaffrey,
      mod= mod,
      this_y = this_y,
      marginal_not_efficient_COV_estimator = marginal_not_efficient_COV_estimator,
      robust_std.errors_marginal =  robust_std.errors_marginal,
      which_marginal = which_marginal,
      X_A_work = X_A_work,
      SIGMA = SIGMA)
    
  }
   
  ## conditional 
  if(c("conditional") %in% perspective) {
    
  # If which_conditional is not specified, test all coefficients
    if(is.null(which_conditional)){
      which_conditional <- 1:ncol(K_inv.Ct.R_inv_C_block$C_block)}
    
    testvec_conditional_list <- testvec_conditional(
    conditional_COV_estimator = conditional_COV_estimator,
    bayesian = bayesian,
    which_conditional = which_conditional,
    C = K_inv.Ct.R_inv_C_block$C_block,
    K_inv_list = K_inv.Ct.R_inv_C_block$K_inv_list,
    Ct.R_inv_list = K_inv.Ct.R_inv_C_block$Ct.R_inv_list,
    Ct.R_inv.C_list = K_inv.Ct.R_inv_C_block$Ct.R_inv.C_list,
    R_list = R_R_inv_list$R)
    
  }
  
  
  #-------------------------------------------------------------------------------
  #### Compute selective inference outputs
  #-------------------------------------------------------------------------------
  
  # Separately for marginal/conditional (COV_estimator, bayesian/not_bayesian, which),
  # compute selective-inference p-values and selective-inference CIs
  # using the selective_inference function (from one of the other scripts)
  
  # For a non-symmetric selection event, the centered
  # distribution should be used as the MC distribution.
  # For a symmetric selection event, the proposal
  # distribution is used unless specified otherwise.
  if(is.null(MC_distribution)){
    if(symmetric_selection == T) {MC_distribution <- "proposal"}
    else if(symmetric == F){MC_distribution <- "centered"}
}
  
  # Create a list with elements corresponding to the perspective, i.e.,
  # one list element $marginal and one $conditional
  selinf_results_list <- setNames(vector("list",length(perspective)),perspective)
  
  ## marginal
  if("marginal" %in% perspective){
    selinf_results_list[["marginal"]] <-
      lapply(marginal_COV_estimator, function(iestimator){
          lapply(which_marginal, function(iwhich){
            
            # Name of the coefficient
            coef_name <- paste("coefficient", iwhich)
            
            ## Input for selective_inference          
            # Which test vector vt:
            vt = testvec_marginal_list$vt_marginal_list[[iestimator]][[coef_name]]
            # Which variance kappa_square of the test statistic
            kappa_square = testvec_marginal_list$kappa_square_marginal_list[[iestimator]][[coef_name]]
            # Which dirV
            dirV = testvec_marginal_list$dirV_marginal_list[[iestimator]][[coef_name]]
            
            
            result <- selective_inference(
                n = n,
                vt = vt,
                kappa_square = kappa_square, 
                dirV = dirV,
                this_y = this_y,
                B_MCsamples = B_MCsamples,
                min_B_MCsamples = min_B_MCsamples,
                checkFun = checkFun,
                MC_distribution = MC_distribution,
                symmetric_selection = symmetric_selection,
                alpha = alpha,
                allow_smaller_alpha = allow_smaller_alpha,
                rho_H0 = rho_H0,
                checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
                checkFun_checkConv = checkFun_checkConv
            )
            # Generating a dynamic name for the result
             result_name <- paste(iestimator,coef_name, sep = "_")
           
             # Assign result with the generated name to the list
             return(setNames(list(result), result_name))
      
        })
      })
    }

  ## marginal not efficient
  if("marginal_not_efficient" %in% perspective){
    selinf_results_list[["marginal_not_efficient"]] <-
      lapply(marginal_not_efficient_COV_estimator, function(iestimator){
        lapply(which_marginal, function(iwhich){
          
          # Name of the coefficient
          coef_name <- paste("coefficient", iwhich)
          
          ## Input for selective_inference          
          # Which test vector vt:
          vt = testvec_marginal_not_efficient_list$vt_marginal_not_efficient_list[[iestimator]][[coef_name]]
          # Which variance kappa_square of the test statistic
          kappa_square = testvec_marginal_not_efficient_list$kappa_square_marginal_not_efficient_list[[iestimator]][[coef_name]]
          # Which dirV
          dirV = testvec_marginal_not_efficient_list$dirV_marginal_not_efficient_list[[iestimator]][[coef_name]]
        
        
          result <- selective_inference(
            n = n,
            vt = vt,
            kappa_square = kappa_square, 
            dirV = dirV,
            this_y = this_y,
            B_MCsamples = B_MCsamples,
            min_B_MCsamples = min_B_MCsamples,
            checkFun = checkFun,
            MC_distribution = MC_distribution,
            symmetric_selection = symmetric_selection,
            alpha = alpha,
            allow_smaller_alpha = allow_smaller_alpha,
            rho_H0 = rho_H0,
            checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
            checkFun_checkConv = checkFun_checkConv
        )
        # Generating a dynamic name for the result
        result_name <- paste(iestimator,coef_name, sep = "_")
        
        # Assign result with the generated name to the list
        return(setNames(list(result), result_name))
        
      })
    })
}


## conditional
if("conditional" %in% perspective){
  selinf_results_list[["conditional"]] <-
      lapply(conditional_COV_estimator, function(iestimator){
       lapply(bayesian, function(ibayesian){
        lapply(which_conditional, function(iwhich){
          
           # Name of the coefficient
           coef_name <- paste("coefficient", iwhich)
          
     ## Input for selective_inference
           # Which test vector
           vt = 
           testvec_conditional_list$vt_conditional_list[[iestimator]][[coef_name]]
           # Which variance kappa_square of the test statistic
           kappa_square = 
           testvec_conditional_list$kappa_square_conditional_list[[iestimator]][[ibayesian]][[coef_name]]
           # Which dirV
           dirV = 
           testvec_conditional_list$dirV_conditional_list[[iestimator]][[ibayesian]][[coef_name]]
           
           result <- selective_inference(
             n = n, 
             vt = vt,
             kappa_square = kappa_square, 
             dirV = dirV,
             this_y = this_y,
             B_MCsamples = B_MCsamples,
             min_B_MCsamples = min_B_MCsamples,
             checkFun = checkFun,
             MC_distribution = MC_distribution,
             symmetric_selection = symmetric_selection,
             alpha = alpha,
             allow_smaller_alpha = allow_smaller_alpha,
             rho_H0 = rho_H0,
             checkConvergence_for_check_Fun = checkConvergence_for_check_Fun,
             checkFun_checkConv = checkFun_checkConv)
           
           # Generating a dynamic name for the result
           result_name <- paste(iestimator,ibayesian,coef_name, sep = "_")
           
           # Assign result with the generated name to the list
           return(setNames(list(result), result_name))
        })
       })
      })
}
    
  
  
  #-------------------------------------------------------------------------------
  #### Return output
  #-------------------------------------------------------------------------------
  
return(selinf_results_list)

}
