################################################################################
# Main function (mocasin-function)
################################################################################

source("install.packages.R")
source("helper_functions.R")
source("CR-Standard_errors.R") # not important
source("CovarianceMatrix-Function.R")
source("parallel_selective-inference-Function.R") # or selective-inference-Function.R
source("Testvector-Functions.R")

mocasin <- function(mod, # fitted model (selected model fit)
                    balanced = TRUE, # equal number of level-1 units per cluster?
                    checkFun = NULL, # function to check agreement with the initial selection event
                    checkConvergence_for_check_Fun = F, # if TRUE, convergence is checked within checkFun
                    # if convergence fails, resampling is performed
                    # in that case, checkFun_checkConv must also be provided
                    checkFun_checkConv = NULL, # checkFun with additional convergence check
                    # must return two outputs: (1) agreement with initial selection event
                    # and (2) whether convergence is acceptable
                    this_y, # vector of observed y values
                    B_MCsamples = 1000, # number of Monte Carlo samples
                    MC_distribution = NULL, # MC samples drawn either from
                                            # proposal distribution ("proposal")
                                            # or from N(0, kappa^2) ("centered")
                    bayesian = c("not_bayesian",
                                 "bayesian"), # whether to use Bayesian covariance matrix
                    which_marginal = NULL, # selection of fixed effects (marginal perspective)
                    which_conditional = NULL, # selection of fixed + random effects (conditional perspective)
                    rho_H0 = 0, # null hypothesis value for all tested coefficients
                    perspective = c("marginal","conditional"), # inference perspective:
                                                               # "conditional", "marginal" (efficient),
                                                               # or "marginal_not_efficient", or all combined
                    robust_std.errors_marginal = F, # whether to use robust standard errors
                                                    # in marginal perspective
                                                    # options: FALSE, "CR0", "CR2". Must be FALSE because it does not work

                    
                    R_diag = TRUE, # whether R = sigma^2 * I (homoskedastic + independent errors)
                    
                    marginal_COV_estimator = c("ModelEstimate", 
                                               "ICM",
                                               "supplied"), # estimator for covariance matrix (marginal)
                    marginal_not_efficient_COV_estimator = c("ModelEstimate", 
                                                             "ICM",
                                                             "supplied"), # covariance estimator for not-efficient marginal case
                    
                    diagF_not_efficient_SIGMA_or_R = c("SIGMA_replaced", "R_replaced"), 
                    # if R_diag = FALSE, proposed covariance matrix extimator: should it replace full SIGMA or only R?
                    
                    conditional_COV_estimator = c("ModelEstimate", 
                                                  "ICM",
                                                  "VarY",
                                                  "supplied"), # estimator for covariance matrix (conditional)
                    
                    HC_type = "HC0", # heteroskedasticity-consistent type (HC0, HC1, HC2, HC3)
                    cadjust = F, # cluster adjustment (TRUE/FALSE)
                    
                    symmetric_selection = T, # symmetric selection event?
                    
                    alpha = 0.05, # significance level
                    allow_smaller_alpha = T, # allow computation of larger confidence intervals?
                    
                    multiple = 7, # controls sequence width in p-value/CI function
                                  # used to generate candidate rho values for CI computation
                    
                    
                    Sigma_supplied = NULL,
                    Z_supplied = NULL,
                    G_supplied = NULL,
                    Lambda_supplied = NULL,
                    Lambdat_supplied = NULL,
                    G_work_supplied = NULL, # supplied G for conditional perspective
                    sigma_square_epsilon_supplied = NULL,
                    R_supplied = NULL){

#-----------------------------------------------------------------------------
#### Extract important parameters from the fitted model `mod`
#-----------------------------------------------------------------------------

# Extract y vector if needed
# if (is.null(this_y)) { this_y <- mod@resp }

# Extract number of observations
n <- length(this_y)

# Extract number of clusters (level-2 units)
number_of_cluster <- ngrps(mod)[[1]]

# Extract working design matrices X_A (fixed effects) and Z_B (random effects)
# Note: In the marginal perspective, no random effects are selected,
# so Z_B_work typically coincides with Z_supplied.
# However, we do NOT pass it as Z_supplied, in case we want to
# simulate violations of this assumption later.
X_A_work <- getME(mod, "X")
Z_B_work <- getME(mod, "Z")

# Extract the estimated working Cholesky factor Lambda from the model.
# This is required for:
# - all estimators in the conditional perspective
# - the ModelEstimate in the marginal perspective
# Therefore, we compute it unconditionally.
#
# Relationship:
# G_est = sigma^2 * Lambda %*% t(Lambda)
# See: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf (Formula 4)
Lambda_work_est  <- getME(mod, "Lambda")
Lambdat_work_est <- getME(mod, "Lambdat") # = t(Lambda_work_est)

# Compute estimated random effects covariance matrix G
G_work_est <- sigma(mod)^2 * Lambda_work_est %*% Lambdat_work_est

  
#-----------------------------------------------------------------------------
#### Construct or estimate the (inverse) covariance matrices
#-----------------------------------------------------------------------------
  
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
  
  
  
  
## Calculation of R,R_invers for the different estimators --------------
  
  fit_ICM <- NULL
  Sigma_ICM <- NULL
  
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

# Currently, working with non-diagonal R is only supported when using
# a supplied covariance estimator in both marginal and conditional perspectives
if (!any(c("ICM","ModelEstimate","VarY") %in% c(marginal_COV_estimator, conditional_COV_estimator))) {
  
  # Case: No supported estimator is used and R is not supplied
  if (is.null(R_supplied) && any(c("marginal","conditional") %in% perspective)) { 
    stop("R_diag = FALSE indicates that a non-diagonal level-1 residual covariance matrix should be used. 
Currently, this is only possible if a non-diagonal level-1 residual covariance matrix is supplied 
for the marginal and/or conditional perspective.")
  }
  
  # Case: Conflicting estimators are specified together with non-diagonal R
  else if (any(c("ICM","ModelEstimate","VarY") %in% c(marginal_COV_estimator, conditional_COV_estimator))) {
    stop("R_diag = FALSE indicates that a non-diagonal level-1 residual covariance matrix should be used. 
Currently, this is only supported for diagonal level-1 residual covariance matrices, or if a non-diagonal 
matrix is explicitly supplied. This error also occurs if 'supplied' is combined with any other estimator.")
  }  
      
  ## marginal not efficient perspective 
    if("marginal_not_efficient" %in% perspective){
      
      # Elemente without "supplied"
      marginal_not_efficient_COV_estimator_no_supplied <- marginal_not_efficient_COV_estimator[!grepl("supplied", marginal_not_efficient_COV_estimator)]
      
      # combinations for not supplied
      combos <- as.vector(outer(marginal_not_efficient_COV_estimator_no_supplied, diagF_not_efficient_SIGMA_or_R, paste, sep = "_"))
      
      # add supplied 
      marginal_not_efficient_COV_estimator  <- c(combos, marginal_not_efficient_COV_estimator[grepl("supplied", marginal_not_efficient_COV_estimator)])
      
      
      if( (is.null(R_supplied)|is.null(Z_supplied)|is.null(G_supplied)) && 
           is.null(Sigma_supplied) && 
           "supplied" %in% marginal_not_efficient_COV_estimator){
        
       stop("diagF_not_efficient_SIGMA_or_R indicates that R, Z, G, or Sigma is supplied, 
             but the corresponding object (R_supplied, Z_supplied, G_supplied, or Sigma_supplied) is NULL.")
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
    
  if(!balanced){stop("Estimation of clustered errors assumes a balanced design.")}
    
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
#### Testvector v  and vt, Variance of test statistic and dirV
  #-----------------------------------------------------------------------------

  
## marginal efficient 
  if(c("marginal") %in% perspective) {
    # if which_marginal not specified, test all coefficients                  
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
     # if which_marginal not specified, test all coefficients    
    if(is.null(which_marginal)){which_marginal <- 1:ncol(X_A_work)}
    
    testvec_marginal_not_efficient_list <- testvec_marginal_not_efficient(
      marginal_not_efficient_COV_estimator = marginal_not_efficient_COV_estimator,
      which_marginal = which_marginal,
      X_A_work = X_A_work,
      SIGMA = SIGMA)
    
  }
   
## conditional 
  if(c("conditional") %in% perspective) {
    
  #  # if which_conditional not specified, test all coefficients    
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
#### Computation of selective inference outputs
#-------------------------------------------------------------------------------

# Separately for marginal / conditional perspective (COV_estimator, bayesian/not_bayesian, which)
# Compute selective inference p-values and confidence intervals
# using the selective_inference function 

# Choice of MC distribution:
# - For non-symmetric selection events: use the "centered" distribution
# - For symmetric selection events: use the "proposal" distribution (default)
if (is.null(MC_distribution)) {
  if (symmetric_selection == TRUE) {
    MC_distribution <- "proposal"
  } else if (symmetric_selection == FALSE) {
    MC_distribution <- "centered"
  }
}

# Create a list to store results for each perspective
# The list will contain elements named "marginal" and/or "conditional"
selinf_results_list <- setNames(vector("list", length(perspective)), perspective)
  
## marginal
if("marginal" %in% perspective){
  selinf_results_list[["marginal"]] <-
    lapply(marginal_COV_estimator, function(iestimator){
        lapply(which_marginal, function(iwhich){
          
          # names of coefficients
          coef_name <- paste("coefficient", iwhich)
          
    ## Input for selective_inference          
          # which vt:
          vt = testvec_marginal_list$vt_marginal_list[[iestimator]][[coef_name]]
          # which variance kappa_square
          kappa_square = testvec_marginal_list$kappa_square_marginal_list[[iestimator]][[coef_name]]
          # which dirV
          dirV = testvec_marginal_list$dirV_marginal_list[[iestimator]][[coef_name]]
          
          
          result <- selective_inference(
              n = n,
              vt = vt,
              kappa_square = kappa_square, 
              dirV = dirV,
              this_y = this_y,
              B_MCsamples = B_MCsamples,
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
        
        # names of coefficients
        coef_name <- paste("coefficient", iwhich)
        
        ## Input for selective_inference          
        # which vt:
        vt = testvec_marginal_not_efficient_list$vt_marginal_not_efficient_list[[iestimator]][[coef_name]]
        # which variance kappa_square
        kappa_square = testvec_marginal_not_efficient_list$kappa_square_marginal_not_efficient_list[[iestimator]][[coef_name]]
        # which dirV
        dirV = testvec_marginal_not_efficient_list$dirV_marginal_not_efficient_list[[iestimator]][[coef_name]]
        
        
        result <- selective_inference(
          n = n,
          vt = vt,
          kappa_square = kappa_square, 
          dirV = dirV,
          this_y = this_y,
          B_MCsamples = B_MCsamples,
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
          
           # names of coefficients
           coef_name <- paste("coefficient", iwhich)
          
    ## Input for selective_inference
           # which vt
           vt = 
           testvec_conditional_list$vt_conditional_list[[iestimator]][[coef_name]]
           # which variance kappa_square
           kappa_square = 
           testvec_conditional_list$kappa_square_conditional_list[[iestimator]][[ibayesian]][[coef_name]]
           # which dirV
           dirV = 
           testvec_conditional_list$dirV_conditional_list[[iestimator]][[ibayesian]][[coef_name]]
           
           result <- selective_inference(
             n = n, 
             vt = vt,
             kappa_square = kappa_square, 
             dirV = dirV,
             this_y = this_y,
             B_MCsamples = B_MCsamples,
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


