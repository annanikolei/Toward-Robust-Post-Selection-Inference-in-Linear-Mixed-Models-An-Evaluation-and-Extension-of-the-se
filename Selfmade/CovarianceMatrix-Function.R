################################################################################
# Covariance matrix estimation functions
################################################################################

#-------------------------------------------------------------------------------
# Function to compute R and R^{-1} for supplied, ModelEstimate, and ICM
#-------------------------------------------------------------------------------
## Note: also used for the VarY estimator in the conditional perspective

R_R_invers <- function(mod, # fitted model
                       this_y = NULL, # vector of observed y values
                       n, # length of this_y
                       marginal_COV_estimator = c("ModelEstimate", 
                                                  "ICM",
                                                  "supplied"), # covariance estimator
                       conditional_COV_estimator = c("ModelEstimate", 
                                                     "ICM",
                                                     "VarY",
                                                     "supplied"), # covariance estimator
                       fit_ICM = NULL,
                       Sigma_supplied = NULL,
                       sigma_square_epsilon_supplied = NULL){
  
  # Prepare list for return
  COV_estimator <- unique(c(marginal_COV_estimator, conditional_COV_estimator))
  R_return_list <- R_inv_return_list <- setNames(vector("list",
                                                        length(COV_estimator)),
                                                 COV_estimator)
  
  ## SUPPLIED
  # If "supplied" is used for marginal and Sigma_supplied is missing OR
  # if "supplied" is used for conditional, compute R from supplied sigma^2
  if((c("supplied") %in% marginal_COV_estimator && is.null(Sigma_supplied))|
     c("supplied") %in% conditional_COV_estimator){
    
    if(!is.null(sigma_square_epsilon_supplied)){
      
      R_supplied <- Matrix(diag(sigma_square_epsilon_supplied,n),sparse=T)
      # store in list
      R_return_list$supplied <- R_supplied
      
      R_supplied_invers <- Matrix(diag(1/sigma_square_epsilon_supplied,n),
                                  sparse=T) # inverse diagonal matrix
      # store in list
      R_inv_return_list$supplied <- R_supplied_invers
    } 
    else if(is.null(sigma_square_epsilon_supplied)){
      stop("If COV_estimator = 'supplied', either Sigma_supplied (marginal) 
            or sigma_square_epsilon_supplied (marginal, conditional) 
            must be provided")
    }
  }
  
  ## MODELESTIMATE
  if(c("ModelEstimate") %in% COV_estimator){
    # Extract estimated sigma^2_epsilon (level-1 residual variance)
    sigma_square_epsilon <- sigma(mod)^2
    
    R_est <- Matrix(diag(sigma_square_epsilon,n), sparse = T )
    # store in list
    R_return_list$ModelEstimate <- R_est
    
    R_est_invers <- Matrix(diag((1/sigma_square_epsilon),n),
                           sparse = T) # inverse diagonal matrix
    # store in list
    R_inv_return_list$ModelEstimate <- R_est_invers
  }
  
  ## ICM
  # Estimate minimal intercept model for later ICM estimator computation
  if(c("ICM") %in%  COV_estimator){
    sigma_square_epsilon_ICM <- getME(fit_ICM, "sigma")^2
    
    R_ICM <- Matrix(diag(sigma_square_epsilon_ICM,n), sparse = T)
    # store in list
    R_return_list$ICM<- R_ICM
    
    R_ICM_invers <- Matrix(diag((1/sigma_square_epsilon_ICM),n),
                           sparse = T) # inverse diagonal matrix
    # store in list
    R_inv_return_list$ICM <- R_ICM_invers
  }
  
  ## VarY
  ## only works for diagonal covariance matrix
  if(c("VarY") %in%  conditional_COV_estimator){
    # plug-in estimator for sigma^2_epsilon (Tibshirani et al. 2018)
    sigma_square_epsilon_VarY <- as.numeric(var(this_y)*(n-1)/n)
    
    R_VarY <- Matrix(diag(sigma_square_epsilon_VarY,n), sparse = T)
    # store in list
    R_return_list$VarY<- R_VarY
    
    R_VarY_invers <- Matrix(diag((1/sigma_square_epsilon_VarY),n),
                            sparse = T) # inverse diagonal matrix
    # store in list
    R_inv_return_list$VarY <- R_VarY_invers
  }
  
  R_R_inv_return_list <- setNames(vector("list",2),c("R","R_inv"))
  R_R_inv_return_list$R <- R_return_list
  R_R_inv_return_list$R_inv <- R_inv_return_list
  return(R_R_inv_return_list)
}


#-------------------------------------------------------------------------------
# Function to compute R and R^{-1} for non-diagonal R
#-------------------------------------------------------------------------------

R_R_invers_nondiagonal <- function(
    marginal_COV_estimator = c("ModelEstimate", 
                               "ICM",
                               "supplied"), # covariance estimator
    conditional_COV_estimator = c("ModelEstimate", 
                                  "ICM",
                                  "VarY",
                                  "supplied"),# covariance estimator
    R_supplied){
  
  # Prepare list for return
  COV_estimator <- unique(c(marginal_COV_estimator, conditional_COV_estimator))
  R_return_list <- R_inv_return_list <- setNames(vector("list",
                                                        length(COV_estimator)),
                                                 COV_estimator)
  
  ## SUPPLIED
  # If "supplied" is used (marginal or conditional),
  # compute R and R^{-1} from the supplied matrix
  if(c("supplied") %in% marginal_COV_estimator |
     c("supplied") %in% conditional_COV_estimator){
    
    R_supplied <- as.matrix(R_supplied)
    
    # store in list
    R_return_list$supplied <- R_supplied
    
    R_supplied_invers <- solve(R_supplied)
    
    # store in list
    R_inv_return_list$supplied <- R_supplied_invers
    
  }
  
  ## MODELESTIMATE
  
  ## ICM
  
  R_R_inv_return_list <- setNames(vector("list",2),c("R","R_inv"))
  R_R_inv_return_list$R <- R_return_list
  R_R_inv_return_list$R_inv <- R_inv_return_list
  return(R_R_inv_return_list)
}


#-------------------------------------------------------------------------------
# Function to compute Sigma for non-diagonal R in the not efficient case
#-------------------------------------------------------------------------------

SIGMA_function_diagF <- function(
    mod,
    marginal_not_efficient_COV_estimator,
    G_work_est = NULL,
    Z_B_work = NULL,
    R_supplied = NULL,
    Z_supplied = NULL,
    G_supplied = NULL,
    Sigma_supplied = NULL,
    balanced = T,
    HC_type, 
    cadjust, 
    n,
    number_of_cluster){
  
  # Prepare list for return
  COV_estimator <- marginal_not_efficient_COV_estimator
  SIGMA_return_list <- setNames(vector("list", length(COV_estimator)),COV_estimator)
  
  ## SUPPLIED
  if("supplied" %in% marginal_not_efficient_COV_estimator){
    
    if(!is.null(Sigma_supplied) && !is.null(R_supplied)){
      stop("Either Sigma_supplied or R_supplied must be provided, not both.")
    }
    
    if(!is.null(Sigma_supplied)) {
      SIGMA_return_list$supplied <- Sigma_supplied
    }
    
    else {
      SIGMA_return_list$supplied <- Z_supplied %*% G_supplied %*% t(Z_supplied)+ R_supplied
    }
    
  }
  
  ## MODELESTIMATE
  if("ModelEstimate_R_replaced" %in% marginal_not_efficient_COV_estimator |
     "ModelEstimate_SIGMA_replaced" %in% marginal_not_efficient_COV_estimator ){
    
    # Estimate covariance matrix via clustered residuals
    KOVMatrix_ModelEstimate <- make_clustered_errors_KOVbdiag( balanced = balanced,
                                                                mod = mod, 
                                                                HC_type = HC_type, 
                                                                cadjust = cadjust, 
                                                                clustersize = n/number_of_cluster,
                                                                number_of_cluster = number_of_cluster )
    
    if("ModelEstimate_R_replaced" %in% marginal_not_efficient_COV_estimator) {
      SIGMA_return_list$ModelEstimate_R_replaced <- Z_B_work %*% G_work_est %*% t(Z_B_work) + KOVMatrix_ModelEstimate
    }
    
    if("ModelEstimate_SIGMA_replaced" %in% marginal_not_efficient_COV_estimator) {
      SIGMA_return_list$ModelEstimate_SIGMA_replaced <- KOVMatrix_ModelEstimate
    }
    
  }
  
  ## ICM
  if("ICM_R_replaced" %in% marginal_not_efficient_COV_estimator|
     "ICM_SIGMA_replaced" %in% marginal_not_efficient_COV_estimator){
    
    KOVMatrix_ICM <- make_clustered_errors_KOVbdiag( balanced = balanced,
                                                      mod = minMod(mod), 
                                                      HC_type = HC_type, 
                                                      cadjust = cadjust, 
                                                      clustersize = n/number_of_cluster,
                                                      number_of_cluster = number_of_cluster )
    
    if("ICM_R_replaced" %in% marginal_not_efficient_COV_estimator) {
      SIGMA_return_list$ICM_R_replaced <- Z_B_work %*% G_work_est %*% t(Z_B_work) + KOVMatrix_ICM
    }
    
    if("ICM_SIGMA_replaced" %in% marginal_not_efficient_COV_estimator) {
      SIGMA_return_list$ICM_SIGMA_replaced <- KOVMatrix_ICM
    }
    
  }
  
  return(SIGMA_return_list)
}

#-------------------------------------------------------------------------------
# Function to compute Sigma for the marginal not efficient perspective when
# R = sigma^2_epsilon * I_n (diagonal residual covariance matrix)
#-------------------------------------------------------------------------------

SIGMA_function_diagT <- function( 
    marginal_not_efficient_COV_estimator =  c("ModelEstimate", "ICM","supplied"),
    Sigma_supplied = NULL,
    Z_supplied = NULL,
    G_supplied = NULL,
    Lambda_supplied = NULL,
    Lambdat_supplied = NULL,
    R_supplied = NULL,
    R_est = NULL,
    Z_B_work = NULL,
    G_work_est = NULL,
    fit_ICM = NULL,
    R_ICM = NULL,
    n = NULL){ 
  
  # Create empty list for SIGMA so that each element can store
  # one covariance matrix estimator
  SIGMA <- setNames(vector("list", length(marginal_not_efficient_COV_estimator)),
                    marginal_not_efficient_COV_estimator)
  
  ## SUPPLIED
  ## Case: marginal_not_efficient_COV_estimator == "supplied",
  ## i.e., Sigma is provided
  if(c("supplied") %in% marginal_not_efficient_COV_estimator){
    
    # If Sigma_supplied is provided, assign it directly
    if(!is.null(Sigma_supplied)){SIGMA$supplied <- Sigma_supplied}
    
    # If Sigma_supplied is not provided but R, Z, Lambda/Lambdat
    # or G_supplied are provided, compute Sigma and store it
    else if (is.null(Sigma_supplied) && !is.null(Z_supplied)){
      
      # If G_supplied is not provided but Lambda_supplied is,
      # compute G from Lambda; otherwise use G_supplied
      if(!is.null(G_supplied)){G_supplied <- G_supplied}
      
      else if (is.null(G_supplied) && !is.null(Lambda_supplied)) {
        if(is.null(Lambdat_supplied)){Lambda_supplied <- t(Lambda_supplied)}
        G_supplied <- Lambda_supplied %*% Lambdat_supplied}
      
      else if(is.null(G_supplied) && is.null(Lambda_supplied)){
        stop("If marginal_not_efficient_COV_estimator = 'supplied', either 
              Lambda_supplied or G_supplied must be provided for the marginal 
              not efficient perspective")}
      
      # Compute Sigma = Z G Z' + R
      SIGMA$supplied <- Z_supplied %*% G_supplied %*% t(Z_supplied) + R_supplied
    }
    
    # Otherwise: required inputs for 'supplied' are missing
    else {stop("If marginal_not_efficient_COV_estimator = 'supplied', either 
                Sigma_supplied or (G_supplied / Lambda_supplied) and Z_supplied 
                must be provided")}
  }
  
  ## MODELESTIMATE
  ## Case: Sigma is estimated via plug-in estimator (model-based)
  if(c("ModelEstimate") %in% marginal_not_efficient_COV_estimator){
    
    # Compute plug-in estimator for Sigma
    SIGMA$ModelEstimate <- Z_B_work %*% G_work_est %*% t(Z_B_work) + R_est
    
  }
  
  ## ICM
  ## Case: Sigma is estimated via conservative minMod (ICM estimator)
  if(c("ICM") %in% marginal_not_efficient_COV_estimator){
    
    # Compute ICM estimator for Sigma
    model <- fit_ICM
    SIGMA$ICM <- Sigma_ICM(model = model,
                           n = n,
                           R_ICM = R_ICM,
                           Z = Z_B_work)
    
  }
  
  return(SIGMA)
  
}

#-------------------------------------------------------------------------------
# Function to compute Sigma^{-1} for the marginal perspective
#-------------------------------------------------------------------------------

SIGMA_invers_function <- function(
    marginal_COV_estimator =  c("ModelEstimate", "ICM","supplied"),
    Sigma_supplied = NULL,
    Z_supplied = NULL,
    G_supplied = NULL,
    Lambda_supplied = NULL,
    Lambdat_supplied = NULL,
    R_supplied = NULL,
    R_supplied_invers = NULL,
    R_est = NULL,
    R_est_invers = NULL,
    Z_B_work = NULL,
    G_work_est = NULL,
    fit_ICM = NULL,
    R_ICM = NULL,
    R_ICM_invers = NULL,
    n = NULL){ 
  
  
  # Create empty list for SIGMA_invers so that each element can store
  # one covariance matrix estimator
  SIGMA_invers <- setNames(vector("list", length(marginal_COV_estimator)),
                           marginal_COV_estimator)
  
  ## SUPPLIED
  ## Case: marginal_COV_estimator == "supplied", i.e., Sigma is provided
  if(c("supplied") %in% marginal_COV_estimator){
    
    # If Sigma_supplied is provided, compute and store its inverse
    if(!is.null(Sigma_supplied)){
      SIGMA_invers$supplied <- solve(Sigma_supplied)}
    
    # If Sigma_supplied is not provided but R, Z, Lambda/Lambdat
    # or G_supplied are provided, compute Sigma and its inverse
    else if (is.null(Sigma_supplied) && !is.null(Z_supplied)){
      
      # If G_supplied is not provided but Lambda_supplied is,
      # compute G from Lambda; otherwise use G_supplied
      if(!is.null(G_supplied)){G_supplied <- G_supplied}
      
      else if (is.null(G_supplied) && !is.null(Lambda_supplied)) {
        if(is.null(Lambdat_supplied)){Lambda_supplied <- t(Lambda_supplied)}
        G_supplied <- Lambda_supplied %*% Lambdat_supplied}
      
      else if(is.null(G_supplied) && is.null(Lambda_supplied)){
        stop("If marginal_COV_estimator = 'supplied', either Lambda_supplied 
              or G_supplied must be provided for the marginal perspective")}
      
      # Use Woodbury matrix identity for efficient inversion
      SIGMA_invers$supplied <- tryCatch({
        woodbury_identity(R_inv = R_supplied_invers,
                          Z = Z_supplied,
                          G = G_supplied)
      }, error = function(e) {
        solve(Z_supplied %*% G_supplied %*% t(Z_supplied) + R_supplied)
      })
      
      
    }
    
    # Otherwise: required inputs for 'supplied' are missing
    else {stop("If marginal_COV_estimator = 'supplied', either 
                  Sigma_supplied or (G_supplied / Lambda_supplied)
                  and Z_supplied must be provided")}
  }
  
  ## MODELESTIMATE
  ## Case: Sigma is estimated via plug-in estimator (model-based)
  if(c("ModelEstimate") %in% marginal_COV_estimator){
    # Compute plug-in estimator for SIGMA^{-1} using Woodbury identity
    
    SIGMA_invers$ModelEstimate <- tryCatch({
      woodbury_identity(R_inv = R_est_invers,
                        Z = Z_B_work,
                        G = G_work_est)
    }, error = function(e) {
      solve(Z_B_work %*% G_work_est %*% t(Z_B_work) + R_est)
    })
    # = solve( Z_B %*% G_work %*% t(Z_B) + R_est) 
  }
  
  ## ICM
  ## Case: Sigma is estimated via conservative minMod estimator (ICM)
  if(c("ICM") %in% marginal_COV_estimator){
    
    # Compute ICM estimator for SIGMA^{-1}
    model <- fit_ICM
    SIGMA_invers$ICM <- Sigma_invers_ICM(model = model,
                                         n = n,
                                         R_ICM = R_ICM,
                                         R_ICM_invers = R_ICM_invers,
                                         Z = Z_B_work)
    
    
  }
  return(SIGMA_invers)
}


#-------------------------------------------------------------------------------
# Function to compute K^{-1}, C'R^{-1}, and C'R^{-1}C for the conditional perspective
#-------------------------------------------------------------------------------

K_inv.Ct.R_inv_list_function <- function(
    conditional_COV_estimator = c("ModelEstimate","ICM","VarY","supplied"),
    X_A_work,
    Z_B_work,
    G_work_est,
    R_supplied_invers,
    G_work_supplied,
    R_est_invers,
    R_ICM_invers,
    R_VarY_invers
){
  
  # Initialize lists
  K_inv_list <- setNames(vector("list",length(conditional_COV_estimator)),
                         conditional_COV_estimator)
  
  Ct.R_inv_list <- setNames(vector("list",length(conditional_COV_estimator)),
                            conditional_COV_estimator)
  
  Ct.R_inv.C_list <- setNames(vector("list",length(conditional_COV_estimator)),
                              conditional_COV_estimator)
  
  # Construct block matrix C = [X | Z]
  C_block <- cbind(X_A_work,Z_B_work)
  C <- C_block
  Ct <- t(C)
  
  # Compute inverse of G_work_est and corresponding Q matrix
  if(any(c("ModelEstimate","ICM","VarY") %in% conditional_COV_estimator)){
    
    G_work_est_invers <- solve(G_work_est) 
    
    Q_est <- bdiag(matrix(0, ncol=ncol(X_A_work), nrow=ncol(X_A_work)),
                   G_work_est_invers)
  }
  
  ## SUPPLIED
  if(c("supplied") %in% conditional_COV_estimator){
    
    if(!is.null(G_work_supplied)){
      
      G_work_supplied_invers <- solve(G_work_supplied)
      
      Q_supplied <- bdiag(matrix(0, ncol=ncol(X_A_work), nrow=ncol(X_A_work)),
                          G_work_supplied_invers)
      
      Ct.R_inv_supplied <- Ct %*% R_supplied_invers
      Ct.R_inv_list$supplied <- Ct.R_inv_supplied
      
      Ct.R_inv.C_supplied <- Ct.R_inv_supplied %*% C
      Ct.R_inv.C_list$supplied <- Ct.R_inv.C_supplied
      
      # Woodbury not applicable because Q is not invertible
      K_inv_list$supplied <- solve(Ct.R_inv.C_supplied + Q_supplied)
    }
    
    else {
      stop("If conditional_COV_estimator = 'supplied', G_work_supplied must be provided.")
    }
  }
  
  ## MODELESTIMATE
  if(c("ModelEstimate") %in% conditional_COV_estimator){
    
    Ct.R_inv_ModelEstimate <- Ct %*% R_est_invers 
    Ct.R_inv_list$ModelEstimate <- Ct.R_inv_ModelEstimate
    
    Ct.R_inv.C_ModelEstimate <- Ct.R_inv_ModelEstimate %*% C
    Ct.R_inv.C_list$ModelEstimate <- Ct.R_inv.C_ModelEstimate
    
    K_inv_list$ModelEstimate <- solve(Ct.R_inv.C_ModelEstimate + Q_est)
  }
  
  ## ICM
  if(c("ICM") %in% conditional_COV_estimator){
    
    Ct.R_inv_ICM <- Ct %*% R_ICM_invers
    Ct.R_inv_list$ICM <- Ct.R_inv_ICM
    
    Ct.R_inv.C_ICM <- Ct.R_inv_ICM %*% C
    Ct.R_inv.C_list$ICM <- Ct.R_inv.C_ICM
    
    K_inv_list$ICM <- solve(Ct.R_inv.C_ICM + Q_est)
  }
  
  ## VarY
  if(c("VarY") %in% conditional_COV_estimator){
    
    Ct.R_inv_VarY <- Ct %*% R_VarY_invers
    Ct.R_inv_list$VarY <- Ct.R_inv_VarY
    
    Ct.R_inv.C_VarY <- Ct.R_inv_VarY %*% C
    Ct.R_inv.C_list$VarY <- Ct.R_inv.C_VarY
    
    K_inv_list$VarY <- solve(Ct.R_inv.C_VarY + Q_est)
  }
  
  # Return object
  return_K_inv.Ct.R_inv_list <- setNames(vector("list",4),
                                         c("K_inv_list","Ct.R_inv_list","Ct.R_inv.C_list","C_block"))
  
  return_K_inv.Ct.R_inv_list$K_inv_list <- K_inv_list
  return_K_inv.Ct.R_inv_list$Ct.R_inv_list <- Ct.R_inv_list
  return_K_inv.Ct.R_inv_list$Ct.R_inv.C_list <- Ct.R_inv.C_list
  return_K_inv.Ct.R_inv_list$C_block <- C_block
  
  return(return_K_inv.Ct.R_inv_list)
}

