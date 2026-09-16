################################################################################
# Estimation functions for covariance matrices
################################################################################

#-------------------------------------------------------------------------------
# Function to compute R, R_inverse for supplied, ModelEstimate and ICM
#-------------------------------------------------------------------------------
## Note: for the VarY estimator under the conditional perspective


R_R_invers <- function(mod, # fit of the selected model
                       this_y = NULL, # vector of observed y values
                       n, # length of this_y
                       marginal_COV_estimator = c("ModelEstimate", 
                                                  "ICM",
                                                  "supplied"), # COV estimator
                       conditional_COV_estimator = c("ModelEstimate", 
                                                     "ICM",
                                                     "VarY",
                                                     "supplied"),# COV estimator
                       fit_ICM = NULL,
                       Sigma_supplied = NULL,
                       sigma_square_epsilon_supplied = NULL){
  
  # Prepare list for the return value
  COV_estimator <- unique(c(marginal_COV_estimator, conditional_COV_estimator))
  R_return_list <- R_inv_return_list <- setNames(vector("list",
                                                        length(COV_estimator)),
                                                 COV_estimator)
  
  
  
  ## SUPPLIED
  # If "supplied" is requested for marginal and no Sigma_supplied is given, or
  # if "supplied" is requested for conditional, compute R_inv_supplied from
  # the supplied sigma_square_error
  if((c("supplied") %in% marginal_COV_estimator && is.null(Sigma_supplied))|
     c("supplied") %in% conditional_COV_estimator){
    
    if(!is.null(sigma_square_epsilon_supplied)){
      
      R_supplied <- Matrix(diag(sigma_square_epsilon_supplied,n),sparse=T)
      # store in the return list
      R_return_list$supplied <- R_supplied
      
      R_supplied_invers <- Matrix(diag(1/sigma_square_epsilon_supplied,n),
                                  sparse=T) # inverse diagonal matrix                             #?# Do I also need to use Rmpfr for R_supplied_invers?
      # store in the return list
      R_inv_return_list$supplied <- R_supplied_invers
    } 
    else if(is.null(sigma_square_epsilon_supplied)){
      stop("Because COV_estimator = 'supplied', either Sigma_supplied (marginal)
            or sigma_square_epsilon_supplied (marginal, conditional) must be
            provided")
    }
  }
  
  ## MODELESTIMATE
  if(c("ModelEstimate") %in% COV_estimator){
    # Extract the estimated sigma_square_epsilon (level-1 residual variance)
    sigma_square_epsilon <- sigma(mod)^2
    
    R_est <- Matrix(diag(sigma_square_epsilon,n), sparse = T )
    # store in the return list
    R_return_list$ModelEstimate <- R_est
    
    
    R_est_invers <- Matrix(diag((1/sigma_square_epsilon),n),
                           sparse = T) # inverse diagonal matrix                            #?# Rmpfr?
    # store in the return list
    R_inv_return_list$ModelEstimate <- R_est_invers
  }
  
  ## ICM
  # Estimate the minimal intercept model, needed later to compute the
  # ICM estimator. We do this here already so it doesn't need to be
  # loaded twice (once for the marginal and once for the conditional
  # perspective)
  if(c("ICM") %in%  COV_estimator){
    sigma_square_epsilon_ICM <- getME(fit_ICM, "sigma")^2
    
    R_ICM <- Matrix(diag(sigma_square_epsilon_ICM,n), sparse = T)
    # store in the return list
    R_return_list$ICM<- R_ICM
    
    R_ICM_invers <- Matrix(diag((1/sigma_square_epsilon_ICM),n),
                           sparse = T) # inverse diagonal matrix                          #?# Rmpfr?
    # store in the return list
    R_inv_return_list$ICM <- R_ICM_invers
  }
  
  ## VarY
  ## only works for a diagonal COV matrix
  if(c("VarY") %in%  conditional_COV_estimator){
    # plug-in estimator for sigma_square_epsilon by Tibshirani et al. 2018
    sigma_square_epsilon_VarY <- as.numeric(var(this_y)*(n-1)/n)
    
    R_VarY <- Matrix(diag(sigma_square_epsilon_VarY,n), sparse = T)
    # store in the return list
    R_return_list$VarY<- R_VarY
    
    R_VarY_invers <- Matrix(diag((1/sigma_square_epsilon_VarY),n),
                            sparse = T) # inverse diagonal matrix
    # store in the return list
    R_inv_return_list$VarY <- R_VarY_invers
  }
  
  R_R_inv_return_list <- setNames(vector("list",2),c("R","R_inv"))
  R_R_inv_return_list$R <- R_return_list
  R_R_inv_return_list$R_inv <- R_inv_return_list
  return(R_R_inv_return_list)
  
}


#-------------------------------------------------------------------------------
# Function to compute R, R_inverse for a non-diagonal R
#-------------------------------------------------------------------------------


R_R_invers_nondiagonal <- function(
    marginal_COV_estimator = c("ModelEstimate", 
                               "ICM",
                               "supplied"), # COV estimator
    conditional_COV_estimator = c("ModelEstimate", 
                                  "ICM",
                                  "VarY",
                                  "supplied"),# COV estimator
    R_supplied){
  
  # Prepare list for the return value
  COV_estimator <- unique(c(marginal_COV_estimator, conditional_COV_estimator))
  R_return_list <- R_inv_return_list <- setNames(vector("list",
                                                        length(COV_estimator)),
                                                 COV_estimator)
  
  
  
  ## SUPPLIED
  # If "supplied" is requested for marginal, or
  # if "supplied" is requested for conditional, compute R and R_inverse
  # from the supplied R
  if(c("supplied") %in% marginal_COV_estimator |
     c("supplied") %in% conditional_COV_estimator){
    
    R_supplied <- as.matrix(R_supplied)
    
    # store in the return list
    R_return_list$supplied <- R_supplied
    
    R_supplied_invers <- solve(R_supplied)
    
    # store in the return list
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
# Function to compute Sigma for a non-diagonal R under the "not efficient"
# perspective
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
  
  # Prepare list for the return value
  COV_estimator <- marginal_not_efficient_COV_estimator
  SIGMA_return_list <- setNames(vector("list", length(COV_estimator)),COV_estimator)
  
  ## SUPPLIED
  if("supplied" %in% marginal_not_efficient_COV_estimator){
    
    if(!is.null(Sigma_supplied) && !is.null(R_supplied)){stop("Provide either
                                                              Sigma_supplied or R_supplied.")}
    
    if(!is.null(Sigma_supplied)) {SIGMA_return_list$supplied <- Sigma_supplied}
    
    else {SIGMA_return_list$supplied <- Z_supplied %*% G_supplied %*% t(Z_supplied)+ R_supplied}
    
  }
  
  
  ## MODELESTIMATE
  if("ModelEstimate_R_replaced" %in% marginal_not_efficient_COV_estimator |
     "ModelEstimate_SIGMA_replaced" %in% marginal_not_efficient_COV_estimator ){
    
    # Fit the fixed-effects model
    KOVMatrix_ModelEstimate <- make_clustered_errors_KOVbdiag( balanced = balanced,
                                                                mod = mod, 
                                                                HC_type = HC_type, 
                                                                cadjust = cadjust, 
                                                                clustersize = n/number_of_cluster,
                                                                number_of_cluster = number_of_cluster )
    
    if("ModelEstimate_R_replaced" %in% marginal_not_efficient_COV_estimator) SIGMA_return_list$ModelEstimate_R_replaced <- Z_B_work %*% G_work_est %*% t(Z_B_work) + KOVMatrix_ModelEstimate
    
    if("ModelEstimate_SIGMA_replaced" %in% marginal_not_efficient_COV_estimator) SIGMA_return_list$ModelEstimate_SIGMA_replaced <- KOVMatrix_ModelEstimate
    
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
    
    if("ICM_R_replaced" %in% marginal_not_efficient_COV_estimator) SIGMA_return_list$ICM_R_replaced <- Z_B_work %*% G_work_est %*% t(Z_B_work) + KOVMatrix_ICM
    
    if("ICM_SIGMA_replaced" %in% marginal_not_efficient_COV_estimator) SIGMA_return_list$ICM_SIGMA_replaced <- KOVMatrix_ICM
    
  }
  
  
  return(SIGMA_return_list)
  
}

#-------------------------------------------------------------------------------
# Function to compute Sigma for the marginal "not efficient" perspective, when
# R = sigma^2eps * In
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
  
  # Create an empty list for SIGMA, so that each list element can hold
  # one covariance-matrix estimator
  SIGMA <- setNames(vector("list", length(marginal_not_efficient_COV_estimator)),
                    marginal_not_efficient_COV_estimator)
  
  ## SUPPLIED
  ## If marginal_not_efficient_COV_estimator == "supplied", i.e. if Sigma is provided
  if(c("supplied") %in% marginal_not_efficient_COV_estimator){
    
    # If Sigma_supplied is given, write it into the SIGMA matrix
    if(!is.null(Sigma_supplied)){SIGMA$supplied <- Sigma_supplied}
    
    # If Sigma_supplied is not given but R, Z, Lambda, Lambdat_supplied
    # or G_supplied are, then compute Sigma_supplied and enter it into
    # SIGMA
    else if (is.null(Sigma_supplied) && !is.null(Z_supplied)){
      
      # If G_supplied is not given but Lambda_supplied is, use it to
      # compute the estimate of Sigma; otherwise use G_supplied
      if(!is.null(G_supplied)){G_supplied <- G_supplied}
      
      else if (is.null(G_supplied) && !is.null(Lambda_supplied)) {
        if(is.null(Lambdat_supplied)){Lambda_supplied <- t(Lambda_supplied)}
        G_supplied <- Lambda_supplied %*% Lambdat_supplied}
      
      else if(is.null(G_supplied) && is.null(Lambda_supplied)){
        stop("Because marginal_not_efficient_COV_estimator = 'supplied', either
              Lambda_supplied or G_supplied must be provided for the marginal
              not-efficient perspective")}
      
      # Woodbury matrix identity for a more efficient inverse computation
      SIGMA$supplied <- Z_supplied %*% G_supplied %*% t(Z_supplied) + R_supplied
    }
    
    # Else: since "supplied" was requested, either Sigma_supplied must be
    # non-NULL, or at least Z_supplied and G_supplied (resp. Lambda) and
    # R_supplied must be provided
    else {stop("Because marginal_not_efficient_COV_estimator = 'supplied', for the
                marginal not-efficient perspective either Sigma_supplied or
                (G_supplied resp. Lambda_supplied) and Z_supplied must be provided")}
  }
  
  ## MODELESTIMATE
  ## If marginal_not_efficient_COV_estimator == "ModelEstimate", i.e. if Sigma
  ## should be estimated via the plug-in estimator Sigma_hat, i.e. ModelEstimate
  if(c("ModelEstimate") %in% marginal_not_efficient_COV_estimator){
    # Compute the plug-in estimator for SIGMA
    
    SIGMA$ModelEstimate <- Z_B_work %*% G_work_est %*% t(Z_B_work) + R_est
    
  }
  
  ## ICM
  ## If marginal_not_efficient_COV_estimator == "ICM", i.e. if Sigma should be
  ## estimated via the conservative minMod estimator, i.e. ICM
  if(c("ICM") %in% marginal_not_efficient_COV_estimator){
    
    # Compute the ICM estimator for SIGMA
    model <- fit_ICM
    SIGMA$ICM <- Sigma_ICM(model = model,
                           n = n,
                           R_ICM = R_ICM,
                           Z = Z_B_work)
    
    
  }
  
  return(SIGMA)
  
}

#-------------------------------------------------------------------------------
# Function to compute Sigma_inverse for the marginal perspective
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
  
  
  # Create an empty list for SIGMA_invers, so that each list element can
  # hold one covariance-matrix estimator
  SIGMA_invers <- setNames(vector("list", length(marginal_COV_estimator)),
                           marginal_COV_estimator)
  
  ## SUPPLIED
  ## If marginal_COV_estimator == "supplied", i.e. if Sigma is provided
  if(c("supplied") %in% marginal_COV_estimator){
    
    # If Sigma_supplied is given, write its inverse into the SIGMA_invers matrix
    if(!is.null(Sigma_supplied)){
      SIGMA_invers$supplied <- solve(Sigma_supplied)}
    
    # If Sigma_supplied is not given but R, Z, Lambda, Lambdat_supplied
    # or G_supplied are, then compute Sigma_supplied and enter its inverse
    # into SIGMA_invers
    else if (is.null(Sigma_supplied) && !is.null(Z_supplied)){
      
      # If G_supplied is not given but Lambda_supplied is, use it to
      # compute the estimate of Sigma_invers; otherwise use G_supplied
      if(!is.null(G_supplied)){G_supplied <- G_supplied}
      
      else if (is.null(G_supplied) && !is.null(Lambda_supplied)) {
        if(is.null(Lambdat_supplied)){Lambda_supplied <- t(Lambda_supplied)}
        G_supplied <- Lambda_supplied %*% Lambdat_supplied}
      
      else if(is.null(G_supplied) && is.null(Lambda_supplied)){
        stop("Because marginal_COV_estimator = 'supplied', either Lambda_supplied
                  or G_supplied must be provided for the marginal perspective")}
      
      # Woodbury matrix identity for a more efficient inverse computation
      SIGMA_invers$supplied <- tryCatch({
        woodbury_identity(R_inv = R_supplied_invers,
                          Z = Z_supplied,
                          G = G_supplied)
      }, error = function(e) {
        solve(Z_supplied %*% G_supplied %*% t(Z_supplied) + R_supplied)
      })
      
      
    }
    
    # Else: since "supplied" was requested, either Sigma_supplied must be
    # non-NULL, or at least Z_supplied and G_supplied (resp. Lambda) and
    # R_supplied must be provided
    else {stop("Because marginal_COV_estimator = 'supplied', for the marginal
                  perspective either Sigma_supplied or (G_supplied resp.
                  Lambda_supplied) and Z_supplied must be provided")}
  }
  
  ## MODELESTIMATE
  ## If marginal_COV_estimator == "ModelEstimate", i.e. if Sigma should be
  ## estimated via the plug-in estimator Sigma_hat, i.e. ModelEstimate
  if(c("ModelEstimate") %in% marginal_COV_estimator){
    # Compute the plug-in estimator for SIGMA_invers using the Woodbury matrix identity
    
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
  ## If marginal_COV_estimator == "ICM", i.e. if Sigma should be estimated via
  ## the conservative minMod estimator, i.e. ICM
  if(c("ICM") %in% marginal_COV_estimator){
    
    # Compute the ICM estimator for SIGMA_invers
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
# Function to compute R_list, K_inv_list, Ct.R_inv_list, Ct.R_inv.C_list
# for the conditional perspective
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
  
  
  # Create an empty list for K_invers := (t(C)%*%R_inv%*%C+Q)_inv, so that
  # each list element can hold one estimator for K_inv
  K_inv_list <- setNames(vector("list",length(conditional_COV_estimator)),
                         conditional_COV_estimator)
  
  # Create empty lists for Ct.R_inv.C := t(C)%*%R_inv%*%C and
  # Ct.R_inv := t(C)%*%R_inv
  Ct.R_inv_list <- setNames(vector("list",length(conditional_COV_estimator)),
                            conditional_COV_estimator)
  Ct.R_inv.C_list <- setNames(vector("list",length(conditional_COV_estimator)),
                              conditional_COV_estimator)
  
  # Compute the block matrix C, needed to compute K_inv
  C_block <- cbind(X_A_work,Z_B_work)
  C <- C_block
  Ct <- t(C)
  
  # Compute the inverse of G_work_est (when it was taken as an estimator
  # from mod), as well as Q_est
  if(any(c("ModelEstimate","ICM","VarY") %in% conditional_COV_estimator)){
    
    # Compute the inverse of G_work_est
    G_work_est_invers <- solve(G_work_est) 
    
    # Compute Q_est
    Q_est <- bdiag(matrix(0, ncol=ncol(X_A_work), nrow=ncol(X_A_work)),
                   G_work_est_invers)
  }
  
  ## SUPPLIED
  # If conditional_COV_estimator == "supplied", use G_work_supplied
  # [i.e. the covariance matrix of the random effects that one would obtain
  # if the working model were fit to the population (i.e. when the true
  # population model actually has random slopes, but is only estimated
  # with random intercepts, for example)]
  # and R_supplied_invers to compute K_invers
  if(c("supplied") %in% conditional_COV_estimator){
    
    if(!is.null(G_work_supplied)){
      # Compute the block matrix Q: first compute G_work_invers,
      # then form the block-diagonal matrix
      G_work_supplied_invers <- solve(G_work_supplied)
      Q_supplied <- bdiag(matrix(0, ncol=ncol(X_A_work), nrow=ncol(X_A_work)),
                          G_work_supplied_invers)
      
      # Compute the list entries of Ct.R_inv_list and Ct.R_inv.C_list
      Ct.R_inv_supplied <- Ct %*% R_supplied_invers
      Ct.R_inv_list$supplied <- Ct.R_inv_supplied
      
      Ct.R_inv.C_supplied <- Ct.R_inv_supplied %*% C
      Ct.R_inv.C_list$supplied <- Ct.R_inv.C_supplied
      
      # Compute the corresponding entry of K_inv                                    #!# Woodbury identity doesn't work because Q is not invertible
      K_inv_list$supplied <- solve(Ct.R_inv.C_supplied + Q_supplied)
    }
    
    
    # Else: since "supplied" was requested, G_work_supplied must be provided.
    else {stop("Because conditional_COV_estimator = 'supplied',
                 G_work_supplied must be provided")}
    
  }
  
  ## MODELESTIMATE
  # If conditional_COV_estimator == "ModelEstimate", use Q_est
  # and R_est_invers to compute K_invers
  if(c("ModelEstimate") %in% conditional_COV_estimator){
    
    # Compute the list entries of Ct.R_inv_list and Ct.R_inv.C_list
    Ct.R_inv_ModelEstimate <- Ct %*% R_est_invers 
    Ct.R_inv_list$ModelEstimate <- Ct.R_inv_ModelEstimate
    
    Ct.R_inv.C_ModelEstimate <- Ct.R_inv_ModelEstimate %*% C
    Ct.R_inv.C_list$ModelEstimate <- Ct.R_inv.C_ModelEstimate
    
    # Compute the corresponding entry of K_inv
    K_inv_list$ModelEstimate <- solve(Ct.R_inv.C_ModelEstimate + Q_est)}
  
  ## ICM
  # If conditional_COV_estimator == "ICM", use Q_est
  # and R_ICM_invers to compute K_invers
  # !!! this should only be used when no random effects are to be
  # selected
  if(c("ICM") %in% conditional_COV_estimator){
    
    
    # Compute the list entries of Ct.R_inv_list and Ct.R_inv.C_list
    Ct.R_inv_ICM <- Ct %*% R_ICM_invers
    Ct.R_inv_list$ICM <- Ct.R_inv_ICM
    
    Ct.R_inv.C_ICM <- Ct.R_inv_ICM %*% C
    Ct.R_inv.C_list$ICM <- Ct.R_inv.C_ICM
    
    # Compute the corresponding entry of K_inv                                    #!# Woodbury identity doesn't work because Q is not invertible
    K_inv_list$ICM <- solve(Ct.R_inv.C_ICM + Q_est)}
  
  ## VarY
  # If conditional_COV_estimator == "VarY", use Q_est
  # and R_VarY_invers to compute K_invers
  if(c("VarY") %in% conditional_COV_estimator){
    
    # Compute the list entries of Ct.R_inv_list and Ct.R_inv.C_list
    #!# Woodbury identity doesn't work because Q is not invertible
    Ct.R_inv_VarY <- Ct %*% R_VarY_invers
    Ct.R_inv_list$VarY <- Ct.R_inv_VarY
    
    Ct.R_inv.C_VarY <- Ct.R_inv_VarY %*% C
    Ct.R_inv.C_list$VarY <- Ct.R_inv.C_VarY
    
    # Compute the corresponding entry of K_inv
    K_inv_list$VarY <- solve(Ct.R_inv.C_VarY + Q_est)}
  
  ## Return objects
  return_K_inv.Ct.R_inv_list <- setNames(vector("list",4),
                                         c("K_inv_list","Ct.R_inv_list","Ct.R_inv.C_list","C_block"))
  
  return_K_inv.Ct.R_inv_list$K_inv_list <- K_inv_list
  return_K_inv.Ct.R_inv_list$Ct.R_inv_list <- Ct.R_inv_list
  return_K_inv.Ct.R_inv_list$Ct.R_inv.C_list <- Ct.R_inv.C_list
  return_K_inv.Ct.R_inv_list$C_block <- C_block
  
  
  return(return_K_inv.Ct.R_inv_list)
}
