################################################################################
# Test-vector functions
################################################################################


#-------------------------------------------------------------------------------
# marginal test-vector function
#-------------------------------------------------------------------------------

testvec_marginal <- function(
    marginal_COV_estimator,
    which_marginal,
    X_A_work,
    SIGMA_invers,
    robust_std.errors_marginal,
    mod,
    this_y,
    splitRatio = 1
){
  
  # Create lists so the parameters can be computed per marginal_COV_estimator
  
  vt_marginal_list <- kappa_square_marginal_list <- dirV_marginal_list <-
    setNames(vector("list",length(marginal_COV_estimator)),
             marginal_COV_estimator) 
  
  # Since it's possible that "which" specifies more than one coefficient,
  # we also create sub-lists within the lists, so that each sub-list
  # refers to one coefficient
  
  # Create sub-lists
  for(iestimator in marginal_COV_estimator) { 
    vt_marginal_list[[iestimator]] <- 
      kappa_square_marginal_list[[iestimator]] <- 
      dirV_marginal_list[[iestimator]] <-
      setNames(vector("list", length(which_marginal)),
               paste("coefficient", which_marginal))
  }
  
  
  ## Test vector vt and variance kappa_square
  # Note: we only compute the efficient test statistic
  # vt_eff := e_j' * (X_A' * Sigma_inv * X_A)_inv' * X_A' * Sigma_inv'
  X_A_workt <- t(X_A_work)
  
  for(iestimator in marginal_COV_estimator) {
    # the matrix product (X_A' * Sigma_inv * X_A)_inv can be reused for
    # both the test-vector definition and the variance
    X_At.Sigma_inv.X_A_inv <- solve( 
      X_A_workt %*% 
        SIGMA_invers[[iestimator]] %*%
        X_A_work)
    X_A.X_At.Sigma_inv.X_A_inv <- X_A_work %*% X_At.Sigma_inv.X_A_inv 
    
    # jointly for all coefficients in which_marginal:
    # determine the test vector vt 
    vt_marginal_all_coeffs <- t(X_At.Sigma_inv.X_A_inv) %*%
      X_A_workt %*% t(SIGMA_invers[[iestimator]])
    
    # robust covariance matrix for all coefficients
    if(length(robust_std.errors_marginal) != 1){ 
      stop("robust_std.errors_marginal may only contain one of the
           following: FALSE, 'CR0' or 'CR2'.")}
    
    if(robust_std.errors_marginal != FALSE){
      
      CR_kappa_square_all_coeffs <-
        CR_std.errors(mod = mod,
                      this_y = this_y,
                      X = X_A_work,
                      Sigma = solve(SIGMA_invers[[iestimator]]),                                             
                      bread = X_At.Sigma_inv.X_A_inv,
                      vt_marginal_all_coeffs = vt_marginal_all_coeffs,
                      robust_std.errors_marginal = robust_std.errors_marginal)
      
    }
    
    # write the entry for each coefficient into the corresponding list element
    for (iwhich in which_marginal){
      # coefficient name for the list entry
      coef_name <- paste("coefficient", iwhich)
      
      # test vectors vt
      vt_marginal_list[[iestimator]][[coef_name]] <- 
        splitRatio * matrix(vt_marginal_all_coeffs[iwhich,],nrow = 1)
      
      # variance kappa_square
      # not robust: k^2 := e_j'(X_A'*Sigma_inv*X_A)_inv e_j
      if(robust_std.errors_marginal == FALSE){
        
        kappa_square_marginal_list[[iestimator]][[coef_name]] <- 
          splitRatio^2 *X_At.Sigma_inv.X_A_inv[iwhich,iwhich] 
        
      }
      #robust
      if(robust_std.errors_marginal != FALSE){
        
        kappa_square_marginal_list[[iestimator]][[coef_name]] <- 
          splitRatio^2 * CR_kappa_square_all_coeffs[iwhich,iwhich] 
        
      }
      
      # dirV := Sigma * v_eff/(v_eff' * Sigma * v_eff) (Rügamer 3.1.2)
      #      = Sigma * v_eff/kappa_square_marginal 
      #      = (X_A * X_At.Sigma_inv.X_A_inv * e_j) /kappa_square_marginal
      dirV_marginal_list[[iestimator]][[coef_name]] <- 
        splitRatio * matrix(X_A.X_At.Sigma_inv.X_A_inv[,iwhich]/
                 (splitRatio^2 *X_At.Sigma_inv.X_A_inv[iwhich,iwhich]), ncol = 1)
    }
  }
  
  result_list <- setNames(vector("list",3),c("vt_marginal_list",
                                             "kappa_square_marginal_list",
                                             "dirV_marginal_list"))
  result_list$vt_marginal_list <- vt_marginal_list
  result_list$kappa_square_marginal_list <- kappa_square_marginal_list
  result_list$dirV_marginal_list <- dirV_marginal_list
  
  return(result_list)
  
}

#-------------------------------------------------------------------------------
# marginal "not efficient" test-vector function
#-------------------------------------------------------------------------------

testvec_marginal_not_efficient <- function(
    AdjMcCaffrey=TRUE,
    mod = NULL,
    this_y = NULL,
    robust_std.errors_marginal,
    marginal_not_efficient_COV_estimator,
    which_marginal,
    X_A_work,
    SIGMA,
    splitRatio = 1 ){
  
  # Create lists so the parameters can be computed per
  # marginal_not_efficient_COV_estimator
  
  vt_marginal_not_efficient_list <- kappa_square_marginal_not_efficient_list <- dirV_marginal_not_efficient_list <-
    setNames(vector("list",length(marginal_not_efficient_COV_estimator)),
             marginal_not_efficient_COV_estimator) 
  
  # Since it's possible that "which" specifies more than one coefficient,
  # we also create sub-lists within the lists, so that each sub-list
  # refers to one coefficient
  
  # Create sub-lists
  for(iestimator in marginal_not_efficient_COV_estimator) { 
    vt_marginal_not_efficient_list[[iestimator]] <- 
      kappa_square_marginal_not_efficient_list[[iestimator]] <- 
      dirV_marginal_not_efficient_list[[iestimator]] <-
      setNames(vector("list", length(which_marginal)),
               paste("coefficient", which_marginal))
  }
  
  
  ## Test vector vt and variance kappa_square
  # Note: here we compute the "not efficient" test statistic
  # vt := e_j' * (X_A'*X_A)_inv' * X_A' 
  X_A_workt <- t(X_A_work)
  X_At.X_A_inv <- solve(X_A_workt %*% X_A_work)
  
  # jointly for all coefficients in which_marginal:
  # determine the test vector vt 
  vt_marginal_not_efficient_all_coeffs <- t(X_At.X_A_inv) %*% X_A_workt           
  
  for(iestimator in marginal_not_efficient_COV_estimator) {
    # the matrix product (X_A'*Sigma_inv*X_A)_inv can be reused for
    # both the test-vector definition and the variance
    X_At.Sigma.X_A <- X_A_workt %*% SIGMA[[iestimator]] %*% X_A_work
    kappasquare_marginal_not_efficient_all_coeffs <- X_At.X_A_inv %*% X_At.Sigma.X_A %*% X_At.X_A_inv
    
    # robust covariance matrix for all coefficients
    if(length(robust_std.errors_marginal) != 1){ 
      stop("robust_std.errors_marginal may only contain one of the
           following: FALSE, 'CR0' or 'CR2'.")}
    
    if(robust_std.errors_marginal != FALSE){
      
      CR_kappa_square_all_coeffs <-
        CR_std.errors_OLS(AdjMcCaffrey=AdjMcCaffrey,
                          mod = mod,
                          this_y = this_y,
                          X = X_A_work,
                          Sigma = SIGMA[[iestimator]],
                          bread =   X_At.X_A_inv,
                          vt_marginal_all_coeffs = vt_marginal_not_efficient_all_coeffs,
                          robust_std.errors_marginal = robust_std.errors_marginal)
      
    }
    
    # write the entry for each coefficient into the corresponding list element
    for (iwhich in which_marginal){
      # coefficient name for the list entry
      coef_name <- paste("coefficient", iwhich)
      
      # test vectors vt
      vt_marginal_not_efficient_list[[iestimator]][[coef_name]] <- 
        splitRatio * matrix(vt_marginal_not_efficient_all_coeffs[iwhich,],nrow = 1)
      
      # variance kappa_square
      # not robust:
      if(robust_std.errors_marginal == FALSE){
        
        kappa_square_marginal_not_efficient_list[[iestimator]][[coef_name]] <- 
          splitRatio^2 *kappasquare_marginal_not_efficient_all_coeffs[iwhich,iwhich] 
        
      }
      #robust
      if(robust_std.errors_marginal != FALSE){
        
        kappa_square_marginal_not_efficient_list[[iestimator]][[coef_name]] <- 
          splitRatio^2 * CR_kappa_square_all_coeffs[iwhich,iwhich] 
        
      }
      
      
      # dirV := Sigma * v/kappa_square_marginal 
      
      dirV_marginal_not_efficient_list[[iestimator]][[coef_name]] <- 
        SIGMA[[iestimator]]%*%t(vt_marginal_not_efficient_list[[iestimator]][[coef_name]])/
        (splitRatio^2 *kappasquare_marginal_not_efficient_all_coeffs[iwhich,iwhich])
      
    }
  }
  
  result_list <- setNames(vector("list",3),c("vt_marginal_not_efficient_list",
                                             "kappa_square_marginal_not_efficient_list",
                                             "dirV_marginal_not_efficient_list"))
  result_list$vt_marginal_not_efficient_list <- vt_marginal_not_efficient_list
  result_list$kappa_square_marginal_not_efficient_list <- kappa_square_marginal_not_efficient_list
  result_list$dirV_marginal_not_efficient_list <- dirV_marginal_not_efficient_list
  
  return(result_list)
  
}

#-------------------------------------------------------------------------------
# conditional test-vector function
#-------------------------------------------------------------------------------


testvec_conditional <- function(
    conditional_COV_estimator,
    bayesian,
    which_conditional,
    C,
    K_inv_list,
    Ct.R_inv_list,
    Ct.R_inv.C_list,
    R_list,
    splitRatio = 1
    
){
  # Create lists so the parameters can be computed per conditional_COV_estimator
  vt_conditional_list <- kappa_square_conditional_list <- 
    dirV_conditional_list <-
    setNames(vector("list",length(conditional_COV_estimator)),
             conditional_COV_estimator) 
  
  # For the variance kappa_square, create a sub-list that distinguishes
  # between bayesian and not_bayesian
  for(iestimator in conditional_COV_estimator) { 
    kappa_square_conditional_list[[iestimator]] <- 
      dirV_conditional_list[[iestimator]] <-                                    
      setNames(vector("list",length(bayesian)),paste(bayesian))
  }
  
  # Since it's possible that "which" specifies more than one coefficient,
  # we also create sub-lists within the lists, so that each sub-list
  # refers to one coefficient
  # Create sub-lists
  for(iestimator in conditional_COV_estimator) { 
    vt_conditional_list[[iestimator]] <- 
      setNames(vector("list",length(which_conditional)),
               paste("coefficient", which_conditional))
    
    for(ibayesian in bayesian){
      kappa_square_conditional_list[[iestimator]][[ibayesian]] <-
        dirV_conditional_list[[iestimator]][[ibayesian]] <-
        setNames(vector("list",length(which_conditional)),
                 paste("coefficient", which_conditional))}
  }
  
  
  ## Test vector vt, variance kappa_square, and dirV
  for(iestimator in conditional_COV_estimator) {
    
    # R
    R <- R_list[[iestimator]]
    
    # jointly for all coefficients in which_conditional:
    # determine the test matrix V for the test vector vt 
    V_conditional <- K_inv_list[[iestimator]] %*% 
      Ct.R_inv_list[[iestimator]]
    
    
    # determine the variance kappa_square: K_inv*C'*R_inv*C*K_inv
    # for not_bayesian
    kappa_square_conditional_all_coeffs <- K_inv_list[[iestimator]] %*% 
      Ct.R_inv.C_list[[iestimator]] %*%
      K_inv_list[[iestimator]]
    
    # determine the variance kappa_square: K_inv 
    # for bayesian 
    bayesian_kappa_square_conditional_all_coeffs <- K_inv_list[[iestimator]]
    
    
    
    
    # write the entry for each coefficient into the corresponding list element
    for (iwhich in which_conditional){
      # coefficient name for the list entry
      coef_name <- paste("coefficient", iwhich)
      
      # test vectors vt
      vt_conditional_list[[iestimator]][[coef_name]] <-
        splitRatio * matrix(V_conditional[iwhich,], nrow = 1)
      
      # variances kappa_square and dirV 
      if("not_bayesian" %in% bayesian){
        # variance kappa_square := e_j'*K_inv*C'*R_inv*C*K_inv*e_j
        kappa_square_conditional_list[[iestimator]][["not_bayesian"]][[coef_name]] <- 
          splitRatio^2 *kappa_square_conditional_all_coeffs[iwhich,iwhich]                                    #?# this is a bit inefficient since it always computes the full matrix V rather than only the rows relevant for "which" 
        
        #dirV := R*v/kappa_square, where v = t(vt)
        dirV_conditional_list[[iestimator]][["not_bayesian"]][[coef_name]] <- 
          R %*% t(vt_conditional_list[[iestimator]][[coef_name]])/
          kappa_square_conditional_list[[iestimator]][["not_bayesian"]][[coef_name]]} 
      
      
      if("bayesian" %in% bayesian){
        # variance kappa_square := e_j'*K_inv*e_j
        kappa_square_conditional_list[[iestimator]][["bayesian"]][[coef_name]] <- 
          splitRatio^2 *bayesian_kappa_square_conditional_all_coeffs[iwhich,iwhich]                                  
        
        #dirV := R*v/kappa_square, where v = t(vt)
        dirV_conditional_list[[iestimator]][["bayesian"]][[coef_name]] <- 
          R %*% t(vt_conditional_list[[iestimator]][[coef_name]])/
          kappa_square_conditional_list[[iestimator]][["bayesian"]][[coef_name]]}
      
      if("bayesian_corrected" %in% bayesian){
        # variance kappa_square := e_j'*K_inv*e_j
        kappa_square_conditional_list[[iestimator]][["bayesian_corrected"]][[coef_name]] <- 
          splitRatio^2 *bayesian_kappa_square_conditional_all_coeffs[iwhich,iwhich]                                    
        #dirV := R*v/kappa_square_freq, where v = t(vt)
        dirV_conditional_list[[iestimator]][["bayesian_corrected"]][[coef_name]] <- 
          R %*% t(vt_conditional_list[[iestimator]][[coef_name]])/
          (splitRatio^2 *kappa_square_conditional_all_coeffs[iwhich,iwhich])  }
      
    }
  }
  result_list <- setNames(vector("list",3),c("vt_conditional_list",
                                             "kappa_square_conditional_list",
                                             "dirV_conditional_list"))
  result_list$vt_conditional_list <- vt_conditional_list
  result_list$kappa_square_conditional_list <- kappa_square_conditional_list
  result_list$dirV_conditional_list <- dirV_conditional_list
  
  return(result_list)
  
}
