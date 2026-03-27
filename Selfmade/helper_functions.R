################################################################################
# Helper functions
################################################################################

# Copied from the cAIC4 package
cnms2formula <-
  function(cnms) {
    # Builds a random effects formula from the component names (cnms),
    # which can be extracted from an lmerMod object via @cnms or getME(., "cnms").
    #
    # Args:
    #   cnms: List from an lmerMod object (@cnms or getME(., "cnms"))
    #
    # Returns:
    #   reFormula: Random effects part of an lmer formula
    
    len      <- unlist(lapply(cnms, length))
    cnms     <- cnms[which(len != 0)]
    charForm <- character(length(cnms))
    
    for (i in 1:length(cnms)) {
      if (cnms[[i]][1] == "(Intercept)") {
        cnms[[i]][1] <- "1"
      } else {
        tpv <- cnms[[i]]
        cnms[[i]] <- append("", tpv)
        cnms[[i]][1] <- "-1"
      }
      charForm[i] <- paste("(", paste(cnms[[i]], collapse = " + "),
                           " | ", names(cnms)[i], ")", sep = "")
    }
    
    reFormula <- paste(charForm, collapse = " + ")
    return(reFormula)
  }

#-------------------------------------------------------------------------------
## Efficient matrix inversion (Woodbury identity)
#-------------------------------------------------------------------------------

# Woodbury identity for efficient inversion:
# (R + Z G Z')^{-1} = R^{-1} - R^{-1} Z (G^{-1} + Z' R^{-1} Z)^{-1} Z' R^{-1}
#
# Useful when:
# - R^{-1} is easy to compute (e.g., R is diagonal)
# - G has lower dimension than R
# Example: R = diag(sigma^2) (n x n), G = random effects covariance matrix

woodbury_identity <- function(R_inv, Z, G) {
  # Transpose of Z
  Zt <- t(Z)
  
  # Precompute repeated matrix product for efficiency
  R_inv.Z <- R_inv %*% Z  
  
  # Compute inverse using Woodbury identity
  Inverse <- R_inv -
    R_inv.Z %*% solve(
      solve(G) +
        Zt %*% R_inv.Z) %*%
    Zt %*% R_inv
  
  return(Inverse)
}

#-------------------------------------------------------------------------------
## ICM estimator
#-------------------------------------------------------------------------------

# Fit intercept-only model with fixed random effects structure
minMod <- function(mod) {
  
  lhs <- formula(mod)[[2]] # left-hand side (response)
  reFormula <- cnms2formula(mod@cnms) # random effects part
  
  # Fit model with fixed effects reduced to intercept only
  update(mod, formula = reformulate(c("1", reFormula), lhs)) 
}

# For marginal model:
# Estimate Sigma^{-1} using the intercept-only model with fixed random effects
Sigma_invers_ICM <- function(model,
                             n = NROW(model@resp$y),
                             sigma_epsilon_ICM = NULL,
                             R_ICM = NULL,
                             R_ICM_invers = NULL,
                             Z = getME(model,"Z")) {
  
  # Extract Lambda and its transpose
  Lambda_ICM  <- getME(model,"Lambda")
  Lambdat_ICM <- getME(model,"Lambdat")
  
  # Estimate residual variance if not provided
  if (is.null(sigma_epsilon_ICM)) {
    sigma_epsilon_ICM <- getME(model, "sigma")
  }
  
  # Compute R^{-1} if not provided
  if (is.null(R_ICM_invers)) {
    R_ICM_invers <- diag(1 / sigma_epsilon_ICM^2, n)
  }
  
  # Compute R if not provided
  if (is.null(R_ICM)) {
    R_ICM <- diag(sigma_epsilon_ICM^2, n)
  }
  
  # Compute random effects covariance matrix G
  G <- sigma_epsilon_ICM^2 * Lambda_ICM %*% Lambdat_ICM
  
  # Try Woodbury identity, fallback to direct inversion if needed
  Inverse_ICM <- tryCatch({
    woodbury_identity(R_inv = R_ICM_invers,
                      Z = Z,
                      G = G)
  }, error = function(e) {
    solve(Z %*% G %*% t(Z) + R_ICM)
  })
  
  return(Inverse_ICM)
}

# Compute Sigma (not inverse) for ICM estimator
Sigma_ICM <- function(model,
                      n = NROW(model@resp$y),
                      sigma_epsilon_ICM = NULL,
                      R_ICM = NULL,
                      R_ICM_invers = NULL,
                      Z = getME(model,"Z")) {
  
  # Extract Lambda and its transpose
  Lambda_ICM  <- getME(model,"Lambda")
  Lambdat_ICM <- getME(model,"Lambdat")
  
  # Estimate residual variance if not provided
  if (is.null(sigma_epsilon_ICM)) {
    sigma_epsilon_ICM <- getME(model, "sigma")
  }
  
  # Compute R if not provided
  if (is.null(R_ICM)) {
    R_ICM <- diag(sigma_epsilon_ICM^2, n)
  }
  
  # Compute G
  G <- sigma_epsilon_ICM^2 * Lambda_ICM %*% Lambdat_ICM
  
  # Compute Sigma = Z G Z' + R
  SigmaMatrix_ICM <- Z %*% G %*% t(Z) + R_ICM
  
  return(SigmaMatrix_ICM)
}

#-------------------------------------------------------------------------------
## Computation of C_eta
#-------------------------------------------------------------------------------

Ceta <- function(wbs,
                 wbs_geq_tobs,
                 rho,
                 kappa_square,
                 Tbs_survived,
                 Tbs_surv_geq_tobs) {
  
  sum(
    exp(Rmpfr::mpfr(Tbs_surv_geq_tobs * rho / kappa_square, 200)) * wbs_geq_tobs
  ) /
  sum(
    exp(Rmpfr::mpfr(Tbs_survived * rho / kappa_square, 200)) * wbs
  ) 
}

#-------------------------------------------------------------------------------
## Construct block-diagonal covariance matrix from residuals
#-------------------------------------------------------------------------------

make_clustered_errors_KOVbdiag <- function(balanced = TRUE,
                                            mod, 
                                            HC_type = "HC0", # HC0, HC1, HC2, HC3
                                            cadjust = FALSE, # apply cluster adjustment?
                                            clustersize, 
                                            number_of_cluster) {
  
  # Only implemented for balanced designs
  if (!balanced) {
    stop("transform_residuals_to_KOVbdiag currently only supports balanced designs.")
  }
  
  # Number of observations
  n_obs <- number_of_cluster * clustersize
  
  # Design matrix
  X <- getME(mod,"X")
  
  # Fit linear model (fixed effects only) and extract residuals
  lm_mod  <- lm(formula(mod, fixed = TRUE), data = mod@frame)
  resid_lm_mod <- resid(lm_mod)
  
  # Number of predictors (excluding intercept)
  p_no_intercept <- sum(names(coef(lm_mod)) != "(Intercept)")
  
  # Residual cross-product depending on HC type
  if (HC_type %in% c("HC0")) {
    ResidualCrossprod <- resid_lm_mod %*% t(resid_lm_mod)
  }
  
  if (HC_type %in% c("HC1")) {
    ResidualCrossprod <- (n_obs - 1) / (n_obs - p_no_intercept - 1) *
      resid_lm_mod %*% t(resid_lm_mod)
  }
  
  # Hat matrix for HC2/HC3
  if (HC_type %in% c("HC2","HC3")) {
    hat_Matrix <- X %*% solve(t(X) %*% X) %*% t(X)
  }
  
  if (HC_type %in% c("HC2")) {
    ResidualCrossprod <- 1 / (1 - hat_Matrix) *
      resid_lm_mod %*% t(resid_lm_mod)
  }
  
  if (HC_type %in% c("HC3")) {
    ResidualCrossprod <- 1 / (1 - hat_Matrix)^2 *
      resid_lm_mod %*% t(resid_lm_mod)
  }
  
  # Create list of cluster-specific blocks
  Res_blocks <- lapply(seq_len(number_of_cluster), function(icluster) {
    rows <- ((icluster - 1) * clustersize + 1):(icluster * clustersize)
    ResidualCrossprod[rows, rows]
  })
  
  # Cluster adjustment (Bell & McCaffrey / Cameron & Miller)
  if ((!cadjust & HC_type != "HC3") | (cadjust & HC_type == "HC3")) {
    cluster_adjustment <- 1
  }
  if (!cadjust & HC_type == "HC3") {
    cluster_adjustment <- (number_of_cluster - 1) / number_of_cluster
  }
  if (cadjust & HC_type != "HC3") {
    cluster_adjustment <- number_of_cluster / (number_of_cluster - 1)
  }
  
  # Construct block-diagonal covariance matrix
  KOVbdiag <- cluster_adjustment * bdiag(Res_blocks)
  
  return(KOVbdiag)
}
