################################################################################
# Helper functions
################################################################################

# copied from cAIC4 package
cnms2formula <-
  function(cnms) {
    # A function that builds a random effects formula from the ?component names?,
    # a list that can be extracted from an lmerMod object by .@cnms or
    # getME(., "cnms").
    #
    # Args:
    #   cnms     = List from an lmerMod object by .@cnms or getME(., "cnms").
    #
    # Returns:
    #   reFormula = random effects part of a lmerMod formula
    #
    len      <- unlist(lapply(cnms, length))
    cnms     <- cnms[which(len != 0)]
    charForm <- character(length(cnms))
    
    for(i in 1:length(cnms)) {
      if (cnms[[i]][1] == "(Intercept)") {
        cnms[[i]][1] <- "1"
      } else {
        tpv <- cnms[[i]]
        cnms[[i]] <- append("",tpv)
        cnms[[i]][1] <- "-1"
      }
      charForm[i] <- paste("(", paste(cnms[[i]], collapse = " + "),
                           " | ",names(cnms)[i], ")", sep = "")
    }
    
    reFormula <- paste(charForm, collapse = " + ")
    
    return(reFormula)
  }

#-------------------------------------------------------------------------------
## More efficient matrix inversion
#-------------------------------------------------------------------------------

# Woodbury identity: computes the inverse of certain matrices more efficiently:
# (R+ZGZ')^{-1}=R^{-1}-R^{-1}Z(G^{-1}+Z'R^{-1}Z)^{-1}Z'R^{-1}
# https://de.wikipedia.org/wiki/Woodbury-Matrix-Identit%C3%A4t
# For example, this is more efficient when R_inv can be computed easily
# (as is the case when R is diagonal) and G has lower dimensionality than R.
# This is the case, for example, when R = diag(sigma_square_epsilon)
# (n_obs x n_obs) and G is the random effects covariance matrix when there is
# only a random intercept (m = nl2 x m = nl2).

woodbury_identity <- function(R_inv, Z, G) {
  # Transpose the Z matrix
  Zt <- t(Z)
  # Compute matrix products that occur repeatedly in advance for efficiency
  R_inv.Z <- R_inv %*% Z  
  
  # Compute the inverse
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

# Estimates the (minMod) intercept model with a fixed random effects structure
minMod <- function(mod)
{
  
  lhs <- formula(mod)[[2]] #left hand side of formula
  reFormula <- cnms2formula(mod@cnms) # random effects part of formula
  #modIC <-
  #fit model with fixed effects part set to zero
  update(mod, formula = reformulate(c("1", reFormula), lhs)) 
  # return(list(sigma = getME(modIC, "sigma"),
  #             tau = getME(modIC, "theta"))
  # )
  
}

# Only relevant for the marginal perspective:
# Here, we implement what the vcov_RI function from the selfmade package
# appears to be intended to do, namely estimating Sigma^{-1} in the intercept
# model with the fixed random effects structure
Sigma_invers_ICM <- function(model, # typically the minMod(mod) model 
                             n = NROW(model@resp$y),
                             sigma_epsilon_ICM = NULL,
                             R_ICM = NULL,
                             R_ICM_invers = NULL,
                             Z = getME(model,"Z") # by default with Z_B, because
                             # when using ICM,
                             # the random effects are fixed,
                             # and therefore should not
                             # have been selected.
){
  # Compute Lambda_ICM and Lambdat_ICM
  Lambda_ICM = getME(model,"Lambda")
  Lambdat_ICM = getME(model,"Lambdat")
  
  # If sigma_epsilon is not provided, compute it
  if(is.null(sigma_epsilon_ICM)){sigma_epsilon_ICM <- getME(model, "sigma")}
  
  # Compute R_ICM_invers
  if(is.null(R_ICM_invers)){
    
    R_ICM_invers <- diag(1/sigma_epsilon_ICM^2, n)}     
  
  # Compute R
  if(is.null(R_ICM)){
    R_ICM <- diag(sigma_epsilon_ICM^2, n)}
  
  # Compute the inverse of the Sigma matrix using the more efficient Woodbury
  # identity
  # Compute G according to Equation (4) in
  # https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  # because Lambda was obtained using the getME function
  
  G = sigma_epsilon_ICM^2 * Lambda_ICM %*% Lambdat_ICM 
  
  Inverse_ICM <- tryCatch({
    woodbury_identity(R_inv = R_ICM_invers,
                      Z = Z,
                      G = G)
  }, error = function(e) {
    
    solve(Z %*% G %*% t(Z) + R_ICM)
    
  })
  
  return(Inverse_ICM)                                                            
  
}


Sigma_ICM <- function(model, # typically the minMod(mod) model 
                             n = NROW(model@resp$y),
                             sigma_epsilon_ICM = NULL,
                             R_ICM = NULL,
                             R_ICM_invers = NULL,
                             Z = getME(model,"Z") # by default with Z_B, because
                             # when using ICM,
                             # the random effects are fixed,
                             # and therefore should not
                             # have been selected.
){
  # Compute Lambda_ICM and Lambdat_ICM
  Lambda_ICM = getME(model,"Lambda")
  Lambdat_ICM = getME(model,"Lambdat")
  
  # If sigma_epsilon is not provided, compute it
  if(is.null(sigma_epsilon_ICM)){sigma_epsilon_ICM <- getME(model, "sigma")}
  
  # Compute R
  if(is.null(R_ICM)) {R_ICM <- diag(sigma_epsilon_ICM^2, n)}
  
  # Compute the inverse of the Sigma matrix using the more efficient Woodbury
  # identity
  # Compute G according to Equation (4) in
  # https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  # because Lambda was obtained using the getME function
  
  G = sigma_epsilon_ICM^2 * Lambda_ICM %*% Lambdat_ICM 
  
  SigmaMatrix_ICM <- Z %*% G %*% t(Z) + R_ICM
  
  return(SigmaMatrix_ICM)                                                             
  
}

#-------------------------------------------------------------------------------
## Computation of Ceta
#-------------------------------------------------------------------------------
Ceta <- function(wbs,
                 wbs_geq_tobs,
                 rho,
                 kappa_square,
                 Tbs_survived,
                 Tbs_surv_geq_tobs){
  sum(
    exp(Rmpfr::mpfr(Tbs_surv_geq_tobs * rho / (kappa_square),200)) * wbs_geq_tobs
  )/
    sum(
      exp(Rmpfr::mpfr(Tbs_survived * rho / (kappa_square),200)) * wbs
    ) 
}


#-------------------------------------------------------------------------------
## Construct block-diagonal matrix of residual variances from a given residual vector
#-------------------------------------------------------------------------------

make_clustered_errors_KOVbdiag <- function(balanced = T,
                                            mod, 
                                            HC_type = "HC0", # or "HC1","HC2","HC3"
                                            cadjust = F, # logical. Should a cluster adjustment be applied?
                                            clustersize, 
                                            number_of_cluster){
  
  if(!balanced){ stop("transform_residuals_to_KOVbdiag currently only supports balanced designs.")}
  
  # Number of observations
  n_obs <- number_of_cluster * clustersize
  X <- getME(mod,"X")
  

  # Fit an lm model using the fixed-effects formula and extract the residual vector
  lm_mod  <- lm(formula(mod, fixed = TRUE), data = mod@frame)
  resid_lm_mod <- resid(lm_mod)
  
  # Number of predictors
  p_no_intercept <- sum(names(coef(lm_mod)) != "(Intercept)")
  
  # 
  # Calculate residual cross-product
  if(HC_type %in% c("HC0")){ ResidualCrossprod <- resid_lm_mod %*% t(resid_lm_mod) }
  if(HC_type %in% c("HC1")){ ResidualCrossprod <- (n_obs - 1)/(n_obs - p_no_intercept - 1) * resid_lm_mod %*% t(resid_lm_mod) }
   
  # Hat values 
  if(HC_type %in% c("HC2","HC3")){ hat_Matrix <-  X %*% solve(t(X) %*% X) %*% t(X) }
  
  if(HC_type %in% c("HC2")){ 
    # Element-wise matrix multiplication
    ResidualCrossprod <- 1/(1 - hat_Matrix) * resid_lm_mod %*% t(resid_lm_mod) }
  
  if(HC_type %in% c("HC3")){
    # Element-wise matrix multiplication
    ResidualCrossprod <- 1/(1 - hat_Matrix)^2 * resid_lm_mod %*% t(resid_lm_mod) }
  
  
  # Create a list of submatrices
  Res_blocks <- lapply(seq_len(number_of_cluster), function(icluster) {
      
    rows <- ((icluster - 1) * clustersize + 1):(icluster * clustersize)
    
    ResidualCrossprod[rows, rows]
  })
  
  
  # Assemble the block-diagonal matrix
  ## "inverse" cluster adjustment that Bell & McCaffrey (2002) and hence also
  ## Cameron & Miller (2005, Eq. 25) recommend for HC3 (but not HC2)
  ## -> canceled out again if cadjust = TRUE
  if( (!cadjust & HC_type!= "HC3") | (cadjust & HC_type == "HC3") ){ cluster_adjustment <- 1 }
  if(!cadjust & HC_type == "HC3" ){ cluster_adjustment <- (number_of_cluster - 1)/number_of_cluster }
  if(cadjust  & HC_type!= "HC3"){ cluster_adjustment <-  number_of_cluster/(number_of_cluster - 1) }
  
  KOVbdiag = cluster_adjustment * as.matrix(bdiag(Res_blocks))
  
  return(KOVbdiag)
  
}









