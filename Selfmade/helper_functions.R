################################################################################
# Hilfsfunktionen 
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
## Effizientere Inversenberechnung ##
#-------------------------------------------------------------------------------

# Woodbury Identität: berechnet Inverse von bestimmten Matrizen effizienter:
# (R+ZGZ')^{-1}=R^{-1}-R^{-1}Z(G^{-1}+Z'R^{-1}Z)^{-1}Z'R^{-1}
# https://de.wikipedia.org/wiki/Woodbury-Matrix-Identit%C3%A4t
# z.B. ist es effizienter, wenn R_invers sich leicht berechnen lässt (was der Fall
# ist, wenn R diagonal) und G ne geringere Dimension hat als R. Das wäre alles
# z.B. dann der Fall, wenn R=diag(sigma_square_epsilon) (n_obs x n_obs) und G die
#  random effects Kovarianzmatrix wenn nur random intercept (m=nl2 x m=nl2).

woodbury_identity <- function(R_inv, Z, G) {
  # Z-Matrix transponieren
  Zt <- t(Z)
  # Matrixprodukte die mehrmals vorkommen, werden vorher berechnet (Effizienz)
  R_inv.Z <- R_inv %*% Z  
  
  # Berechnung der Inverse
  Inverse <- R_inv -
    R_inv.Z %*% solve(
      solve(G) +
        Zt %*% R_inv.Z) %*%
    Zt %*% R_inv
  return(Inverse)
}

#-------------------------------------------------------------------------------
## ICM-estimator ##
#-------------------------------------------------------------------------------

#Schätzt das (minMod) intercept model mit fixed random effects structure 
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

# nur fürs marginale Modell interessant: 
# Wir machen hier das, was die vcov_RI-Funktion aus dem selfmade package
# aus meiner Sicht eigentlich beabsichtigt hat, nämlich die Schätzung von dem 
# Sigma invers in dem Intercept-Modell mit der fixed random effects structure
Sigma_invers_ICM <- function(model, # typischerweise das minMod(mod)-Modell 
                             n = NROW(model@resp$y),
                             sigma_epsilon_ICM = NULL,
                             R_ICM = NULL,
                             R_ICM_invers = NULL,
                             Z = getME(model,"Z") # per default mit Z_B, denn
                             # wenn man ICM benutzt, 
                             # sind random effects fest,
                             # sollten also nicht
                             # selektiert worden sein.
){
  # Berechne Lambda_ICM und Lambdat_ICM
  Lambda_ICM = getME(model,"Lambda")
  Lambdat_ICM = getME(model,"Lambdat")
  
  # falls kein sigma_epsilon gegeben, berechne es
  if(is.null(sigma_epsilon_ICM)){sigma_epsilon_ICM <- getME(model, "sigma")}
  
  # zur Berechnung von R_ICM_invers
  if(is.null(R_ICM_invers)){
    
    R_ICM_invers <- diag(1/sigma_epsilon_ICM^2, n)}     
  
  # Berechnung R
  if(is.null(R_ICM)){
    R_ICM <- diag(sigma_epsilon_ICM^2, n)}
  
  # Berechne inverse Sigma Matrix mit effizienterer Woodbury-Identität
  # Berechnung von G wieder nach Formel (4) aus
  # https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  # weil Lambda mit getMe-Befehl gewählt wurde 
  
  G = sigma_epsilon_ICM^2 * Lambda_ICM %*% Lambdat_ICM 
  
  Inverse_ICM <- tryCatch({
    woodbury_identity(R_inv = R_ICM_invers,
                      Z = Z,
                      G = G)
  }, error = function(e) {
    
    solve(Z %*% G %*% t(Z) + R_ICM)
    
  })
  
  return(Inverse_ICM)                                                              #?# ich hoffe, dass das so die Intention war bei dem ICM-estimator, also dass auch da Sigma_invers berechnet werden soll von dem minMod
  
}


Sigma_ICM <- function(model, # typischerweise das minMod(mod)-Modell 
                             n = NROW(model@resp$y),
                             sigma_epsilon_ICM = NULL,
                             R_ICM = NULL,
                             R_ICM_invers = NULL,
                             Z = getME(model,"Z") # per default mit Z_B, denn
                             # wenn man ICM benutzt, 
                             # sind random effects fest,
                             # sollten also nicht
                             # selektiert worden sein.
){
  # Berechne Lambda_ICM und Lambdat_ICM
  Lambda_ICM = getME(model,"Lambda")
  Lambdat_ICM = getME(model,"Lambdat")
  
  # falls kein sigma_epsilon gegeben, berechne es
  if(is.null(sigma_epsilon_ICM)){sigma_epsilon_ICM <- getME(model, "sigma")}
  
  # Berechnung R
  if(is.null(R_ICM)) {R_ICM <- diag(sigma_epsilon_ICM^2, n)}
  
  # Berechne inverse Sigma Matrix mit effizienterer Woodbury-Identität
  # Berechnung von G wieder nach Formel (4) aus
  # https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
  # weil Lambda mit getMe-Befehl gewählt wurde 
  
  G = sigma_epsilon_ICM^2 * Lambda_ICM %*% Lambdat_ICM 
  
  SigmaMatrix_ICM <- Z %*% G %*% t(Z) + R_ICM
  
  return(SigmaMatrix_ICM)                                                              #?# ich hoffe, dass das so die Intention war bei dem ICM-estimator, also dass auch da Sigma_invers berechnet werden soll von dem minMod
  
}

#-------------------------------------------------------------------------------
## Berechnung von Ceta ##
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
## Blockdiagonalmatrix der Residualvarianzen aus angegebenem Residualvektor erstellen
#-------------------------------------------------------------------------------

make_clustered_errors_KOVbdiag <- function(balanced = T,
                                            mod, 
                                            HC_type = "HC0", # or "HC1","HC2","HC3"
                                            cadjust = F, # logical. Should a cluster adjustment be applied?
                                            clustersize, 
                                            number_of_cluster){
  
  if(!balanced){ stop("transform_residuals_to_KOVbdiag momentan nur für balancierte Designs.")}
  
  # number of observations
  n_obs <- number_of_cluster * clustersize
  X <- getME(mod,"X")
  

  # fit lm with fixed effects formula and extract residual vector
  lm_mod  <- lm(formula(mod, fixed = TRUE), data = mod@frame)
  resid_lm_mod <- resid(lm_mod)
  
  # number of predictors 
  p_no_intercept <- sum(names(coef(lm_mod)) != "(Intercept)")
  
  # 
  # calculate residual crossproduct
  if(HC_type %in% c("HC0")){ ResidualCrossprod <- resid_lm_mod %*% t(resid_lm_mod) }
  if(HC_type %in% c("HC1")){ ResidualCrossprod <- (n_obs - 1)/(n_obs - p_no_intercept - 1) * resid_lm_mod %*% t(resid_lm_mod) }
   
  # Hat values 
  if(HC_type %in% c("HC2","HC3")){ hat_Matrix <-  X %*% solve(t(X) %*% X) %*% t(X) }
  
  if(HC_type %in% c("HC2")){ 
    # elementweise Matrixmultiplikation
    ResidualCrossprod <- 1/(1 - hat_Matrix) * resid_lm_mod %*% t(resid_lm_mod) }
  
  if(HC_type %in% c("HC3")){
    # elementweise Matrixmultiplikation
    ResidualCrossprod <- 1/(1 - hat_Matrix)^2 * resid_lm_mod %*% t(resid_lm_mod) }
  
  
  # Liste mit Submatrizen erzeugen
  Res_blocks <- lapply(seq_len(number_of_cluster), function(icluster) {
      
    rows <- ((icluster - 1) * clustersize + 1):(icluster * clustersize)
    
    ResidualCrossprod[rows, rows]
  })
  
  
  # Blockdiagonalmatrix zusammensetzen
  ## "inverse" cluster adjustment that Bell & McCaffrey (2002) and hence also
  ## Cameron & Miller (2005, Eq. 25) recommend for HC3 (but not HC2)
  ## -> canceled out again if cadjust = TRUE
  if( (!cadjust & HC_type!= "HC3") | (cadjust & HC_type == "HC3") ){ cluster_adjustment <- 1 }
  if(!cadjust & HC_type == "HC3" ){ cluster_adjustment <- (number_of_cluster - 1)/number_of_cluster }
  if(cadjust  & HC_type!= "HC3"){ cluster_adjustment <-  number_of_cluster/(number_of_cluster - 1) }
  
  KOVbdiag = cluster_adjustment * as.matrix(bdiag(Res_blocks))
  
  return(KOVbdiag)
  
}










