################################################################################
# Testvektor-Funktionen 
################################################################################


#-------------------------------------------------------------------------------
# marginale Testvektor-Funktion
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
  
  # Erstelle Listen, um pro marginal_COV_estimator die Parameter berechnen zu 
  # können
  
  vt_marginal_list <- kappa_square_marginal_list <- dirV_marginal_list <-
    setNames(vector("list",length(marginal_COV_estimator)),
             marginal_COV_estimator) 
  
  # Da es möglich ist, dass in which mehr als ein Koeffizient angegeben wird
  # erstellen wir noch Unterlisten in den Listen, sodass sich jede auf einen
  # Koeffizienten bezieht
  
  # Erstelle Unterlisten
  for(iestimator in marginal_COV_estimator) { 
    vt_marginal_list[[iestimator]] <- 
      kappa_square_marginal_list[[iestimator]] <- 
      dirV_marginal_list[[iestimator]] <-
      setNames(vector("list", length(which_marginal)),
               paste("coefficient", which_marginal))
  }
  
  
  ## Testvektor vt und Varianz kappa_square
  # Bemerkung: Wir berechnen nur die efficient Teststatistik
  # vt_eff:= e_j'*(X_A'*Sigma_inv*X_A)_inv' * X_A' * Sigma_inv'
  X_A_workt <- t(X_A_work)
  
  for(iestimator in marginal_COV_estimator) {
    # das Matrixprodukt (X_A'*Sigma_inv*X_A)_inv kann man sowohl für
    # die Testvektordefinition als auch für die Varianz gebrauchen
    X_At.Sigma_inv.X_A_inv <- solve( 
      X_A_workt %*% 
        SIGMA_invers[[iestimator]] %*%
        X_A_work)
    X_A.X_At.Sigma_inv.X_A_inv <- X_A_work %*% X_At.Sigma_inv.X_A_inv 
    
    # für alle Koeffizienten in which_marginal gemeinsam:
    # Testvektor vt bestimmen 
    vt_marginal_all_coeffs <- t(X_At.Sigma_inv.X_A_inv) %*%
      X_A_workt %*% t(SIGMA_invers[[iestimator]])
    
    # robust covariance matrix for all coefficients
    if(length(robust_std.errors_marginal) != 1){ 
      stop("robust_std.errors_marginal darf nur eines der folgenden
           Elemente enthalten: FALSE, 'CR0' oder 'CR2'.")}
    
    if(robust_std.errors_marginal != FALSE){
      
      CR_kappa_square_all_coeffs <-
        CR_std.errors(mod = mod,
                      this_y = this_y,
                      X = X_A_work,
                      Sigma = solve(SIGMA_invers[[iestimator]]),                                                #?# ja das könnte man noch effizienter machen
                      bread = X_At.Sigma_inv.X_A_inv,
                      vt_marginal_all_coeffs = vt_marginal_all_coeffs,
                      robust_std.errors_marginal = robust_std.errors_marginal)
      
    }
    
    # Eintrag pro Koeffizient in das entsprechende Listenelement
    for (iwhich in which_marginal){
      # Name des Koeffizienten zur Eintragung in Liste
      coef_name <- paste("coefficient", iwhich)
      
      # Testvektoren vt
      vt_marginal_list[[iestimator]][[coef_name]] <- 
        splitRatio * matrix(vt_marginal_all_coeffs[iwhich,],nrow = 1)
      
      # Varianz kappa_square
      # not robust: k^2:= e_j'(X_A'*Sigma_inv*X_A)_inv e_j
      if(robust_std.errors_marginal == FALSE){
        
        kappa_square_marginal_list[[iestimator]][[coef_name]] <- 
          splitRatio^2 *X_At.Sigma_inv.X_A_inv[iwhich,iwhich] 
        
      }
      #robust
      if(robust_std.errors_marginal != FALSE){
        
        kappa_square_marginal_list[[iestimator]][[coef_name]] <- 
          splitRatio^2 * CR_kappa_square_all_coeffs[iwhich,iwhich] 
        
      }
      
      # dirV:= Sigma * v_eff/(v_eff' * Sigma * v_eff) (Rügamer 3.1.2)
      #      = Sigma * v_eff/kappa_square_marginal 
      #      = (X_A * X_At.Sigma_inv.X_A_inv * e_j) /kappa_square_marginal
      dirV_marginal_list[[iestimator]][[coef_name]] <- 
        splitRatio * matrix(X_A.X_At.Sigma_inv.X_A_inv[,iwhich]/
                 kappa_square_marginal_list[[iestimator]][[coef_name]], ncol = 1)
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
# marginale not efficient Testvektor-Funktion
#-------------------------------------------------------------------------------

testvec_marginal_not_efficient <- function(
    robust_std.errors_marginal = F,
    marginal_not_efficient_COV_estimator,
    which_marginal,
    X_A_work,
    SIGMA,
    splitRatio = 1 ){
  
  # Erstelle Listen, um pro marginal_not_efficient_COV_estimator die Parameter berechnen zu 
  # können
  
  vt_marginal_not_efficient_list <- kappa_square_marginal_not_efficient_list <- dirV_marginal_not_efficient_list <-
    setNames(vector("list",length(marginal_not_efficient_COV_estimator)),
             marginal_not_efficient_COV_estimator) 
  
  # Da es möglich ist, dass in which mehr als ein Koeffizient angegeben wird
  # erstellen wir noch Unterlisten in den Listen, sodass sich jede auf einen
  # Koeffizienten bezieht
  
  # Erstelle Unterlisten
  for(iestimator in marginal_not_efficient_COV_estimator) { 
    vt_marginal_not_efficient_list[[iestimator]] <- 
      kappa_square_marginal_not_efficient_list[[iestimator]] <- 
      dirV_marginal_not_efficient_list[[iestimator]] <-
      setNames(vector("list", length(which_marginal)),
               paste("coefficient", which_marginal))
  }
  
  
  ## Testvektor vt und Varianz kappa_square
  # Bemerkung: Wir berechnen  die not efficient Teststatistik
  # vt:= e_j'*(X_A'*X_A)_inv' * X_A' 
  X_A_workt <- t(X_A_work)
  X_At.X_A_inv <- solve(X_A_workt %*% X_A_work)
  
  # für alle Koeffizienten in which_marginal gemeinsam:
  # Testvektor vt bestimmen 
  vt_marginal_not_efficient_all_coeffs <- t(X_At.X_A_inv) %*% X_A_workt            #!# das t() könnte man weglassen
  
  for(iestimator in marginal_not_efficient_COV_estimator) {
    # das Matrixprodukt (X_A'*Sigma_inv*X_A)_inv kann man sowohl für
    # die Testvektordefinition als auch für die Varianz gebrauchen
    X_At.Sigma.X_A <- X_A_workt %*% SIGMA[[iestimator]] %*% X_A_work
    kappasquare_marginal_not_efficient_all_coeffs <- X_At.X_A_inv %*% X_At.Sigma.X_A %*% X_At.X_A_inv
    
    # robust covariance matrix for all coefficients
    if(length(robust_std.errors_marginal) != 1){ 
      stop("robust_std.errors_marginal nicht möglich für marginal_not_efficient.")}
    
    if(robust_std.errors_marginal != FALSE){warning("robust_std.errors_marginal nicht 
                                                   möglich für marginal_not_efficient. Normale Standardfehler werden berichtet.")}
    
    
    # Eintrag pro Koeffizient in das entsprechende Listenelement
    for (iwhich in which_marginal){
      # Name des Koeffizienten zur Eintragung in Liste
      coef_name <- paste("coefficient", iwhich)
      
      # Testvektoren vt
      vt_marginal_not_efficient_list[[iestimator]][[coef_name]] <- 
        splitRatio * matrix(vt_marginal_not_efficient_all_coeffs[iwhich,],nrow = 1)
      
      # Varianz kappa_square
      # not efficient:
      kappa_square_marginal_not_efficient_list[[iestimator]][[coef_name]] <- 
        splitRatio^2 *kappasquare_marginal_not_efficient_all_coeffs[iwhich,iwhich] 
      
      # dirV:= Sigma * v/kappa_square_marginal 
      
      dirV_marginal_not_efficient_list[[iestimator]][[coef_name]] <- 
        SIGMA[[iestimator]]%*%t(vt_marginal_not_efficient_list[[iestimator]][[coef_name]])/
        kappa_square_marginal_not_efficient_list[[iestimator]][[coef_name]]
      
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
# konditionale Testvektor-Funktion
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
  # Erstelle Listen, um pro conditional_COV_estimator die Parameter berechnen  
  # zu können
  vt_conditional_list <- kappa_square_conditional_list <- 
    dirV_conditional_list <-
    setNames(vector("list",length(conditional_COV_estimator)),
             conditional_COV_estimator) 
  
  # Erstelle für die Varianz kappa_square eine Unterliste, die nach bayesian
  # und not_bayesian unterscheidet
  for(iestimator in conditional_COV_estimator) { 
    kappa_square_conditional_list[[iestimator]] <- 
      dirV_conditional_list[[iestimator]] <-                                        #?# hier bin ich unser bei dem dirV im paper unterscheiden die nicht nach bayesian dafür aber im package schon, deswegen unterscheide ich jetzt auch mal
      setNames(vector("list",length(bayesian)),paste(bayesian))
  }
  
  # Da es möglich ist, dass in which mehr als ein Koeffizient angegeben wird
  # erstellen wir noch Unterlisten in den Listen, sodass sich jede auf einen
  # Koeffizienten bezieht
  # Erstelle Unterlisten
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
  
  
  ## Testvektor vt und Varianz kappa_square und dirV
  for(iestimator in conditional_COV_estimator) {
    
    # R
    R <- R_list[[iestimator]]
    
    # für alle Koeffizienten in which_conditional gemeinsam:
    # Testmatrix V für Testvektor vt bestimmen 
    V_conditional <- K_inv_list[[iestimator]] %*% 
      Ct.R_inv_list[[iestimator]]
    
    
    # Varianz kappa_square bestimmen: K_inv*C'*R_inv*C*K_inv
    # für not_bayesian
    kappa_square_conditional_all_coeffs <- K_inv_list[[iestimator]] %*% 
      Ct.R_inv.C_list[[iestimator]] %*%
      K_inv_list[[iestimator]]
    
    # Varianz kappa_square bestimmen: K_inv 
    # für bayesian 
    bayesian_kappa_square_conditional_all_coeffs <- K_inv_list[[iestimator]]
    
    
    
    
    # Eintrag pro Koeffizient in das entsprechende Listenelement
    for (iwhich in which_conditional){
      # Name des Koeffizienten zur Eintragung in Liste
      coef_name <- paste("coefficient", iwhich)
      
      # Testvektoren vt
      vt_conditional_list[[iestimator]][[coef_name]] <-
        splitRatio * matrix(V_conditional[iwhich,], nrow = 1)
      
      #Varianzen kappa_square und dirV 
      if("not_bayesian" %in% bayesian){
        # Varianz kappa_square:= e_j'*K_inv*C'*R_inv*C*K_inv*e_j
        kappa_square_conditional_list[[iestimator]][["not_bayesian"]][[coef_name]] <- 
          splitRatio^2 *kappa_square_conditional_all_coeffs[iwhich,iwhich]                                    #?# hmm so ist das natürlich ein wenig ineffizient weil der ja immer die gesamte Matrix V berechnet und nicht direkt nur die für which relevanten Zeilen 
        
        #dirV:= R*v/kappa_square wobei v=t(vt)
        dirV_conditional_list[[iestimator]][["not_bayesian"]][[coef_name]] <- 
          R %*% t(vt_conditional_list[[iestimator]][[coef_name]])/
          kappa_square_conditional_list[[iestimator]][["not_bayesian"]][[coef_name]]} 
      
      
      if("bayesian" %in% bayesian){
        # Varianz kappa_square:= e_j'*K_inv*e_j
        kappa_square_conditional_list[[iestimator]][["bayesian"]][[coef_name]] <- 
          splitRatio^2 *bayesian_kappa_square_conditional_all_coeffs[iwhich,iwhich]                                     #?# hmm so ist das natürlich ein wenig ineffizient weil der ja immer die gesamte Matrix V berechnet und nicht direkt nur die für which relevanten Zeilen 
        
        #dirV:= R*v/kappa_square wobei v=t(vt)
        dirV_conditional_list[[iestimator]][["bayesian"]][[coef_name]] <- 
          R %*% t(vt_conditional_list[[iestimator]][[coef_name]])/
          kappa_square_conditional_list[[iestimator]][["bayesian"]][[coef_name]]}
      
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


