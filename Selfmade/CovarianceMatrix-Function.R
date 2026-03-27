################################################################################
# Schätzung der Kovarianzmatrizen-Funktionen
################################################################################

#-------------------------------------------------------------------------------
# Funktion zur Berechnung von R, R_invers für supplied, ModelEstimate und ICM
#-------------------------------------------------------------------------------
## Bemerkung: für den VarY-estimator bei der conditional perspective


R_R_invers <- function(mod, # fit des selektierten Modell
                       this_y = NULL, # Vektor an vorliegenden y-Werten
                       n, # Länge von this_y
                       marginal_COV_estimator = c("ModelEstimate", 
                                                  "ICM",
                                                  "supplied"), # KOV-estimator
                       conditional_COV_estimator = c("ModelEstimate", 
                                                     "ICM",
                                                     "VarY",
                                                     "supplied"),# KOV-estimator
                       fit_ICM = NULL,
                       Sigma_supplied = NULL,
                       sigma_square_epsilon_supplied = NULL){
  
  # Liste für return vorbereiten
  COV_estimator <- unique(c(marginal_COV_estimator, conditional_COV_estimator))
  R_return_list <- R_inv_return_list <- setNames(vector("list",
                                                        length(COV_estimator)),
                                                 COV_estimator)
  
  
  
  ## SUPPLIED
  # Falls supplied bei marginal und kein Sigma_supplied angegeben ist oder
  # falls supplied bei conditional, berechne aus dem supplied sigma-square_error
  # das R_inv_supplied
  if((c("supplied") %in% marginal_COV_estimator && is.null(Sigma_supplied))|
     c("supplied") %in% conditional_COV_estimator){
    
    if(!is.null(sigma_square_epsilon_supplied)){
      
      R_supplied <- Matrix(diag(sigma_square_epsilon_supplied,n),sparse=T)
      # in Listen übergeben 
      R_return_list$supplied <- R_supplied
      
      R_supplied_invers <- Matrix(diag(1/sigma_square_epsilon_supplied,n),
                                  sparse=T) # Inverse Diagonalmatrix                             #?# Muss ich beim R_supplied_invers auch mit Rmpfr arbeiten?
      # in Listen übergeben 
      R_inv_return_list$supplied <- R_supplied_invers
    } 
    else if(is.null(sigma_square_epsilon_supplied)){
      stop("Da COV_estimator = 'supplied' muss entweder Sigma_supplied(marginal) 
            oder sigma_square_epsilon_supplied (marginal, conditional) 
            angegeben sein")
    }
  }
  
  ## MODELESTIMATE
  if(c("ModelEstimate") %in% COV_estimator){
    # Extrahiere estimated sigma_square_epsilon (Level-1-Residualvarianz)
    sigma_square_epsilon <- sigma(mod)^2
    
    R_est <- Matrix(diag(sigma_square_epsilon,n), sparse = T )
    # in Listen übergeben 
    R_return_list$ModelEstimate <- R_est
    
    
    R_est_invers <- Matrix(diag((1/sigma_square_epsilon),n),
                           sparse = T) # Inverse Diagonalmatrix                            #?# Rmpfr?
    # in Listen übergeben 
    R_inv_return_list$ModelEstimate <- R_est_invers
  }
  
  ## ICM
  # Schätzung des minimalen Intercept Modells zur späteren Berechnung des
  # ICM-Schätzers. Wir machen das aber bereits hier, damit es nicht bei marginal
  # und conditional perspective, also ggf. zweimal geladen werden muss
  if(c("ICM") %in%  COV_estimator){
    sigma_square_epsilon_ICM <- getME(fit_ICM, "sigma")^2
    
    R_ICM <- Matrix(diag(sigma_square_epsilon_ICM,n), sparse = T)
    # in Listen übergeben 
    R_return_list$ICM<- R_ICM
    
    R_ICM_invers <- Matrix(diag((1/sigma_square_epsilon_ICM),n),
                           sparse = T) # Inverse Diagonalmat                          #?# Rmpfr?
    # in Listen übergeben 
    R_inv_return_list$ICM <- R_ICM_invers
  }
  
  ## VarY
  ## geht nur bei diagonaler KOV-Matrix
  if(c("VarY") %in%  conditional_COV_estimator){
    # plugin für sigma_square_epsilon by Tibshirani et al. 2018
    sigma_square_epsilon_VarY <- as.numeric(var(this_y)*(n-1)/n)
    
    R_VarY <- Matrix(diag(sigma_square_epsilon_VarY,n), sparse = T)
    # in Listen übergeben 
    R_return_list$VarY<- R_VarY
    
    R_VarY_invers <- Matrix(diag((1/sigma_square_epsilon_VarY),n),
                            sparse = T) # Inverse Diagnlmtrx
    # in Listen übergeben 
    R_inv_return_list$VarY <- R_VarY_invers
  }
  
  R_R_inv_return_list <- setNames(vector("list",2),c("R","R_inv"))
  R_R_inv_return_list$R <- R_return_list
  R_R_inv_return_list$R_inv <- R_inv_return_list
  return(R_R_inv_return_list)
  
}


#-------------------------------------------------------------------------------
# Funktion zur Berechnung von R, R_invers für nicht-diagonales R
#-------------------------------------------------------------------------------


R_R_invers_nondiagonal <- function(
    marginal_COV_estimator = c("ModelEstimate", 
                               "ICM",
                               "supplied"), # KOV-estimator
    conditional_COV_estimator = c("ModelEstimate", 
                                  "ICM",
                                  "VarY",
                                  "supplied"),# KOV-estimator
    R_supplied){
  
  # Liste für return vorbereiten
  COV_estimator <- unique(c(marginal_COV_estimator, conditional_COV_estimator))
  R_return_list <- R_inv_return_list <- setNames(vector("list",
                                                        length(COV_estimator)),
                                                 COV_estimator)
  
  
  
  ## SUPPLIED
  # Falls supplied bei marginal oder
  # falls supplied bei conditional, berechne aus dem supplied R und Rinvers
  if(c("supplied") %in% marginal_COV_estimator |
     c("supplied") %in% conditional_COV_estimator){
    
    R_supplied <- as.matrix(R_supplied)
    
    # in Listen übergeben 
    R_return_list$supplied <- R_supplied
    
    R_supplied_invers <- solve(R_supplied)
    
    # in Listen übergeben 
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
# Funktion zur Berechnung von Sigma für nicht-diagonales R für not efficient
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
  
  # Liste für return vorbereiten
  COV_estimator <- marginal_not_efficient_COV_estimator
  SIGMA_return_list <- setNames(vector("list", length(COV_estimator)),COV_estimator)
  
  ## SUPPLIED
  if("supplied" %in% marginal_not_efficient_COV_estimator){
    
    if(!is.null(Sigma_supplied) && !is.null(R_supplied)){stop("Entweder Sigma_supplied oder
                                                              R_supplied angeben.")}
    
    if(!is.null(Sigma_supplied)) {SIGMA_return_list$supplied <- Sigma_supplied}
    
    else {SIGMA_return_list$supplied <- Z_supplied %*% G_supplied %*% t(Z_supplied)+ R_supplied}
    
  }
  
  
  ## MODELESTIMATE
  if("ModelEstimate_R_replaced" %in% marginal_not_efficient_COV_estimator |
     "ModelEstimate_SIGMA_replaced" %in% marginal_not_efficient_COV_estimator ){
    
    # Fixed Effects Modell schätzen
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
# Funktion zur Berechnung von Sigma für die marginale not efficient Perspektive, wenn
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
  
  # Leere Liste fürSIGMA erstellen, sodass jedes Listenelement einen
  #  Kovarianzmatrix-Schätzer enthalten kann
  SIGMA <- setNames(vector("list", length(marginal_not_efficient_COV_estimator)),
                    marginal_not_efficient_COV_estimator)
  
  ## SUPPLIED
  ## Falls marginal_not_efficient_COV_estimator == "supplied",d.h. falls Sigma angegeben ist
  if(c("supplied") %in% marginal_not_efficient_COV_estimator){
    
    # Falls Sigma_supplied angegeben, schreibe es in SIGMA-Matrix
    if(!is.null(Sigma_supplied)){SIGMA$supplied <- Sigma_supplied}
    
    # Falls Sigma_supplied nicht angegeben ist aber R,Z,Lambda,Lamdat_supplied
    #  oder G_supplied,dann berechne Sigma_supplied und trage dies in
    #  SIGMA ein
    else if (is.null(Sigma_supplied) && !is.null(Z_supplied)){
      
      # Falls G_supplied nicht angegeben ist aber Lambda_supplied, berechne 
      # damit die Schätzung von Sigma, ansonsten mit G_supplied
      if(!is.null(G_supplied)){G_supplied <- G_supplied}
      
      else if (is.null(G_supplied) && !is.null(Lambda_supplied)) {
        if(is.null(Lambdat_supplied)){Lambda_supplied <- t(Lambda_supplied)}
        G_supplied <- Lambda_supplied %*% Lambdat_supplied}
      
      else if(is.null(G_supplied) && is.null(Lambda_supplied)){
        stop("Da marginal_not_efficient_COV_estimator = 'supplied' muss für marginale not efficient 
              Perspektive entweder Lambda_supplied oder G_supplied angegeben sein")}
      
      # Woodbury-Matrix-Identität zur effizienteren Berechnung der Inverse
      SIGMA$supplied <- Z_supplied %*% G_supplied %*% t(Z_supplied) + R_supplied
    }
    
    # Else: Da "supplied angegeben wurde, muss auch Sigma_supplied nichtNULL 
    # sein oder zumindest muss dann Z_supplied und G_supplied (bzw.Lambda) und 
    # R_supplied angegeben sein 
    else {stop("Da Da marginal_not_efficient_COV_estimator  = 'supplied' muss für 
                marginale not efficient Perspektive entweder Sigma_supplied oder
                (G_supplied bzw.Lambda_supplied) und Z_supplied angegeben sein")}
  }
  
  ## MODELESTIMATE
  ## Falls marginal_not_efficient_COV_estimator == "ModelEstimate", d.h. falls Sigma durch 
  ## den plug-in Schätzer Sigma_dach also ModelEstimate geschätzt werden soll
  if(c("ModelEstimate") %in% marginal_not_efficient_COV_estimator){
    # Berechne den plug-in Schätzer für SIGMA 
    
    SIGMA$ModelEstimate <- Z_B_work %*% G_work_est %*% t(Z_B_work) + R_est
    
  }
  
  ## ICM
  ## Falls marginal_not_efficient_COV_estimator == "ICM", d.h. falls Sigma durch 
  ## den konservativen minMod Schätzer also ICM geschätzt werden soll
  if(c("ICM") %in% marginal_not_efficient_COV_estimator){
    
    # Berechne den ICM Schätzer für SIGMA
    model <- fit_ICM
    SIGMA$ICM <- Sigma_ICM(model = model,
                           n = n,
                           R_ICM = R_ICM,
                           Z = Z_B_work)
    
    
  }
  
  return(SIGMA)
  
}

#-------------------------------------------------------------------------------
# Funktion zur Berechnung von Sigma_invers für die marginale Perspektive
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
  
  
  # Leere Liste für SIGMA_invers erstellen, sodass jedes Listenelement einen
  #  Kovarianzmatrix-Schätzer enthalten kann
  SIGMA_invers <- setNames(vector("list", length(marginal_COV_estimator)),
                           marginal_COV_estimator)
  
  ## SUPPLIED
  ## Falls marginal_COV_estimator == "supplied",d.h. falls Sigma angegeben ist
  if(c("supplied") %in% marginal_COV_estimator){
    
    # Falls Sigma_supplied angegeben, schreibe Inverse in SIGMA_invers-Matrix
    if(!is.null(Sigma_supplied)){
      SIGMA_invers$supplied <- solve(Sigma_supplied)}
    
    # Falls Sigma_supplied nicht angegeben ist aber R,Z,Lambda,Lamdat_supplied
    #  oder G_supplied,dann berechne Sigma_supplied und trage diese Inverse in
    #  SIGMA_invers ein
    else if (is.null(Sigma_supplied) && !is.null(Z_supplied)){
      
      # Falls G_supplied nicht angegeben ist aber Lambda_supplied, berechne 
      # damit die Schätzung von Sigma_invers, ansonsten mit G_supplied
      if(!is.null(G_supplied)){G_supplied <- G_supplied}
      
      else if (is.null(G_supplied) && !is.null(Lambda_supplied)) {
        if(is.null(Lambdat_supplied)){Lambda_supplied <- t(Lambda_supplied)}
        G_supplied <- Lambda_supplied %*% Lambdat_supplied}
      
      else if(is.null(G_supplied) && is.null(Lambda_supplied)){
        stop("Da marginal_COV_estimator = 'supplied' muss für marginale 
                  Perspektive 
                  entweder Lambda_supplied oder G_supplied angegeben sein")}
      
      # Woodbury-Matrix-Identität zur effizienteren Berechnung der Inverse
      SIGMA_invers$supplied <- tryCatch({
        woodbury_identity(R_inv = R_supplied_invers,
                          Z = Z_supplied,
                          G = G_supplied)
      }, error = function(e) {
        solve(Z_supplied %*% G_supplied %*% t(Z_supplied) + R_supplied)
      })
      
      
    }
    
    # Else: Da "supplied angegeben wurde, muss auch Sigma_supplied nichtNULL 
    # sein oder zumindest muss dann Z_supplied und G_supplied (bzw.Lambda) und 
    # R_supplied angegeben sein 
    else {stop("Da marginal_COV_estimator = 'supplied' muss für marginale 
                  Perspektive 
                  entweder Sigma_supplied oder (G_supplied bzw.Lambda_supplied)
                  und Z_supplied angegeben sein")}
  }
  
  ## MODELESTIMATE
  ## Falls marginal_COV_estimator == "ModelEstimate", d.h. falls Sigma durch 
  ## den plug-in Schätzer Sigma_dach also ModelEstimate geschätzt werden soll
  if(c("ModelEstimate") %in% marginal_COV_estimator){
    # Berechne den plug-in Schätzer für SIGMA_invers mit Woodbury-Matrix
    
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
  ## Falls marginal_COV_estimator == "ICM", d.h. falls Sigma durch 
  ## den konservativen minMod Schätzer also ICM geschätzt werden soll
  if(c("ICM") %in% marginal_COV_estimator){
    
    # Berechne den ICM Schätzer für SIGMA_invers
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
# Funktion zur Berechnung von R_list,K_inv_list, Ct.R_inv_list, Ct.R_inv.C_list
# für die konditionale Perspektive
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
  
  
  # Leere Liste für K_invers:=(t(C)%*%R_inv%*%C+Q)_inv erstellen, sodass jedes 
  # Listenelement einen Schätzer für K_inv enthalten kann
  K_inv_list <- setNames(vector("list",length(conditional_COV_estimator)),
                         conditional_COV_estimator)
  
  # Leere Liste für Ct.R_inv.C:=t(C)%*%R_inv%*%C  und Ct.R_inv:=t(C)%*%R_inv
  # erstellen
  Ct.R_inv_list <- setNames(vector("list",length(conditional_COV_estimator)),
                            conditional_COV_estimator)
  Ct.R_inv.C_list <- setNames(vector("list",length(conditional_COV_estimator)),
                              conditional_COV_estimator)
  
  # Berechne Blockmatrix C zur Berechnung von K_inv
  C_block <- cbind(X_A_work,Z_B_work)
  C <- C_block
  Ct <- t(C)
  
  # Berechnung der Inverse von G_work_est, wenn dieses als Schätzer aus 
  # mod übernommen wurde. Sowie Berechnung von Q_est
  if(any(c("ModelEstimate","ICM","VarY") %in% conditional_COV_estimator)){
    
    # Berechnung der Inverse von G_work_est
    G_work_est_invers <- solve(G_work_est) 
    
    # Berechnung von Q_est
    Q_est <- bdiag(matrix(0, ncol=ncol(X_A_work), nrow=ncol(X_A_work)),
                   G_work_est_invers)
  }
  
  ## SUPPLIED
  # Falls conditional_COV_estimator == "supplied", nutze G_work_supplied
  # [also
  # die Kovarianzmatrix der random effects, die man erhalten würde, wenn man
  # mit der Population das working model schätzen würde (also wenn das
  # eigentlich ein Populationsmodell mit random slopes ist, aber man das nur
  # mit random intercepts schätzt z.B)]
  # und R_supplied_invers zur Berechnung von K_invers
  if(c("supplied") %in% conditional_COV_estimator){
    
    if(!is.null(G_work_supplied)){
      # Berechnung der Block-Matrix Q: zunächst Berechnung von G_work_invers,
      # danach die Bildung der Blockdiagonalmatrix
      G_work_supplied_invers <- solve(G_work_supplied)
      Q_supplied <- bdiag(matrix(0, ncol=ncol(X_A_work), nrow=ncol(X_A_work)),
                          G_work_supplied_invers)
      
      # Berechnung der Listeneinträge von Ct.R_inv_list und Ct.R_inv.C_list
      Ct.R_inv_supplied <- Ct %*% R_supplied_invers
      Ct.R_inv_list$supplied <- Ct.R_inv_supplied
      
      Ct.R_inv.C_supplied <- Ct.R_inv_supplied %*% C
      Ct.R_inv.C_list$supplied <- Ct.R_inv.C_supplied
      
      # Berechnung des enstprechenden Eintrags von K_inv                                    #!# Woodbury-Identität geht nicht weil Q nicht invertierbar ist 
      K_inv_list$supplied <- solve(Ct.R_inv.C_supplied + Q_supplied)
    }
    
    
    # Else: Da "supplied angegeben wurde, muss G_work_supplied
    # angegeben sein.
    else {stop("Da conditional_COV_estimator = 'supplied' muss
                 G_work_supplied angegeben sein")}
    
  }
  
  ## MODELESTIMATE
  # Falls conditional_COV_estimator == "ModelEstimate", nutze Q_est
  # und R_est_invers zur Berechnung von K_invers
  if(c("ModelEstimate") %in% conditional_COV_estimator){
    
    # Berechnung der Listeneinträge von Ct.R_inv_list und Ct.R_inv.C_list
    Ct.R_inv_ModelEstimate <- Ct %*% R_est_invers 
    Ct.R_inv_list$ModelEstimate <- Ct.R_inv_ModelEstimate
    
    Ct.R_inv.C_ModelEstimate <- Ct.R_inv_ModelEstimate %*% C
    Ct.R_inv.C_list$ModelEstimate <- Ct.R_inv.C_ModelEstimate
    
    # Berechnung des enstprechenden Eintrags von K_inv
    K_inv_list$ModelEstimate <- solve(Ct.R_inv.C_ModelEstimate + Q_est)}
  
  ## ICM
  # Falls conditional_COV_estimator == "ICM", nutze Q_est
  # und R_ICM_invers zur Berechnung von K_invers
  # !!! das sollte aber nur verwendet werden, wenn keine random effects
  # selektiert werden sollen
  if(c("ICM") %in% conditional_COV_estimator){
    
    
    # Berechnung der Listeneinträge von Ct.R_inv_list und Ct.R_inv.C_list
    Ct.R_inv_ICM <- Ct %*% R_ICM_invers
    Ct.R_inv_list$ICM <- Ct.R_inv_ICM
    
    Ct.R_inv.C_ICM <- Ct.R_inv_ICM %*% C
    Ct.R_inv.C_list$ICM <- Ct.R_inv.C_ICM
    
    # Berechnung des enstprechenden Eintrags von K_inv                                    #!# Woodbury-Identität geht nicht weil Q nicht invertierbar ist 
    K_inv_list$ICM <- solve(Ct.R_inv.C_ICM + Q_est)}
  
  ## VarY
  # Falls conditional_COV_estimator == "ICM", nutze Q_est
  # und R_ICM_invers zur Berechnung von K_invers
  if(c("VarY") %in% conditional_COV_estimator){
    
    # Berechnung der Listeneinträge von Ct.R_inv_list und Ct.R_inv.C_list
    #!# Woodbury-Identität geht nicht weil Q nicht invertierbar ist 
    Ct.R_inv_VarY <- Ct %*% R_VarY_invers
    Ct.R_inv_list$VarY <- Ct.R_inv_VarY
    
    Ct.R_inv.C_VarY <- Ct.R_inv_VarY %*% C
    Ct.R_inv.C_list$VarY <- Ct.R_inv.C_VarY
    
    # Berechnung des enstprechenden Eintrags von K_inv
    K_inv_list$VarY <- solve(Ct.R_inv.C_VarY + Q_est)}
  
  ## Rückgabeobjekte
  return_K_inv.Ct.R_inv_list <- setNames(vector("list",4),
                                         c("K_inv_list","Ct.R_inv_list","Ct.R_inv.C_list","C_block"))
  
  return_K_inv.Ct.R_inv_list$K_inv_list <- K_inv_list
  return_K_inv.Ct.R_inv_list$Ct.R_inv_list <- Ct.R_inv_list
  return_K_inv.Ct.R_inv_list$Ct.R_inv.C_list <- Ct.R_inv.C_list
  return_K_inv.Ct.R_inv_list$C_block <- C_block
  
  
  return(return_K_inv.Ct.R_inv_list)
}

