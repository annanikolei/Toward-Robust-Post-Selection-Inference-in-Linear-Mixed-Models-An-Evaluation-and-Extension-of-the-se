################################################################################
# CR-Standardfehler
## copied from Huang; https://francish.net/post/robust-standard-errors/
################################################################################
library(CR2)
library(clubSandwich)

CR_std.errors <- function(mod,
                          this_y,
                          X,
                          Sigma,
                          bread,
                          vt_marginal_all_coeffs,
                          robust_std.errors_marginal
){
  
  # Extahiere die Cluster-Variable
  cluster <- getME(mod, "flist")
  
  # Falls nur eine nesting-Variable verwendet wurde, können die robust std.err
  # berechnet werden -----------------------------------------------------------
  
  if(length(cluster)>1){ 
  stop("Berechnung von robust standard errors nicht möglich, da mehr als eine
  cluster-Variable verwendet wurde.") 
      }
  
  if(length(cluster)==1){
    
    cluster <- cluster[[1]]
    
    # residuals
    beta_fixef <- vt_marginal_all_coeffs %*% this_y
    rr <- this_y - X %*% beta_fixef
    
    Xl <- Sigmal <- Rl <- list()
    cnames <- names(table(cluster))
    
    ## to get Xc
    Xlist <- function(x){
      sel <- which(cluster == x)
      Xl[[x]] <- X[sel,]
    }
    
    Xl <- lapply(cnames, Xlist)
    ## lists become numeric? just add the names back later...
    
    Sigmalist <- function(x){
      sel <- which(cluster == x)
      Sigmal[[x]] <- solve(Sigma[sel, sel]) #can use chol2inv but using solve for transparency
    }
    
    Sigmal <- lapply(cnames, Sigmalist)
    names(Xl) <- names(Sigmal) <- cnames ## adding names back?
    
    ## residuals in a list -- for the meat
    Rl <- split(rr, cluster) #has names; this is easier for vectors
    
  ## CR0 std.err ---------------------------------------------------------------
    if(robust_std.errors_marginal == "CR0"){
      
      meat <- function(x){
        t(Rl[[x]]) %*% Sigmal[[x]] %*% Xl[[x]]  #compute part of the the meat matrix
      }
      
      u <- t(sapply(cnames, meat))
      u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
      colnames(u) <- colnames(u[[1]])
      
      
      meat <- t(u) %*% u #the full meat
      clvc.cr <- bread %*% meat %*% bread #variance covariance matrix
      
    }
    
  ## CR2 std.err ---------------------------------------------------------------
    if(robust_std.errors_marginal == "CR2"){
      
      tXs <- function(s) {
        sel <- which(cluster == s)
        Sigmac <- Sigma[sel, sel] #V matrix per cluster
        Ijj <- diag(nrow(Xl[[s]])) #Identity matrix
        Hjj <- Xl[[s]] %*% bread %*% t(Xl[[s]]) %*% Sigmal[[s]] 
        IHjj <- Ijj - Hjj #identity less Hat matrix
        #MatSqrtInverse(Ijj - Hjj) #early adjustment / valid
        V3 <- chol(Sigmac) #based on MBB; #see footnote 5 in HWZ; upper triangle
        Bi <- V3 %*% IHjj %*% Sigmac %*% t(V3)
        t(V3) %*% MatSqrtInverse(Bi) %*% V3 #need the CR2 package installed
      } 
      
      tX <- lapply(cnames, tXs) #these are the adjustment matrices
      names(tX) <- cnames
      
      cc2 <- function(x){
        Rl[[x]] %*% tX[[x]] %*% Sigmal[[x]] %*% Xl[[x]] # r x Adj x Vinv x Xc
      }
      
      u <- t(sapply(cnames, cc2)) 
      u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
      colnames(u) <- colnames(u[[1]])
      
      meat <- t(u) %*% u #the full meat matrix
      clvc.cr <- bread %*% meat %*% bread #variance covariance matrix
      
    }
    
    return(clvc.cr)
    
  }
  
}

CR_std.errors_OLS <- function(AdjMcCaffrey=TRUE,
                              mod,
                              this_y,
                              X,
                              cluster,
                              Sigma,
                              bread,
                              vt_marginal_all_coeffs,
                              robust_std.errors_marginal
){
  
  robust_std.errors <- robust_std.errors_marginal
  
  # Extract Cluster-Variable
  cluster <- getME(mod, "flist")
  
  if(length(cluster)!=1){ 
    stop("Berechnung von robust standard errors nicht möglich, da mehr als eine
  cluster-Variable verwendet wurde.") 
  }
  
  cluster <- cluster[[1]]
  
  
  # residuals -----------------------------------------------------------------
  beta_fixef <- vt_marginal_all_coeffs %*% this_y
  rr <- this_y - X %*% beta_fixef
  
  cnames <- names(table(cluster))
  
  ## to get Xc
  Xl <- list()
  Xlist <- function(x){
    sel <- which(cluster == x)
    Xl[[x]] <- X[sel, , drop = FALSE]
  }
  Xl <- lapply(cnames, Xlist)
  names(Xl) <- cnames
  
  ## residuals in a list -- for the meat
  Rl <- split(rr, cluster) # has names; easier for vectors
  
  ## CR0 std.err ---------------------------------------------------------------
  if(robust_std.errors == "CR0"){
    
    cc0 <- function(x){
      t(Rl[[x]]) %*% Xl[[x]]   # ê_c' X_c  (row vector, p-dim)
    }
    
    u <- t(sapply(cnames, cc0))
    #u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
    colnames(u) <- colnames(u[[1]])
    
    meat <- t(u) %*% u
    clvc.cr <- bread %*% meat %*% bread
    
  }
  
  ## CR2 std.err ---------------------------------------------------------------
  
  if(robust_std.errors == "CR2"){
   
    if(AdjMcCaffrey){ 
    
    tXs <- function(s){
      sel <- which(cluster == s)
      Xc  <- Xl[[s]]
      Sigmac <- Sigma[sel, sel, drop = FALSE]  # V_c = Sigma_c, the working (estimated) marginal covariance
      # -- used here even though beta_fixef is OLS, so that kappa^2_CR
      # rests on the same covariance assumption as the model-based
      # kappa^2 in Eq.~\eqref{eq:OLS_NV}, rather than assuming V_c = I_c
      Ijj <- diag(nrow(Xc))            # I_c
      Hjj <- Xc %*% bread %*% t(Xc)    # H_c = X_c (X'X)^{-1} X_c'  (unweighted OLS hat matrix -- bread stays OLS)
      IHjj <- Ijj - Hjj                # I_c - H_c
      U3 <- chol(Sigmac)               # Cholesky root of V_c = Sigma_c; based on MBB, see footnote 5 in HWZ; upper triangle
      Bi <- U3 %*% IHjj %*% Sigmac %*% t(U3)
      t(U3) %*% MatSqrtInverse(Bi) %*% U3   # need the CR2 package installed
    }
    
    tX <- lapply(cnames, tXs)
    names(tX) <- cnames
    
    cc2 <- function(x){
      Rl[[x]] %*% tX[[x]] %*% Xl[[x]]  # ê_c' Adj_c X_c
    }
    
    u <- t(sapply(cnames, cc2))
    u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
    colnames(u) <- colnames(u[[1]])
    
    meat <- t(u) %*% u
    clvc.cr <- bread %*% meat %*% bread
    
    }
    else if(!AdjMcCaffrey){

        tXs <- function(s){
          sel <- which(cluster == s)
          Xc  <- Xl[[s]]
          Ijj <- diag(nrow(Xc))            # I_c
          Hjj <- Xc %*% bread %*% t(Xc)    # H_c = X_c (X'X)^{-1} X_c'
          IHjj <- Ijj - Hjj                # I_c - H_c
          MatSqrtInverse(IHjj)             # [I_c - H_c]^{-1/2}  (needs clubSandwich)
        }

        tX <- lapply(cnames, tXs)
        names(tX) <- cnames

        cc2 <- function(x){
          Rl[[x]] %*% tX[[x]] %*% Xl[[x]]  # ê_c' [I_c-H_c]^{-1/2} X_c
        }

        u <- t(sapply(cnames, cc2))
        # u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
        colnames(u) <- colnames(u[[1]])

        meat <- t(u) %*% u
        clvc.cr <- bread %*% meat %*% bread

  }
    
  }
    
  return(clvc.cr)
}

# CR_std.errors_OLS <- function(mod,
#                               this_y,
#                               X,
#                               cluster,
#                               bread,
#                               vt_marginal_all_coeffs,
#                               robust_std.errors_marginal
# ){
#   
#   robust_std.errors <- robust_std.errors_marginal
#   
#   # ExtractCluster-Variable
#   cluster <- getME(mod, "flist")
#   
#   if(length(cluster)!=1){ 
#     stop("Berechnung von robust standard errors nicht möglich, da mehr als eine
#   cluster-Variable verwendet wurde.") 
#   }
#   
#  
#    cluster <- cluster[[1]]
#   
#   
#   # residuals -----------------------------------------------------------------
#   beta_fixef <- vt_marginal_all_coeffs %*% this_y
#   rr <- this_y - X %*% beta_fixef
#   
#   cnames <- names(table(cluster))
#   
#   ## to get Xc
#   Xl <- list()
#   Xlist <- function(x){
#     sel <- which(cluster == x)
#     Xl[[x]] <- X[sel, , drop = FALSE]
#   }
#   Xl <- lapply(cnames, Xlist)
#   names(Xl) <- cnames
#   
#   ## residuals in a list -- for the meat
#   Rl <- split(rr, cluster) # has names; easier for vectors
#   
#   ## CR0 std.err ---------------------------------------------------------------
#   if(robust_std.errors == "CR0"){
#     
#     cc0 <- function(x){
#       t(Rl[[x]]) %*% Xl[[x]]   # ê_c' X_c  (row vector, p-dim)
#     }
#     
#     u <- t(sapply(cnames, cc0))
#     #u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
#     colnames(u) <- colnames(u[[1]])
#     
#     meat <- t(u) %*% u
#     clvc.cr <- bread %*% meat %*% bread
#     
#   }
#   
#   ## CR2 std.err ---------------------------------------------------------------
#   if(robust_std.errors == "CR2"){
#     
#     tXs <- function(s){
#       sel <- which(cluster == s)
#       Xc  <- Xl[[s]]
#       Ijj <- diag(nrow(Xc))            # I_c
#       Hjj <- Xc %*% bread %*% t(Xc)    # H_c = X_c (X'X)^{-1} X_c'
#       IHjj <- Ijj - Hjj                # I_c - H_c
#       MatSqrtInverse(IHjj)             # [I_c - H_c]^{-1/2}  (needs clubSandwich)
#     }
#     
#     tX <- lapply(cnames, tXs)
#     names(tX) <- cnames
#     
#     cc2 <- function(x){
#       Rl[[x]] %*% tX[[x]] %*% Xl[[x]]  # ê_c' [I_c-H_c]^{-1/2} X_c
#     }
#     
#     u <- t(sapply(cnames, cc2))
#     # u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
#     colnames(u) <- colnames(u[[1]])
#     
#     meat <- t(u) %*% u
#     clvc.cr <- bread %*% meat %*% bread
#     
#   }
#   
#   return(clvc.cr)
# }

