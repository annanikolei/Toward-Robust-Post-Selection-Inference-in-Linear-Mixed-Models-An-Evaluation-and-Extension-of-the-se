################################################################################
# CR-Standardfehler
## https://francish.net/post/robust-standard-errors/
################################################################################
install.packages("CR2")
library(CR2)

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


