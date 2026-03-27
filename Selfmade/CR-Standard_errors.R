################################################################################
# Cluster-robust standard errors
# Source: https://francish.net/post/robust-standard-errors/
################################################################################

install.packages("CR2")
library(CR2)

CR_std.errors <- function(mod,
                          this_y,
                          X,
                          Sigma,
                          bread,
                          vt_marginal_all_coeffs,
                          robust_std.errors_marginal) {
  
  # Extract clustering variable(s)
  cluster <- getME(mod, "flist")
  
  # If more than one clustering variable is present, computation is not supported
  if (length(cluster) > 1) { 
    stop("Computation of robust standard errors is not possible when more than one clustering variable is used.") 
  }
  
  # Case: exactly one clustering variable
  if (length(cluster) == 1) {
    
    cluster <- cluster[[1]]
    
    # Compute residuals
    beta_fixef <- vt_marginal_all_coeffs %*% this_y
    rr <- this_y - X %*% beta_fixef
    
    # Initialize lists
    Xl <- Sigmal <- Rl <- list()
    cnames <- names(table(cluster))  # cluster labels
    
    ## Build list of design matrices per cluster
    Xlist <- function(x) {
      sel <- which(cluster == x)
      Xl[[x]] <- X[sel, ]
    }
    Xl <- lapply(cnames, Xlist)
    
    ## Build list of inverse covariance matrices per cluster
    Sigmalist <- function(x) {
      sel <- which(cluster == x)
      Sigmal[[x]] <- solve(Sigma[sel, sel])  # could use chol2inv for efficiency
    }
    Sigmal <- lapply(cnames, Sigmalist)
    
    # Restore names
    names(Xl) <- names(Sigmal) <- cnames
    
    ## Residuals split by cluster (used for "meat" matrix)
    Rl <- split(rr, cluster)
    
    #---------------------------------------------------------------------------
    ## CR0 standard errors
    #---------------------------------------------------------------------------
    if (robust_std.errors_marginal == "CR0") {
      
      # Compute cluster-wise contributions to the meat matrix
      meat <- function(x) {
        t(Rl[[x]]) %*% Sigmal[[x]] %*% Xl[[x]]
      }
      
      u <- t(sapply(cnames, meat))
      u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
      colnames(u) <- colnames(u[[1]])
      
      # Full meat matrix
      meat <- t(u) %*% u
      
      # Cluster-robust variance-covariance matrix
      clvc.cr <- bread %*% meat %*% bread
    }
    
    #---------------------------------------------------------------------------
    ## CR2 standard errors
    #---------------------------------------------------------------------------
    if (robust_std.errors_marginal == "CR2") {
      
      # Compute adjustment matrices per cluster
      tXs <- function(s) {
        sel <- which(cluster == s)
        
        Sigmac <- Sigma[sel, sel]            # cluster-specific covariance matrix
        Ijj <- diag(nrow(Xl[[s]]))           # identity matrix
        Hjj <- Xl[[s]] %*% bread %*% t(Xl[[s]]) %*% Sigmal[[s]]  # hat matrix
        IHjj <- Ijj - Hjj                    # I - H
        
        # Adjustment based on CR2 methodology
        V3 <- chol(Sigmac)                   # Cholesky factor
        Bi <- V3 %*% IHjj %*% Sigmac %*% t(V3)
        
        # Requires CR2 package
        t(V3) %*% MatSqrtInverse(Bi) %*% V3
      } 
      
      # Compute adjustment matrices for all clusters
      tX <- lapply(cnames, tXs)
      names(tX) <- cnames
      
      # Compute cluster-wise contributions
      cc2 <- function(x) {
        Rl[[x]] %*% tX[[x]] %*% Sigmal[[x]] %*% Xl[[x]]
      }
      
      u <- t(sapply(cnames, cc2)) 
      u <- do.call(rbind, lapply(u, function(mat) as.numeric(mat)))
      colnames(u) <- colnames(u[[1]])
      
      # Full meat matrix
      meat <- t(u) %*% u
      
      # Cluster-robust variance-covariance matrix
      clvc.cr <- bread %*% meat %*% bread
    }
    
    return(clvc.cr)
  }
}
