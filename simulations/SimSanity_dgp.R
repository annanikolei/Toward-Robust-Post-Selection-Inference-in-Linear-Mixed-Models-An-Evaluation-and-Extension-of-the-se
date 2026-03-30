################################################################################
# Sanity Check Simulation - data generating process 
################################################################################

## similar to candidate models in Simulation Study 2 but with random X##

SanCheck_dgp  <- function(n_obs = 40^2, # number of observations
                          m = 40, # number of level-2 units
                          R_squared = 0.6, # pseudo coefficient of determination
                          D = matrix(c(2,0.2,
                                       0.2,1.5), nrow=2, byrow=T)
                          
){
  
  
  ## balanced design
  l2.id <- rep(1:m, each = n_obs/m)
  
  ## number of fixed effects (-intercept)
  p.l1 <- 2 # X1.l1,..., X4.l1
  
  ## number of random effects
  q <- 2 # u0,u2
  
  
  ## draw random effects and X from multivariate normal distribution -------------
  
  # covariance matrix of predictor variables X1.l1 and X2.l1 
  # (var(X1.l1) = var(X2.l1) = 1)
  vcovMatrix_X <- diag(p.l1)
  
  # random effects covariance matrix per level2-unit
  G <- bdiag(replicate(m, D, simplify=F))
  
  ## draw (X1.l1ij,X2.l1ij,X.l2j) and (u0j,u1j) from (joint) multivariate normal distribution 
  # X level-1 
  X.l1 <- rmvnorm(n = n_obs, 
                  mean = rep(2, p.l1), # E(X) = 2
                  sigma = vcovMatrix_X,
                  checkSymmetry = F)
  X.l1 <- cbind(1, X.l1) # add intercept-column
  colnames(X.l1) <- c("intercept","X1","X2")
  
  # X 
  X <- X.l1
  
  
  # u
  u_matrix <- rmvnorm(n = m,
                      mean = rep(0, q), # E(u) = 0
                      sigma = D,
                      checkSymmetry = F)
  
  u <- as.vector(t(u_matrix)) 
  
  
  ## beta
  beta <- c(1,4,0) # c00,c10,c20
  
  ## use random effects structure to generate Z ----------------------------------
  Z <- XzuZ(X = X,
            level2.index = l2.id,
            random.cols = c(1,3)) # intercept, X2
  
  
  ## get level-1 residuals -------------------------------------------------------
  
  # R^2 = Var(X*beta)/Var(Y)
  # Var(X*beta)= Var(-4*X1)=16*1
  
  # Var(uk*Xk) = E(xk)^2 * tauk^2 + tauk^2 * Var(Xk) = (4+1)*tauk^2 
  # KOV(u0,uk*Xk) = E(xk) * tau0k = 2* tau0k
  # KOV(uk*Xk, uj*Xj) = E(Xk)*E(Xj) * taujk = 4*taujk
  # KOV(Xk, uk*Xk)= 0
  
  # Var(y)= Var(-4*X1+u0+u2*X2+e)
  #       = 16*1 +
  #         tau0^2 + 2*(2*tau02) +
  #         5* tau2^2 +
  #         sigma_epsilon^2
  
  # sigma_square_epsilon = Var(X * beta)/R_squared - Var(Y\epsilon)
  
  tau0 <- sqrt(D[1,1])
  tau02 <- D[1,2]
  tau2 <- sqrt(D[2,2])
  
  sigma_square_epsilon <- (16*1)/R_squared - (16*1 + 
                                                   tau0^2 + 2*(2*tau02) +
                                                   5* tau2^2)
  epsilon <- rnorm(n = n_obs,
                   mean = 0,
                   sd = sqrt(sigma_square_epsilon))
  
  
  ## generate y ------------------------------------------------------------------
  y <-  X %*% beta + Z %*% u + epsilon
  y <- as.vector(y)
  
  
  ## output ----------------------------------------------------------------------
  
  return(list(y = y,
              X = X,
              Z = Z,
              G = G,
              sigma_square_epsilon = sigma_square_epsilon,
              l2.id = l2.id))
  }


