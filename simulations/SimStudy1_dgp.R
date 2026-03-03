################################################################################
# Simulation Study 1: Data generating process
################################################################################


######################## dgp function ########################################## 

SimStudy1_dgp<- function( n_obs, # number of observations
                          ICC_pseudo # pseudo intraclass-correlation
){
  
  ## set number of level2-units
  m <- 30
  
  ## balanced design
  l2.id <- rep(1:m, each = n_obs/m)
  
  ## number of fixed effects (-intercept)
  p.l1 <- 2 # X1.l1 and X2.l1
  p.l2 <- 1 # P
  
  ## number of random effects
  q <- 2 # u0 and u1
  
  
  ## draw random effects and X from multivariate normal distribution -------------
  
  # covariance matrix of predictor variables X1.l1 and X2.l1 
  # (var(X1.l1) = var(X2.l1) = 1)
  vcovMatrix_X <- diag(p.l1)
  
  # random effects covariance matrix per level2-unit
  tau0_square <- 1.3
  tau1_square <- 1
  tau01 <- 0.3
  D <- matrix(c(tau0_square, tau01, tau01, tau1_square), nrow = q, byrow = T)
  G <- bdiag(replicate(m, D, simplify=F))
  
  ## draw (X1.l1ij,X2.l1ij,X.l2j) and (u0j,u1j) from (joint) multivariate normal distribution 
  # X level-1 
  X.l1 <- rmvnorm(n = n_obs, 
                  mean = rep(0, p.l1), # E(X) = 0
                  sigma = vcovMatrix_X,
                  checkSymmetry = F)
  X.l1 <- cbind(1, X.l1) # add intercept-column
  colnames(X.l1) <- c("intercept","X1.l1","X2.l1")
  
  # X level-2
  X.l2_per_l2unit <- rnorm(n = m,
                           mean = 0,
                           sd = 1)
  X.l2 <- vector(length = n_obs)
  for(i in 1:m){
    for(j in ((i - 1) * (n_obs / m) + 1):(i * (n_obs / m))){
      X.l2[j] <- X.l2_per_l2unit[i]
    }
  }
  X.l2 <- X.l2/sd(X.l2)
  
  # X 
  X <- cbind(X.l1,X.l2)
  
  
  # u
  u_matrix <- rmvnorm(n = m,
                      mean = rep(0, q), # E(u) = 0
                      sigma = D,
                      checkSymmetry = F)
  
  u <- as.vector(t(u_matrix)) # (u01,u11,...,u0j,u1j,u0j+1,u1j+1,...,u0m,u1m)
  
  
  ## use random effects structure to generate Z ----------------------------------
  Z <- XzuZ(X = X.l1,
            level2.index = l2.id,
            random.cols = c(1,2))
  
  
  ## get level-1 residuals -------------------------------------------------------
  
  # ICC = be_var/(be_var + sigma_square_epsilon)
  # <=> sigma_square_epsilon = be_var * (1-ICC)/ICC
  # be_var = tau_0**2 + var(X1.l1) * tau_1**2 = tau_0**2 + tau_1**2
  sigma_square_epsilon <- (tau0_square + tau1_square) * (1-ICC_pseudo)/ICC_pseudo
  epsilon <- rnorm(n = n_obs,
                   mean = 0,
                   sd = sqrt(sigma_square_epsilon))
  
  
  ## generate y ------------------------------------------------------------------
  # y = c00 + 0 * X.l2_j + 0 * X1.l1_ij + 0 * X2.l1_ij + u0_j + u1_j * X1.l1_ij
  #     + epsilon_ij
  #   = X %*% beta + Z %*% u + epsilon
  # c00 := 0
  y <-  0 + Z %*% u + epsilon
  y <- as.vector(y)
  
  
  ## output ----------------------------------------------------------------------
  
  return(list(y = y,
              X = X,
              Z = Z,
              G = G,
              sigma_square_epsilon = sigma_square_epsilon,
              l2.id = l2.id)
  )
  
  
}
