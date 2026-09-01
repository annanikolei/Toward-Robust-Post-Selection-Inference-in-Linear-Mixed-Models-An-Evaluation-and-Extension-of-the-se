################################################################################
# Simulation Study 2: Data generating process (quadratic dgp)
################################################################################

library(Matrix)
library(whitening)

######################## dgp function ########################################## 

# generate X seperately so that X and Z do not need to be constructed ----------
# multiple times --------------------------------------------------------------

SimStudy2_dgp_X <- function( n_obs = 40^2, # number of observations
                             m = 40 # number of level2-units 
){
  
  ## create level-2 id variable in balanced design
  l2.id <- rep(1:m, each = n_obs/m)
  
  ## number of fixed (linear) effects ( = number of predictor variables)
  p <- 4 # X1, X2, X3, X4
  
  ## Generation of X -----------------------------------------------------------
  
  X.l1_per_l2_half <- mvtnorm::rmvnorm(n = n_obs/m/2, 
                                       mean = rep(0, p), 
                                       sigma = diag(p),
                                       checkSymmetry = F)
  
  # symmetric duplication around mean( X.l1_per_l2_half) 
  X.l1_per_l2 <- rbind(X.l1_per_l2_half, 
                       -X.l1_per_l2_half + mean( X.l1_per_l2_half) ) 
  
  # whitening
  X.l1_per_l2.white <- whiten(X.l1_per_l2,center=T) 
  
  # complete X-matrix for all level-2 units
  X.l1 <- do.call(rbind, replicate(m, X.l1_per_l2.white, simplify = FALSE))
  
  # E(X1) = ... = E(X4) = 2 and Var(X1) = ... = Var(X4) = 1
  X.l1 <- X.l1/apply(X.l1,2, sd)
  X.l1 <- X.l1 + 2 - colMeans(X.l1)
  
  # add column of 1s for intercept 
  X <- cbind(1, X.l1) # add intercept-column
  colnames(X) <- c("intercept","X1","X2","X3","X4")
  
  # square all predictor variables 
  X.quadratic <- cbind(X, X.l1^2)
  colnames(X.quadratic) <- c("intercept","X1","X2","X3","X4",
                             "X1_squared","X2_squared","X3_squared","X4_squared")
  X.quadratic <- X.quadratic[,c("intercept",
                                "X1","X1_squared",
                                "X2","X2_squared",
                                "X3","X3_squared",
                                "X4","X4_squared")]
  
  ## Generation of Z -----------------------------------------------------------
  
  Z <- XzuZ(X = X,
            level2.index = l2.id,
            random.cols = c(1,3,5)) # u0, u3(X2), u7(X4)
  
  
  ## output
  return(list( l2.id = l2.id,
               n_obs = n_obs,
               m = m,
               X = X,
               X.quadratic = X.quadratic,
               Z = Z ))
  
  
}



## generate y ------------------------------------------------------------------

SimStudy2_dgp_y <- function( l2.id, # id variable for level-2 units
                             n_obs, # number of observations
                             m, # number of level-2 units
                             X, # X-matrix 
                             X.quadratic, # X-matrix (including squared variables)
                             Z, # design matrix corresponding to X
                             R_squared = 0.1, # pseudo coefficient of determination
                             D = matrix(c(2,0.2,0.1,
                                          0.2,1.5,0.2,
                                          0.1,0.2,1.5),
                                        nrow=3, byrow=T) # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             
                             
){
  
  ## number of random effects
  q <- ncol(D)
  
  ## draw random effects from multivariate normal distribution -----------------
  
  # random effects covariance matrix 
  G <- bdiag(replicate(m, D, simplify = F))
  
  ## draw  (u0j,u3j,u7j) from multivariate normal distribution 
  # u
  u_matrix <- mvtnorm::rmvnorm(n = m,
                               mean = rep(0, q), # E(u) = 0
                               sigma = D,
                               checkSymmetry = F)
  
  u <- as.vector(t(u_matrix)) 
  # (u01,u31,u71,...,u0j,u3j,u7j,u0j+1,u3j+1,u7j+1,...,u0m,u3m,u7m)
  
  
  ## get level-1 residuals -----------------------------------------------------
  # R_squared := Var(X * beta)/Var(Y) 
  # <=> Var(Y\epsilon) + sigma_square_epsilon = Var(X * beta)/R_squared
  # <=> sigma_square_epsilon = Var(X * beta)/R_squared - Var(Y\epsilon)
  #     = Var(X * beta)/R-squared - ( Var(X * beta) + tau0^2 + 
  #                                   [E(X2)^2* tau3^2 + tau3^2 *Var(X2)] + 
  #                                   [E(X4)^2* tau7^2 + tau7^2 *Var(X4)] +
  #                                   2* [E(X2) * tau03] + 2* [E(X4) * tau07] +
  #                                   2* [E(X2) * E(X4) * tau37]  )
  #
  # Var(X * beta) = var( X.quadratic %*% beta ) since X is fixed
  # X includes all quadratic predictors (i.e. X = X.quadratic)
  # KOV(X1,X2) = 0
  
  tau0_squared = diag(D)[1]
  tau3_squared = diag(D)[2]
  tau7_squared = diag(D)[3]
  tau03 = D[1,2]
  tau07 = D[1,3]
  tau37 = D[2,3]
  
  beta <- c(1,0,-0.5,0,0.5,-4,1,4,-1) # (c00,c10,c20,c30(rs),c40,c50,c60,c70(rs),c80)
  
  # var(X beta)
  Var_Xbeta <- var(X.quadratic %*% beta)
  
  sigma_squared_epsilon <- Var_Xbeta/R_squared - 
    (Var_Xbeta + tau0_squared + 
       (2^2* tau3_squared + tau3_squared *1) + 
       (2^2* tau7_squared + tau7_squared *1) +
       2* (2 * tau03) + 2* (2 * tau07) +
       2* (2 * 2 * tau37) )
  
  epsilon <- rnorm(n = n_obs,
                   mean = 0,
                   sd = sqrt(sigma_squared_epsilon))
  
  
  ## generate y ----------------------------------------------------------------
  # y = X.quadratic %*% beta + Z %*% u + epsilon
  y <-  X.quadratic %*% beta + Z %*% u + epsilon
  y <- as.vector(y)
  
  ## output --------------------------------------------------------------------
  
  return(list(y = y,
              X = X,
              Z = Z,
              G = G,
              sigma_squared_epsilon = sigma_squared_epsilon,
              l2.id = l2.id)
  )
  
  
}


