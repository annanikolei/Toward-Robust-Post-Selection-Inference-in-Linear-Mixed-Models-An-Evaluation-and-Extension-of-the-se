################################################################################
# Simulation Study 3-Discussion: Conditional Fixed RE (DGP)
################################################################################

library(Matrix)
library(whitening)

######################## dgp function ########################################## 

SimStudy3_dgp_u <- function(m = 40,
                            D = matrix(c(2,0.2,0.1,
                                         0.2,1.5,0.2,
                                         0.1,0.2,1.5),
                                       nrow=3, byrow=T),
                            seed){
  
  q <- ncol(D) # number of random effects
  
  set.seed(seed)
  u_matrix <- mvtnorm::rmvnorm(n = m,
                               mean = rep(0, q), # E(u) = 0
                               sigma = D,
                               checkSymmetry = F)
  
  u <- as.vector(t(u_matrix)) 
  

  mean_u <- colMeans(u_matrix)
  
  return(list(D=D,
              u =u,
              mean_u=mean_u)
  )
}


# u_fixed <- SimStudy3_dgp_u(m = 40,
#                 D = matrix(c(2,0.2,0.1,
#                              0.2,1.5,0.2,
#                              0.1,0.2,1.5),
#                            nrow=3, byrow=T),
#                 seed = 170) #check


## generate y ------------------------------------------------------------------

SimStudy3_dgp_y <- function( l2.id, # id variable for level-2 units
                             n_obs, # number of observations
                             m, # number of level-2 units
                             X, # X-matrix 
                             X.quadratic, # X-matrix (including squared variables)
                             Z, # design matrix corresponding to X
                             R_squared = 0.1, # pseudo coefficient of determination
                             u, # vector of fixed random effects
                             D # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                             
                             
){
  
  
  
  # random effects covariance matrix 
  G <- bdiag(replicate(m, D, simplify = F))
  

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
