################################################################################
# Random X Simulation Study 3: Data generating process (quadratic dgp)
################################################################################

library(Matrix)

######################## dgp function ########################################## 

SimStudy3_dgp<- function( n_obs = 40^2, # number of observations
                          m = 40, # number of level2-units
                          R_squared = 0.3, # pseudo coefficient of determination
                          D = matrix(c(2,0.2,0.1,
                                       0.2,1.5,0.2,
                                       0.1,0.2,1.5),
                                     nrow=3, byrow=T) # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
){
  
  ## balanced design
  l2.id <- rep(1:m, each = n_obs/m)
  
  ## number of fixed (linear) effects ( = number of predictor variables)
  p <- 4 # X1,X2,X3,X4
  
  ## number of random effects
  q <- ncol(D)
  
  
  ## draw random effects and X from multivariate normal distribution -----------
  
  # covariance matrix of predictor variables X1 and X2 
  # COV(Xi,Xj) = 0 
  # var(X1) = ... = var(X4) = 1
  vcovMatrix_X <- diag(p) 
  
  # random effects covariance matrix 
  G <- bdiag(replicate(m, D, simplify = F))
  
  ## draw (X1ij,...,X4ij) and (u0j,u3j,u7j) from (joint) multivariate normal 
  ## distribution 
  # X level-1 
  X.l1 <- mvtnorm::rmvnorm(n = n_obs, 
                           mean = rep(2, p), # E(X1) = ... = E(X4) = 2
                           sigma = vcovMatrix_X,
                           checkSymmetry = F)
  X <- cbind(1, X.l1) # add intercept-column
  colnames(X) <- c("intercept","X1","X2","X3","X4")
  
  X.quadratic <- cbind(X, X.l1^2)
  colnames(X.quadratic) <- c("intercept","X1","X2","X3","X4",
                             "X1_squared","X2_squared","X3_squared","X4_squared")
  X.quadratic <- X.quadratic[,c("intercept",
                                "X1","X1_squared",
                                "X2","X2_squared",
                                "X3","X3_squared",
                                "X4","X4_squared")]
  
  # u
  u_matrix <- mvtnorm::rmvnorm(n = m,
                               mean = rep(0, q), # E(u) = 0
                               sigma = D,
                               checkSymmetry = F)
  
  u <- as.vector(t(u_matrix)) 
  # (u01,u31,u71,...,u0j,u3j,u7j,u0j+1,u3j+1,u7j+1,...,u0m,u3m,u7m)
  
  
  ## use random effects structure to generate Z --------------------------------
  Z <- XzuZ(X = X,
            level2.index = l2.id,
            random.cols = c(1,3,5)) # u0,u3(X2),u7(X4)
  
  
  ## get level-1 residuals -----------------------------------------------------
  # R_squared := Var(X * beta)/Var(Y) 
  # <=> Var(Y\sigma_square_epsilon) + sigma_square_epsilon = Var(X * beta)/R_squared
  # <=> sigma_square_epsilon = Var(X * beta)/R_squared - Var(Y\sigma_square_epsilon)
  #     = Var(X * beta)/R-squared - ( Var(X * beta) + tau0^2 + 
  #                                   [E(X2)^2* tau3^2 + tau3^2 *Var(X2)] + 
  #                                   [E(X4)^2* tau7^2 + tau7^2 *Var(X4)] +
  #                                   2* [E(X2) * tau03] + 2* [E(X4) * tau07] +
  #                                   2* [E(X2) * E(X4) * tau37]  )
  #
  # Var(X * beta) = c20^2 * [4*E(X1)^2 * Var(X1) + 2* Var(X1)^2] + 
  #                 c40^2 * [4*E(X2)^2 * Var(X2) + 2* Var(X2)^2] +
  #                 c50^2 * Var(X3) + c70^2 * Var(X4) +
  #                 c60^2 * [4*E(X3)^2 * Var(X3) + 2* Var(X3)^2] + 
  #                 c80^2 * [4*E(X4)^2 * Var(X4) + 2* Var(X4)^2] +
  #                 2 *c50 *c60 * [2*E(X3)*Var(X3)] +
  #                 2 *c70 *c80 * [2*E(X4)*Var(X4)] 
  # X includes all quadratic predictors (i.e. X = X.quadratic)
  
  tau0_squared = diag(D)[1]
  tau3_squared = diag(D)[2]
  tau7_squared = diag(D)[3]
  tau03 = D[1,2]
  tau07 = D[1,3]
  tau37 = D[2,3]
  
  beta <- c(0,0,-0.5,0,0.5,-4,1,4,-1) # (c00,c10,c20,c30,c40,c50,c60,c70,c80)
  c20 = beta[3]
  c40 = beta[5]
  c50 = beta[6]
  c60 = beta[7]
  c70 = beta[8]
  c80 = beta[9]
  
  # E(X1) = ... = E(X4) = 2; Var(X1) = ... = Var(X4) = 1
  Var_Xbeta <- c20^2 * (4*2^2 * 1 + 2* 1^2) + 
    c40^2 * (4*2^2 * 1 + 2* 1^2) +
    c50^2 * 1 + c70^2 * 1 +
    c60^2 * (4*2^2 * 1 + 2* 1^2) + 
    c80^2 * (4*2^2 * 1 + 2* 1^2) +
    2 *c50 *c60 * (2*2*1) +
    2 *c70 *c80 * (2*2*1)
  
  if(R_squared != "great"){
  sigma_squared_epsilon <- Var_Xbeta/R_squared - 
    (Var_Xbeta + tau0_squared + (2^2* tau3_squared + tau3_squared *1) + 
       (2^2* tau7_squared + tau7_squared *1) +
       2* (2 * tau03) + 2* (2 * tau07) +
       2* (2 * 2 * tau37) )}
  
  if(R_squared == "great"){
  sigma_squared_epsilon <-  Var_Xbeta/4 }
  
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
