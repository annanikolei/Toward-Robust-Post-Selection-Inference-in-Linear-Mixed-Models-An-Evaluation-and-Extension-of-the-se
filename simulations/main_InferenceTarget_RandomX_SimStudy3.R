# ==============================================================================
# Simulation Study 3 RandomX -- Monte Carlo estimation of the OLS and GLS inference
# targets for the linear working model (X1,X2,X3,X4), based on the DGM of
# Simulation Study 3, but with X redrawn (iid per observation) in every
# iteration and b redrawn independently as well.
# ==================

# Start this script from the repository root.
setwd("Selfmade")
source("selfmade**.R")

setwd("../simulations")

## Install and load all necessary packages
source("Project1_install_packages.R")

# Load SimStudy3-RandomX files
source("SimStudy3_helper.R")
source("RandomX_SimStudy3_dgp.R")
source("SimStudy3_SelectionProcess.R")
source("RandomX_SimStudy3_one_simulation.R")
source("RandomX_SimStudy3_simulation_study.R")


# ------------------------------------------------------------------------------
# one_simulation
# ------------------------------------------------------------------------------

one_simulation_RandomX <- function(n_obs = 40^2,
                                   m = 40,
                                   R_squared = 0.3,
                                   D = matrix(c(2,0.2,0.1,
                                                0.2,1.5,0.2,
                                                0.1,0.2,1.5),
                                              nrow = 3, byrow=T),
                                   perspective = "marginal"){
  
  # data.frame 
  dgp_output <- SimStudy3_dgp(n_obs = n_obs,
                              m = m, 
                              R_squared = R_squared, 
                              D = D )
  
  dat <- data.frame(y = dgp_output$y, dgp_output$X, l2.id = dgp_output$l2.id)
  
  mod <- lmer(y ~ X1 + X2 + X3 + X4 +(1|l2.id), data = dat, REML=FALSE)
  
  checkFun <- function(yb){TRUE}
  
  selfmade_output <-
    mocasin(mod = mod, 
            checkFun = checkFun,
            checkConvergence_for_check_Fun = FALSE,
            checkFun_checkConv = NULL,
            this_y = dat$y, 
            B_MCsamples = 1, # number of MC-samples
            MC_distribution = "centered",
            which_marginal = 4:5, 
            rho_H0 = 0, 
            perspective = perspective, 
            robust_std.errors_marginal = FALSE,
            marginal_COV_estimator = c("ModelEstimate","supplied"), 
            marginal_not_efficient_COV_estimator = c("ModelEstimate","supplied"),
            diagF_not_efficient_SIGMA_or_R = c("SIGMA_replaced", "R_replaced"),
            symmetric_selection = T, 
            alpha = 0.05, 
            allow_smaller_alpha = T, 
            multiple = 7, 
            R_diag = TRUE, 
            Sigma_supplied = NULL,
            Z_supplied = dgp_output$Z,
            G_supplied = dgp_output$G,
            Lambda_supplied = NULL,
            Lambdat_supplied = NULL,
            G_work_supplied = NULL, 
            sigma_square_epsilon_supplied = 
              as.numeric(dgp_output$sigma_squared_epsilon),
            R_supplied = Matrix(diag(as.numeric(dgp_output$sigma_squared_epsilon),n_obs), sparse = T)) 
  
  
  if(perspective == "marginal"){
  df <- data.frame(
    perspective = rep("marginal(eff)",4),
    estimator = c("ModelEstimate", "ModelEstimate", "supplied", "supplied"),
    coefficient = c(4, 5, 4, 5),
    t_obs = c(
      selfmade_output$marginal[[1]][[1]][[1]]$t_obs,
      selfmade_output$marginal[[1]][[2]][[1]]$t_obs,
      selfmade_output$marginal[[2]][[1]][[1]]$t_obs,
      selfmade_output$marginal[[2]][[2]][[1]]$t_obs
    )
  )}
  
  if(perspective == "marginal_not_efficient"){
    df <- data.frame(
      perspective = rep("marginal(not_eff)",4),
      estimator = c("ModelEstimate", "ModelEstimate", "supplied", "supplied"),
      coefficient = c(4, 5, 4, 5),
      t_obs = c(
        selfmade_output$marginal_not_efficient[[1]][[1]][[1]]$t_obs,
        selfmade_output$marginal_not_efficient[[1]][[2]][[1]]$t_obs,
        selfmade_output$marginal_not_efficient[[2]][[1]][[1]]$t_obs,
        selfmade_output$marginal_not_efficient[[2]][[2]][[1]]$t_obs
      )
    )}
  
  return(df)
}

setwd("Results_InferenceTarget_RandomX_SimStudy3")


results_10000 <- vector("list", 10000)

for (i in 1:10000) {
  results_10000[[i]] <- one_simulation_RandomX()
  results_10000[[i]]$simulation <- i
}

results_10000_marginaleff <- do.call(rbind, results_10000)
saveRDS(results_10000_marginaleff, "10000_marginaleff.RDS")


results_10000_noteff <- vector("list", 10000)

for (i in 1:10000) {
  results_10000_noteff[[i]] <- one_simulation_RandomX(perspective = "marginal_not_efficient")
  results_10000_noteff[[i]]$simulation <- i
}

results_10000_marginalnoteff<- do.call(rbind, results_10000_noteff)
saveRDS(results_10000_marginalnoteff, "10000_marginalnoteff.RDS")



results_10000_marginaleff <- readRDS("10000_marginaleff.RDS")
library(dplyr)

mean_t_obs_eff <- results_10000_marginaleff %>%
  group_by(estimator, coefficient) %>%
  summarise(
    mean_t_obs = mean(t_obs, na.rm = TRUE),
    .groups = "drop"
  )

sd_t_obs_eff <- results_10000_marginaleff %>%
  group_by(estimator, coefficient) %>%
  summarise(
    sd_t_obs = sd(t_obs, na.rm = TRUE),
    .groups = "drop"
  )

results_10000_marginalnoteff <- readRDS("10000_marginalnoteff.RDS")

mean_t_obs_noteff <- results_10000_marginalnoteff %>%
  group_by(estimator, coefficient) %>%
  summarise(
    mean_t_obs = mean(t_obs, na.rm = TRUE),
    .groups = "drop"
  )

sd_t_obs_noteff <- results_10000_marginalnoteff %>%
  group_by(estimator, coefficient) %>%
  summarise(
    sd_t_obs = sd(t_obs, na.rm = TRUE),
    .groups = "drop"
  )

mean_t_obs_eff
sd_t_obs_eff
mean_t_obs_noteff
sd_t_obs_noteff/sd_t_obs_eff
