################################################################################
# Simulation Study 3 under Random X: Results
################################################################################

message("Start this script from the repository root.")

## Install and load necessary packages
source("simulations/Project1_install_packages.R")

setwd("Selfmade")
source("selfmade**.R")

setwd("../simulations")

## Load Simulation Study 3 files
source("SimStudy3_helper.R")
source("RandomX_SimStudy3_dgp.R")
source("SimStudy3_SelectionProcess.R")
source("RandomX_SimStudy3_one_simulation.R")
source("RandomX_SimStudy3_simulation_study.R")



n_obs = 40^2 ; m = 40

RUN_SIMULATIONS <- FALSE 
results_dir <- "RandomX_SimStudy3_Results"

if (RUN_SIMULATIONS && !dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

if (!dir.exists(results_dir)) {
  stop("Results folder not found. Extract RandomX_SimStudy3_Results.zip first.")
}

setwd(results_dir)

if (RUN_SIMULATIONS) {
########################
# Run Simulations
########################
N <- 1000


##################################### R_squared = 0.3 ##########################
#######################
# CR2
#######################
R_squared = 0.3
#### X1234
save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_CR2_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag=TRUE,
                                 AdjMcCaffrey = TRUE,
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}


save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_marginal(noteff)CR2_AdjMcCaffreyFALSE_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag=TRUE,
                                 AdjMcCaffrey = FALSE,
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}

########################
# no_HC_correction
########################

## all perspectives no_HC_correction
#### X1234
save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}

########################
# HC3 correction
########################

#### X1234
N <- 1000

save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = FALSE, # (which) robust std.err in marginal perspective
                                 R_diag =F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}

# HC3 + standard
SimulationStudy3_no_HC_X1234 <- readRDS(paste0("RandomX_SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_",R_squared,"_1000MC.RDS"))
SimulationStudy3_HC3_X1234 <- readRDS(paste0("RandomX_SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_",R_squared,"_1000MC.RDS"))
SimulationStudy3 <- c(SimulationStudy3_no_HC_X1234, 
                      SimulationStudy3_HC3_X1234)
saveRDS(SimulationStudy3,paste0("RandomX_SimStudy3_workingModel=X1234_R2=",R_squared,"_1000MC.RDS"))


##################################### R_squared = "great" ######################

#######################
# CR2
#######################
#### X1234
R_squared = "great" # this is the condition with the small level-1 resiudal variance

save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_CR2_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 AdjMcCaffrey= TRUE,
                                 R_diag=TRUE,
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}

save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_marginal(noteff)CR2_AdjMcCaffreyFALSE_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = "CR2", # (which) robust std.err in marginal perspective
                                 R_diag=TRUE,
                                 AdjMcCaffrey = FALSE,
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}

########################
# no_HC_correction
########################

## all perspectives no_HC_correction
#### X1234
save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal","marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = F, # (which) robust std.err in marginal perspective
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}





########################
# HC3 correction
########################

#### X1234
N <- 1000

save_file <- paste0(
  "RandomX_SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_", R_squared, "_1000MC.RDS"
)


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue With Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      #suppressMessages(
      SimStudy3_simulation_study(n_iter = 1, #number of replications per condition
                                 n_obs = n_obs, # number of observations
                                 m = m, # number of level-2 units
                                 R_squared = R_squared, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),nrow=3, byrow=T), # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 perspective = c("marginal_not_efficient"), # marginal and/or conditional (for robust std. errors we only need marginal)
                                 robust_std.errors_marginal = FALSE, # (which) robust std.err in marginal perspective
                                 R_diag =F,
                                 HC_type = "HC3", # for marginal_diagF_not_efficient: (or "HC1","HC2","HC3")
                                 cadjust = T, # for marginal_diagF_not_efficient: logical. Should a cluster adjustment be applied?
                                 B_MCsamples = 1000, # number of MC samples in selective inference process
                                 remove_singular_fits = F,
                                 interesting_vars = c("X3","X4"), # for which variables do we want to do post-selection inference
                                 interesting_working_model_vars = c("X1","X2","X3","X4"),
                                 checkConvergence_for_check_Fun = T # convergence check for MC sampling?
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}

# HC3 + standard
SimulationStudy3_no_HC_X1234 <- readRDS(paste0("RandomX_SimStudy3_workingModel=X1234_no_HC_correction_allperspectives_R2_",R_squared,"_1000MC.RDS"))
SimulationStudy3_HC3_X1234 <- readRDS(paste0("RandomX_SimStudy3_workingModel=X1234_HC3_marginal(not_eff)_R2_",R_squared,"_1000MC.RDS"))
SimulationStudy3 <- c(SimulationStudy3_no_HC_X1234, 
                      SimulationStudy3_HC3_X1234)
saveRDS(SimulationStudy3,paste0("RandomX_SimStudy3_workingModel=X1234_R2=",R_squared,"_1000MC.RDS"))

}
################################################################################
## QQ-Plots, Convergence ##
################################################################################

R_squared_values <- c("0.3", "great")



# ============================================================
# Read simulation results for both R² conditions
# ============================================================

# -----------------------------
# Standard
# -----------------------------

SimulationStudy3_all <- lapply(R_squared_values, function(R_squared) {
  
  SimulationStudy3 <- readRDS(
    paste0(
      "RandomX_SimStudy3_workingModel=X1234_R2=",
      R_squared,
      "_1000MC.RDS"
    )
  )
  
  results <- do.call(
    rbind,
    lapply(
      SimulationStudy3,
      function(x) x$post_selection_pvalues
    )
  )
  
  results$R_squared <- R_squared
  
  results
})

SimStudy3_Results <- do.call(rbind, SimulationStudy3_all)


# -----------------------------
# CR2
# -----------------------------

SimulationStudy3_CR2_all <- lapply(R_squared_values, function(R_squared) {
  
  SimulationStudy3_CR2 <- readRDS(
    paste0(
      "RandomX_SimStudy3_workingModel=X1234_CR2_R2_",
      R_squared,
      "_1000MC.RDS"
    )
  )
  
  results <- do.call(
    rbind,
    lapply(
      SimulationStudy3_CR2,
      function(x) x$post_selection_pvalues
    )
  )
  
  results$R_squared <- R_squared
  
  results
})

SimStudy3_CR2_Results <- do.call(
  rbind,
  SimulationStudy3_CR2_all
)


SimStudy3_CR2_Results$perspective <- ifelse(
  SimStudy3_CR2_Results$perspective %in%
    c("marginal", "marginal_not_efficient", "naive"),
  paste0(
    SimStudy3_CR2_Results$perspective,
    "_CR2"
  ),
  SimStudy3_CR2_Results$perspective
)


# -----------------------------
# CR2 -- McCaffrey TRUE
# -----------------------------

SimulationStudy3_CR2_MCCaffreyFALSE_all <- lapply(
  R_squared_values,
  function(R_squared) {
    
    SimulationStudy3_CR2_MCCaffreyFALSE <- readRDS(
      paste0(
        "RandomX_SimStudy3_workingModel=X1234_marginal(noteff)CR2_AdjMcCaffreyFALSE_R2_",
        R_squared,
        "_1000MC.RDS"
      )
    )
    
    results <- do.call(
      rbind,
      lapply(
        SimulationStudy3_CR2_MCCaffreyFALSE,
        function(x) x$post_selection_pvalues
      )
    )
    
    results$R_squared <- R_squared
    
    results
  }
)

SimStudy3_CR2_Results_MCCaffreyFALSE <- do.call(
  rbind,
  SimulationStudy3_CR2_MCCaffreyFALSE_all
)


SimStudy3_CR2_Results_MCCaffreyFALSE$perspective <- ifelse(
  SimStudy3_CR2_Results_MCCaffreyFALSE$perspective %in%
    c("marginal", "marginal_not_efficient", "naive"),
  paste0(
    SimStudy3_CR2_Results_MCCaffreyFALSE$perspective,
    "_CR2_notMCCaffrey"
  ),
  SimStudy3_CR2_Results_MCCaffreyFALSE$perspective
)


# ============================================================
# Combine all results
# ============================================================

SimStudy3_Results <- rbind(
  SimStudy3_CR2_Results,
  SimStudy3_Results,
  SimStudy3_CR2_Results_MCCaffreyFALSE
)


SimStudy3_Results$pval_type <- ifelse(
  SimStudy3_Results$estimator %in% c("naive"),
  "Naive (CR2)",
  ifelse(
    SimStudy3_Results$estimator %in%
      c("fullCandidate", "fullCandidate_supplied"),
    "Standard (no selection)",
    "Selfmade"
  )
)


SimStudy3_Results$estimator <- ifelse(
  SimStudy3_Results$estimator == "naive",
  "ModelEstimate",
  SimStudy3_Results$estimator
)


SimStudy3_Results$perspective <- ifelse(
  SimStudy3_Results$estimator %in%
    c(
      "ModelEstimate_SIGMA_replaced",
      "ModelEstimate_R_replaced",
      "ICM_SIGMA_replaced",
      "ICM_R_replaced"
    ),
  "marginal_not_efficient_HC3",
  SimStudy3_Results$perspective
)


# ============================================================
# Perspectives
# ============================================================

new_persp <- c(
  "marginal",
  "marginal_not_efficient",
  "marginal_CR2",
  "marginal_not_efficient_CR2",
  "marginal_not_efficient_CR2_notMCCaffrey",
  "marginal_not_efficient_HC3",
  "conditional",
  "fullCandidate"
)


# -----------------------------
# Expand naive CR2 rows
# -----------------------------

naive_rows <- SimStudy3_Results[
  SimStudy3_Results$perspective == "naive_CR2",
]

naive_expanded <- naive_rows[
  rep(seq_len(nrow(naive_rows)), each = length(new_persp)),
]

naive_expanded$perspective <- rep(
  new_persp,
  times = nrow(naive_rows)
)


SimStudy3_Results <- SimStudy3_Results[
  SimStudy3_Results$perspective != "naive_CR2",
]

SimStudy3_Results <- rbind(
  naive_expanded,
  SimStudy3_Results
)


# -----------------------------
# Supplied
# -----------------------------

# supplied is the same if HC3/CR2 estimator is used => duplicate this

# Supplied:
# marginal_not_efficient -> marginal_not_efficient_HC3

supplied_rows <- SimStudy3_Results[
  SimStudy3_Results$estimator == "supplied" &
    SimStudy3_Results$perspective == "marginal_not_efficient",
]


supplied_expanded <- supplied_rows
supplied_expanded$perspective <- "marginal_not_efficient_HC3"


# Supplied:
# marginal_not_efficient -> marginal_not_efficient_CR2

supplied_expanded_CR2_noteff <- supplied_rows
supplied_expanded_CR2_noteff$perspective <-
  "marginal_not_efficient_CR2"


# Supplied:
# marginal_not_efficient -> marginal_not_efficient_CR2_notMCCaffrey

supplied_expanded_CR2_notMCCaffrey_noteff <- supplied_rows
supplied_expanded_CR2_notMCCaffrey_noteff$perspective <-
  "marginal_not_efficient_CR2_notMCCaffrey"


# Supplied:
# marginal -> marginal_CR2

supplied_rows_eff <- SimStudy3_Results[
  SimStudy3_Results$estimator == "supplied" &
    SimStudy3_Results$perspective == "marginal",
]


supplied_expanded_CR2_eff <- supplied_rows_eff
supplied_expanded_CR2_eff$perspective <- "marginal_CR2"


SimStudy3_Results <- rbind(
  SimStudy3_Results,
  supplied_expanded,
  supplied_expanded_CR2_noteff,
  supplied_expanded_CR2_notMCCaffrey_noteff,
  supplied_expanded_CR2_eff
)


# ============================================================
# Labels
# ============================================================

SimStudy3_Results$coefficient <-
  ifelse(
    SimStudy3_Results$coefficient == "X3",
    "L1-noise~x[3]",
    "L1-noise~(RS)~x[4]"
  )


library(ggplot2)
library(dplyr)
library(cowplot)


# ============================================================
# Prepare data
# ============================================================

results_all_conditions_df <- SimStudy3_Results %>%
  mutate(
    label_facet = paste0(
      "Selected Fixed = ",
      working_model
    )
  )


results_all_conditions_df$label_facet <- gsub(
  "Selected Fixed = X",
  "Sel. FE = X",
  results_all_conditions_df$label_facet
)


results_all_conditions_df$label_facet <- gsub(
  "X(?=[0-9])",
  "",
  results_all_conditions_df$label_facet,
  perl = TRUE
)


results_all_conditions_df$label_facet <- sub(
  "Sel\\. FE = ",
  "Sel. FE = X",
  results_all_conditions_df$label_facet
)


# ============================================================
# Perspective factor
# ============================================================

results_all_conditions_df$perspective <-
  factor(results_all_conditions_df$perspective)


results_all_conditions_df$perspective <- factor(
  results_all_conditions_df$perspective,
  
  levels = c(
    "conditional_uncorrected",
    "conditional_corrected",
    "marginal",
    "marginal_not_efficient",
    "marginal_CR2",
    "marginal_not_efficient_CR2",
    "marginal_not_efficient_CR2_notMCCaffrey",
    "marginal_not_efficient_HC3",
    "fullCandidate"
  ),
  
  labels = c(
    "Conditional\n(uncorrected)",
    "Conditional\n(corrected)",
    "Marginal\n(efficient)",
    "Marginal\n(not eff)",
    "Marginal\n(eff) CR2",
    "Marginal\n(not eff) CR2_MCCaffrey",
    "Marginal\n(not eff) CR2",
    "Marginal\n(not eff) HC3",
    "fullCandidate"
  )
)


# ============================================================
# Estimator factor
# ============================================================

results_all_conditions_df$estimator <-
  factor(results_all_conditions_df$estimator)


results_all_conditions_df$estimator <- factor(
  results_all_conditions_df$estimator,
  
  levels = c(
    "fullCandidate_supplied",
    "supplied",
    "ModelEstimate",
    "ICM",
    "VarY",
    "ModelEstimate_SIGMA_replaced",
    "ModelEstimate_R_replaced",
    "ICM_SIGMA_replaced",
    "ICM_R_replaced",
    "fullCandidate"
  ),
  
  labels = c(
    "fullCandidate_supplied",
    "supplied",
    "ModelEstimate",
    "ICM",
    "VarY",
    "ModelEstimate_SIGMA",
    "ModelEstimate_R",
    "ICM_SIGMA",
    "ICM_R",
    "fullCandidate"
  )
)


results_all_conditions_df$pval_type <-
  relevel(
    as.factor(results_all_conditions_df$pval_type),
    "Selfmade"
  )


results_all_conditions_df$coefficient <- factor(
  results_all_conditions_df$coefficient,
  
  levels = c(
    "L1-noise~x[3]",
    "L1-noise~(RS)~x[4]"
  )
)


# ============================================================
# Manual color palette
# ============================================================

estimator_colors <- c(
  "fullCandidate_supplied" = "black",
  "supplied" = "#009E73",
  "ModelEstimate" = "#8C510A",
  "ICM" = "#D55E00",
  "VarY" = "#0072B2",
  "ModelEstimate_SIGMA" = "#56B4E9",
  "ModelEstimate_R" = "#E69F00",
  "ICM_SIGMA" = "#CC79A7",
  "ICM_R" = "#6A3D9A",
  "fullCandidate" = "white"
)


# ============================================================
# Data for Simulation Study 3
# ============================================================

results_all_conditions_df_SimStudy3 <-
  results_all_conditions_df[
    results_all_conditions_df$perspective %in%
      c(
        "Marginal\n(efficient)",
        "Marginal\n(not eff)",
        "Marginal\n(eff) CR2",
        "Marginal\n(not eff) CR2",
        "Marginal\n(not eff) HC3",
        "Conditional\n(corrected)"
      ),
  ]


# ============================================================
# R² factor
# ============================================================

results_all_conditions_df_SimStudy3$R_squared <-
  factor(
    results_all_conditions_df_SimStudy3$R_squared,
    levels = c("0.3", "great"),
    labels = c("R² = 0.3", "R² = great")
  )


# ============================================================
# Plot
# ============================================================

plot_SimStudy3 <- ggplot(
  results_all_conditions_df_SimStudy3,
  
  aes(
    sample = pval_two,
    color = estimator,
    linetype = pval_type
  )
) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "grey"
  ) +
  
  geom_qq(
    distribution = stats::qunif,
    size = 1.2,
    geom = "line"
  ) +
  
  facet_grid(
    perspective ~
      R_squared +
      label_facet +
      coefficient,
    
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  ) +
  
  theme_minimal(base_size = 20) +
  
  theme(
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    
    strip.text = element_text(
      color = "black",
      size = 24,
      lineheight = 0.8
    ),
    
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    axis.text = element_text(size = 19),
    axis.title = element_text(size = 24)
  ) +
  
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
  ) +
  
  scale_color_manual(
    values = estimator_colors,
    
    labels = function(x) parse(
      text = c(
        supplied = "supplied",
        ModelEstimate = "ModelEstimate",
        ICM = "ICM",
        VarY = "VarY",
        ModelEstimate_SIGMA = "ModelEstimate[Sigma]",
        ICM_SIGMA = "ICM[Sigma]",
        ModelEstimate_R = "ModelEstimate[R]",
        ICM_R = "ICM[R]"
      )[x]
    )
  )


# ============================================================
# Show plot
# ============================================================

plot_SimStudy3


# ============================================================
# Save plot
# ============================================================

ggsave(
  "plot_RandomX_SimStudy3_R2_both.png",
  plot = plot_SimStudy3,
  dpi = 300,
  width = 21,
  height = 21,
  units = "in"
)

