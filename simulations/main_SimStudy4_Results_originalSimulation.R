# Start this script from the repository root.
setwd("Selfmade")
source("selfmade**.R")

setwd("../simulations")

## Install and load all necessary packages
source("Project1_install_packages.R")

## Load Simulation Study 4 files
source("SimStudy4_helper.R")
source("SimStudy4_dgp.R")
source("SimStudy4_SelectionProcess.R")
source("SimStudy4_one_simulation.R")
source("SimStudy4_simulation_study.R")


R_squared = 0.3
n_obs = 40^2 ; m = 40

# load fixed predictor matriced
dgp_X_output <- readRDS("dgp_X_output_Sim4_n_obs=1600_m=40.RDS")

# load list_supplied_D_work
load(paste0("SimulationResults_supplied_D/SIM4_list_supplied_D_work_fullrandom=(Intercept)X2X4_n=400_Rsquared=",R_squared,"_nrsamples=10000.RData"))



################################################################################
## Run the Simulation ##
################################################################################
setwd("SimStudy4_Results")

# Number of planned iterations
N <- 1000


########################
# CR2 standard errors
########################

save_file <- paste0(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal_CR2_n=40_R2_", R_squared, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Weiter bei Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy4_simulation_study(
        n_iter = 1, 
        n_obs = n_obs,
        m = m,
        R_squared = R_squared,
        dgp_X_output = dgp_X_output,
        D = matrix(c(2,0.2,0.1,
                     0.2,1.5,0.2,
                     0.1,0.2,1.5), nrow=3, byrow=TRUE),
        perspective = c("marginal_not_efficient","marginal"),
        R_diag = TRUE,
        MC_distribution = "centered", # or proposal, does not matter
        HC_type = "HC3",
        cadjust = TRUE,
        robust_std.errors_marginal = "CR2",
        B_MCsamples = 1000,
        remove_singular_fits = FALSE,
        list_supplied_D_work = list_supplied_D_work,
        interesting_working_model_ranef = c("(Intercept)","X2"),
        checkConvergence_for_check_Fun = TRUE
      )
    },
    error = function(e) {
      warning("error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}

########################
# normal standard errors
########################

save_file <- paste0(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal(not_eff)_RdiagTRUE_n=40_R2_", R_squared, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Weiter bei Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy4_simulation_study(
        n_iter = 1, 
        n_obs = n_obs,
        m = m,
        R_squared = R_squared,
        dgp_X_output = dgp_X_output,
        D = matrix(c(2,0.2,0.1,
                     0.2,1.5,0.2,
                     0.1,0.2,1.5), nrow=3, byrow=TRUE),
        perspective = c("marginal_not_efficient"),
        R_diag = TRUE,
        MC_distribution = "centered", # or proposal, does not matter
        HC_type = "HC3",
        cadjust = TRUE,
        robust_std.errors_marginal = FALSE,
        B_MCsamples = 1000,
        remove_singular_fits = FALSE,
        list_supplied_D_work = list_supplied_D_work,
        interesting_working_model_ranef = c("(Intercept)","X2"),
        checkConvergence_for_check_Fun = TRUE
      )
    },
    error = function(e) {
      warning("error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}


save_file <- paste0(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal(eff)_RdiagTRUE_n=40_R2_", R_squared, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Weiter bei Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy4_simulation_study(
        n_iter = 1, 
        n_obs = n_obs,
        m = m,
        R_squared = R_squared,
        dgp_X_output = dgp_X_output,
        D = matrix(c(2,0.2,0.1,
                     0.2,1.5,0.2,
                     0.1,0.2,1.5), nrow=3, byrow=TRUE),
        perspective = c("marginal"),
        R_diag = TRUE,
        MC_distribution = "centered", # or proposal, does not matter
        HC_type = "HC3",
        cadjust = TRUE,
        robust_std.errors_marginal = FALSE,
        B_MCsamples = 1000,
        remove_singular_fits = FALSE,
        list_supplied_D_work = list_supplied_D_work,
        interesting_working_model_ranef = c("(Intercept)","X2"),
        checkConvergence_for_check_Fun = TRUE
      )
    },
    error = function(e) {
      warning("error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}


########################
# HC3 correction
########################

save_file <- paste0(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal_not_efficient_RdiagFALSE_n=40_R2_", R_squared, "_1000MC.RDS"
)

if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Weiter bei Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy4_simulation_study(
        n_iter = 1, 
        n_obs = n_obs,
        m = m,
        R_squared = R_squared,
        dgp_X_output = dgp_X_output,
        D = matrix(c(2,0.2,0.1,
                     0.2,1.5,0.2,
                     0.1,0.2,1.5), nrow=3, byrow=TRUE),
        perspective = c("marginal_not_efficient"),
        R_diag = FALSE,
        MC_distribution = "centered", # or proposal, does not matter
        HC_type = "HC3",
        cadjust = TRUE,
        robust_std.errors_marginal = FALSE,
        B_MCsamples = 1000,
        remove_singular_fits = FALSE,
        list_supplied_D_work = list_supplied_D_work,
        interesting_working_model_ranef = c("(Intercept)","X2"),
        checkConvergence_for_check_Fun = TRUE
      )
    },
    error = function(e) {
      warning("error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
 
  saveRDS(results_list, save_file)
}

################################################################################
## QQ-Plots, Convergence ##
################################################################################
setwd("SimStudy4_Results")

# ============================================================
# Read simulation results
# ============================================================

SimulationStudy4_marginal_not_efficient_RdiagTRUE <- readRDS(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal(not_eff)_RdiagTRUE_n=40_R2_0.3_1000MC.RDS"
)

SimulationStudy4_marginal_not_efficient_RdiagFALSE <- readRDS(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal_not_efficient_RdiagFALSE_n=40_R2_0.3_1000MC.RDS"
)

SimulationStudy4_marginal_eff <- readRDS(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal(eff)_RdiagTRUE_n=40_R2_0.3_1000MC.RDS"
)


SimulationStudy4 <- c(
  SimulationStudy4_marginal_not_efficient_RdiagTRUE,
  SimulationStudy4_marginal_not_efficient_RdiagFALSE,
  SimulationStudy4_marginal_eff
)


SimStudy4_Results <- do.call(
  rbind,
  lapply(
    SimulationStudy4,
    function(x) x$post_selection_pvalues
  )
)


# ============================================================
# CR2 results
# ============================================================

SimulationStudy4_marginalCR2 <- readRDS(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal_CR2_n=40_R2_0.3_1000MC.RDS"
)


SimStudy4_marginalCR2_Results <- do.call(
  rbind,
  lapply(
    SimulationStudy4_marginalCR2,
    function(x) x$post_selection_pvalues
  )
)


# Add "_CR2" to CR2 perspectives
SimStudy4_marginalCR2_Results$perspective <- ifelse(
  SimStudy4_marginalCR2_Results$perspective %in%
    c(
      "marginal",
      "marginal_not_efficient",
      "naive"
    ),
  paste0(
    SimStudy4_marginalCR2_Results$perspective,
    "_CR2"
  ),
  SimStudy4_marginalCR2_Results$perspective
)


# ICM should not be used when random effects are selected
# => also exclude ICM from CR2
SimStudy4_marginalCR2_Results <-
  SimStudy4_marginalCR2_Results[
    SimStudy4_marginalCR2_Results$estimator != "ICM",
  ]


# ============================================================
# Combine standard and CR2 results
# ============================================================

SimStudy4_Results <- rbind(
  SimStudy4_Results,
  SimStudy4_marginalCR2_Results
)


# ============================================================
# P-value type
# ============================================================

SimStudy4_Results$pval_type <- ifelse(
  SimStudy4_Results$perspective == "naive_CR2",
  "Naive (CR2)",
  ifelse(
    SimStudy4_Results$estimator == "naive",
    "Naive",
    ifelse(
      SimStudy4_Results$estimator %in%
        c(
          "fullCandidate",
          "fullCandidate_supplied"
        ),
      "Standard (no selection)",
      "Selfmade"
    )
  )
)


# ============================================================
# Rename naive estimator
# ============================================================

SimStudy4_Results$estimator <- ifelse(
  SimStudy4_Results$estimator == "naive",
  "ModelEstimate",
  SimStudy4_Results$estimator
)


# ============================================================
# HC3 estimators
# ============================================================

SimStudy4_Results$perspective <- ifelse(
  SimStudy4_Results$estimator %in%
    c(
      "ModelEstimate_SIGMA_replaced",
      "ModelEstimate_R_replaced",
      "ICM_SIGMA_replaced",
      "ICM_R_replaced"
    ),
  "marginal_not_efficient_HC3",
  SimStudy4_Results$perspective
)


# ============================================================
# Expand naive results to ALL perspectives
# ============================================================

new_persp <- c(
  "marginal",
  "marginal_not_efficient",
  "marginal_CR2",
  "marginal_not_efficient_CR2",
  "marginal_not_efficient_HC3",
  "conditional",
  "fullCandidate"
)


# ------------------------------------------------------------
# Standard naive
# ------------------------------------------------------------

naive_rows <- SimStudy4_Results[
  SimStudy4_Results$perspective == "naive",
]


naive_expanded <- naive_rows[
  rep(
    seq_len(nrow(naive_rows)),
    each = length(new_persp)
  ),
]


naive_expanded$perspective <- rep(
  new_persp,
  times = nrow(naive_rows)
)


naive_expanded$pval_type <- "Naive"


# ------------------------------------------------------------
# Naive CR2
# ------------------------------------------------------------

naive_CR2_rows <- SimStudy4_Results[
  SimStudy4_Results$perspective == "naive_CR2",
]


naive_CR2_expanded <- naive_CR2_rows[
  rep(
    seq_len(nrow(naive_CR2_rows)),
    each = length(new_persp)
  ),
]


naive_CR2_expanded$perspective <- rep(
  new_persp,
  times = nrow(naive_CR2_rows)
)


naive_CR2_expanded$pval_type <- "Naive (CR2)"


# ------------------------------------------------------------
# Remove original naive rows
# ------------------------------------------------------------

SimStudy4_Results <- SimStudy4_Results[
  !SimStudy4_Results$perspective %in%
    c(
      "naive",
      "naive_CR2"
    ),
]


# ------------------------------------------------------------
# Add expanded naive results
# ------------------------------------------------------------

SimStudy4_Results <- rbind(
  SimStudy4_Results,
  naive_expanded,
  naive_CR2_expanded
)


# ============================================================
# Supplied results
# ============================================================

# Supplied results are the same if HC3 or CR2 estimator is used
# => duplicate these


# ------------------------------------------------------------
# Supplied: marginal_not_efficient -> HC3
# ------------------------------------------------------------

supplied_rows <- SimStudy4_Results[
  SimStudy4_Results$estimator == "supplied" &
    SimStudy4_Results$perspective == "marginal_not_efficient",
]


supplied_expanded_HC3 <- supplied_rows

supplied_expanded_HC3$perspective <-
  "marginal_not_efficient_HC3"


# ------------------------------------------------------------
# Supplied: marginal_not_efficient -> CR2
# ------------------------------------------------------------

supplied_expanded_CR2_noteff <- supplied_rows

supplied_expanded_CR2_noteff$perspective <-
  "marginal_not_efficient_CR2"


# ------------------------------------------------------------
# Supplied: marginal -> CR2
# ------------------------------------------------------------

supplied_rows_eff <- SimStudy4_Results[
  SimStudy4_Results$estimator == "supplied" &
    SimStudy4_Results$perspective == "marginal",
]


supplied_expanded_CR2_eff <- supplied_rows_eff

supplied_expanded_CR2_eff$perspective <-
  "marginal_CR2"


# ------------------------------------------------------------
# Add supplied copies
# ------------------------------------------------------------

SimStudy4_Results <- rbind(
  SimStudy4_Results,
  supplied_expanded_HC3,
  supplied_expanded_CR2_noteff,
  supplied_expanded_CR2_eff
)


# ============================================================
# Coefficient labels
# ============================================================

SimStudy4_Results$coefficient <-
  ifelse(
    SimStudy4_Results$coefficient == "X3",
    "L1-noise~x[3]",
    "L1-noise~(RS)~x[4]"
  )


# ============================================================
# Libraries
# ============================================================

library(ggplot2)
library(dplyr)


# ============================================================
# Prepare data
# ============================================================

results_all_conditions_df <- SimStudy4_Results %>%
  mutate(
    label_facet = paste0(
      "Selected Random Effects = ",
      working_model
    )
  )


results_all_conditions_df$label_facet <- gsub(
  "Selected Random Effects =X",
  "Sel. RE = X",
  results_all_conditions_df$label_facet
)


results_all_conditions_df$label_facet <- gsub(
  "X(?=[0-9])",
  "",
  results_all_conditions_df$label_facet,
  perl = TRUE
)


results_all_conditions_df$label_facet <- sub(
  "Sel\\. RE = ",
  "Sel. RE = X",
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
    "marginal",
    "marginal_not_efficient",
    "marginal_CR2",
    "marginal_not_efficient_CR2",
    "marginal_not_efficient_HC3",
    "fullCandidate"
  ),
  
  labels = c(
    "Marginal\n(eff)",
    "Marginal\n(not eff)",
    "Marginal\n(eff) CR2",
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


# ============================================================
# P-value type factor
# ============================================================

results_all_conditions_df$pval_type <-
  relevel(
    as.factor(results_all_conditions_df$pval_type),
    "Selfmade"
  )


# ============================================================
# Coefficient factor
# ============================================================

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
# Data for Simulation Study 4 plot
# ============================================================

results_all_conditions_df_SimStudy4 <-
  results_all_conditions_df[
    results_all_conditions_df$perspective %in%
      c(
        "Marginal\n(eff)",
        "Marginal\n(not eff)",
        "Marginal\n(eff) CR2",
        "Marginal\n(not eff) CR2",
        "Marginal\n(not eff) HC3"
      ),
  ]


# ============================================================
# Plot
# ============================================================

plot_SimStudy4 <- ggplot(
  
  results_all_conditions_df_SimStudy4,
  
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
  
  theme_minimal(base_size = 21) +
  
  theme(
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    
    strip.text = element_text(
      color = "black",
      size = 21,
      lineheight = 0.9
    ),
    
    title = element_text(size = 21), 
    legend.text = element_text(size = 21),
    legend.title = element_text(size = 21),
    axis.text = element_text(size = 16)
  ) +
  
  geom_qq(
    distribution = stats::qunif,
    size = 1.2,
    geom = "line"
  ) +
  
  facet_grid(
   coefficient  ~perspective,
    
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  ) +
  
  labs(
    title = "Selected Random Effects = (Intercept) X2",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
  ) +
  
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  
  scale_color_manual(
    values = estimator_colors,
    
    labels = function(x) parse(
      text = c(
        fullCandidate_supplied = "fullCandidate (supplied)",
        supplied = "supplied",
        ModelEstimate = "ModelEstimate",
        ICM = "ICM",
        VarY = "VarY",
        ModelEstimate_SIGMA = "ModelEstimate[Sigma]",
        ModelEstimate_R = "ModelEstimate[R]",
        ICM_SIGMA = "ICM[Sigma]",
        ICM_R = "ICM[R]"
      )[x]
    )
  )


# ============================================================
# Display and save
# ============================================================

plot_SimStudy4

results_all_conditions_df_SimStudy4[
  which(is.na(results_all_conditions_df_SimStudy4$pval_two)),
] # Identify results for which no MC samples survived the inference-region check


ggsave(
  "plot_SimStudy4.png",
  plot = plot_SimStudy4,
  dpi = 300,
  width = 18,
  height = 7.5,
  units = "in"
)


# ============================================================
# Convergence issues
# ============================================================

# Initial model fitting

SimStudy4_convergence <- do.call(
  rbind,
  lapply(
    SimulationStudy4,
    function(x) x$convergence_singular_info
  )
)


colMeans(
  SimStudy4_convergence[
    sapply(SimStudy4_convergence, is.numeric)
  ],
  na.rm = F
)


29.037 / (84.4986666667 + 29.037)

0.0003333333/84.4986666667 


# ============================================================
# Monte Carlo resampling process
# ============================================================

plot_df <- results_all_conditions_df %>%
  filter(
    perspective %in% c(
      "Marginal\n(eff)",
      "Marginal\n(not eff)",
      "Marginal\n(eff) CR2",
      "Marginal\n(not eff) CR2",
      "Marginal\n(not eff) HC3"
    ),
    pval_type != "Naive",
    pval_type != "Naive (CR2)"
  ) %>%
  mutate(
    estimator = factor(estimator),
    pval_type = factor(pval_type),
    perspective = factor(perspective),
    label_facet = factor(label_facet),
    coefficient = factor(coefficient)
  )


plot_MCResampling_SimStudy4 <-
  ggplot(
    plot_df,
    
    aes(
      x = "",
      y = total_resamplings_due_to_convergence,
      color = estimator,
      group = interaction(
        pval_type,
        estimator
      )
    )
  ) +
  
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.alpha = 0.25,
    fill = NA
  ) +
  
  facet_grid(
    coefficient  ~perspective,
    
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )  +
  
  labs(
    x = "",
    y = "Number of MC Resamplings due to convergence",
    color = "Estimator"
  ) +
  
  theme_minimal(base_size = 20) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    
    strip.text = element_text(
      color = "black",
      size = 19
    ),
    
    legend.position = "right",
    panel.spacing.y = unit(1.5, "lines"),
    panel.spacing.x = unit(1, "lines")
  ) +
  
  scale_color_manual(
    values = estimator_colors,
    
    labels = function(x) parse(
      text = c(
        fullCandidate_supplied = "fullCandidate (supplied)",
        supplied = "supplied",
        ModelEstimate = "ModelEstimate",
        ICM = "ICM",
        VarY = "VarY",
        ModelEstimate_SIGMA = "ModelEstimate[Sigma]",
        ModelEstimate_R = "ModelEstimate[R]",
        ICM_SIGMA = "ICM[Sigma]",
        ICM_R = "ICM[R]",
        fullCandidate = "fullCandidate"
      )[x]
    )
  )


plot_MCResampling_SimStudy4


ggsave(
  "plot_MCResampling_SimStudy4.png",
  plot = plot_MCResampling_SimStudy4,
  dpi = 300,
  width = 15,
  height = 7,
  units = "in"
)


