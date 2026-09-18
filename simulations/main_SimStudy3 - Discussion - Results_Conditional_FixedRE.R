################################################################################
# Simulation Study 3-Discussion: Conditional Fixed RE
################################################################################

## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

## Load Simulation Study files
setwd("~/Downloads/selective-inference-lmm-simulation-main/simulations")
source("SimStudy3_helper.R")
source("SimStudy3_Discussion_dgp_Conditional_FixedRE.R")
source("SimStudy3_SelectionProcess.R")
source("SimStudy3_Discussion_one_simulation_Conditional_FixedRE.R")
source("SimStudy3_Discussion_simulation_study_Conditional_FixedRE.R")

R_squared = 0.3
n_obs = 40^2 ; m = 40

# load fixed predictor matriced
dgp_X_output <- readRDS("dgp_X_output_n_obs=1600_m=40.RDS")

# load list_supplied_D_work
load(paste0("SimulationResults_supplied_D/SIM3_list_supplied_D_work_full=X1X2X3X4_interesting=X3X4_n=400_Rsquared=",R_squared,"_nrsamples=10000.RData"))

u_fixed <- SimStudy3_dgp_u(m = m,
                           D = matrix(c(2,0.2,0.1,
                                        0.2,1.5,0.2,
                                        0.1,0.2,1.5),
                                          nrow=3, byrow=T),
                           seed = 170)
########################
# Run Simulations
########################
setwd("Conditional_FixedRE_SimStudy3_Results")


# conditional X124, target = rho
N <- 1000

save_file <- paste0(
  "FixedRE_rho_SimStudy3_workingModel=X124_conditional(corrected)_R2_0.3_1000MC.RDS"
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
                                 dgp_X_output = dgp_X_output, # SimStudy3_dgp_X( n_obs = n_obs, m = m)
                                 R_squared = 0.3, # pseudo coefficient of determination
                                 D = matrix(c(2,0.2,0.1,
                                              0.2,1.5,0.2,
                                              0.1,0.2,1.5),
                                            nrow = 3, byrow=T) , # 3x3 covariance matrix of random effects (u0,u3,u7) per level2 unit
                                 list_supplied_D_work = list_supplied_D_work,
                                 u_fixed = u_fixed,
                                 perspective = "conditional", 
                                 bayesian= c("bayesian_corrected"),
                                 B_MCsamples = 1000, # number of MC-samples
                                 remove_singular_fits = F, #remove singular fit 
                                 interesting_var = "X4", # "X3" or "X4"
                                 interesting_working_model_vars = c("X1","X2","X4"), # for which working model do we want to make post-selective inference, if NULL then no specific working model should be considered
                                 checkConvergence_for_check_Fun = T,
                                 target = "rho" )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    })
  
  
  results_list[[i]] <- res
  
  
  saveRDS(results_list, save_file)
}


################################################################################
## QQ-Plots, Convergence ##
################################################################################
setwd("Conditional_FixedRE_SimStudy3_Results")

SimulationStudy3 <- readRDS(paste0(
  "FixedRE_rho_SimStudy3_workingModel=X124_conditional(corrected)_R2_0.3_1000MC.RDS"
))


estimates_df <- do.call(
    rbind,
    lapply(SimulationStudy3, function(x) x$estimate)
  )

colMeans(estimates_df["estimate"]) - u_fixed$mean_u[3]
  


SimStudy3_Results <- do.call(rbind, lapply(SimulationStudy3 ,
                                           function(x) x$post_selection_pvalues)
)




SimStudy3_Results$coefficient <-
  ifelse(SimStudy3_Results$coefficient == "X4",
         "L1-noise~(RS)~x[4]","L1-noise~x[3]")


library(ggplot2)
library(dplyr)


results_all_conditions_df <- SimStudy3_Results  %>%
  mutate(
    label_facet = paste0(
      #"R² = ", R_squared_pseudo,
      "Selected Fixed = ", working_model
      
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

results_all_conditions_df$perspective <- factor(results_all_conditions_df$perspective)
results_all_conditions_df$perspective <- factor(
  results_all_conditions_df$perspective,
  levels = c(
    "conditional"
  ),
  labels = c(
    "Conditional\n(corrected)"
  )
)


results_all_conditions_df$estimator <- factor(results_all_conditions_df$estimator)
results_all_conditions_df$estimator <- factor(
  results_all_conditions_df$estimator,
  levels = c("supplied_bayesian_corrected",
             "ModelEstimate_bayesian_corrected", 
             "ICM_bayesian_corrected",
             "VarY_bayesian_corrected"),
  labels = c("supplied","ModelEstimate", "ICM","VarY"))

results_all_conditions_df$pval_type<-"Selfmade"
results_all_conditions_df$coefficient <- factor(
  results_all_conditions_df$coefficient,
  levels = c( "L1-noise~x[3]",
              "L1-noise~(RS)~x[4]")
)


# -----------------------------
# Manual color palette
# -----------------------------

library(ggplot2)
library(cowplot)

# -----------------------------
# Estimator colors
# -----------------------------

estimator_colors <- c(
  "fullCandidate_supplied"              = "black",
  "supplied"                            = "#009E73",  
  "ModelEstimate"                       = "#8C510A",  
  "ICM"                                 = "#D55E00",  
  "VarY"                                = "#0072B2", 
  "ModelEstimate_SIGMA"                 = "#56B4E9",  
  "ModelEstimate_R"                     = "#E69F00", 
  "ICM_SIGMA"                           = "#CC79A7",  
  "ICM_R"                               = "#6A3D9A",  
  "fullCandidate"                       = "white"
  
)

# -----------------------------
# Data
# -----------------------------

results_all_conditions_df_SimStudy3 <- 
  results_all_conditions_df

# -----------------------------
# Base plot
# -----------------------------

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
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  ) +
  theme_minimal(base_size = 18) +
  theme(
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    strip.text = element_text(
      color = "black",
      size = 19,
      lineheight = 0.8
    ),
    legend.text = element_text(size = 19),
    legend.title = element_text(size = 19),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 19)
  )+
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
  ) +
  scale_color_manual(
    values = estimator_colors,
    labels = function(x) parse(text = c(
      supplied                         = "supplied",
      ModelEstimate                    = "ModelEstimate",
      ICM                               = "ICM",
      VarY                              ="VarY",
      ModelEstimate_SIGMA               = "ModelEstimate[Sigma]",
      ICM_SIGMA                         = "ICM[Sigma]",
      ModelEstimate_R                   = "ModelEstimate[R]",
      ICM_R                             = "ICM[R]"
    )[x])
  )

plot_SimStudy3

ggsave("plot_SimStudy3.png", 
       plot = plot_SimStudy3, 
       dpi = 300, width = 16, height = 12, units = "in")



SimStudy3_convergence <- do.call(rbind, lapply(c(SimulationStudy3,SimulationStudy3_marginalCR2) ,
                                               function(x) x$convergence_singular_info)
)



colMeans(SimStudy3_convergence[sapply(SimStudy3_convergence, is.numeric)], na.rm = F)

0.0004666667/( 0.0004666667+ 189.8884666667)  # less than 0.001%

library(dplyr)

View(SimStudy3_convergence %>%
       group_by(working_model) %>%
       summarise(across(where(is.numeric), ~ mean(.x, na.rm = FALSE))))

0.0004/(0.0004+ 200.9648) 
0.0010/(0.0010+
          365.8078)# less than 0.001%
0.0004/
  70.8144 # less than 0.001%

## monte carlo resampling process

plot_df <- results_all_conditions_df_SimStudy3 %>%
  filter(
    pval_type != "Naive"
  ) %>%
  mutate(
    estimator = factor(estimator),
    pval_type = factor(pval_type),
    perspective = factor(perspective),
    label_facet = factor(label_facet),
    coefficient = factor(coefficient)
  )

plot_MCResampling_SimStudy3 <-
  ggplot(
    plot_df,
    aes(
      x ="",
      y = total_resamplings_due_to_convergence,
      color = estimator,
      group = interaction(pval_type, estimator)
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.alpha = 0.25,
    fill = NA
  ) +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )+
  labs(
    x = "",
    y = "Number of MC Resamplings due to convergence",
    color = "Estimator"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    strip.text = element_text(
      color = "black",
      size = 19,
      lineheight = 0.8
    ),
    legend.text = element_text(size = 19),
    legend.title = element_text(size = 19),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 19)
  )+
  scale_color_manual(
    values = estimator_colors,
    labels = function(x) parse(text = c(
      supplied                         = "supplied",
      ModelEstimate                    = "ModelEstimate",
      ICM                               = "ICM",
      VarY                              ="VarY",
      ModelEstimate_SIGMA               = "ModelEstimate[Sigma]",
      ICM_SIGMA                         = "ICM[Sigma]",
      ModelEstimate_R                   = "ModelEstimate[R]",
      ICM_R                             = "ICM[R]"
    )[x])
  )


plot_MCResampling_SimStudy3
ggsave("plot_MCResampling_SimStudy3.png", 
       plot = plot_MCResampling_SimStudy3, 
       dpi = 300, width = 16, height = 12, units = "in")



# bayesian conditional (uncorrected)
results_all_conditions_df_SimStudy3_Appendix <- 
  results_all_conditions_df[
    results_all_conditions_df$perspective %in% 
      c(
        "Conditional\n(uncorrected)"
      ),
  ]


# -----------------------------
# Base plot
# -----------------------------

plot_SimStudy3_Appendix <- ggplot(
  results_all_conditions_df_SimStudy3_Appendix,
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
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  ) +
  theme_minimal(base_size = 18) +
  theme(
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    strip.text = element_text(
      color = "black",
      size = 19,
      lineheight = 0.8
    ),
    legend.text = element_text(size = 19),
    legend.title = element_text(size = 19),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 19)
  )+
  labs(
    title = "",
    x = "Expected p-values",
    y = "Observed p-values",
    color = "Estimator",
    linetype = "p-value type"
  ) +
  scale_color_manual(
    values = estimator_colors,
    labels = function(x) parse(text = c(
      supplied                         = "supplied",
      ModelEstimate                    = "ModelEstimate",
      ICM                               = "ICM",
      VarY                              ="VarY",
      ModelEstimate_SIGMA               = "ModelEstimate[Sigma]",
      ICM_SIGMA                         = "ICM[Sigma]",
      ModelEstimate_R                   = "ModelEstimate[R]",
      ICM_R                             = "ICM[R]"
    )[x])
  )

plot_SimStudy3_Appendix

ggsave("plot_SimStudy3_Appendix.png", 
       plot = plot_SimStudy3_Appendix, 
       dpi = 300, width = 16, height = 4, units = "in")




plot_df <- results_all_conditions_df_SimStudy3_Appendix %>%
  filter(
    pval_type != "Naive"
  ) %>%
  mutate(
    estimator = factor(estimator),
    pval_type = factor(pval_type),
    perspective = factor(perspective),
    label_facet = factor(label_facet),
    coefficient = factor(coefficient)
  )

plot_MCResampling_SimStudy3_Appendix <-
  ggplot(
    plot_df,
    aes(
      x ="",
      y = total_resamplings_due_to_convergence,
      color = estimator,
      group = interaction(pval_type, estimator)
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.alpha = 0.25,
    fill = NA
  ) +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )+
  labs(
    x = "",
    y = "Number of MC Resamplings due to convergence",
    color = "Estimator"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    strip.background = element_rect(
      fill = "lightgrey",
      color = NA
    ),
    strip.text = element_text(
      color = "black",
      size = 19,
      lineheight = 0.8
    ),
    legend.text = element_text(size = 19),
    legend.title = element_text(size = 19),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 19)
  )+
  scale_color_manual(
    values = estimator_colors,
    labels = function(x) parse(text = c(
      supplied                         = "supplied",
      ModelEstimate                    = "ModelEstimate",
      ICM                               = "ICM",
      VarY                              ="VarY",
      ModelEstimate_SIGMA               = "ModelEstimate[Sigma]",
      ICM_SIGMA                         = "ICM[Sigma]",
      ModelEstimate_R                   = "ModelEstimate[R]",
      ICM_R                             = "ICM[R]"
    )[x])
  )


plot_MCResampling_SimStudy3_Appendix
ggsave("plot_MCResampling_SimStudy3_Appendix.png", 
       plot = plot_MCResampling_SimStudy3_Appendix, 
       dpi = 300, width =16, height = 6.5, units = "in")


