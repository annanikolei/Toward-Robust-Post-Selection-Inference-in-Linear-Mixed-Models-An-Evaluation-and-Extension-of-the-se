## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

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
load(paste0("~/SimulationResults_supplied_D/SIM4_list_supplied_D_work_fullrandom=(Intercept)X2X4_n=400_Rsquared=",R_squared,"_nrsamples=10000.RData"))



################################################################################
## Run the Simulation ##
################################################################################
setwd("~/SimStudy4_Results")

# Number of planned iterations
N <- 1000

########################
# normal standard errors
########################

save_file <- paste0(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_marginal_not_efficient_RdiagTRUE_n=40_R2_", R_squared, "_1000MC.RDS"
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

save_file <- paste0(
  "SelFun=anova_SimulationStudy4_random=InterceptX2_conditional_RdiagTRUE_n=40_R2_", R_squared, "_1000MC.RDS"
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
        perspective = c("conditional"),
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
# HC3 standard errors
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
setwd("~/SimStudy4_Results")

SimulationStudy4_marginal_not_efficient_RdiagTRUE <- readRDS("SelFun=anova_SimulationStudy4_random=InterceptX2_marginal(not_eff)_RdiagTRUE_n=40_R2_0.3_1000MC.RDS")
SimulationStudy4_marginal_not_efficient_RdiagFALSE <- readRDS("SelFun=anova_SimulationStudy4_random=InterceptX2_marginal_not_efficient_RdiagFALSE_n=40_R2_0.3_1000MC.RDS")
SimulationStudy4_marginal_eff <- readRDS("SelFun=anova_SimulationStudy4_random=InterceptX2_marginal(eff)_RdiagTRUE_n=40_R2_0.3_1000MC.RDS")
SimulationStudy4_conditional <- readRDS("SelFun=anova_SimulationStudy4_random=InterceptX2_conditional_RdiagTRUE_n=40_R2_0.3_1000MC.RDS")
SimulationStudy4 <- c(SimulationStudy4_marginal_not_efficient_RdiagTRUE,
                      SimulationStudy4_marginal_not_efficient_RdiagFALSE, 
                      SimulationStudy4_marginal_eff, 
                      SimulationStudy4_conditional)
SimStudy4_Results <- do.call(rbind, lapply(SimulationStudy4 ,
                                           function(x) x$post_selection_pvalues)
)


SimStudy4_Results$pval_type <- ifelse(SimStudy4_Results$estimator%in% c("naive"),"Naive", ifelse(SimStudy4_Results$estimator%in% c("fullCandidate","fullCandidate_supplied"),"Standard (no selection)","Selfmade"))
SimStudy4_Results$estimator <- ifelse(SimStudy4_Results$estimator=="naive","ModelEstimate",SimStudy4_Results$estimator)
SimStudy4_Results$perspective <- ifelse(SimStudy4_Results$estimator %in% c("ModelEstimate_SIGMA_replaced","ModelEstimate_R_replaced",
                                                                           "ICM_SIGMA_replaced","ICM_R_replaced"),"marginal_not_efficient_HC3", SimStudy4_Results$perspective  )                                                                        



new_persp <- c(
  "marginal",
  "marginal_not_efficient",
  "marginal_not_efficient_HC3",
  "conditional",
  "fullCandidate"
)
naive_rows <- SimStudy4_Results[SimStudy4_Results$perspective == "naive", ]
naive_expanded <- naive_rows[rep(seq_len(nrow(naive_rows)), each = length(new_persp)), ]
naive_expanded$perspective <- rep(new_persp, times = nrow(naive_rows))
SimStudy4_Results <- SimStudy4_Results[SimStudy4_Results$perspective != "naive", ]
SimStudy4_Results <- rbind(SimStudy4_Results, naive_expanded)

# supplied is the same if HC3 estimator is used => duplicate this
supplied_rows <- SimStudy4_Results[SimStudy4_Results$estimator == "supplied" & SimStudy4_Results$perspective=="marginal_not_efficient", ]
supplied_expanded <- supplied_rows[rep(seq_len(nrow(supplied_rows)), each = length("marginal_not_efficient_HC3")), ]
supplied_expanded$perspective <- rep("marginal_not_efficient_HC3", times = nrow(supplied_rows))

SimStudy4_Results <- rbind(SimStudy4_Results, supplied_expanded )


SimStudy4_Results$coefficient <-
  ifelse(SimStudy4_Results$coefficient == "X3",
         "L1-noise~x[3]",
         "L1-noise~(random~slope)~x[4]")

  
library(ggplot2)
library(dplyr)


results_all_conditions_df <- SimStudy4_Results  %>%
  mutate(
    label_facet = paste0(
      #"R² = ", R_squared_pseudo,
      "Selected Random Effects = ", working_model
      
    )
  )

results_all_conditions_df$perspective <- factor(results_all_conditions_df$perspective)
results_all_conditions_df$perspective <- factor(
  results_all_conditions_df$perspective,
  levels = c("conditional", "marginal", "marginal_not_efficient","marginal_not_efficient_HC3","fullCandidate"),
  labels = c("Conditional", "Marginal (eff)", "Marginal (not eff)","Marginal (not eff) HC3","fullCandidate")
)

results_all_conditions_df$estimator <- factor(results_all_conditions_df$estimator)
results_all_conditions_df$estimator <- factor(
  results_all_conditions_df$estimator,
  levels = c("fullCandidate_supplied","supplied","ModelEstimate", "ICM","VarY",
             "ModelEstimate_SIGMA_replaced","ModelEstimate_R_replaced",
             "ICM_SIGMA_replaced","ICM_R_replaced",
             "fullCandidate"),
  labels = c("fullCandidate_supplied","supplied","ModelEstimate", "ICM","VarY",
             "ModelEstimate_SIGMA","ModelEstimate_R",
             "ICM_SIGMA","ICM_R",
             "fullCandidate")
)

results_all_conditions_df$pval_type<-relevel(as.factor(results_all_conditions_df$pval_type),"Selfmade")
results_all_conditions_df$coefficient <- factor(
  results_all_conditions_df$coefficient,
  levels = c("L1-noise~x[3]", "L1-noise~(random~slope)~x[4]")
)

# -----------------------------
# Manual color palette
# -----------------------------

estimator_colors <- c(
  "fullCandidate_supplied" ="black",
  "supplied"            = "#009E73",  
  "ModelEstimate"       =  "#8C510A",  
  "ICM"                 = "#D55E00",  
  "VarY"                =  "#0072B2", 
  "ModelEstimate_SIGMA" = "#56B4E9",  
  "ModelEstimate_R"     = "#E69F00", 
  "ICM_SIGMA"           = "#CC79A7",  
  "ICM_R"               = "#6A3D9A",  
  "fullCandidate"           = "white"
  
)

# -----------------------------
# Plot
# -----------------------------


results_all_conditions_df_SimStudy4 <- 
  results_all_conditions_df[results_all_conditions_df$perspective %in% c("Marginal (not eff)","Marginal (not eff) HC3"),]



plot_SimStudy4 <- ggplot(results_all_conditions_df_SimStudy4, 
                         aes(
                           sample = pval_two,
                           color  = estimator,
                           linetype = pval_type)) +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  theme_minimal(base_size = 18) +  # Basis-Schriftgröße
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA),  
    strip.text = element_text(color = "black"),
    legend.position = "bottom",       # <--- Hier Legende unten setzen
    legend.box = "vertical"
  ) +
  geom_qq(  
    distribution = stats::qunif, size = 0.8, geom = "line") +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )+
  labs(
    title = "",
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
    labels = function(x) parse(text = c(
      fullCandidate_supplied = "fullCandidate (supplied)",
      supplied                 = "supplied",
      ModelEstimate               = "ModelEstimate",
      ICM                         = "ICM",
      VarY                        = "VarY",
      ModelEstimate_SIGMA         = "ModelEstimate[Sigma]",
      ModelEstimate_R             = "ModelEstimate[R]",
      ICM_SIGMA                   = "ICM[Sigma]",
      ICM_R                       = "ICM[R]"
    )[x])
  )

plot_SimStudy4

ggsave("plot_SimStudy4.png", 
       plot = plot_SimStudy4, 
       dpi = 900, width = 9, height = 9, units = "in")


results_all_conditions_df_SimStudy4_Appendix <- 
  results_all_conditions_df[results_all_conditions_df$perspective %in% c("Conditional","Marginal (eff)"),]

plot_SimStudy4_Appendix <- ggplot(results_all_conditions_df_SimStudy4_Appendix, 
                                  aes(
                                    sample = pval_two,
                                    color  = estimator,
                                    linetype = pval_type)) +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  theme_minimal(base_size = 18) +  # Basis-Schriftgröße
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA),  
    strip.text = element_text(color = "black"),
    legend.position = "bottom",       # <--- Hier Legende unten setzen
    legend.box = "vertical"
  ) +
  geom_qq(  
    distribution = stats::qunif, size = 0.8, geom = "line") +
  facet_grid(
    perspective ~ label_facet + coefficient,
    labeller = labeller(
      coefficient = label_parsed,
      label_facet = label_wrap_gen(width = 35)
    )
  )+
  labs(
    title = "",
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
    labels = function(x) parse(text = c(
      fullCandidate_supplied = "fullCandidate (supplied)",
      supplied                 = "supplied",
      ModelEstimate               = "ModelEstimate",
      ICM                         = "ICM",
      VarY                        = "VarY",
      ModelEstimate_SIGMA         = "ModelEstimate[Sigma]",
      ModelEstimate_R             = "ModelEstimate[R]",
      ICM_SIGMA                   = "ICM[Sigma]",
      ICM_R                       = "ICM[R]",
      fullCandidate               = "fullCandidate"
    )[x])
  )
plot_SimStudy4_Appendix

ggsave("plot_SimStudy4_Appendix.png", 
       plot = plot_SimStudy4_Appendix, 
       dpi = 900, width = 9, height = 7, units = "in")




### convergence issues
# initial model fitting

SimStudy4_convergence <- do.call(rbind, lapply(SimulationStudy4 ,
                                               function(x) x$convergence_singular_info)
)



colMeans(SimStudy4_convergence[sapply(SimStudy4_convergence, is.numeric)], na.rm = F)

29.5290/( 85.5325 +29.5290)  
0.0005/85.5325   

## monte carlo resampling process

plot_df <- results_all_conditions_df %>%
  filter(
    perspective %in% c("Marginal (not eff) HC3","Marginal (not eff)"), 
    pval_type != "Naive"
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
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    panel.spacing.y = unit(1.5, "lines"),  
    panel.spacing.x = unit(1, "lines")     
  ) +
  scale_color_manual(
    values = estimator_colors,
    labels = function(x) parse(text = c(
      fullCandidate_supplied = "fullCandidate (supplied)",
      supplied                 = "supplied",
      ModelEstimate               = "ModelEstimate",
      ICM                         = "ICM",
      VarY                        = "VarY",
      ModelEstimate_SIGMA         = "ModelEstimate[Sigma]",
      ModelEstimate_R             = "ModelEstimate[R]",
      ICM_SIGMA                   = "ICM[Sigma]",
      ICM_R                       = "ICM[R]",
      fullCandidate               = "fullCandidate"
    )[x])
  )


plot_MCResampling_SimStudy4
ggsave("plot_MCResampling_SimStudy4.png", 
       plot = plot_MCResampling_SimStudy4, 
       dpi = 900, width = 9, height = 7.5, units = "in")





plot_df <- results_all_conditions_df %>%
  filter(
    perspective %in% c("Conditional","Marginal (eff)"), 
    pval_type != "Naive"
  ) %>%
  mutate(
    estimator = factor(estimator),
    pval_type = factor(pval_type),
    perspective = factor(perspective),
    label_facet = factor(label_facet),
    coefficient = factor(coefficient)
  )

plot_MCResampling_SimStudy4_Appendix <-
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
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    panel.spacing.y = unit(1.5, "lines"),  
    panel.spacing.x = unit(1, "lines")     
  ) +
  scale_color_manual(
    values = estimator_colors,
    labels = function(x) parse(text = c(
      fullCandidate_supplied = "fullCandidate (supplied)",
      supplied                 = "supplied",
      ModelEstimate               = "ModelEstimate",
      ICM                         = "ICM",
      VarY                        = "VarY",
      ModelEstimate_SIGMA         = "ModelEstimate[Sigma]",
      ModelEstimate_R             = "ModelEstimate[R]",
      ICM_SIGMA                   = "ICM[Sigma]",
      ICM_R                       = "ICM[R]",
      fullCandidate               = "fullCandidate"
    )[x])
  )


plot_MCResampling_SimStudy4_Appendix
ggsave("plot_MCResampling_SimStudy4_Appendix.png", 
       plot = plot_MCResampling_SimStudy4_Appendix, 
       dpi = 900, width = 9, height = 7.5, units = "in")
