################################################################################
# Simulation Study 1: Results
################################################################################

## Install and load all necessary packages
source("Project1_install_packages.R")

## Load selfmade implementation
source("selfmade**.R")

## Load Simulation Study 1 files
source("SimStudy1_helper.R")
source("SimStudy1_dgp.R")
source("SimStudy1_one_simulation.R")
source("SimStudy1_simulation_study.R")


## Run Simulation

N <- 3000  # Number of Iterations


for(n_obs in c(990,4830)){
  
  for(ICC_pseudo in c(0.05,0.25)){

save_file <- paste0("Simulation_Study1_all_n=", n_obs,
                                 "_ICC=", ICC_pseudo, ".RDS")


if (file.exists(save_file)) {
  results_list <- readRDS(save_file)
  start_iter <- length(results_list) + 1
  message("Continue with Iteration: ", start_iter)
} else {
  results_list <- list()
  start_iter <- 1
}

for (i in start_iter:N) {
  message("Iteration ", i, "/", N)
  
  res <- tryCatch(
    {
      SimStudy1_simulation_study(
        perspective = c("marginal_not_efficient", "marginal", "conditional"),
        n_iter = 1,
        n_obs = n_obs,
        ICC_pseudo = ICC_pseudo,
        checkConvergence_for_check_Fun = TRUE,
        ignore_singular_fits = TRUE,
        seed = i
      )
    },
    error = function(e) {
      warning("Error in Iteration ", i, ": ", e$message)
      return(NULL)
    }
  )
  
  results_list[[i]] <- res
  
  # nach jeder Iteration sichern
  saveRDS(results_list, save_file)
  
   }
  }
}



## Create QQ-Plots

## Create QQ-plot --------------------------------------------------------------
# 0. load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)


for(n_obs in c(990,4830)){
  
  for(ICC_pseudo in c(0.05,0.25)){
    
    Simulation_Study1 <- readRDS(paste0("Simulation_Study1_all_n=", n_obs,
                                        "_ICC=", ICC_pseudo, ".RDS"))   

# 1. Extract p-values and combine them into a long dataframe
# Preparation: Combine p-values, n_obs, ICC_pseudo in a list

df_input <- tibble(
  pval = sapply(Simulation_Study1, function(x) x$pvalues_and_MCresamplings),
  n_obs = n_obs,
  ICC_pseudo = ICC_pseudo
)

# 1. Extract p-values and convert to long dataframe format
pval_df <- pmap_dfr(df_input, function(pval, n_obs, ICC_pseudo) {
  as.data.frame(pval) %>%
    mutate(
      term = rownames(pval),
      n_obs = n_obs,
      ICC_pseudo = ICC_pseudo
    )
}, .id = "i_iter")


# 2. Long format for ggplot
pval_long <- pval_df %>%
  pivot_longer(
    cols = ends_with("_pvals"),
    names_to = "perspective",
    values_to = "p_value"
  )

DEFT = function(ICC,n){
  
  sqrt(1+(n/30 -1)*ICC)
}

pval_long$DEFT <- DEFT(pval_long$ICC_pseudo,pval_long$n_obs)

# 3. Rename for better readability 
pval_long$perspective <- recode(pval_long$perspective,
                                "naive_pvals" = "Naive",
                                "selfmade_marginal_pvals" = "Marginal (efficient)",
                                "selfmade_marginal_not_efficient_pvals" = "Marginal (not efficient)",
                                "selfmade_conditional_bayes_pvals" = "Conditional (Bayes)",
                                "selfmade_conditional_not_bayes_pvals" = "Conditional (frequentist)"
)

pval_long <- pval_long %>%
  mutate(
    n_obs_label = paste0("n == ", n_obs),
    ICC_label = paste0("ICC == ", ICC_pseudo)
  )

pval_long$term <- recode(pval_long$term, 
                         "X1.l1" = "L1-noise~(random~slope)",
                         "X2.l1" = "L1-noise",
                         "X.l2" = "L2-noise")

# Relevel für n_obs
pval_long$n_obs_label <- factor(
  pval_long$n_obs_label,
  levels = paste0("n == ", sort(unique(pval_long$n_obs)))
)
pval_long$n_obs_label <- gsub("^n", "N", pval_long$n_obs_label)


# Relevel für ICC
pval_long$ICC_label <- factor(
  pval_long$ICC_label,
  levels = paste0("ICC == ", sort(unique(pval_long$ICC_pseudo)))
)



# Round DEFT 
pval_long$DEFT_label <- paste0("DEFT == ", round(pval_long$DEFT, 2))



pval_long_clean <- pval_long %>%
  filter(
    !is.na(DEFT_label),
    !is.na(n_obs_label),
    !is.na(ICC_label),
    !is.na(term)
  )

plot_SimStudy1 <- ggplot(data = pval_long_clean,  
                         aes(sample = p_value,  
                             colour = perspective
                         )) +
  theme_minimal(base_size = 18) +  # Basis-Schriftgröße
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA),  
    strip.text = element_text(color = "black"),
    legend.position = "bottom"       # <--- Hier Legende unten setzen
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey") + 
  geom_qq(  
    distribution = stats::qunif, size = 0.8, geom = "line") +
  xlab("Expected p-values") + ylab("Observed p-values") +  
  facet_grid( DEFT_label + n_obs_label + ICC_label ~ term, labeller = label_parsed) +  
  labs(colour = "",
       title = "", 
       subtitle = "")


ggsave("plot_SimStudy1.png", 
       plot = paste0("plotSimStudy1_n=", n_obs,
                     "_ICC=", ICC_pseudo, ".RDS"), 
       dpi = 900, width = 12.2, height = 10, units = "in")


  }
   }

#-------------------------------------------------------------------------------
# convergence issues with initial model fitting

all_tables <- list()

for(n_obs in c(990, 4830)){
  
  for(ICC_pseudo in c(0.05, 0.25)){
    
    file_name <- paste0("Simulation_Study1_all_n=", 
                        n_obs,
                        "_ICC=", 
                        ICC_pseudo, 
                        ".RDS")
    
    if(!file.exists(file_name)) next
    
    Simulation_Study1 <- readRDS(file_name)
    
    # Remove NULL entries (failed iterations)
    Simulation_Study1 <- Simulation_Study1[!sapply(Simulation_Study1, is.null)]
    
    retry_counts <- sapply(Simulation_Study1, function(x) x$retry_counts)
    
    # Store table with descriptive name
    table_name <- paste0("n=", n_obs, "_ICC=", ICC_pseudo)
    all_tables[[table_name]] <- table(retry_counts)
  }
}

all_tables


#-------------------------------------------------------------------------------
# convergence issues with Monte Carlo sampling approach 

for(n_obs in c(990, 4830)){
  
  for(ICC_pseudo in c(0.05, 0.25)){
    
    file_name <- paste0("Simulation_Study1_all_n=", 
                        n_obs,
                        "_ICC=", 
                        ICC_pseudo, 
                        ".RDS")

Simulation_Study1_MCresamplings <- 
  as.data.frame(do.call(rbind,sapply(file_name, function(x) x$pvalues_and_MCresamplings)))


df_long_MCresamplings <- Simulation_Study1_MCresamplings %>%
  tibble::rownames_to_column("term") %>%
  pivot_longer(
    cols = -term,
    names_to = "name",
    values_to = "value"
  ) %>%
  mutate(
    type = if_else(str_detect(name, "pvals$"), "pvals", "resamplings"),
    perspective = name %>%
      str_remove("_pvals$") %>%
      str_remove("_total_resamplings_due_to_convergence$")
  ) %>%
  select(term, perspective, type, value) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  )
df_long_MCresamplings$N <- n_obs
df_long_MCresamplings$ICC <- ICC_pseudo
df_long_MCresamplings$DEFT <- DEFT(ICC_pseudo,n_obs)
df_long_MCresamplings <-
  df_long_MCresamplings[df_long_MCresamplings$perspective %in% c("selfmade_marginal_not_efficient",
                                                                                 "selfmade_marginal",
                                                                                 "selfmade_conditional_bayes"),]
 
plot_df <- rbind(df_long_MCresamplings) %>% mutate(
                   perspective = factor(perspective),
                   DEFT = factor(DEFT),
                   ICC = factor(ICC),
                   N=factor(N)
                 )

plot_df$perspective <- recode(plot_df$perspective,
                              "selfmade_marginal" = "Marginal (efficient)",
                              "selfmade_marginal_not_efficient" = "Marginal (not efficient)",
                              "selfmade_conditional_bayes" = "Conditional (Bayes)",
                              "selfmade_conditional_not_bayes" = "Conditional (frequentist)"
)



plot_df$term <- sub("\\.[0-9999999]+$", "", plot_df$term)
plot_df$term <- sub("(l[0-9])([0-9]+)$", "\\1", plot_df$term)

plot_df$term <- factor(recode(plot_df$term, 
                              "X1.l1" = "L1-noise~(random~slope)",
                              "X2.l1" = "L1-noise",
                              "X.l2" = "L2-noise"))

plot_df$N <- factor(
  plot_df$N,
  levels = sort(unique(plot_df$N)),
  labels = paste0("N == ", sort(unique(plot_df$N)))
)



# Relevel für ICC
plot_df$ICC <- factor(
  plot_df$ICC,
  levels = sort(unique(plot_df$ICC)),
  labels = paste0("ICC == ", sort(unique(plot_df$ICC)))
)



# Runde DEFT auf z.B. 2 Nachkommastellen für schönere Labels
plot_df$DEFT <- factor(
  plot_df$DEFT,
  levels = sort(unique(plot_df$DEFT)),
  labels = paste0("DEFT ==", sort(unique(plot_df$DEFT)))
)


boxplot <-
  ggplot(
    plot_df,
    aes(
      x = "",
      y = resamplings,
      colour = perspective
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.alpha = 0.25,
    fill = NA
  ) +
  facet_grid(DEFT + N + ICC ~ term, labeller = label_parsed) +
  labs(
    colour = "",
    x = "",
    y = "Number of MC Resamplings due to convergence"
  )+
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(color = "black", size = 12),
    legend.position = "bottom",
    panel.spacing.y = unit(1.5, "lines"),  # ⬅️ Abstand zwischen Zeilen erhöhen
    panel.spacing.x = unit(1, "lines")     # optional: Abstand zwischen Spalten
  ) +scale_color_manual(values = c(
    "#F8766D",
    "#7CAE00",
    "#00BFC4",
    "#00B0F6"
  ))
boxplot



ggsave("boxplot_SimStudy1.png", 
       plot = paste0("boxplotSimStudy1_n=", n_obs,
                     "_ICC=", ICC_pseudo, ".RDS"), 
       dpi = 900, width = 12.2, height = 10, units = "in")

  }
  
  }
