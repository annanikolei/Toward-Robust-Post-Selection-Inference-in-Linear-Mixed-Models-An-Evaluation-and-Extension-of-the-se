################################################################################
# Simulation Study 4 - Selection process functions
################################################################################

library(lmerTest)

## Function: Fit Fullest candidate model ---------------------------------------

modFun <- function(dat){
  
  full_model <-
    lmer( y ~ 1 + X1 + X2 + X3 + X4 + (1 + X2 + X4 | l2.id), 
          REML = FALSE, data = dat)
  
  return(full_model)
}

## selection function ----------------------------------------------------------

SelFun <- function(dat, full_model, alpha = 0.000001) {
  
  library(lme4)
  
  # Volles Modell
  full_model <- full_model
  
  # Modell ohne X2-Slope
  no_X4 <- lmer(y ~ 1 + X1 + X2 + X3 + X4 + 
                  (1 + X2 | l2.id),
                REML = FALSE, data = dat)
  
  # Nur Random Intercept
  intercept_only <- lmer(y ~ 1 + X1 + X2 + X3 + X4 + 
                           (1 | l2.id),
                         REML = FALSE, data = dat)
  
  current_model <- full_model
  
  # Test X4-Slope
  test_X4 <- anova(current_model, no_X4)
  if(test_X4$`Pr(>Chisq)`[2] > alpha) {
    current_model <- no_X4
  }
  
  # Test X2-Slope (gegen aktuelles Modell!)
  test_X2 <- anova(current_model, intercept_only)
  if(test_X2$`Pr(>Chisq)`[2] > alpha) {
    current_model <- intercept_only
  }
  
  return(current_model)
}





## extract selection event -----------------------------------------------------

extractSelFun_random <- function(this_mod){
  
  return(names(ranef(this_mod)$l2.id))
}


