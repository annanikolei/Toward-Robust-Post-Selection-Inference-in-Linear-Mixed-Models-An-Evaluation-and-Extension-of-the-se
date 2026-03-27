################################################################################
# Simulation Study 3 - Selection process functions
################################################################################

library(lmerTest)

## Function: Fit Fullest candidate model ---------------------------------------

modFun <- function(dat){
  
  full_model <-
    lmer( y ~ 1 + X1 + X2 + X3 + X4 + (1| l2.id), 
          REML = FALSE, data = dat)
  
  return(full_model)
}

## selection function ----------------------------------------------------------



SelFun <- function(model) {
  
  dat <- model.frame(model)
  
  # Fix: overwrite the model call to carry the actual data frame
  model@call$data <- dat
  
  attr( lmerTest:::step.lmerModLmerTest(
    model,
    alpha.fixed = 0.05,
    reduce.random = FALSE), "model")
  
}






## extract selection event -----------------------------------------------------

extractSelFun_fixed <- function(this_mod){
  
  return(setdiff(names(fixef(this_mod)), "(Intercept)"))
}


