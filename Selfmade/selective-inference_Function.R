################################################################################
# selective-inference functions
################################################################################

library(performance)

# Description: In this script, MC samples are drawn and used to compute Ceta,
#              which is then used to determine selective-inference p-values and
#              selective-inference confidence intervals.
#
# The first function (name: pval_CI; analogous to selinf in the selfmade package)
# computes p-values and confidence intervals based on the given MC samples using
# the approximate Ceta formula (7) [thesis p. 18].
#
# The second function (name: selective_inference; analogous to pval_vT_cov
# function in the selfmade package) includes the Ceta function and first draws
# MC samples from the proposal distribution.

#-------------------------------------------------------------------------------
# pval_CI
#-------------------------------------------------------------------------------

pval_CI <- function(Tbs_survived, # all MC test statistics that lie in the inference region
                    t_obs, # observed test statistic
                    wbs, # weights for approximating Ceta: formula (7) in thesis
                    kappa_square, # variance of the test statistic
                    rho_H0 = 0, # null hypothesis value of the parameter
                    symmetric_selection = T, # whether selection is symmetric
                    alpha = 0.05, # significance level
                    allow_smaller_alpha = T, # allow computation of wider CIs if needed
                    multiple = 7, # controls range of rho_sequence
                    total_resamplings_due_to_convergence
){
  
  
  # number of MC samples that lie in the inference region
  number_survived_MC <- length(Tbs_survived)
  
  ## If computation of p-values and CI bounds is not possible
  # If none of the MC test statistics lie in the inference region,
  # neither p-values nor confidence intervals can be computed
  
  # check whether at least 1 MC sample lies in the inference region
  if(number_survived_MC == 0) {
    warning("All MC samples lie outside the inference region, therefore 
    no p-values or confidence intervals could be computed. Reduce 
    rejection rate in MC sampling using the following options:
    (1) Increase B_MC_samples
    (2) Use a different MC_distribution than 'centered'")
    return(data.frame(number_survived_MC = 0,
                      pval_right = NA,
                      pval_two = NA,
                      cil = NA, 
                      ciu = NA,
                      cil.smalleralpha=NA,
                      ciu.smalleralpha=NA,
                      total_resamplings_due_to_convergence))
  }
  
  # which Tbs in the inference region are greater (or equal) to t_obs
  which.Tbs_surv_geq_tobs <- Tbs_survived >= rep(t_obs, length(Tbs_survived))
  Tbs_surv_geq_tobs <- Tbs_survived[which.Tbs_surv_geq_tobs]
  
  if(symmetric_selection == FALSE){
    # which absolute values of Tbs in the inference region are greater (or equal)
    # to the absolute value of t_obs
    which_abs.Tbs_surv_geq_tobs <- abs(Tbs_survived) >= rep(abs(t_obs), 
                                                            length(Tbs_survived))
    abs.Tbs_surv_geq_tobs <- Tbs_survived[which_abs.Tbs_surv_geq_tobs]}
  
  
  # Computation of p-values--------------------------------------------------------
  
  # p-value:
  # Usually computed using the Ceta formula [formula (7) in thesis], but if
  # testing against rho_H0 = 0, the formula simplifies since exp(0) = 1
  
  if(rho_H0 == 0){
    
    # denominator of the Ceta formula
    sum_wbs <- sum(wbs)
    
    # right-sided p-value
    pval_right <- as.numeric(sum(wbs[which.Tbs_surv_geq_tobs]) / sum_wbs)
    
    # two-sided p-value
    # defined here as P(|T| >= |t_obs| | H0), since the original distribution
    # of the test statistic is symmetric (and often remains so after selection)
    if(symmetric_selection == TRUE){
      pval_two <- 2* min(pval_right, 1 - pval_right)
    }
    else if(symmetric_selection == FALSE){
      pval_two <- as.numeric(sum(wbs[which_abs.Tbs_surv_geq_tobs]) / sum_wbs)
    }
    
  }
  
  else if(rho_H0 != 0){ # for rho_H0 != 0, the simplified Ceta formula cannot be used
    
    # right-sided p-value
    pval_right <- as.numeric(
      Ceta(wbs = wbs,
           wbs_geq_tobs = wbs[which.Tbs_surv_geq_tobs],
           rho = rho_H0,
           kappa_square = kappa_square,
           Tbs_survived = Tbs_survived,
           Tbs_surv_geq_tobs = Tbs_surv_geq_tobs))
    
    # two-sided p-value
    if(symmetric_selection == TRUE){
      pval_two <- min(pval_right, 1 - pval_right)
    }
    if(symmetric_selection == FALSE){
      pval_two <-  as.numeric(
        Ceta(wbs = wbs,
             wbs_geq_tobs = wbs[which_abs.Tbs_surv_geq_tobs],
             rho = rho_H0,
             kappa_square = kappa_square,
             Tbs_survived = Tbs_survived,
             Tbs_surv_geq_tobs = abs.Tbs_surv_geq_tobs))
    }
    
  }
  
  
  
  # Computation of confidence interval bounds-------------------------------------
  
  ## If computation of (1-alpha)% CI bounds is not possible
  # This happens if all surviving MC Tbs are either greater or smaller than t_obs,
  # since a rho must be found such that Ceta or 1-Ceta equals alpha/2 or 1-alpha/2.
  if(symmetric_selection == TRUE){
    if(all(Tbs_survived > rep(t_obs,number_survived_MC)) |
       all(Tbs_survived < rep(t_obs,number_survived_MC)) ){
      return(data.frame(number_survived_MC = number_survived_MC,
                        pval_right = pval_right,
                        pval_two = pval_two,
                        cil = NA, 
                        ciu = NA,
                        cil.smalleralpha=NA,
                        ciu.smalleralpha=NA,
                        total_resamplings_due_to_convergence))
      
    }
    
    ## If CI computation is possible, the Ceta formula must be inverted w.r.t. rho
    
    # generate a sequence of possible rho values
    rho_sequence <- seq(min(Tbs_survived ) - multiple*sqrt(kappa_square),
                        max(Tbs_survived ) + multiple*sqrt(kappa_square),
                        l = multiple * 1000)
    
    # function to compute Ceta values for given rho
    Ceta_function <- function(rho){Ceta(wbs = wbs,
                                        wbs_geq_tobs = wbs[which.Tbs_surv_geq_tobs],
                                        rho = rho,
                                        kappa_square,
                                        Tbs_survived,
                                        Tbs_surv_geq_tobs)}
    
    # compute lower and upper CI bounds: find rho such that
    # Ceta(rho_low) = alpha/2 and Ceta(rho_up) = 1 - alpha/2
    
    ## lower bound 
    ci_low <- NA
    cil.smalleralpha <- NA
    
    ci_low <- try(uniroot(f = function(rho){
      as.numeric(Ceta_function(rho)) - alpha/2}, 
      interval = range(rho_sequence),
      extendInt = "upX",
      maxiter = 100000)$root, silent = TRUE)
    
    # if this fails, compute alternative CI bound
    # NOTE: The logic for alternative CI bounds is currently not correct
    # and should be revisited. 
    if(allow_smaller_alpha == T & class(ci_low) =="try-error"){ 
      if(t_obs >= rho_H0){
        cil.smalleralpha <- 
          max(rho_sequence[rho_sequence > rho_H0 & rho_sequence < t_obs]) 
      }
      else if(t_obs < rho_H0){
        cil.smalleralpha <- 
          max(rho_sequence[rho_sequence < rho_H0 & rho_sequence < t_obs]) 
      }
      
      warning("Lower interval limit is not the ",
              alpha/2*100, "%-quantile, but the (",
              signif(as.numeric(Ceta_function(cil.smalleralpha)*100),3), 
              ")%-quantile.") 
    }
    
    ## upper bound
    ci_up <- NA
    ciu.smalleralpha <- NA
    
    ci_up <- try(uniroot(f = function(rho){
      as.numeric(Ceta_function(rho)) - (1-alpha/2)}, 
      interval = range(rho_sequence),
      extendInt = "upX",
      maxiter = 100000)$root, silent = TRUE)
    
    # if this fails, compute alternative CI bound
    if(allow_smaller_alpha == T & class(ci_low) =="try-error"){ 
      ci_up <- NA
      
      if(t_obs >= rho_H0){
        ciu.smalleralpha <- 
          min(rho_sequence[rho_sequence > rho_H0 & rho_sequence > t_obs]) 
      }
      else if(t_obs < rho_H0){
        ciu.smalleralpha <- 
          min(rho_sequence[rho_sequence < rho_H0 & rho_sequence > t_obs]) 
      }
      
      warning("Upper interval limit is not the ",
              (1-alpha/2)*100, "%-quantile, but the (",
              signif(100-(as.numeric(Ceta_function(ciu.smalleralpha))*100),3), 
              ")%-quantile.") 
    }
    
    return(data.frame(number_survived_MC = number_survived_MC,
                      pval_right = pval_right,
                      pval_two = pval_two,
                      cil = ci_low, 
                      ciu = ci_up,
                      cil.smalleralpha = cil.smalleralpha,
                      ciu.smalleralpha = ciu.smalleralpha,
                      total_resamplings_due_to_convergence))
  }       
  
  
  
}


#-------------------------------------------------------------------------------
# selective_inference
#-------------------------------------------------------------------------------
if (!("future" %in% rownames(installed.packages()))) {
  install.packages("future")
}
if (!("future.apply" %in% rownames(installed.packages()))) {
  install.packages("future.apply")
}
library(future.apply)       
plan(multisession) #?# nutzt alle verfügbaren Kerne. Ist das dumm?

selective_inference <- function(n,
                                vt,
                                kappa_square,
                                dirV,
                                this_y,
                                rho_H0 = 0,
                                B_MCsamples = 1000,
                                checkFun,
                                MC_distribution = "proposal",
                                symmetric_selection = T,
                                alpha = 0.05, 
                                allow_smaller_alpha = T, 
                                multiple = 7,
                                checkConvergence_for_check_Fun,
                                checkFun_checkConv,
                                output = "pvalCI", # or "selective_inference"
                                definition_sampFun = NULL, # function(num) random generation
                                definition_dsampFun = NULL # function(tsta,log) # density
){
  
  # Berechnung der beobachteten Teststatistik t_obs
  t_obs <- as.numeric(vt %*% this_y) 
  
  # Berechnung der orthogonalen Projektion := (I-dirV) * vt *y
  orth_proj <- (diag(n) - dirV %*% vt) %*% this_y
  
  # MC-samples aus zentrierter Normalverteilung N(0,kappa_square) bzw. aus 
  # proposal distribution N(t_obs,kappa_square) oder proposal distribution, die
  # sich aus mehreren NVs zusammensetzt, ziehen und dann Übereinstimmung mit dem
  # initialen selection event überprüfen#-----------------------------------------
  
  # MC_distribution festlegen:                                                             #?# Einbauen von check, ob proposal distribution gut funktioniert
  if(MC_distribution == "proposal"){
    sampFun = function(num) rnorm(num,
                                  mean = t_obs,
                                  sd = sqrt(kappa_square))
    dsampFun = function(tsta,log) dnorm(tsta,
                                        mean = t_obs,
                                        sd = sqrt(kappa_square),
                                        log = log)
  }
  
  if(MC_distribution == "definition"){
    sampFun = definition_sampFun
    dsampFun = definition_dsampFun
  }
  
  if(MC_distribution == "centered"){
    sampFun = function(num) rnorm(num, 
                                  mean = 0, 
                                  sd = sqrt(kappa_square))
  }
  
  # Hier könnte eine andere Möglichkeit für sampFun stehen                                 #?# eventuell hier das noch mit der mixture distribution?
  
  if(checkConvergence_for_check_Fun == F){
    
    if(is.null(checkFun)){  
      stop("checkFun muss angegeben sein, wenn convergence bei MCsampling nicht
         überprüft werden soll!") }
    
    # MC-samples ziehen
    Tbs <- sampFun(num = B_MCsamples) 
    
    # Transformation der Tbs zu den y_bs. Pro b: y_b = T_b * dirV + P_ortho
    ybs <- lapply(Tbs, function(Tb) as.numeric(Tb * dirV + orth_proj)) 
    
    # Übereinstimmung mit initialem selection event überprüfen
    # # unparallelisierte Version
    # indicator_inference_region<- c() # leerer Vektor, in den eingetragen wird, 
    #                                  # ob Übereinstimmung
    #                                  # mit selection event S(thisy) oder 
    #                                  # nicht und somit, ob T_b
    #                                  # in Inferenzregion liegt
    # 
    # for(iyb in 1:length(ybs)){
    # indicator_inference_region[iyb] <- checkFun(ybs[[iyb]])}                                #?# hier parallelisieren?  #?# das liegt vermutlich an der Komplexität von simulation study 2 bzw. an checkFun aber man kriegt da immer singular fit
    
    
    # Übereinstimmung mit initialem selection event überprüfen
    # parallelisierte Version
    #plan(multisession)
    indicator_inference_region <- future_sapply(ybs, function(yb) {
      checkFun(yb)
    }, future.seed = TRUE)
    #plan(sequential)
    
    # Extrahiere alle MC-samples der MC-Teststatistiken Tbs, die in der inference 
    # region liegen, die selection also überlebt haben
    Tbs_survived <- Tbs[indicator_inference_region]
    
    total_resamplings_due_to_convergence <- 0 # since no resampling is done
  }
  
  if (checkConvergence_for_check_Fun == TRUE) {
    
    if (is.null(checkFun_checkConv)) {  
      stop("checkFun_checkConv muss angegeben sein, wenn convergence bei 
         MCsampling überprüft werden soll!") }
    
    #plan(multisession)
    result_list <- future_lapply(seq_len(B_MCsamples), function(b) {
      convergence_okay <- FALSE
      attempts <- 0
      
      while (!convergence_okay) {
        attempts <- attempts + 1
        Tb <- sampFun(num = 1)
        yb <- as.numeric(Tb * dirV + orth_proj)
        selection_convergence_check <- checkFun_checkConv(yb)
        convergence_okay <- selection_convergence_check$convergence_okay
      }
      
      indicator_inference_region <- selection_convergence_check$in_inference_region
      
      list(
        Tb = if (indicator_inference_region) Tb else NA,
        attempts = attempts
      )
    },
    future.seed = TRUE)
   # plan(sequential)
    
    # Extrahiere die Tb-Werte und Anzahl der Versuche
    Tbs_survived <- sapply(result_list, function(res) res$Tb)
    attempts_per_sample <- sapply(result_list, function(res) res$attempts)
    
    # Entferne NA-Einträge
    Tbs_survived <- Tbs_survived[!is.na(Tbs_survived)]
    
    # Anzahl der *zusätzlichen* Samplings wegen convergence issues
    total_resamplings_due_to_convergence <- sum(attempts_per_sample) - B_MCsamples
    
  }
  
  
  # Berechnung der Gewichte für die importance selection--------------------------
  
  # Muss natürlich nur gemacht werden, wenn man nicht die die zentrierte 
  # N(0,kappa_square) als MC_distribution genommen hat, denn fürs zentrierte
  # etzen wir die Gewichte einfach auf 1
  
  if(output == "pvalCI"){ 
    
    if(MC_distribution == "centered"){ 
      wbs <- Rmpfr::mpfr(rep(1, length(Tbs_survived)),2)                                   #?# ich habe das hier auch einfach mal mit dem mpfr gemacht, damit das von der Struktur her immer noch gleich ist
    }
    
    else if(MC_distribution != "centered"){
      
      # berechne Gewichte wb := f_N(0,kappa_square)(Tb)/f_Nproposal(T_b)
      # mit Rmpfr, falls der Nenner zu klein wird und wb sonst als inf betrachtet 
      # würde. Dazu auch die log-exp-Transformation.
      
      # f_N(0,kappa_square)(Tbs)
      log_noms <- Rmpfr::mpfr(dnorm(Tbs_survived, 
                                    mean = 0,
                                    sd = sqrt(kappa_square),
                                    log = TRUE), 200)
      # f_Nproposal(T_bs)
      log_denoms <- Rmpfr::mpfr(dsampFun(tsta = Tbs_survived, log = TRUE), 200)
      
      # Berechne log(wb) = log(nom/denom) = log(nom) - log(denom)
      log_wbs <- log_noms - log_denoms 
      
      # Berechne wbs 
      wbs <- exp(log_wbs)
    }
  }
  
  #-------------------------------------------------------------------------------
  
  if(output == "pvalCI"){
    # Ausgabe der selective-inference p-Werte und KIs mithilfe der pval_CI-Funktion
    return(
      pval_CI(
        Tbs_survived = Tbs_survived,
        t_obs = t_obs,
        wbs = wbs,
        kappa_square = kappa_square,
        rho_H0 = rho_H0,
        symmetric_selection = symmetric_selection,
        alpha = alpha,
        allow_smaller_alpha = allow_smaller_alpha,
        multiple = multiple,
        total_resamplings_due_to_convergence
      )
    )
  }
  
  if(output == "selective_inference"){
    
    return(
      list(
        Tbs_survived = Tbs_survived,
        number_survived_MC = length(Tbs_survived),
        t_obs = t_obs,
        total_resamplings_due_to_convergence =  total_resamplings_due_to_convergence 
      )
    )
  }
  
  
}




