################################################################################
# Selective-inference functions
################################################################################

library(performance)

# Description: This script draws the Monte Carlo samples and uses them to
#              compute Ceta, selective-inference p-values, and
#              selective-inference confidence intervals.
#
# The first function (pval_CI; analogous to selinf in the selfmade package)
# computes p-values and confidence intervals from the supplied Monte Carlo
# samples based on the approximate Ceta formula (7) [MA, p. 18].
#
# The second function (selective_inference; analogous to the pval_vT_cov
# function in the selfmade package) contains the Ceta function and first draws
# the Monte Carlo samples from the proposal distribution.

#-------------------------------------------------------------------------------
# pval_CI
#-------------------------------------------------------------------------------

pval_CI <- function(Tbs_survived, # All MC test statistics that fall within the
                    # inference region
                    t_obs, # Observed test statistic
                    wbs, # Weights for approximating Ceta: formula (7) in MA
                    kappa_square, # Variance of the test statistic
                    rho_H0 = 0, # Null-hypothesis value of the tested parameter
                    symmetric_selection = T, # Whether the selection event is symmetric
                    alpha = 0.05, # Significance level
                    allow_smaller_alpha = T, # Whether to compute alternative CI bounds
                    multiple = 7, # Controls the width of the rho sequence
                    total_resamplings_due_to_convergence
){
  
  
  # Number of MC samples that fall within the inference region
  number_survived_MC <- length(Tbs_survived)
  
  ## If p-values and CI bounds cannot be determined
  # If none of the MC test statistics falls within the inference region,
  # neither a p-value nor a CI can be computed.
  
  # Check whether at least one MC sample falls within the inference region;
  # otherwise, return NA.
  if(number_survived_MC == 0) {
    warning("All MC samples fall outside the inference region, so no p-values
    or CI bounds could be determined. Reduce the rejection rate during MC
    sampling using one of the following options:
    (1) Increase B_MC_samples
    (2) Use an MC_distribution other than 'centered'")
    return(data.frame(number_survived_MC = 0,
                      t_obs = t_obs,
                      kappa_square = kappa_square,
                      pval_right = NA,
                      pval_two = NA,
                      cil = NA, 
                      ciu = NA,
                      cil.smalleralpha=NA,
                      ciu.smalleralpha=NA,
                      total_resamplings_due_to_convergence))
  }
  
  # Identify the survived Tbs that are greater than or equal to t_obs
  which.Tbs_surv_geq_tobs <- Tbs_survived >= rep(t_obs, length(Tbs_survived))
  Tbs_surv_geq_tobs <- Tbs_survived[which.Tbs_surv_geq_tobs]
  
  if(symmetric_selection == FALSE){
    # Identify the survived Tbs whose absolute values are greater than or equal
    # to the absolute value of t_obs
    which_abs.Tbs_surv_geq_tobs <- abs(Tbs_survived) >= rep(abs(t_obs), 
                                                            length(Tbs_survived))
    abs.Tbs_surv_geq_tobs <- Tbs_survived[which_abs.Tbs_surv_geq_tobs]}
  
  
  # Computation of p-values-----------------------------------------------------
  
  # p-value:
  # The Ceta formula (formula (7) in MA) would normally be used. However, when
  # testing rho_H0 = 0, the formula simplifies because exp(0) = 1.
  
  if(rho_H0 == 0){
    
    # Denominator of the Ceta formula
    sum_wbs <- sum(wbs)
    
    # One-sided (right-tailed) p-value
    pval_right <- as.numeric(sum(wbs[which.Tbs_surv_geq_tobs]) / sum_wbs)
    
    # Two-sided p-value
    # We define the two-sided p-value differently from the selfmade package,
    # namely as P(|T| >= |t_obs| | H0), because the original distribution of
    # the test statistic was symmetric and remains symmetric for most selection
    # events. If the distribution is not exhaustively sampled, this still suffers
    # from the same issue as min(pr, pl), but this issue does not arise under
    # exhaustive sampling.
    if(symmetric_selection == TRUE){
      pval_two <- 2* min(pval_right, 1 - pval_right)
    }
    else if(symmetric_selection == FALSE){
      pval_two <- as.numeric(sum(wbs[which_abs.Tbs_surv_geq_tobs]) / sum_wbs)
      # NOTE: This may be problematic if the distribution is not sampled
      # exhaustively.
    }
    
  }
  
  else if(rho_H0 != 0){ # For rho_H0 != 0, the simplified Ceta formula cannot
    # be used, so the general formula is required.
    
    # One-sided (right-tailed) p-value
    pval_right <- as.numeric(
      Ceta(wbs = wbs,
           wbs_geq_tobs = wbs[which.Tbs_surv_geq_tobs],
           rho = rho_H0,
           kappa_square = kappa_square,
           Tbs_survived = Tbs_survived,
           Tbs_surv_geq_tobs = Tbs_surv_geq_tobs))
    
    # Two-sided p-value
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
  
  
  
  # Computation of CI bounds----------------------------------------------------
  
  ## If the (1-alpha)% CI bounds cannot be determined
  # This occurs if all survived MC test statistics are either larger or smaller
  # than t_obs. A value of rho must be found such that Ceta or 1-Ceta equals
  # 1-alpha/2 or alpha/2, respectively. If all MC test statistics are larger
  # than t_obs, Ceta = 1; if all are smaller than t_obs, Ceta = 0.
  
  if(all(Tbs_survived > rep(t_obs,number_survived_MC)) |
     all(Tbs_survived < rep(t_obs,number_survived_MC)) ){
    return(data.frame(number_survived_MC = number_survived_MC,
                      t_obs = t_obs,
                      kappa_square = kappa_square,
                      pval_right = pval_right,
                      pval_two = pval_two,
                      cil = NA, 
                      ciu = NA,
                      cil.smalleralpha=NA,
                      ciu.smalleralpha=NA,
                      total_resamplings_due_to_convergence))
    
  }
  
  ## If the CI bounds can be computed, the Ceta formula must be inverted with
  # respect to rho.
  
  # Generate a sequence of possible rho values
  rho_sequence <- seq(min(Tbs_survived ) - multiple*sqrt(kappa_square),
                      max(Tbs_survived ) + multiple*sqrt(kappa_square),
                      l = multiple * 1000)
  
  # Function for computing Ceta values for the rho sequence
  # using the observed MC samples
  Ceta_function <- function(rho){Ceta(wbs = wbs,
                                      wbs_geq_tobs = wbs[which.Tbs_surv_geq_tobs],
                                      rho = rho,
                                      kappa_square,
                                      Tbs_survived,
                                      Tbs_surv_geq_tobs)}
  
  # Compute the lower and upper CI bounds by finding rho_low and rho_up such
  # that Ceta(rho_low) = alpha/2 <=> Ceta(rho_low) - alpha/2 = 0
  # and Ceta(rho_up) = 1 - alpha/2 <=> Ceta(rho_up) - (1 - alpha/2) = 0.
  # NOTE: For non-symmetric selection events, this should be checked carefully.
  
  ## Lower bound
  # Initialize the alpha-level CI bound as NA to ensure that a value is
  # returned even if the root-finding procedure fails.
  ci_low <- NA
  cil.smalleralpha <- NA
  
  # Compute the alpha-level CI bound
  ci_low <- try(uniroot(f = function(rho){
    as.numeric(Ceta_function(rho)) - alpha/2}, 
    interval = range(rho_sequence),
    extendInt = "upX",
    maxiter = 100000)$root, silent = TRUE)
  
  # If the alpha-level CI bound cannot be obtained, compute an alternative
  # interval bound.
  # NOTE: This was adopted from the selfmade package, but it is unclear whether
  # this procedure actually corresponds to a smaller-alpha confidence bound.
  if(allow_smaller_alpha == T & class(ci_low) =="try-error"){ 
    # NOTE: This procedure may not be theoretically correct. It could perhaps
    # be modified to choose the value closer to the minimum rather than the
    # maximum, but the issue is left unchanged for now because it is expected
    # to occur relatively rarely.
    if(t_obs >= rho_H0){
      cil.smalleralpha <- 
        max(rho_sequence[rho_sequence > rho_H0 & rho_sequence < t_obs]) 
    }
    else if(t_obs < rho_H0){
      cil.smalleralpha <- 
        max(rho_sequence[rho_sequence < rho_H0 & rho_sequence < t_obs]) 
    }
    
    # Issue a warning and report which quantile was selected instead.
    warning("Lower interval limit is not the ",
            alpha/2*100, "%-quantile, but the (",
            signif(as.numeric(Ceta_function(cil.smalleralpha)*100),3), 
            ")%-quantile.") 
  }
  
  ## Upper bound
  # Initialize the alpha-level CI bound as NA to ensure that a value is
  # returned even if the root-finding procedure fails.
  ci_up <- NA
  ciu.smalleralpha <- NA
  
  # Compute the alpha-level CI bound
  ci_up <- try(uniroot(f = function(rho){
    as.numeric(Ceta_function(rho)) - (1-alpha/2)}, 
    interval = range(rho_sequence),
    extendInt = "upX",
    maxiter = 100000)$root, silent = TRUE)
  
  # If the alpha-level CI bound cannot be obtained, compute an alternative
  # interval bound.
  # NOTE: This was adopted from the selfmade package, but it is unclear whether
  # this procedure actually corresponds to a smaller-alpha confidence bound.
  if(allow_smaller_alpha == T & class(ci_low) =="try-error"){ 
    # Set the alpha-level CI bound to NA.
    ci_up <- NA
    
    # Compute the alternative CI bound
    if(t_obs >= rho_H0){
      ciu.smalleralpha <- 
        min(rho_sequence[rho_sequence > rho_H0 & rho_sequence > t_obs]) 
    }
    else if(t_obs < rho_H0){
      ciu.smalleralpha <- 
        min(rho_sequence[rho_sequence < rho_H0 & rho_sequence > t_obs]) 
    }
    
    # Issue a warning and report which quantile was selected instead.
    warning("Upper interval limit is not the ",
            (1-alpha/2)*100, "%-quantile, but the (",
            signif(100-(as.numeric(Ceta_function(ciu.smalleralpha))*100),3), 
            ")%-quantile.") 
  }
    
  return(data.frame(number_survived_MC = number_survived_MC,  
                    t_obs = t_obs,
                    kappa_square = kappa_square,
                    pval_right = pval_right,
                    pval_two = pval_two,
                    cil = ci_low, 
                    ciu = ci_up,
                    cil.smalleralpha = cil.smalleralpha,
                    ciu.smalleralpha = ciu.smalleralpha,
                    total_resamplings_due_to_convergence))
}


#-------------------------------------------------------------------------------
# selective_inference
#-------------------------------------------------------------------------------
# if (!("future" %in% rownames(installed.packages()))) {
#   install.packages("future")
# }
# if (!("future.apply" %in% rownames(installed.packages()))) {
#   install.packages("future.apply")
# }
# library(future.apply)   
# plan(multisession)

if (!("parallel" %in% rownames(installed.packages()))) {
  install.packages("parallel")
}
library(parallel) 

selective_inference <- function(n,
                                vt,
                                kappa_square,
                                dirV,
                                this_y,
                                rho_H0 = 0,
                                B_MCsamples = 1000,
                                min_B_MCsamples = 0,
                                checkFun,
                                MC_distribution = "proposal",
                                symmetric_selection = T,
                                alpha = 0.05, 
                                allow_smaller_alpha = T, 
                                multiple = 7,
                                checkConvergence_for_check_Fun,
                                checkFun_checkConv,
                                output = "pvalCI", # or "selective_inference"
                                definition_sampFun = NULL, # function(num) for random generation
                                definition_dsampFun = NULL # function(tsta, log) for the density
){
  
  # Compute the observed test statistic t_obs
  t_obs <- as.numeric(vt %*% this_y) 
  
  # Compute the orthogonal projection := (I - dirV * vt) * y
  orth_proj <- (diag(n) - dirV %*% vt) %*% this_y
  
  # Draw MC samples from a centered normal distribution N(0, kappa_square),
  # from the proposal distribution N(t_obs, kappa_square), or from a proposal
  # distribution composed of multiple normal distributions, and check whether
  # the samples satisfy the original selection event.
  
  # Specify the MC distribution:
  # NOTE: Add a check for whether the proposal distribution is sufficiently
  # efficient.
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
  
  
  # Another option for sampFun could be implemented here.
  # NOTE: Consider implementing a mixture distribution here.
  
  
  if(checkConvergence_for_check_Fun == FALSE){
    
    if(is.null(checkFun)){  
      stop("checkFun must be provided when convergence is not checked during
         MC sampling!") 
    }
    
    Tbs_survived <- c()
    MC_resamplings <- 0
    
    while(length(Tbs_survived) < min_B_MCsamples|MC_resamplings==0){
      
      # Determine how many samples are still missing
      n_missing <- min_B_MCsamples - length(Tbs_survived)
      
      # Optionally draw additional samples for greater efficiency
      n_to_sample <- max(10*n_missing, B_MCsamples)
      
      # Draw MC samples
      Tbs <- sampFun(num = n_to_sample) 
      
      # Transformation
      # ybs <- lapply(Tbs, function(Tb) as.numeric(Tb * dirV + orth_proj)) 
      # 
      # 
      # # Parallel check
      # indicator_inference_region <- unlist(mclapply(ybs, function(yb) {
      #   checkFun(yb)
      # }, mc.set.seed = TRUE, mc.cores = 30))
      
      
      indicator_inference_region <- unlist(
        mclapply(
          Tbs,
          function(Tb) {
            yb <- Tb * dirV + orth_proj
            checkFun(as.numeric(yb))
          },
          mc.cores = 10
        )
      )
      
      # Extract the samples that fall within the inference region
      Tbs_new_survived <- Tbs[indicator_inference_region]
      
      # Append the survived samples
      Tbs_survived <- c(Tbs_survived, Tbs_new_survived)
      
      MC_resamplings <- MC_resamplings + 1
      
      # Optional safeguard against an infinite loop
      if(MC_resamplings >1000000){
        warning("Too many resampling attempts; the inference region may be very small.")
        break
      }
    }
    
    total_resamplings_due_to_convergence <- 0 # No resampling is performed
  }
  
  # if (checkConvergence_for_check_Fun == TRUE) {
  #   
  #   if (is.null(checkFun_checkConv)) {  
  #     stop("checkFun_checkConv must be provided when convergence is checked during
  #        MC sampling!") }
  #   
  #   
  #   result_list <- mclapply(seq_len(B_MCsamples), function(b) {
  #     convergence_okay <- FALSE
  #     attempts <- 0
  #     
  #     while (!convergence_okay) {
  #       attempts <- attempts + 1
  #       Tb <- sampFun(num = 1)
  #       yb <- as.numeric(Tb * dirV + orth_proj)
  #       selection_convergence_check <- checkFun_checkConv(yb)
  #       convergence_okay <- selection_convergence_check$convergence_okay
  #     }
  #     
  #     indicator_inference_region <- selection_convergence_check$in_inference_region
  #     
  #     list(
  #       Tb = if (indicator_inference_region) Tb else NA,
  #       attempts = attempts
  #     )
  #   }, mc.set.seed = TRUE, mc.cores = 10)
  #   
  #   # Extract the Tb values and number of attempts
  #   Tbs_survived <- sapply(result_list, function(res) res$Tb)
  #   attempts_per_sample <- sapply(result_list, function(res) res$attempts)
  #   
  #   # Remove NA entries
  #   Tbs_survived <- Tbs_survived[!is.na(Tbs_survived)]
  #   
  #   # Number of *additional* samples drawn due to convergence issues
  #   total_resamplings_due_to_convergence <- sum(attempts_per_sample) - B_MCsamples
  #   
  # }
  
  if (checkConvergence_for_check_Fun == TRUE) {
    
    if (is.null(checkFun_checkConv)) {  
      stop("checkFun_checkConv must be provided when convergence is checked during
         MC sampling!") 
    }
    
    Tbs_survived <- c()
    total_resamplings_due_to_convergence <- 0
    MC_batches <- 0
    
    while(length(Tbs_survived) < min_B_MCsamples| MC_batches == 0){
      
      result_list <- mclapply(seq_len(B_MCsamples), function(b) {
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
      }, mc.set.seed = TRUE, mc.cores = 10)
      
      # Extract results
      Tbs_batch <- sapply(result_list, function(res) res$Tb)
      attempts_per_sample <- sapply(result_list, function(res) res$attempts)
      
      # Retain only survived samples
      Tbs_survived <- c(Tbs_survived, Tbs_batch[!is.na(Tbs_batch)])
      
      # Count additional resampling attempts
      total_resamplings_due_to_convergence <- 
        total_resamplings_due_to_convergence + (sum(attempts_per_sample) - B_MCsamples)
      
      MC_batches <- MC_batches + 1
      
      # Optional safeguard against an infinite loop
      if(MC_batches >1000000){
        warning("Too many batches; the inference region may be very small.")
        break
      }
    }
  }
  
  
  # Compute weights for importance sampling-------------------------------------
  
  # This is only necessary when the centered N(0, kappa_square) is not used as
  # the MC distribution. For the centered distribution, all weights are set to 1.
  if(output == "pvalCI"){ 
    
    if(MC_distribution == "centered"){ 
      wbs <- Rmpfr::mpfr(rep(1, length(Tbs_survived)),2)
      # NOTE: mpfr is used here so that the data structure remains consistent
      # across the different sampling distributions.
    }
    
    else if(MC_distribution != "centered"){
      
      # Compute weights wb := f_N(0,kappa_square)(Tb) / f_Nproposal(Tb).
      # Rmpfr is used in case the denominator becomes very small and the weight
      # would otherwise be represented as Inf. The computation is performed on
      # the log scale and then transformed back.
      
      # f_N(0,kappa_square)(Tbs)
      log_noms <- Rmpfr::mpfr(dnorm(Tbs_survived, 
                                    mean = 0,
                                    sd = sqrt(kappa_square),
                                    log = TRUE), 200)
      # f_Nproposal(Tbs)
      log_denoms <- Rmpfr::mpfr(dsampFun(tsta = Tbs_survived, log = TRUE), 200)
      
      # Compute log(wb) = log(numerator/denominator) = log(numerator) -
      # log(denominator)
      log_wbs <- log_noms - log_denoms 
      
      # Compute wbs
      wbs <- exp(log_wbs)
    }
  }
  
  #-------------------------------------------------------------------------------
  
  if(output == "pvalCI"){
    # Return selective-inference p-values and CIs using the pval_CI function
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
        total_resamplings_due_to_convergence = total_resamplings_due_to_convergence 
      )
    )
  }
  
  
  
}

