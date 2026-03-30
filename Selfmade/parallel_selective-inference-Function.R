################################################################################
# selective-inference Funktionen
################################################################################

library(performance)

# Beschreibung: In diesem Skript werden die MC-samples gezogen und daraus Ceta
#               berechnet und somit die selective-inference-p-Werte und die
#               selective-inference-KIs bestimmt. 
#
# Die erste Funktion (Name: pval_CI ;analog zu selinf im selfmade-package)
# berechnet mithilfe der angegebenen MC-samples p-Werte und Konfidenzintervalle
# basierend auf der approximativen Ceta-Formel (7) [MA Seite 18].
#
#
# Die zweite Funktion (Name: selective_inference ;analog zu pval_vT_cov Funktion  
# im selfmade-package) beinhaltet Ceta-Funktion und zieht vorher die MC-samples
# aus der proposal distribution 

#-------------------------------------------------------------------------------
# pval_CI
#-------------------------------------------------------------------------------

pval_CI <- function(Tbs_survived, # alle MC-Teststatistiken, die in inference 
                    # region liegen
                    t_obs, # beobachteter Testwert
                    wbs, # Gewichte für Approximation von Ceta: Formel (7) in MA
                    kappa_square, # Varianz der Teststatistik
                    rho_H0 = 0, # Nullhypothesenwert des zu testenden Gewichts
                    symmetric_selection = T, # vmtl. symmetrische selection?
                    alpha = 0.05, # welches Signifikanzniveau
                    allow_smaller_alpha = T, # Berechnung von größeren KIs?
                    multiple = 7, # reguliert Sequenzbreite bei rho_sequence
                    total_resamplings_due_to_convergence
){
  
  
  # wie viele MC-samples liegen in inference region
  number_survived_MC <- length(Tbs_survived)
  
  ## Wenn Bestimmung von p-Werten und KI-Grenzen nicht möglich ist
  # falls keine der MC-Teststatistiken Tbs in der inference region liegt, kann
  # weder ein p-Wert noch ein KI berechnet werden
  
  # check ob mind. 1 MC-sample in inference region liegt, gebe ansonsten NA aus
  if(number_survived_MC == 0) {
    warning("Alle MC-samples liegen außerhalb der inference region, daher 
    konnten keine p-Werte und keine KI-Grenzen bestimmt werden. Reduziere 
    rejection rate bei MC-sampling durch diese zwei Optionen:
    (1) Erhöhe B_MC_samples
    (2) Nutze andere MC_distribution als 'centered'")
    return(data.frame(number_survived_MC = 0,
                      pval_right = NA,
                      pval_two = NA,
                      cil = NA, 
                      ciu = NA,
                      cil.smalleralpha=NA,
                      ciu.smalleralpha=NA,
                      total_resamplings_due_to_convergence))
  }
  
  # welche Tbs, die in der inference region liegen, sind größer (gleich) t_obs
  which.Tbs_surv_geq_tobs <- Tbs_survived >= rep(t_obs, length(Tbs_survived))
  Tbs_surv_geq_tobs <- Tbs_survived[which.Tbs_surv_geq_tobs]
  
  if(symmetric_selection == FALSE){
    # welche Absolutwerte der Tbs, die in der inference region liegen, sind größer 
    # (gleich) der Absolutbetrag von t_obs
    which_abs.Tbs_surv_geq_tobs <- abs(Tbs_survived) >= rep(abs(t_obs), 
                                                            length(Tbs_survived))
    abs.Tbs_surv_geq_tobs <- Tbs_survived[which_abs.Tbs_surv_geq_tobs]}
  
  
  # Berechnung der p-Werte--------------------------------------------------------
  
  # p-Wert:
  # eigentlich Berechnung über die Ceta-Formel [Formel (7) in MA] aber wenn man
  # Testung gegen rho_H0=0 vornimmt, dann kann man die Formel vereinfachen, denn
  # exp(0)=1
  
  if(rho_H0 == 0){
    
    # Nenner von Ceta-Formel
    sum_wbs <- sum(wbs)
    
    # rechtsseitiger p-Wert
    pval_right <- as.numeric(sum(wbs[which.Tbs_surv_geq_tobs]) / sum_wbs)
    
    # zweiseitiger p-wert
    # wir definieren diesen anders als im selfmade-package und zwar als
    # P(|T|>=|t_obs| | H0) weil die Ursprungsverteilung der Teststatistik ja mal
    # eine symmetrische Verteilung war (und bei den meisten selection events 
    # auch immer noch ist). Wenn nicht ausschöpfend gesampled wird, ist das 
    # dann zwar immer noch das gleiche Problem, wie bei min(pr,pl) aber bei 
    # ausschöpfendem sampling zumindest nicht
    if(symmetric_selection == TRUE){
      pval_two <- 2* min(pval_right, 1 - pval_right)
    }
    else if(symmetric_selection == FALSE){
      pval_two <- as.numeric(sum(wbs[which_abs.Tbs_surv_geq_tobs]) / sum_wbs)               #!# Problem, wenn nicht ausschöpfend aus distribution gesampled wird
    }
    
  }
  
  else if(rho_H0 != 0){ # wenn rho_H0!=0 können wir nicht mehr die vereinfachte
    # Ceta-Formel verwenden, sondern die normale 
    
    # rechtsseitiger p-Wert
    pval_right <- as.numeric(
      Ceta(wbs = wbs,
           wbs_geq_tobs = wbs[which.Tbs_surv_geq_tobs],
           rho = rho_H0,
           kappa_square = kappa_square,
           Tbs_survived = Tbs_survived,
           Tbs_surv_geq_tobs = Tbs_surv_geq_tobs))
    
    # zweiseitiger p-Wert
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
  
  
  
  # Berechnung der KI-Grenzen-----------------------------------------------------
  
  ## Wenn Bestimmung der (1-alpha)% KI-Grenzen nicht möglich ist
  # Ist dann nicht möglich, wenn alle MC-Tbs, die überlebt haben, größer oder 
  # kleiner als t_obs sind, denn es muss ja ein rho gefunden werden, sodass
  # Ceta bzw. 1-Ceta gleich 1-alpha/2 bzw. alpha/2 wird. Aber wenn alle 
  # MC-sample Tbs größer als t_obs sind, dann ist Ceta = 1 und wenn alle kleiner
  # als t_obs sind, dann ist Ceta = 0.
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
    
    ## Wenn Berechnung der KI-Grenzen möglich ist, muss die Ceta-Formel nach rho
    # umgestellt werden
    
    # Erzeugung von Sequenz von möglichen rho-Werten
    rho_sequence <- seq(min(Tbs_survived ) - multiple*sqrt(kappa_square),
                        max(Tbs_survived ) + multiple*sqrt(kappa_square),
                        l = multiple * 1000)
    
    # Funktion zur Berechnung der Ceta-Werte der rho_sequence-Werte
    # unter den spezifischen vorliegenden MC-samples
    Ceta_function <- function(rho){Ceta(wbs = wbs,
                                        wbs_geq_tobs = wbs[which.Tbs_surv_geq_tobs],
                                        rho = rho,
                                        kappa_square,
                                        Tbs_survived,
                                        Tbs_surv_geq_tobs)}
    
    # Berechnung der unteren und oberen KI-Grenzen: suche Werte rho_low und rho_up       #?# wenn man die nicht-symmetrischen selection events hat, muss man auch hier mal gucken, ob das so Sinn ergibt
    # sodass Ceta(rho_low) = alpha/2 <=> Ceta(rho_low) - alpha/2 = 0
    # und Ceta(rho_up) = 1 - alpha/2 <=> Ceta(rho_up) - (1 - alpha/2) = 0
    
    ## untere Grenze 
    # setze Werte der alpha%-KI-Grenze auf NA, damit es auf jeden Fall Werte
    # gibt, die zurückgegeben werden
    ci_low <- NA
    cil.smalleralpha <- NA
    
    # Berechne den alpha%-KI-Wert der Grenze
    ci_low <- try(uniroot(f = function(rho){
      as.numeric(Ceta_function(rho)) - alpha/2}, 
      interval = range(rho_sequence),
      extendInt = "upX",
      maxiter = 100000)$root, silent = TRUE)
    
    # falls das nicht geklappt hat mit dem alpha%-KI, berechne anderes KI                 #?# ich habe das jetzt so übernommen, wie in dem selfmade-package aber ich verstehe nicht, warum die es so berechnen (ist das überhaupt kleineres alpha - nope)
    if(allow_smaller_alpha == T & class(ci_low) =="try-error"){ 
      #!# ja das ist auf jeden Fall falsch, man könnte es halt wenn noch besser zu min als max machen aber es ist eh falsch... ich lasse es jetzt erstmal so weil ich glaube, dass das Problem so gar nicht so häufig auftritt
      # Berechne die alternative KI-Grenze
      if(t_obs >= rho_H0){
        cil.smalleralpha <- 
          max(rho_sequence[rho_sequence > rho_H0 & rho_sequence < t_obs]) 
      }
      else if(t_obs < rho_H0){
        cil.smalleralpha <- 
          max(rho_sequence[rho_sequence < rho_H0 & rho_sequence < t_obs]) 
      }
      
      # Gebe Warnmeldung aus und berichte, welche Qunatil-Grenze gewählt wurde
      warning("Lower interval limit is not the ",
              alpha/2*100, "%-quantile, but the (",
              signif(as.numeric(Ceta_function(cil.smalleralpha)*100),3), 
              ")%-quantile.") 
    }
    
    ## obere Grenze
    # setze Werte der alpha%-KI-Grenze auf NA, damit es auf jeden Fall Werte
    # gibt, die zurückgegeben werden
    ci_up <- NA
    ciu.smalleralpha <- NA
    
    # Berechne den alpha%-KI-Wert der Grenze
    ci_up <- try(uniroot(f = function(rho){
      as.numeric(Ceta_function(rho)) - (1-alpha/2)}, 
      interval = range(rho_sequence),
      extendInt = "upX",
      maxiter = 100000)$root, silent = TRUE)
    
    # falls das nicht geklappt hat mit dem alpha%-KI, berechne anderes KI                 #?# ich habe das jetzt so übernommen, wie in dem selfmade-package aber ich verstehe nicht, warum die es so berechnen
    if(allow_smaller_alpha == T & class(ci_low) =="try-error"){ 
      # setze Wert der alpha%-KI-Grenze auf NA
      ci_up <- NA
      
      # Berechne die alternative KI-Grenze
      if(t_obs >= rho_H0){
        ciu.smalleralpha <- 
          min(rho_sequence[rho_sequence > rho_H0 & rho_sequence > t_obs]) 
      }
      else if(t_obs < rho_H0){
        ciu.smalleralpha <- 
          min(rho_sequence[rho_sequence < rho_H0 & rho_sequence > t_obs]) 
      }
      
      # Gebe Warnmeldung aus und berichte, welche Qunatil-Grenze gewählt wurde
      warning("Upper interval limit is not the ",
              (1-alpha/2)*100, "%-quantile, but the (",
              signif(100-(as.numeric(Ceta_function(ciu.smalleralpha))*100),3), 
              ")%-quantile.") 
    }
    
    return(data.frame(number_survived_MC = number_survived_MC,                             #?# oder soll ich dazu auch noch ausgeben, welches %-KI das jetzt ist?
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
  
  # Compute the observed test statistic t_obs
  t_obs <- as.numeric(vt %*% this_y) 
  
  # Compute the orthogonal projection := (I - dirV) * vt * y
  orth_proj <- (diag(n) - dirV %*% vt) %*% this_y
  
  # Draw MC samples from a centered normal distribution N(0, kappa_square) or from
  # a proposal distribution N(t_obs, kappa_square) or another user-defined proposal,
  # then check agreement with the initial selection event
  #-----------------------------------------
  
  # Define MC_distribution:                                                             #?# consider adding a check whether the proposal distribution performs well
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
  
  
  # Another option for sampFun could be implemented here                                 #?# possibly extend this to a mixture distribution
  
  if(checkConvergence_for_check_Fun == F){
    
    if(is.null(checkFun)){  
      stop("checkFun must be provided if convergence is not checked during MC sampling!") }
    
    # Draw MC samples
    Tbs <- sampFun(num = B_MCsamples) 
    
    # Transform Tbs to y_bs. For each b: y_b = T_b * dirV + P_ortho
    ybs <- lapply(Tbs, function(Tb) as.numeric(Tb * dirV + orth_proj)) 
    
    # Check agreement with the initial selection event
    # # non-parallel version
    # indicator_inference_region<- c() # empty vector to store whether
    #                                  # the sample matches the selection event S(this_y)
    #                                  # and thus whether T_b lies in the inference region
    # 
    # for(iyb in 1:length(ybs)){
    # indicator_inference_region[iyb] <- checkFun(ybs[[iyb]])}                                #?# parallelize here? likely bottleneck due to complexity of simulation study 2 or checkFun (singular fit issues)
    
    
    # Check agreement with the initial selection event
    # parallelized version
    # indicator_inference_region <- future_sapply(ybs, function(yb) {
    #   checkFun(yb)
    # }, future.seed = TRUE)
    
    indicator_inference_region <- unlist(mclapply(ybs, function(yb) {
      checkFun(yb)
    }, mc.set.seed = TRUE, mc.cores = 10))
    
    # Extract MC samples of test statistics Tbs that lie in the inference region
    # (i.e., that survive the selection event)
    Tbs_survived <- Tbs[indicator_inference_region]
    
    total_resamplings_due_to_convergence <- 0 # since no resampling is done
  }
  
  if (checkConvergence_for_check_Fun == TRUE) {
    
    if (is.null(checkFun_checkConv)) {  
      stop("checkFun_checkConv must be provided if convergence during MC sampling is checked!") }
    
    # result_list <- future_lapply(seq_len(B_MCsamples), function(b) {
    #   convergence_okay <- FALSE
    #   attempts <- 0
    #   
    #   while (!convergence_okay) {
    #     attempts <- attempts + 1
    #     Tb <- sampFun(num = 1)
    #     yb <- as.numeric(Tb * dirV + orth_proj)
    #     selection_convergence_check <- checkFun_checkConv(yb)
    #     convergence_okay <- selection_convergence_check$convergence_okay
    #   }
    #   
    #   indicator_inference_region <- selection_convergence_check$in_inference_region
    #   
    #   list(
    #     Tb = if (indicator_inference_region) Tb else NA,
    #     attempts = attempts
    #   )
    # },
    # future.seed = TRUE)
    
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
    
    # Extract Tb values and number of attempts
    Tbs_survived <- sapply(result_list, function(res) res$Tb)
    attempts_per_sample <- sapply(result_list, function(res) res$attempts)
    
    # Remove NA entries
    Tbs_survived <- Tbs_survived[!is.na(Tbs_survived)]
    
    # Number of *additional* samplings due to convergence issues
    total_resamplings_due_to_convergence <- sum(attempts_per_sample) - B_MCsamples
    
  }
  
  
  # Compute importance sampling weights--------------------------
  
  # Only needed if the centered N(0, kappa_square) is not used.
  # For the centered case, weights are simply set to 1
  if(output == "pvalCI"){ 
    
  if(MC_distribution == "centered"){ 
    wbs <- Rmpfr::mpfr(rep(1, length(Tbs_survived)),2)                                   #?# using mpfr here as well to keep structure consistent
  }
  
  else if(MC_distribution != "centered"){
    
    # Compute weights wb := f_N(0,kappa_square)(Tb)/f_Nproposal(T_b)
    # Use Rmpfr in case the denominator becomes too small and would otherwise
    # be treated as infinity. Apply log-exp transformation for stability.
    
    # f_N(0,kappa_square)(Tbs)
    log_noms <- Rmpfr::mpfr(dnorm(Tbs_survived, 
                                  mean = 0,
                                  sd = sqrt(kappa_square),
                                  log = TRUE), 200)
    # f_Nproposal(T_bs)
    log_denoms <- Rmpfr::mpfr(dsampFun(tsta = Tbs_survived, log = TRUE), 200)
    
    # Compute log(wb) = log(nom/denom) = log(nom) - log(denom)
    log_wbs <- log_noms - log_denoms 
    
    # Compute wbs 
    wbs <- exp(log_wbs)
  }
  }
  
  #-------------------------------------------------------------------------------
  
  if(output == "pvalCI"){
  # Return selective inference p-values and confidence intervals using pval_CI
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


