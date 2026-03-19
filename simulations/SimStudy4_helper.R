################################################################################
# SanityCheck: helper functions
################################################################################

## function to get design matrix Z from X --------------------------------------
XzuZ <- function(X,
                 level2.index,
                 random.cols #which cols are random
){
  level2.index<-as.numeric(as.factor((level2.index)))
  X.level2<-cbind(X,level2.index)
  all.level2<-1:max(level2.index)
  if(length(random.cols)==1){X.random=X.level2[,1:2]
  Z <- lapply(all.level2, function(iind){
    X.random[level2.index == iind,]
  })
  Z<-lapply(Z, function(mat) {
    mat <- mat[, -2]  
    return(mat)
  })}
  else{X.random=X.level2[,random.cols]
  Z <- lapply(all.level2, function(iind){
    X.random[level2.index == iind,]
  })}
  
  Z <- do.call(bdiag, Z)
  return(Z)
}

# function to extract data frame of p-values from mocasin-output-object --------
extract_pval_two <- function(data) {
  results <- list()
  
  for (source in c("marginal", "marginal_not_efficient","conditional")) {
    if (!is.null(data[[source]])) {
      entries <- data[[source]]
      
      for (i in seq_along(entries)) {
        entry <- entries[[i]]
        
        if (source == "conditional") {
          entry <- entry[[1]]
        }
        
        
        for(j in 1:length(entry)){
          name <- names(entry[[j]])
          df <- entry[[j]][[1]]
          
          matches <- regmatches(name, 
                                regexec("^(.*?)_coefficient (\\d+)$", name))[[1]]
          type <- matches[2]
          type <- sub("_bayesian$", "", type)
          coefficient <- paste("coefficient", matches[3])
          
          
          pval_two <- df$pval_two[1]
          cil <- df$cil[1]
          ciu <- df$ciu[1]
          cil.smalleralpha <- df$cil.smalleralpha[1]
          ciu.smalleralpha <- df$ciu.smalleralpha[1]
          total_resamplings_due_to_convergence <- df$total_resamplings_due_to_convergence[1]
          
          results[[length(results) + 1]] <- data.frame(
            perspective = source,
            estimator = type,
            coefficient = coefficient,
            pval_two = pval_two,
            cil = cil,
            ciu = ciu,
            cil.smalleralpha = cil.smalleralpha,
            ciu.smalleralpha = ciu.smalleralpha,
            total_resamplings_due_to_convergence,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  do.call(rbind, results)
}


