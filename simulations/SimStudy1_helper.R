################################################################################
# Simulation Study 1: helper functions
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






