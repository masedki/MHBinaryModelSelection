vcramer <- function(x,y)
  {
  if((!is.matrix(x)) ||(!is.vector(y)))
  stop(paste(sQuote("x and y"), "must be a vectors!"))
  
  if(nrow(x) != length(y))
    stop(paste(sQuote("x and y"), "must have the  length!"))
  
  
  vc <- rep(NA, ncol(x))
  for(k in 1:length(vc))
  {
    mat <- table(x[,k], y)
    vc[k] <- sqrt((sum( mat**2/ ((rowSums(mat)%*%t(colSums(mat)))/sum(mat)) ) - sum(mat))/sum(mat))
  }
  return(vc) 
  }