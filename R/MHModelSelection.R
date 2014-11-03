MHModelSelection <- function(data, 
                             maxit)
{
  
  y <- as.vector(data$y)  
  x <- as.matrix(data[2:dim(data)[2]], dim(data)[1], (dim(data)[2]-1))  
  n <- nrow(x)
  p <- ncol(x) 
  #compute v.cramer
  v.cramer <- rep(NA, p)
  for(k in 1:p)
  {
    mat <- table(x[,k], y)
    v.cramer[k] <- sqrt((sum( mat**2/ ((rowSums(mat)%*%t(colSums(mat)))/sum(mat)) ) - sum(mat))/sum(mat))
  }
  #un premier modèle avant d'entrer dans la boucle
  current.model <- rbinom(p, size = 1, prob = rep(1/2, p))
  X <- x[,(1:p)[current.model==1]]
  vars <- paste("X[,",1:ncol(X),"]",sep="")
  fla <- paste("y~", paste(vars, collapse="+"))
  current.fit <- glm(as.formula(fla), family = binomial, data = data.frame(y=y,X=X))
  twice.log.lik <- -current.fit$aic + 2*length(current.fit$coef)
  current.model.bic <-  twice.log.lik - log(n)*length(current.fit$coef)
  
  iter <- 0 
  best.model <- current.model
  best.model.bic <- current.model.bic
  best.fit <- current.fit
  while(iter < maxit)
  {
    iter <- iter + 1
    add <- sample(c(FALSE,TRUE), size=1)
    probatire <- rep(0, p)
    alternative <- rep(0,p)
    candidate.model <- current.model
    if(add)
    {
      probatire[candidate.model==0] <- v.cramer[candidate.model==0]
      alternative[candidate.model==1] <- 1 - v.cramer[candidate.model==1]
    }
    else
    {
      probatire[candidate.model==1] <- 1-v.cramer[candidate.model==1]
      alternative[candidate.model==0] <- v.cramer[candidate.model==0]
    }
    
    probatire <- probatire/sum(probatire)
    drug  <-sample(p, size = 1, prob = probatire)  
    if(add)
      alternative[drug] <- v.cramer[drug]  
    else
      alternative[drug] <- 1 - v.cramer[drug]
    
    alternative<- alternative/sum(alternative)
    candidate.model[drug] <- 1 - candidate.model[drug]
    candidate.model.bic <- log(0)
    normalise <- alternative[drug]/probatire[drug]
    
    #compute candidate.model.bic using the same scheme like the first one
    X <- x[,(1:p)[candidate.model==1]]
    vars <- paste("X[,",1:ncol(X),"]",sep="")
    fla <- paste("y~", paste(vars, collapse="+"))
    current.coef <- rep(0,p)
    current.coef[(1:p)[current.model==1]] <- current.fit$coef[-1]
    current.coef <- c(current.fit$coef[1], current.coef)
    coef.start <- c(current.coef[1],current.coef[-1][(1:p)[candidate.model==1]])
    candidate.fit <- glm(as.formula(fla), family = binomial, data = data.frame(y=y,X=X), start = coef.start)
    twice.log.lik <- -candidate.fit$aic + 2*length(candidate.fit$coef)
    candidate.model.bic <-  twice.log.lik - log(n)*length(candidate.fit$coef)
    
    #compute of the acceptance probability
    rho <- normalise*exp(candidate.model.bic - current.model.bic)
    if(runif(1) < rho)
    {
      current.model <- candidate.model
      current.model.bic <- candidate.model.bic 
      current.fit <- candidate.fit
      if(current.model.bic > best.model.bic)
      {
        print(c("I spent ...",iter," ... iterations to get bic = ...",current.model.bic))  
        best.model <- current.model
        best.model.bic <- current.model.bic
        best.fit <- current.fit
        iter <-0
      }
    }
  }
  best.return <- list(best.model=NULL, best.model.bic=NULL, best.fit=NULL)
  best.return$best.model <- best.model
  best.return$best.model.bic <- best.model.bic
  best.return$best.fit <- best.fit
  return(best.return)
}