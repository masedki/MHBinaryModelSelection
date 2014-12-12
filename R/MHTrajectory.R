MHTrajectory <- function(ae,
                         drug,  
                         vc,
                         maxit,
                         period)
{
  t0 <- proc.time()
  t1 <- proc.time()
  y <- as.vector(ae) 
  x <- as.matrix(drug)  
  v.cramer <- as.vector(vc) 
  n <- nrow(x)
  p <- ncol(x) 
  
  trajectory <- matrix(0, maxit, p)
  acceptedmodel <- rep(0, maxit)
  failedinit <- rep(0, maxit)
  allbic <- rep(0, maxit)
  
  #un premier modèle avant d'entrer dans la boucle
  current.model <- rbinom(p, size = 1, prob = rep(1/2, p))
  X <- x[,(1:p)[current.model==1]]
  vars <- paste("X[,",1:ncol(X),"]",sep="")
  fla <- paste("y~", paste(vars, collapse="+"))
  current.fit <- glm(as.formula(fla), family = binomial, data = data.frame(y=y,X=X))
  twice.log.lik <- -current.fit$aic + 2*length(current.fit$coef)
  current.model.bic <-  twice.log.lik - log(n)*length(current.fit$coef)
  
  iter <- 1 
  best.model <- current.model
  best.model.bic <- current.model.bic
  best.fit <- current.fit
  
  
  trajectory[iter,] <- current.model
  allbic[iter] <- current.model.bic
  acceptedmodel[iter] <- 1
  
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
    candidate.fit <- try(glm(as.formula(fla), family = binomial, data = data.frame(y=y,X=X), start = coef.start), 
                         silent = TRUE)
    if(class(candidate.fit) =="try-error")
    {
      candidate.fit <- glm(as.formula(fla), family = binomial, data = data.frame(y=y,X=X))
      failedinit[iter] <- 1
    }
    
    
    twice.log.lik <- -candidate.fit$aic + 2*length(candidate.fit$coef)
    candidate.model.bic <-  twice.log.lik - log(n)*length(candidate.fit$coef)
    
    trajectory[iter,] <- current.model
    allbic[iter] <- current.model.bic
    #compute of the acceptance probability
    rho <- normalise*exp(candidate.model.bic - current.model.bic)
    if(runif(1) < rho)
    {
      current.model <- candidate.model
      current.model.bic <- candidate.model.bic 
      current.fit <- candidate.fit
      if(current.model.bic > best.model.bic)
      {
        ##print(c("I spent ...",iter," ... iterations to get bic = ...",current.model.bic))  
        best.model <- current.model
        best.model.bic <- current.model.bic
        best.fit <- current.fit
        acceptedmodel[iter] <- 1
        #iter <-0
      }
    }
    
    ## un peu d'affichage
    if((iter %% period) == 0)
    {
      print("<---------------------------------------------------------------------------------->")
      print("<---------------------------------------------------------------------------------->")
      print("<---------------------------------------------------------------------------------->")
      print(paste("we are at iter---->", iter, "........ and bic---->", current.model.bic, sep=""))
      t1 <- proc.time() - t1       
      print(paste("computing time for the last period ---->", t1[3], sep=""))  
      print(paste("failed initialisation---->", sum(failedinit), sep=""))
      
    }
  }
  output <- list(tajectory = NULL, 
                 allbic = NULL, 
                 failedinit = NULL, 
                 acceptedmodel = NULL,
                 best.model = NULL, 
                 best.model.bic = NULL, 
                 best.fit = NULL,
                 time = NULL)
  
  output$trajectory <- trajectory
  output$allbic <- allbic
  output$failedinit <- failedinit
  output$acceptedmodel <- acceptedmodel
  output$best.model <- best.model
  output$best.model.bic <- best.model.bic
  output$best.fit <- best.fit
  output$time <- proc.time() - t0
  return(output)
}