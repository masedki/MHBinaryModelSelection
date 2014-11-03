# cette fonction génère une pseudo-matrice de consommation de médicaments
# pour le moment je génère suivant une seule stratégie, il faut ajouter un arguement "type" pour diverses stratégies
GenerateBinaryData <-function(n, 
                              probs)
  {
 
  d <- length(probs)
  p <- probs
  x <- matrix(rbinom(n*d, size=1, prob=p), byrow = T, n, d)
  coef.s <- seq(2,0.05, len=round(d/2))
  coef.ns <- rep(0, (d - round(d/2)))
  intercept <- 0.1
  reg.coef <- c(intercept, coef.s, coef.ns)
  z <- cbind(rep(1,n), x) %*% reg.coef  # linear combination with intercept
  pr <- 1/(1+exp(-z)) 
  y <- rbinom(n,1,pr)  
  df <- data.frame(y=y,x=x)
  return(df)
}