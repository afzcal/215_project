load('DTIdata.Rda')

# Functions used to implement Storey procedure for finding FDR confidence interval
#################################################################################################
#given a rejection region
WLambda <- function(lambda, p.values){
  sum(p.values > lambda)
}

RGamma <- function(gamma, p.values){
  sum(p.values < gamma)
}

pFDR <- function(p.values, r, gamma){
  # a function that returns estimates of FDR as described in Storey procedure
  # args:
  #   <r> a sequence of lambda values
  #   gamma a rejection region boundary
  #   <p.vals> the sequence of p-values under consideration
  m <- length(p.values)

  FDRs <- sapply(r, function(lambda){
      (WLambda(lambda, p.values) * gamma) / ( (1 - lambda) * max(RGamma(gamma,
        p.values), 1) * (1 - (1 - gamma) ^ m))
  })
  return(FDRs)
}

BootstrapFDR <- function(p.values, n.bs.samples, r, gamma){
    # boostrap to generate an estimate of MSE(lambda) as desbribed in Storey
    #args:
    #  <p.values> the sequence of p-values under consideration
    #  <n.bs.samples> the number of boostrap samples to take
    #  <r> sequence of lambda values
    bs.pFDR <- matrix(rep(NA, n.bs.samples * length(r)), ncol = n.bs.samples)
    
    for (i in 1:n.bs.samples) {
        bs.sample <- sample(p.values, length(p.values), replace = T)
        temp.pFDR <- pFDR(bs.sample, r, gamma)
        bs.pFDR[, i] <- temp.pFDR
    }
    return(bs.pFDR)
}

#Implement Storey procedure
#################################################################################################
#seemse like we want to do this to account for 1 v 2 sided issues
p.vals <- pnorm(DTIdata[, 4])
#p.vals[p.vals > .5] <- 1 - p.vals[p.vals > .5]
gamma <- .001
r <- seq(0, 1 - .025, by = .01)
n.bootstrap <- 1000

# a matrix of estimated positive FRD as a function of lambda and gamma.  Rows
  #correspond to levels of <gammas> and cols to levels of <r> 
sample.pFDR <- pFDR(p.vals, r, gamma)
     	  
# boostrapped estimates of pFDR at various levels of lambda
bs.pFDR <- BootstrapFDR(p.vals, n.bootstrap, r, gamma)

lambda.MSE <- apply(bs.pFDR, 1, function(lambda.pFDR){
  1 / n.bootstrap * sum((lambda.pFDR - min(sample.pFDR)) ^ 2)
})

opt.lambda.indx <- which(lambda.MSE == min(lambda.MSE))
lambda.hat <- r[opt.lambda.indx]
pi.lambda <- WLambda(lambda.hat, p.vals) / ((1 - lambda.hat) * length(p.vals))

# calculate the 1 - <alpha> upper CI for estimate of pFDR
alpha <- .05
pFDR.lambda.hat <- sample.pFDR[opt.lambda.indx]

temp <- BootstrapFDR(p.vals, 10000, lambda.hat, gamma)
upper.bound <- quantile(temp, 1 - alpha)
