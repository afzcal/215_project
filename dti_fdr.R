load('DTIdata.Rda')
p.vals <- pnorm(DTIdata[, 4])

#seemse like we want to do this to account for 1 v 2 sided issues
p.vals[p.vals > .5] <- 1 - p.vals[p.vals > .5]

# Storey procedure for calculating optimal lambda
gamma <- c(.01)
m <- length(p.vals)
r <- seq(0, .475, by = .025)

WLambda <- function(lambda, p.vals){
  sum(p.vals > lambda)
}

RGamma <- function(gamma, p.vals){
  sum(p.vals < gamma)
}

pFDR <- function(r, gamma, p.vals){ 
  # a function that returns estimates of FDR as described in Storey procedure
  # args:
  #   <r> a sequence of lambda values
  #   <gamma> a sequence of rejection region boundaries
  #   <p.vals> the sequence of p-values under consideration

  FDRs <- sapply(r, function(lambda){
      (WLambda(lambda, p.vals) * gamma) / ( (1 - lambda) * max(RGamma(gamma,
        p.vals), 1) * (1 - (1 - gamma) ^ m))
  })
  return(FDRs)
}

# a matrix of estimated positive FRD as a function of lambda and gamma.  Rows 
  #correspond to levels of <gammas> and cols to levels of <r> 
sample.pFDR <- pFDR(r, gamma, p.vals)
     	  

# Generate <B> bootsrap sample of pFDR to get estimate of MSE(lambda) to find
# optimal value for lambda <lambda.hat>
B <- 1000

bs.pFDR <- matrix(rep(NA, B * length(r)), ncol = B)
for (i in 1:B) {
    bs.sample <- sample(p.vals, length(p.vals), replace = T)
    temp.pFDR <- pFDR(r, gamma, bs.sample)
    bs.pFDR[, i] <- temp.pFDR
}

lambda.MSE <- apply(bs.pFDR, 1, function(lambda.pFDR){
  1 / B * sum((lambda.pFDR - min(sample.pFDR)) ^ 2)
})

opt.lambda.indx <- which(lambda.MSE == min(lambda.MSE))
lambda.hat <- r[opt.lambda.indx]

# calculate the 1 - <alpha> upper CI for estimate of pFDR
alpha <- .05
upper.bound <- quantile(bs.pFDR[opt.lambda.indx, ], 1 - alpha)
pFDR.lambda.hat <- sample.pFDR[opt.lambda.indx]
