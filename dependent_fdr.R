library(MASS)
library(foreach)
library(doParallel)

RGamma <- function(gamma, test.stat){
  # function that returns the number of test statistics in a rejection region
  #args:
  #<gamma> rejection region boundary
  #<test.stat> test statistics of interest
  #NOTE: use abs(test.stat) if two sided

  sum(test.stat > gamma)
}

WGamma <- function(gamma, test.stat){
  # function that returns the number of test statistics not contained in rejection region
  #args:
  #<gamma> rejection region boundary
  #<test.stat> test statistic of interest
  #NOTE: use abs(test.stat) for two sided tests

  r <- RGamma(gamma, test.stat)
  return(length(test.stat) - r)
}

EstimateFDR <- function(gamma, gamma.prime, test.stat, n.iters=50){
  #function that estimates pFDR for dependent data in a given rejection region
  #args:
  #<gamma> rejection region boundary under consideration
  #<gamma.prime> rejection region likely to contain mostly null statisitcs, selected by OptGamma
  #<test.stat> test statsitics under consideration
  #<n.iters> the number of times to simulate null data (<B> in Storey)
  #NOTE: this function is written such that the <test.stat> is assumed N(0, 1)
  
  n.test.stat <- length(test.stat)
  test.stats.sim <- sapply(1:n.iters, function(i){
     unif <- runif(n.test.stat)
     test.stat <- qexp(unif, 1)
     return(test.stat)
  })
  
  expected.r.null <- 1 / n.iters * sum(apply(test.stats.sim, 2, function(test.stat){
                       RGamma(gamma, test.stat)
                     }))
  prob.r.positive <- 1 / n.iters * sum(apply(test.stats.sim, 2, function(test.stat){                                                                                                
                       return(RGamma(gamma, test.stat) > 0)
                     }))
  expected.w.null <- n.test.stat - 1 / n.iters * sum(apply(test.stats.sim, 2, function(test.stat){
                       RGamma(gamma.prime, test.stat)
                     }))
  
  pi.hat.null <- WGamma(gamma.prime, test.stat) / expected.w.null

  pFDR.hat <- (pi.hat.null * expected.r.null) / (prob.r.positive * max(RGamma(gamma, test.stat), 1))
  FDR.hat <- (pi.hat.null * expected.r.null) / max(RGamma(gamma, test.stat), 1)

  return(list(pFDR=pFDR.hat, FDR=FDR.hat))
}


OptGamma <- function(gamma, lambda.sequence, test.stat, n.iters=50){
    # Implementation of Storey procedure for choosing Gamma* (optimal lambda) based on MSE when
    #observations are dependent
    #args:
    #<gamma> the rejection region being considered
    #<lambda.sequence> a sequence that determines potential Gamma' values for estimating FDR
    #<test.stat> a matrix of test statistics, NOTE: the procedure assumes there is some independent
    #dimension, here, we assume the columns to be independent
    #<n.iters> the number of times to simulate null data when calling EstimatepFDR
  
  #we are again assuming that our test statistics are N(0,1)
  gamma.primes <- qexp(lambda.sequence, 1)
  pFDR.lambda <- sapply(gamma.primes, function(gam.prime){
    estimate <- EstimateFDR(gamma, gam.prime, test.stat, n.iters=50)
    return(estimate$pFDR)
  })


  bias.lambdas.sq <- sapply(pFDR.lambda, function(pFDR){
    (pFDR - min(pFDR.lambda)) ^ 2
  })

  bias.lambdas.sq.adj <- bias.lambdas.sq / median(bias.lambdas.sq)

  # We are assuming our test statisitics come in matrix w/ independent cols
  var.lambdas <- c()
  for (gam.prime in gamma.primes) {
    n <- ncol(test.stat)
    pFDR.heldout <- foreach (i = 1:n) %dopar% {
      pFDR.heldout <- EstimateFDR(gamma, gam.prime, test.stat[, -i], n.iters=50)
      return(pFDR.heldout$pFDR)
    }
    pFDR.heldout <- unlist(pFDR.heldout)
    pFDR.full <- EstimateFDR(gamma, gam.prime, test.stat, n.iters=50)
    var.lambdas <- c(var.lambdas, ((n - 1) / n) * sum( (pFDR.heldout - pFDR.full$pFDR) ^ 2))
  }

  var.lambdas.adj <- var.lambdas / median(var.lambdas)
  MSE.lambdas <- bias.lambdas.sq.adj + var.lambdas.adj
  opt.indx <- which(MSE.lambdas == min(MSE.lambdas))
  opt.gamma <- gamma.primes[opt.indx]
  
  return(opt.gamma)
}

ncores <- 12
registerDoParallel(ncores)

# A simulation using the prostate data
load('prostatedata.Rda')
test.statistics <- prostatedata
test.statistics <- apply(test.statistics, 1, function(gene){
  ((gene - mean(gene)) / sd(gene)) ^ 2
})

gamma <- 5
lambda.seq <- seq(.1, .9, by=.1)

opt.gamma <- OptGamma(gamma, lambda.seq, test.statistics)
pFDR.est <- EstimateFDR(gamma, lambda.seq, test.statistics)

filename = 'fdr.Rdata'
save(opt.gamma, pFDR.est, file=filename)





