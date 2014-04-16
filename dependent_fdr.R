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

EstimatepFDR <- function(gamma, gamma.prime, test.stat, n.iters){
  #function that estimates pFDR for dependent data in a given rejection region
  #args:
  #<gamma> rejection region boundary under consideration
  #<gamma.prime> rejection region likely to contain mostly null statisitcs, selected by OptGamma
  #<test.stat> test statsitics under consideration
  #<n.iters> the number of times to simulate null data (<B> in Storey)
  #<n.test.stat> the number of test statistics to be simulated in each iteration (<m> in Storey)
  #NOTE: this function is written such that the <test.stat> is assumed N(0, 1)
  
  n.test.stat <- length(test.stat)

  test.stats.sim <- sapply(1:n.iters, function(iter) {
    rnorm(n.test.stat)
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