#setwd("graduate\\stat_215b_wd\\proj")
load("prostatedata.Rda")
X=prostatedata
#dim(X)
#head(X)
#Xt=t(X)
#Xt.sc=scale(Xt)
#X.sc=t(Xt.sc)
#Xp=as.data.frame(X.sc)
#names(Xp)

# 50 
# 52
n0 = 50
n1 = 52
  
X0 = X[,1:50]
X1 = X[,51:102]

#====================================
# 
# functions for getting t-statistics
#
#====================================
# one group of data, observations on columns
# calculate n, xbar(dim p), and sample var (dim p)
group_info = function(X){
  # "patients" on columns
  # Covariates/genes/etc. on rows
  n = ncol(X)
  xbar= rowMeans(X)
  vars=apply(X,1,var)
  
  out=list()
    out[[1]]=n
    out[[2]]=xbar
    out[[3]]=vars
  names(out)=c("n","xbar","svar")  
  
  out
}

# Calculate the t-statistics
t_stats = function( X0, X1){
  a=group_info(X0)
  b=group_info(X1)
  
  n0 = a$n
  n1 = b$n
  dfree=n0+n1-2
  
  s2_pool = 1/(dfree)*( (n0-1)*a$svar + (n1-1)*b$svar )*(1/n0 + 1/n1)
  t.stats = (a$xbar - b$xbar)/sqrt(s2_pool) 
  
  t.stats
}

# convert t-statistics (over df DoF) to z-values
z_vals = function(tstats,df){
  t_two= t.stats^2#-1*t.stats*sign(t.stats) #two-sided
  
  pvals = pf(q=t_two,df1=1,df2=df)
  zvals=qnorm(pvals)
  
  zvals
}

# p-values for t-statistics with df degrees of freedom
p_vals_t = function(tstats, df){
  t_two=-1*t.stats*sign(t.stats) #two-sided
  
  pvals = pt(q=t_two,df=df)
  pvals
}

# p-values for z-statistics
p_vals_z = function(zstats){
  z_two=-1*zstats*sign(zstats) #two-sided
  
  pvals = pnorm(q=z_two)
  pvals
}


#====================================
#
#       relevant calculations
#
#====================================
dfree= n0+n1-2
t.stats = t_stats(X0,X1)
zs = z_vals(t.stats,dfree)
ps.t = p_vals_t(t.stats,dfree)
plot(1:50,sort(ps.t)[1:50]) #not matching efron
# Pretty z-values distribution
hist(zs)

load("prostz.Rda")
Z=prostz
ps.z = p_vals_z(Z)
points(1:50,sort(ps.z)[1:50],pch="x") 
# hmm pretty similar not exact 
# still not matching efron


#====================================
#
#       benj_hoch
#
#====================================
ben_hoch_idx = function(pvals, qval){
  N = length(pvals)
  
  p.ord.idx=order(pvals)
  p.ord = ps.t[p.ord.idx]
  
  if( qval>=1 | qval<=0) stop("q must be in (0,1)")
  
  # Case 1 - one qval
  comp.seq = (1:N)*(qval/N) 
  
  # for scalability could for loop, but not bad here
  comp.test=(p.ord <= comp.seq)
  
  p.bh=max(p.ord[comp.test])
  i.bh = which(p.ord==p.bh)
  
  p.ord.idx[1:i.bh]
}

bh.idx=ben_hoch_idx(ps.t,qval=.1)
