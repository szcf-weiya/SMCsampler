####################
## paras:
##

KM = 4 # number of components in mixture
POP = 1000 # number of particles
R = 100 # number of data points
R1 = 50
MAX = 100 # number of distns ??????
PIN = 1 # invariant dist pi_n (1) or pi_n-1 (0)

## ##################
##
## simulated data points
## ##################

generateData <- function(ntotal, prop, mu, sigma)
{
    part = round(ntotal*prop);
    res = c()
    for (i in 1:length(prop))
    {
        tmp = rnorm(part[i], mu[i], sigma[i])
        res = c(res, tmp)
    }
    return(res)
}

data = generateData(R, rep(1/4, 4), c(-3, 0, 3, 6), rep(0.55, 4))


## ###################
## prior on mean
## 
## mu_j \sim N(xi, kappa^-1)
## ###################

mu.prior.xi = 0.5*(max(data) + min(data))
mu.prior.kappa = 1/(max(data) - min(data))^2
mu.prior.kappa1 = 1/sqrt(mu.prior.kappa)

## ###################
## prior on precision
## 
## precision \sim Gamma(alpha, beta)
## 
## where beta \sim Gamma(g, h)
## ###################

lambda.prior.alpha = 2.0 # fixed?
lambda.prior.beta.g = 0.2
lambda.prior.beta.h = 10*mu.prior.kappa

## ###################
## prior on weight
## 
## weight \sim Dirichlet(delta)
##
## can be simulate from gamma(delta, 1) and standardize
## ###################

w.prior.delta = 1.0 # fixed?

## ##################
## an SMC algorithm with MAX densities
## 
## ##################

beta = numeric(POP)
li = numeric(POP) # log-likelihood
mu = matrix(nrow = KM, ncol = POP)
lambda = matrix(nrow = KM, ncol = POP)
w = matrix(nrow = KM, ncol = POP)
weight = numeric(POP) # unnormalized importance weight
weight1 = numeric(POP)
vec = numeric(POP) # incremental weight
flag.resample = FALSE
mapmean = numeric(KM)
maplambda = numeric(KM)
mapw = numeric(KM)

for (n in 0:(MAX-1))
{
  ## initialize
  if (n == 0)
  {
    for (i in 1:POP)
    {
        beta[i] = lambda.prior.beta.g/lambda.prior.beta.h
        mu[,i] = mu.prior.xi #+TODO 
        lambda[,i] = rgamma(KM, lambda.prior.alpha, beta[i])
        w[,i] = rgamma(KM, w.prior.delta, 1.0)
        w[,i] = w[,i]/sum(w[,i])
        li.tmp = sapply(1:R, function(i) mixden(data[i], mu[,i], lambda[,i], w[,i], KM))
        li[i] = sum(li.tmp)
        weight[i] = 0
    }
  }
  else
  {
    if (n == 1)
    {
      gamma = 15/(100*20)
      gammad = 0
    }
    if (n < 20 && n != 1)
    {
      gammad = gamma
      gamma = gamma + 15/(100*20)
    }
    else if(n < 60 && n >= 20)
    {
      gammad = gamma
      gamma = gamma + 25/(100*40)
    }
    else
    {
      gammad = gamma 
      gamma = gamma + 60/(100*(MAX-60))
    }
    if (n == MAX-1)
    {
      gammad = gamma 
      gamma = 1.0
    }
    ## perform a cycle of moves
    if (n > 60)
    {
      if (n > 80)
      {
        mu.sigma.proposal = (1+sqrt(1/n))*0.03
        lambda.sigma.proposal = (1+sqrt(1/n))*0.18
        w.sigma.proposal = (1+sqrt(1/n))*0.07
      }
      else
      {
        mu.sigma.proposal = (1+sqrt(1/n))*0.08
        lambda.sigma.proposal = (1+sqrt(1/n))*0.32
        w.sigma.proposal = (1+sqrt(1/n))*0.085
      }
    }
    else if (n > 20 && n <= 60)
    {
      mu.sigma.proposal = (1+sqrt(1/n))*0.84
      lambda.sigma.proposal = (1+sqrt(1/n))*0.36
      w.sigma.proposal = (1+sqrt(1/n))*0.15
    }
    else
    {
      mu.sigma.proposal = (1+sqrt(1/n))*3.5
      lambda.sigma.proposal = (1+sqrt(1/n))*0.41
      w.sigma.proposal = (1+sqrt(1/n))*0.4
    }
    ## default AIS = 0
    for (i in 1:POP)
    {
      weight1[i] = weight[i] # previous unnormalized importance weight
      weight[i] = weight[i] + (gamma - gammad)*li[i] #???????????
      vec[i] = weight[i] - weight1[i]
    }
    sum = 0 # ???????????????????
    if (!flag.resample)
    {
      weight1 = selfweight(weight1)
      sum = sum(weight1*exp(vec))
    }
    else
    {
      sum = sum(exp(vec))
      sum = sum/POP # equal weight1
    }
    
    ## nc[1] = sum*nc[1]
    i = resampdec(weight, 0)
    if(i == 1)
    {
      resamplestrat() # TODO!!!!!!!!!!!!!!!!!!!!!!!!!!
      flag.resample = TRUE
    }
    else
    {
      flag.resample = FALSE
    }
    
    vec = numeric(KM)
    vec1 = numeric(KM)
    if (n < 3000)
      sw = 10
    else
      sw = 1
    for (i in 1:POP)
    {
      for (l in 1:sw)
      {
        ## MH on mu
        vec = mu[, i] #+ TODO!!!!!!!! # proposal
        acc.tmp.1 = sapply(1:KM, function(i) norma1(vec[i], mu.prior.xi, mu.prior.kappa))
        acc.tmp.2 = sapply(1:KM, function(i) norma1(mu[,i], mu.prior.xi, mu.prior.kappa))
        acc = sum(acc.tmp.1 - acc.tmp.2)
        
        sum.tmp = sapply(1:R, function(i) mixden(data[i], mu[,i], lambda[,i], w[,i], KM))
        sum = sum(sum.tmp)
        if (PIN == 0)
        {
          acc = acc + gammad*(sum-li[i])
        }
        else
        {
          acc = acc + gamma*(sum-li[i])
        }
        
        u = log(rnorm()) # TO CHECK!!!!!!!!!!!!!
        if (u < acc)
        {
          mu[, i] = vec
        }
        li[i] = sum
        
        ## MH on lambda
        ## #####################
        ## TO DO!!!!!!!!!!!!!!!!!!!
        ## 
        ## ######################
        
        ## MH on weights
        ## #####################
        ## TO DO!!!!!!!!!!!!!!!!!!!
        ## 
        ## ######################
        
      }
    }
  }
  ## default AIS = 0; skip
  ## n > 0 && AIS =1 1
  for (i in 1:POP)
  {
    if (n == 0 && i == 0)
    {
      map = logpost()
      mapmean = mu[,i]
      maplambda = lambda[,i]
      mapw = lambda[,i]
    }
    else
    {
      sum = logpost()
      if (sum > map)
      {
        map = sum
        mapmean = mu[,i]
        maplambda = lambda[,i]
        mapw = w[,i]
      }
    }
  }
  ## skip average unnoramlized log posterior
}

## ###################
## the log of normal density (NOT include constant)
## 
## 
## ###################
norma1 <- function(x, mu, precision)
{
    val = 0.5*log(precision) - 0.5*precision*(x - mu)*(x-mu)
    return(val)
}

## ###################
## mixdensity
##
## ###################
mixden <- function(x, mu, lambda, w, KM)
{
  oneden <- function(x, mu, lambda, w)
  {
    return(exp(norma1(x, mu, lambda))*w)
  }
  res = sapply(1:KM, function(i) oneden(x, mu[i], lambda[i], w[i]))
  return(res)
}


## ###################
## selfweight
## 
## ###################

selfweight <- function(weight)
{
  weight.max = max(weight)
  weight.tmp = exp(weight-weight.max)
  weight.norm = weight.tmp/sum(weight.tmp)
  return(weight.norm)
}

## ###################
## resampdec
## 
## if ess below N/2, then resample
## ###################

resampdec <- function(weight, ess)
{
  weight.max = max(weight)
  nc1 = sum(exp(weight-weight.max))
  ## normalizing constant for weights (log scale)
  nc1 = -weight.max - log(nc1)
  weight.tmp = weight + nc1
  
  weight.tmp.max = max(weight.tmp)
  acc = 0
  for (i in 1:POP)
  {
    acc = acc + exp(2*(weight.tmp - weight.tmp.max))
  }
  acc = -2*weight.tmp.max - log(acc)
  
  criterion = log(R1)
  if (acc < criterion)
    return(1)
  else
    return(0)
}

## #####################
## resamplestrat
## 
## #####################

resamplestrat <- function()
{
  
}

## #####################
## logpost
##
## #####################
logpost <- function()
{
  
}
