####################
## paras:
##

KM = 4 # number of components in mixture
POP = 1000 # number of particles
R = 100 # number of data points
MAX = 100 # number of distns ??????


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
        wei[i] = 0
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
        mu.proposal = (1+sqrt(1/n))*0.03
        lambda.proposal = (1+sqrt(1/n))*0.18
        w.proposal = (1+sqrt(1/n))*0.07
      }
      else
      {
        mu.proposal = (1+sqrt(1/n))*0.08
        lambda.proposal = (1+sqrt(1/n))*0.32
        w.proposal = (1+sqrt(1/n))*0.085
      }
    }
    else if (n > 20 && n <= 60)
    {
      mu.proposal = (1+sqrt(1/n))*0.84
      lambda.proposal = (1+sqrt(1/n))*0.36
      w.proposal = (1+sqrt(1/n))*0.15
    }
    else
    {
      mu.proposal = (1+sqrt(1/n))*3.5
      lambda.proposal = (1+sqrt(1/n))*0.41
      w.proposal = (1+sqrt(1/n))*0.4
    }
      
  }
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

