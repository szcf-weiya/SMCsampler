####################
## paras:
##

KM = 4 # number of components in mixture
POP = 1000 # number of particles
R = 100 # number of data points


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

## ##################
## an SMC algorithm with MAX densities
## 
## ##################

beta = numeric(POP)
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
      
    }
  }
}