
library(gpboost)
library(RandomFields)
# install.packages("remotes")
# library(remotes)
# install_version("RandomFieldsUtils", "1.2.5")
# install_version("RandomFields", "3.3.14")

make_data <- function(n, 
                      likelihood = "bernoulli_logit",
                      sigma2=1, #marginal variance of GP
                      rho=0.1, #range parameter
                      cov_function="exponential",
                      matern_nu=2.5,
                      with_fixed_effects=FALSE,
                      beta,
                      gamma_shape = 1,
                      p) {
  
  #Simulate spatial Gaussian process
  coords <- matrix(runif(n*2),ncol=2)
  if (cov_function == "exponential"){
    RFmodel <- RMexp(var=sigma2, scale=rho)
  } else if(cov_function == "matern") {
    RFmodel <- RMmatern(matern_nu, notinvnu=TRUE, var=sigma2, scale=rho)
  } else {
    stop("cov_function not supported")
  }
  sim <- RFsimulate(RFmodel, x=coords)
  eps <- sim$variable1
  eps <- eps - mean(eps)
  if(with_fixed_effects){
    #   Simulate fixed effects
    X <- cbind(rep(1,n),matrix(runif(n*p)-0.5,ncol=p))
    f <- X %*% beta    
  } else{
    X <- NA
    f <- 0
  }
  b <- f+eps
  #Simulate response variable
  if (likelihood == "bernoulli_probit") {
    probs <- pnorm(b)
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "bernoulli_logit") {
    probs <- 1/(1+exp(-(b)))
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "poisson") {
    mu <- exp(b)
    y <- qpois(runif(n), lambda = mu)
  } else if (likelihood == "gamma") {
    y <- qgamma(runif(n), rate = gamma_shape * exp(-b), shape = gamma_shape) #E[y] = exp(b)
  } else if (likelihood == "gaussian") {
    mu <- b
    y <- rnorm(n,sd=0.05) + mu
  }
  list(y=y, coords=coords, b=b, X=X)
}

