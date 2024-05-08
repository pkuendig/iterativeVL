require(fields)
require(GPvecchia) #V0.1.7

set.seed(1)
################################################################################
# GPVecchia: Estimation
# Procedure from:
# https://github.com/katzfuss-group/GPvecchia-Laplace/blob/master/VL_scripts/Modis_pipeline_VL.R
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
NOAA <- readRDS("../../data_sets/NOAA/NOAA.rds")

####trend estimation############################################################
trend_locations_train <- NOAA$coords_train
trend_y <- NOAA$Y_train
trend_X_train <- NOAA$X_train

beta <- c(1, 0.001, 0.001)
beta_prev <- c(0,0,0)
for(i in 1:100){
  print(i)
  XB = trend_X_train%*% beta
  W = -diag.spam(array(exp(-XB)*trend_y))
  A  = exp(-XB)*trend_y-1
  U  = W%*% XB - A
  beta = solve(t(trend_X_train) %*% W %*% trend_X_train , t(trend_X_train) %*% U)
  if(all(abs(beta-beta_prev)<1e-6)){
    cat("Stopp after ", i, "iterations.")
    break
  }
  beta_prev <- beta
}

print(beta)

#use beta from run
beta <- matrix(beta)
XB <- NOAA$X_train %*% beta

####Shape estimation############################################################
#shape estimation: maximize conditional likelihood

update_a = function(a_init, covparms, vecchia.approx, vecchia.approx.IW, XB){
  a_prev =a_init
  for(i in 1:10){
    posterior = calculate_posterior_VL(NOAA$Y_train,
                                       vecchia.approx,
                                       likelihood_model="gamma",
                                       covparms=c(covparms[1], covparms[2]/sqrt(2*smoothness), smoothness),
                                       likparms = list("alpha"=a_prev),
                                       prior_mean = XB)
    mu = posterior$mean + XB
    llh = function(a) -sum(-a*exp(-mu)*NOAA$Y_train + (a-1)*log(NOAA$Y_train)+a*log(a)-a*mu-log(gamma(a))) # concave in a and XB
    param_est = optim(a_prev, llh, method = "Brent", lower = .01, upper = 1e2)
    a = param_est$par; print(a)
    if(abs(a-a_prev) < 1e-5) {print("convergence criteria met (fitting shape parameter)"); break}
    a_prev = a
  }
  return(a)
}

####Covar param estimation: integrated likelihood###############################

fit_covparms = function(a, covparms_init, vecchia.approx, vecchia.approx.IW, XB){
  vl_likelihood = function(x0){
    theta = exp(x0)
    covparms=c(theta[1], theta[2]/sqrt(2*smoothness), smoothness) # sigma range smoothness
    sprintf("Evaluating covparms = (%.4f %.4f %.4f)", covparms[1], covparms[2],covparms[3])
    default_lh_params = list("alpha"=a, "sigma"=sqrt(.1))
    # Perform inference on latent mean with Vecchia Laplace approximation
    vll = vecchia_laplace_likelihood(NOAA$Y_train,
                                     vecchia.approx,
                                     likelihood_model="gamma",
                                     covparms = covparms,
                                     return_all = FALSE,
                                     likparms = default_lh_params,
                                     prior_mean = XB,
                                     vecchia.approx.IW=vecchia.approx.IW,
                                     y_init = NA )
    sprintf("Likelihood for covparms = (%.4f %.4f %.4f): %.4f",
            covparms[1], covparms[2],covparms[3], vll)
    return(-vll)
  }
  x0 = log(covparms_init)
  vl_likelihood(x0)
  res = optim(x0, vl_likelihood, method = "Nelder-Mead", control = list("trace" = 4, "maxit" = 500, "reltol" = 1e-5))
  print(res$convergence); print(exp(res$par))
  return(exp(res$par))
}

####Do parameter estimation#####################################################

m_rep = 20 #NN fix
smoothness = 1.5 #smoothness, fix

vecchia.approx = vecchia_specify(NOAA$coords_train, m=m_rep, cond.yz = "zy")
vecchia.approx.IW = vecchia_specify(NOAA$coords_train, m=m_rep)

#Iterative method: estimate a, then covparms, then a again
a_prev = 1.8
covparms_prev = c(1, 6.6) # sigma range

iter_count = 1
for(i in 1:5){
  a = update_a(a_prev, covparms_prev, vecchia.approx, vecchia.approx.IW, XB)
  covparms = fit_covparms(a, covparms_prev, vecchia.approx, vecchia.approx.IW, XB)
  print(paste("Found shape parameter a=", a))
  print(paste("Found covariance parameters (sig, range) = ",covparms))
  
  if(abs(a-a_prev)<1e-3 &
     abs(covparms_prev[1] - covparms[1]) < 1e-2 &
     abs(covparms_prev[2] - covparms[2])/covparms_prev[2] < 1e-2){
    print("Convergence criteria met (fitting all parameters)")
    iter_count = i
    break
  }
  a_prev = a
  covparms_prev = covparms
}

cat("alpha: ", a, "\n")
cat("sigma, range: ", covparms, "\n")
cat("Number of iterations: ", iter_count, "\n")

################################################################################
# Prediction
################################################################################

NN <- 20
smoothness = 1.5

estimates <- list()
estimates$alpha <- a
estimates$cov_pars <- c(covparms[1], covparms[2]/sqrt(2*smoothness), smoothness) #sigma, range, smoothness
estimates$betas <- beta

X_testB <- NOAA$X_test %*% beta

vecchia.approx <- vecchia_specify(NOAA$coords_train, m=NN, cond.yz = "zy")
default_lh_params <- list("alpha"= estimates$alpha, "sigma"=sqrt(.1)) #alpha = shape
post_zy <- calculate_posterior_VL(NOAA$Y_train, vecchia.approx, "gamma", estimates$cov_pars, 
                                  likparms = default_lh_params, prior_mean = XB)
vecchia.approx.train.test <- vecchia_specify(NOAA$coords_train, m=NN, locs.pred=NOAA$coords_test)
preds_test <- vecchia_laplace_prediction(post_zy, vecchia.approx.train.test, estimates$cov_pars, pred.mean = X_testB)

#RMSE: Response Prediction
pred_response_mu <- exp(preds_test$mu.pred+1/2*preds_test$var.pred)
RMSE <- sqrt(mean((pred_response_mu-NOAA$Y_test)^2))

#CRPS: Latent Prediction
pred_latent_mu <- preds_test$mu.pred
pred_latent_sd <- sqrt(pmax(preds_test$var.pred, 0)) #as in Modis_pipeline_VL.R

n_samples <- 100 #number of samples
sample_mat = matrix(0, ncol=n_samples, nrow = length(NOAA$Y_test))

set.seed(1)
for(s in 1:n_samples){
  sample_latent_mu = rnorm(n = length(NOAA$Y_test), mean = pred_latent_mu, sd = pred_latent_sd)
  sample_response_mean = exp(sample_latent_mu)
  sample_mat[,s]= rgamma(length(NOAA$Y_test), shape = estimates$alpha, rate = estimates$alpha/sample_response_mean)
}
CRPS <- mean(scoringRules::crps_sample(y = NOAA$Y_test, dat = sample_mat))

#Response variance
pred_response_var <- exp(2*preds_test$mu.pred+2*preds_test$var.pred)/estimates$alpha +
  (exp(preds_test$var.pred)-1)*exp(2*preds_test$mu.pred+preds_test$var.pred)

saveRDS(list(estimates=estimates,
             pred_response_mu=pred_response_mu,
             pred_response_var=pred_response_var,
             RMSE=RMSE,
             CRPS=CRPS), "./data/GPVecchia.rds")
