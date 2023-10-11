require(fields)
require(GPvecchia) #V0.1.4

################################################################################
# Origin: https://github.com/katzfuss-group/GPvecchia-Laplace/blob/master/VL_scripts/Modis_pipeline_VL.R
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
M_MYD05 <- readRDS("../data_sets/MYD05/M_MYD05.rds")

####trend estimation############################################################

trend_locations_train <- M_MYD05$M_locations_train
trend_y <- M_MYD05$M_vapor_train
trend_X_train <- M_MYD05$M_X_train

t_start = Sys.time()
beta = c(1, 0.001, 0.001)
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
t_end = Sys.time()
time_dur = as.double(difftime(t_end, t_start, units = "mins")); print(time_dur)
print(beta)
# [,1]
# [1,]  0.0815648857
# [2,] -0.0000698454
# [3,] -0.0007610925

#use beta from run
beta = matrix(c(0.0815648857, -0.0000698454, -0.0007610925))
XB = M_MYD05$M_X_train %*% beta

####Shape estimation############################################################
#shape estimation: maximize conditional likelihood. Integrated likelihood diverges.

update_a = function(a_init, covparms, vecchia.approx, vecchia.approx.IW, XB){
  a_prev =a_init
  for(i in 1:10){
    t_start = Sys.time()
    posterior = calculate_posterior_VL(M_MYD05$M_vapor_train,
                                       vecchia.approx,
                                       likelihood_model="gamma",
                                       covparms=c(covparms[1], covparms[2]/sqrt(2*smoothness), smoothness),
                                       likparms = list("alpha"=a_prev),
                                       prior_mean = XB)
    mu = posterior$mean + XB
    llh = function(a) -sum(-a*exp(-mu)*M_MYD05$M_vapor_train + (a-1)*log(M_MYD05$M_vapor_train)+a*log(a)-a*mu-log(gamma(a))) # concave in a and XB
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
    vll = vecchia_laplace_likelihood(M_MYD05$M_vapor_train,
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

print("Step 1, generating vecchia approximations")
vecchia.approx = vecchia_specify(M_MYD05$M_locations_train, m=m_rep, cond.yz = "zy")
vecchia.approx.IW = vecchia_specify(M_MYD05$M_locations_train, m=m_rep)

#Iterative method: estimate a, then covparms, then a again
print("Step 2, optimizing parameters")

#Start:
a_prev = 30
covparms_prev = c(.3, 15) #sigma range

t_start = Sys.time()
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
t_end = Sys.time()
time_dur = as.double(difftime(t_end, t_start, units = "mins"))
print(time_dur)

cat("alpha: ", a, "\n")
#alpha: 0.9500872

cat("sigma, range: ", covparms, "\n")
#sigma, range: 0.2557801 86.46232 

cat("Number of iterations: ", iter_count, "\n")
#Number of iterations:  3

################################################################################
# Prediction
################################################################################
require(fields)
require(GPvecchia)

NN <- 20
smoothness = 1.5

estimates_GPVecchia <- list()
estimates_GPVecchia$alpha <- 0.9500872
estimates_GPVecchia$cov_pars <- c(0.2557801, 86.46232/sqrt(2*smoothness), smoothness) #sigma, range, smoothness
estimates_GPVecchia$beta <- matrix(c(0.0815648857, -0.0000698454, -0.0007610925))

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
L_MYD05 <- readRDS("../data_sets/MYD05/L_MYD05.rds")

L_X_train_betas <- c(L_MYD05$L_X_train %*% estimates_GPVecchia$beta)
L_X_test_betas <- c(L_MYD05$L_X_test %*% estimates_GPVecchia$beta)

vecchia.approx <- vecchia_specify(L_MYD05$L_locations_train, m=NN, cond.yz = "zy") #for posterior

default_lh_params <- list("alpha"= estimates_GPVecchia$alpha, "sigma"=sqrt(.1)) #alpha = shape

system.time(post_zy <- calculate_posterior_VL(L_MYD05$L_vapor_train, vecchia.approx, "gamma", estimates_GPVecchia$cov_pars, 
                                  likparms = default_lh_params, prior_mean = L_X_train_betas))[3]
system.time(vecchia.approx.train.test <- vecchia_specify(L_MYD05$L_locations_train, m=NN, locs.pred=L_MYD05$L_locations_test))[3]
system.time(preds_test <- vecchia_laplace_prediction(post_zy, vecchia.approx.train.test, estimates_GPVecchia$cov_pars, pred.mean = L_X_test_betas))[3]

#RMSE: Response Prediction
pred_response_mu <- exp(preds_test$mu.pred+1/2*preds_test$var.pred)
sqrt(mean((pred_response_mu-L_MYD05$L_vapor_test)^2))
#RMSE:0.1205773

#CRPS: Probabilistic Prediction
pred_latent_mu <- preds_test$mu.pred
pred_latent_sd <- sqrt(pmax(preds_test$var.pred, 0)) #as in Modis_pipeline_VL.R

n_samples <- 100 #number of samples
sample_mat = matrix(0, ncol=n_samples, nrow = length(L_MYD05$L_vapor_test))

set.seed(1)
for(s in 1:n_samples){
  sample_latent_mu = rnorm(n = length(L_MYD05$L_vapor_test), mean = pred_latent_mu, sd = pred_latent_sd)
  sample_response_mean = exp(sample_latent_mu)
  sample_mat[,s]= rgamma(length(L_MYD05$L_vapor_test), shape = estimates_GPVecchia$alpha, rate = estimates_GPVecchia$alpha/sample_response_mean)
}

mean(scoringRules::crps_sample(y = L_MYD05$L_vapor_test, dat = sample_mat))
#CRPS:0.1411717

#Response mean
#quilt.plot(L_MYD05$L_locations_test, pred_response_mu, nx=150, ny=300, main="GPVecchia", add.legend = FALSE)
saveRDS(pred_response_mu, "./data/pred_response_mean_GPVecchia.rds")

#Response variance
pred_response_var <- exp(2*preds_test$mu.pred+2*preds_test$var.pred)/estimates_GPVecchia$alpha +
  (exp(preds_test$var.pred)-1)*exp(2*preds_test$mu.pred+preds_test$var.pred)
#quilt.plot(L_MYD05$L_locations_test, pred_response_var, nx=150, ny=300, main="GPVecchia", add.legend = FALSE)
saveRDS(pred_response_var, "./data/pred_response_var_GPVecchia.rds")
