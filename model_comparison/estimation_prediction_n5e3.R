################################################################################
# Estimation and prediction with Laplace model
# Used for comparison with results from NN_estimation_prediction_n5e3.R
################################################################################

library(gpboost)
library(GPvecchia)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data_sets/simulated/make_data.R")

n <- 5000
n_rep <- 100
sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
NN <- 20 
init_cov_pars <- c(1.0000000, 0.1911457) #model_Laplace$get_optim_params()$init_cov_pars

res_cols <- c("sigma2", "rho", "negLL_at_true_covpars", "num_optim_iter", 
              "RMSE_latent_mu", "log_score", "number_of_errors")
results_Laplace <- data.frame(matrix(nrow=n_rep, ncol = length(res_cols)))
colnames(results_Laplace) <- res_cols
results_GPVecchia <- results_iterative_VL <- results_VL <- results_Laplace
results_GPVecchia$number_of_errors <- 0

i <- 1
for(r in 1:n_rep){
  
  ############################################################################  
  # Generate Data
  ############################################################################  
  
  set.seed(i)
  mydata <- make_data(n=2*n,
                      likelihood = "bernoulli_logit",
                      sigma2=sigma2_true,
                      rho=rho_true,
                      cov_function="matern",
                      matern_nu=1.5,
                      with_fixed_effects=FALSE)
  
  i_train <- sample((2*n), n, replace = FALSE)
  i_test <- 1:(2*n)
  i_test <- i_test[!(i_test %in% i_train)]
  
  coords_train <- mydata$coords[i_train,]
  y_train <- mydata$y[i_train]
  b_train <- mydata$b[i_train]
  
  coords_test <- mydata$coords[i_test,]
  y_test <- mydata$y[i_test]
  b_test <- mydata$b[i_test]

  ############################################################################
  # Model Estimation + Prediction on Test Set
  ############################################################################
  
  ##Laplace###################################################################
    
  ##Estimation
  model_Laplace <- GPModel(gp_coords = coords_train,
                           cov_function = "matern",
                           cov_fct_shape=1.5,
                           likelihood="bernoulli_logit",
                           matrix_inversion_method = "cholesky",
                           gp_approx="none")
  
  model_Laplace$set_optim_params(params = list(maxit=1000,
                                          trace=TRUE,
                                          init_cov_pars = init_cov_pars))
  
  model_Laplace$fit(y=y_train)
  
  results_Laplace$num_optim_iter[i] <- model_Laplace$get_num_optim_iter()
  results_Laplace$sigma2[i] <- model_Laplace$get_cov_pars()[1]
  results_Laplace$rho[i] <- model_Laplace$get_cov_pars()[2]
  
  ##Prediction
  prediction_Laplace <- predict(model_Laplace,
                                gp_coords_pred=as.matrix(coords_test),
                                predict_var = TRUE,
                                predict_cov_mat = FALSE,
                                predict_response = FALSE)
  
  results_Laplace$RMSE_latent_mu[i] <- sqrt(mean((prediction_Laplace$mu-b_test)^2))
  results_Laplace$log_score[i] <- -sum(dnorm(b_test, prediction_Laplace$mu, sqrt(prediction_Laplace$var), log = TRUE))
  
  #Likelihood at true parameters
  results_Laplace$negLL_at_true_covpars[i] <- model_Laplace$neg_log_likelihood(cov_pars=true_covpars, y=y_train)
  
  #############################################################################
  i <- i + 1
  
  saveRDS(list(results_Laplace = results_Laplace,
               results_VL = results_VL,
               results_iterative_VL = results_iterative_VL,
               results_GPVecchia = results_GPVecchia), "./data/estimation_predicition_n5e3.rds")
  gc()
}
