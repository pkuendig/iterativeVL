################################################################################
# Prediction accuracy for different NN for Iterative-VL + GPVecchia
################################################################################

library(gpboost)
library(GPvecchia)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data_sets/simulated/make_data.R")

n <- 50000
n_rep <- 50
sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
NN <- c(10,20,40,60)

res_cols <- c("NN","time_negLL", "negLL_at_true_covpars",
              "RMSE_latent_mu", "log_score", "number_of_errors", "number_of_neg_pred_var")
results_iterative_VL <- data.frame(matrix(nrow=length(NN)*n_rep, ncol = length(res_cols)))
colnames(results_iterative_VL) <- res_cols
results_GPVecchia <- results_iterative_VL
results_GPVecchia$number_of_errors <- 0

i <- 1 #counter data set seed
k <- 1 #counter to store results
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
  # For different NN: Prediction on Test Set
  ############################################################################
  for(j in 1:length(NN)){
    
    ##Iterative-VL#########################################################
    # cg_delta_conv = 1e-2
    # num_rand_vec_trace = 50
    # cg_preconditioner_type = "Sigma_inv_plus_BtWB"
    # nsim_var_pred = 2000 
    
    model_iterative_VL <- GPModel(gp_coords = coords_train,
                                  cov_function = "matern",
                                  cov_fct_shape=1.5,
                                  likelihood="bernoulli_logit",
                                  matrix_inversion_method = "iterative",
                                  gp_approx = "vecchia",
                                  vecchia_ordering = "random",
                                  num_neighbors=NN[j])
    
    
    model_iterative_VL$set_optim_params(params = list(trace=TRUE,
                                                      cg_max_num_it=100000,
                                                      cg_max_num_it_tridiag=100000,
                                                      cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                                      seed_rand_vec_trace = k))
    results_iterative_VL$NN[k] <- NN[j]
    
    ##Prediction
    model_iterative_VL$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                                           num_neighbors_pred=NN[j],
                                           nsim_var_pred = 2000)
    
    prediction_iterative_VL <- predict(model_iterative_VL,
                                       y=y_train,
                                       cov_pars=true_covpars,
                                       gp_coords_pred=as.matrix(coords_test),
                                       predict_var = TRUE,
                                       predict_cov_mat = FALSE,
                                       predict_response = FALSE)
    
    results_iterative_VL$RMSE_latent_mu[k] <- sqrt(mean((prediction_iterative_VL$mu-b_test)^2))
    results_iterative_VL$log_score[k] <- -sum(dnorm(b_test, prediction_iterative_VL$mu, sqrt(prediction_iterative_VL$var), log = TRUE))
    
    #Likelihood at true parameters
    results_iterative_VL$time_negLL[k] <- system.time(
      results_iterative_VL$negLL_at_true_covpars[k] <- model_iterative_VL$neg_log_likelihood(cov_pars=true_covpars, y=y_train))[3]
    
    ##GPVecchia##################################################################
    
    vecchia.approx <- vecchia_specify(coords_train, m=NN[j], cond.yz = "zy") #for posterior
    vecchia.approx.IW <- vecchia_specify(coords_train, m=NN[j]) #for integrated likelihood
    
    GPVecchia_mll <- function(cov_pars){
      covparms <- c(exp(cov_pars[1]), exp(cov_pars[2])/sqrt(2*1.5), 1.5) #sigma^2, range, smoothness -> range is divided by sqrt(2*smoothness) such that Matern covariance is identical with GPModel
      
      n_mll <- tryCatch({
        -vecchia_laplace_likelihood(y_train,
                                    vecchia.approx,
                                    likelihood_model="logistic",
                                    covmodel = "matern",
                                    covparms = covparms,
                                    return_all = FALSE,
                                    vecchia.approx.IW=vecchia.approx.IW)
      }, error = function(err) {
        results_GPVecchia$number_of_errors[k] <<- results_GPVecchia$number_of_errors[k] + 1 
        return(Inf)
      })
      
      return(n_mll)
    }
    
    results_GPVecchia$NN[k] <- NN[j]
    
    ##Prediction
    post_zy <- calculate_posterior_VL(y_train, vecchia.approx, "logistic", c(true_covpars[1], true_covpars[2]/sqrt(2*1.5), 1.5)) 
    vecchia.approx.train.test <- vecchia_specify(coords_train, m=NN[j], locs.pred=coords_test)
    preds_test <- vecchia_laplace_prediction(post_zy, vecchia.approx.train.test, c(true_covpars[1], true_covpars[2]/sqrt(2*1.5), 1.5))
    
    results_GPVecchia$RMSE_latent_mu[k] <- sqrt(mean((preds_test$mu.pred-b_test)^2))
    
    if(any(preds_test$var.pred < 0)){
      results_GPVecchia$number_of_neg_pred_var[k] <- sum(preds_test$var.pred < 0)
    }else{
      results_GPVecchia$log_score[k] <- -sum(dnorm(b_test, preds_test$mu.pred, sqrt(preds_test$var.pred), log = TRUE))  
    }
    
    #Likelihood at true parameters
    results_GPVecchia$time_negLL[k] <- system.time(results_GPVecchia$negLL_at_true_covpars[k] <- GPVecchia_mll(log(true_covpars)))[3]
    
    ############################################################################
    k <- k + 1
    saveRDS(list(results_iterative_VL = results_iterative_VL,
                 results_GPVecchia = results_GPVecchia), "./data/NN_prediction_n5e4.rds")
    gc()
  }  
  #############################################################################
  i <- i + 1
}
