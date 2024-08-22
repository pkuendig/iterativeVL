################################################################################
# Prediction accuracy (after parameter estimation) and marginal log-likelihood
# for different numbers of nearest neighbors (NN) for the Vecchia-Laplace 
# approximation relative to Laplace.
# Considered models: Iterative-VL | GPVecchia (v0.1.7) | Laplace (reference)
# Sample size: n=n_p=5'000
################################################################################

library(gpboost)
library(GPvecchia)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/make_data.R")

n <- 5000
n_rep <- 100
sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
NN <- c(10,20,40,60)
init_cov_pars <- c(1.0000000, 0.1911457) #model_iterativeVL$get_optim_params()$init_cov_pars

#Results for Iterative-VL + GPVecchia
res_cols <- c("NN", "sigma2", "rho", "time_negLL", "negLL_at_true_covpars", "num_optim_iter", 
              "RMSE_latent_mu", "log_score", "number_of_errors", "number_of_neg_pred_var")
results_iterativeVL <- data.frame(matrix(nrow=length(NN)*n_rep, ncol = length(res_cols)))
colnames(results_iterativeVL) <- res_cols
results_GPVecchia <- results_iterativeVL
results_GPVecchia$number_of_errors <- 0
results_GPVecchia$number_of_neg_pred_var <- 0

#Results for Laplace
res_cols <- c("sigma2", "rho", "negLL_at_true_covpars", "num_optim_iter", 
              "RMSE_latent_mu", "log_score")
results_Laplace <- data.frame(matrix(nrow=n_rep, ncol = length(res_cols)))
colnames(results_Laplace) <- res_cols

i <- 1 #counter data set seed (for same seed as in Laplace model)
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

  ############################################################################
  # For different NN: Model Estimation + Prediction on Test Set
  ############################################################################
  for(j in 1:length(NN)){
    
    ##Iterative-VL############################################################
    
    ##Estimation
    model_iterativeVL <- GPModel(gp_coords = coords_train,
                                  cov_function = "matern",
                                  cov_fct_shape=1.5,
                                  likelihood="bernoulli_logit",
                                  matrix_inversion_method = "iterative",
                                  gp_approx = "vecchia",
                                  vecchia_ordering = "random",
                                  num_neighbors=NN[j])
    
    model_iterativeVL$set_optim_params(params = list(maxit=1000,
                                                       trace=TRUE,
                                                       cg_max_num_it=10000,
                                                       cg_max_num_it_tridiag=10000,
                                                       init_cov_pars = init_cov_pars,
                                                       cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                                       seed_rand_vec_trace = k))
    
    model_iterativeVL$fit(y=y_train)
    
    results_iterativeVL$NN[k] <- NN[j]
    results_iterativeVL$num_optim_iter[k] <- model_iterativeVL$get_num_optim_iter()
    results_iterativeVL$sigma2[k] <- model_iterativeVL$get_cov_pars()[1]
    results_iterativeVL$rho[k] <- model_iterativeVL$get_cov_pars()[2]
    
    ##Prediction
    model_iterativeVL$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                                           num_neighbors_pred=NN[j],
                                           nsim_var_pred = 2000)
    
    prediction_iterative_VL <- predict(model_iterativeVL,
                                        gp_coords_pred=as.matrix(coords_test),
                                        predict_var = TRUE,
                                        predict_cov_mat = FALSE,
                                        predict_response = FALSE)
    
    results_iterativeVL$RMSE_latent_mu[k] <- sqrt(mean((prediction_iterative_VL$mu-b_test)^2))
    results_iterativeVL$log_score[k] <- -sum(dnorm(b_test, prediction_iterative_VL$mu, sqrt(prediction_iterative_VL$var), log = TRUE))
    
    #Likelihood at true parameters
    results_iterativeVL$time_negLL[k] <- system.time(
      results_iterativeVL$negLL_at_true_covpars[k] <- model_iterativeVL$neg_log_likelihood(cov_pars=true_covpars, y=y_train))[3]
    
    ##GPVecchia##################################################################
    
    ##Estimation
    vecchia.approx <- vecchia_specify(coords_train, m=NN[j], cond.yz = "zy") #for posterior #RF ordering, cond.yz = "zy", ordering = "maxmin"
    vecchia.approx.IW <- vecchia_specify(coords_train, m=NN[j]) #for integrated likelihood #IW ordering, cond.yz = "SGV", ordering = "maxmin"
    
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
    
    opt <- optim(par = log(init_cov_pars), fn = GPVecchia_mll, method = "Nelder-Mead",
                 control = list("trace" = 4, "maxit" = 1000))
    
    final_covparms <- c(exp(opt$par[1]), exp(opt$par[2])/sqrt(2*1.5), 1.5)
    
    results_GPVecchia$num_optim_iter[k] <- opt$counts[1]
    results_GPVecchia$sigma2[k] <- exp(opt$par[1])
    results_GPVecchia$rho[k] <- exp(opt$par[2])
    
    ##Prediction
    post_zy <- calculate_posterior_VL(y_train, vecchia.approx, "logistic", final_covparms)
    vecchia.approx.train.test <- vecchia_specify(coords_train, m=NN[j], locs.pred=coords_test)
    preds_test <- vecchia_laplace_prediction(post_zy, vecchia.approx.train.test, final_covparms)
    
    results_GPVecchia$RMSE_latent_mu[k] <- sqrt(mean((preds_test$mu.pred-b_test)^2))
    
    if(any(preds_test$var.pred < 0)){
      results_GPVecchia$number_of_neg_pred_var[k] <- sum(preds_test$var.pred < 0)
    }else{
      results_GPVecchia$log_score[k] <- -sum(dnorm(b_test, preds_test$mu.pred, sqrt(preds_test$var.pred), log = TRUE))  
    }
    
    #Likelihood at true parameters
    results_GPVecchia$time_negLL[k] <- system.time(results_GPVecchia$negLL_at_true_covpars[k] <- GPVecchia_mll(log(true_covpars)))[3]
    k <- k + 1 #counter to store results
    gc()
  } 
  i <- i + 1
}

############################################################################
# Plotting
############################################################################
library(ggplot2)
library(grid)
library(tidyr)

#Comparison with Laplace - 1:length(NN)
results_iterativeVL$RRMSE_latent_mu <- results_iterativeVL$RMSE_latent_mu / rep(results_Laplace$RMSE_latent_mu, each=length(NN))
results_GPVecchia$RRMSE_latent_mu <- results_GPVecchia$RMSE_latent_mu / rep(results_Laplace$RMSE_latent_mu, each=length(NN))
results_iterativeVL$diff_log_score <- results_iterativeVL$log_score - rep(results_Laplace$log_score, each=length(NN))
results_GPVecchia$diff_log_score <- results_GPVecchia$log_score - rep(results_Laplace$log_score, each=length(NN))
results_iterativeVL$rel_negLL_at_true_covpars <- results_iterativeVL$negLL_at_true_covpars / rep(results_Laplace$negLL_at_true_covpars, each=length(NN))
results_GPVecchia$rel_negLL_at_true_covpars <- results_GPVecchia$negLL_at_true_covpars / rep(results_Laplace$negLL_at_true_covpars, each=length(NN))
results_iterativeVL$model <- "Iterative-VL"
results_GPVecchia$model <- "GPVecchia"

results_iterativeVL$model <- factor(results_iterativeVL$model, levels = c("Iterative-VL", "GPVecchia"), ordered=TRUE)
results <- rbind(results_iterativeVL, results_GPVecchia)

results_long <- gather(results, key="key", value="value", c("rel_negLL_at_true_covpars", "time_negLL", "RRMSE_latent_mu", "diff_log_score"))
results_long$key <- factor(results_long$key, levels = c("rel_negLL_at_true_covpars", "time_negLL", "RRMSE_latent_mu", "diff_log_score"), ordered=TRUE)

#labeling
labels <- list(
  'rel_negLL_at_true_covpars'="relative log-likelihood",
  'time_negLL'="Time (s) for log-likelihood",
  'RRMSE_latent_mu'="RRMSE",
  'diff_log_score'="dLS"
)
my_labeller <- function(variable,value){
  return(labels[value])
}

my_colors <- RColorBrewer::brewer.pal(6,"Set1")[3:4]
options(scipen = 999)

ggplot(results_long, aes(x=NN, y=value, color=model)) + 
  stat_summary(fun = mean, geom="line", linewidth=1) +
  stat_summary(fun = mean,
               geom = "errorbar",
               linewidth=1,
               width = 1.5,
               fun.max = function(x) mean(x) + 2*sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - 2*sd(x) / sqrt(length(x))) + 
  facet_wrap(~key, scales = "free_y", nrow = 2, labeller = my_labeller) + 
  labs(color = "") + scale_color_manual(values = my_colors) +
  theme_bw() + theme(legend.position = "top") + ylab("") + xlab("m")
