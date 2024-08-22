################################################################################
# Estimation and prediction on 100 simulated data sets with n=n_p=2e4
# Considered models: Cholesky-VL | Iterative-VL | GPVeccia (v0.1.7)
################################################################################

library(gpboost)
library(GPvecchia)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/make_data.R")

n <- 20000
n_rep <- 100
sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
NN <- 20 
init_cov_pars <- c(1.0000000, 0.1920107) #model_iterative_VL$get_optim_params()$init_cov_pars

res_cols <- c("time_estimation", "sigma2", "rho", "negLL_at_true_covpars", "num_optim_iter", 
              "RMSE_latent_mu", "log_score", "number_of_errors")
results_VL <- data.frame(matrix(nrow=n_rep, ncol = length(res_cols)))
colnames(results_VL) <- res_cols
results_GPVecchia <- results_iterative_VL <- results_VL
results_GPVecchia$number_of_errors <- 0
results_GPVecchia$number_of_neg_pred_var <- 0

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
  
  ##Cholesky-VL############################################################

  ##Estimation
  model_VL <- GPModel(gp_coords = coords_train,
                      cov_function = "matern",
                      cov_fct_shape=1.5,
                      likelihood="bernoulli_logit",
                      matrix_inversion_method = "cholesky",
                      gp_approx = "vecchia",
                      vecchia_ordering = "random",
                      num_neighbors=NN)
  
  model_VL$set_optim_params(params = list(maxit=1000,
                                           trace=TRUE,
                                           init_cov_pars = init_cov_pars))
  
  results_VL$time_estimation[i] <- system.time(model_VL$fit(y=y_train))[3]
  results_VL$num_optim_iter[i] <- model_VL$get_num_optim_iter()
  results_VL$sigma2[i] <- model_VL$get_cov_pars()[1]
  results_VL$rho[i] <- model_VL$get_cov_pars()[2]
  
  ##Prediction
  model_VL$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                               num_neighbors_pred=NN)
  
  prediction_VL <- predict(model_VL,
                           gp_coords_pred=as.matrix(coords_test),
                           predict_var = TRUE,
                           predict_cov_mat = FALSE,
                           predict_response = FALSE)
  
  results_VL$RMSE_latent_mu[i] <- sqrt(mean((prediction_VL$mu-b_test)^2))
  results_VL$log_score[i] <- -sum(dnorm(b_test, prediction_VL$mu, sqrt(prediction_VL$var), log = TRUE))
  
  #Likelihood at true parameters
  results_VL$negLL_at_true_covpars[i] <- model_VL$neg_log_likelihood(cov_pars=true_covpars, y=y_train)
  
  ##Iterative Veccia-Laplace####################################################
  
  ##Estimation
  model_iterative_VL <- GPModel(gp_coords = coords_train,
                                cov_function = "matern",
                                cov_fct_shape=1.5,
                                likelihood="bernoulli_logit",
                                matrix_inversion_method = "iterative",
                                gp_approx = "vecchia",
                                vecchia_ordering = "random",
                                num_neighbors=NN)
  
  model_iterative_VL$set_optim_params(params = list(maxit=1000,
                                               trace=TRUE,
                                               init_cov_pars = init_cov_pars,
                                               cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                               seed_rand_vec_trace = i))
    
  results_iterative_VL$time_estimation[i] <- system.time(model_iterative_VL$fit(y=y_train))[3]
  results_iterative_VL$num_optim_iter[i] <- model_iterative_VL$get_num_optim_iter()
  results_iterative_VL$sigma2[i] <- model_iterative_VL$get_cov_pars()[1]
  results_iterative_VL$rho[i] <- model_iterative_VL$get_cov_pars()[2]
  
  ##Prediction
  model_iterative_VL$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                                         num_neighbors_pred=NN,
                                         nsim_var_pred = 2000)
  
  prediction_iterative_VL <- predict(model_iterative_VL,
                              gp_coords_pred=as.matrix(coords_test),
                              predict_var = TRUE,
                              predict_cov_mat = FALSE,
                              predict_response = FALSE)
  
  results_iterative_VL$RMSE_latent_mu[i] <- sqrt(mean((prediction_iterative_VL$mu-b_test)^2))
  results_iterative_VL$log_score[i] <- -sum(dnorm(b_test, prediction_iterative_VL$mu, sqrt(prediction_iterative_VL$var), log = TRUE))
  
  #Likelihood at true parameters
  results_iterative_VL$negLL_at_true_covpars[i] <- model_iterative_VL$neg_log_likelihood(cov_pars=true_covpars, y=y_train)
  
  ##GPVecchia##################################################################
  
  ##Estimation
  vecchia.approx <- vecchia_specify(coords_train, m=NN, cond.yz = "zy") #for posterior
  vecchia.approx.IW <- vecchia_specify(coords_train, m=NN) #for integrated likelihood
  
  GPVecchia_mll <- function(cov_pars){
    #sigma^2, range, smoothness 
    #range is divided by sqrt(2*smoothness) such that Matern covariance is identical with GPModel
    covparms <- c(exp(cov_pars[1]), exp(cov_pars[2])/sqrt(2*1.5), 1.5)
    
    n_mll <- tryCatch({
      -vecchia_laplace_likelihood(y_train,
                                   vecchia.approx,
                                   likelihood_model="logistic",
                                   covmodel = "matern",
                                   covparms = covparms,
                                   return_all = FALSE,
                                   vecchia.approx.IW=vecchia.approx.IW)
    }, error = function(err) {
      results_GPVecchia$number_of_errors[i] <<- results_GPVecchia$number_of_errors[i] + 1 
      return(Inf)
    })
    
    cat("Negative marginal LL:",  n_mll, "for covparms:", covparms[1], ",", covparms[2], "\n")
    return(n_mll)
  }
  
  results_GPVecchia$time_estimation[i] <- system.time({
    opt <- optim(par = log(init_cov_pars), fn = GPVecchia_mll, method = "Nelder-Mead",
               control = list("trace" = 4, "maxit" = 1000))})[3]
  
  final_covparms <- c(exp(opt$par[1]), exp(opt$par[2])/sqrt(2*1.5), 1.5)
  
  results_GPVecchia$num_optim_iter[i] <- opt$counts[1]
  results_GPVecchia$sigma2[i] <- exp(opt$par[1])
  results_GPVecchia$rho[i] <- exp(opt$par[2])
  
  ##Prediction
  post_zy <- calculate_posterior_VL(y_train, vecchia.approx, "logistic", final_covparms)
  vecchia.approx.train.test <- vecchia_specify(coords_train, m=NN, locs.pred=coords_test)
  preds_test <- vecchia_laplace_prediction(post_zy, vecchia.approx.train.test, final_covparms)
  
  results_GPVecchia$RMSE_latent_mu[i] <- sqrt(mean((preds_test$mu.pred-b_test)^2))
  
  if(any(preds_test$var.pred < 0)){
    results_GPVecchia$number_of_neg_pred_var[i] <- sum(preds_test$var.pred < 0)
  }else{
    results_GPVecchia$log_score[i] <- -sum(dnorm(b_test, preds_test$mu.pred, sqrt(preds_test$var.pred), log = TRUE))  
  }
  
  #Likelihood at true parameters
  results_GPVecchia$negLL_at_true_covpars[i] <- GPVecchia_mll(log(true_covpars))
  
  #############################################################################
  i <- i + 1
  gc()
}

################################################################################
# Plotting
################################################################################
library(ggplot2)
library(grid)

results_VL$model <- "Cholesky-VL"
results_iterative_VL$model <- "Iterative-VL"
results_GPVecchia$model <- "GPVecchia"

results <- rbind(results_VL,
                 results_iterative_VL,
                 results_GPVecchia[,-c(9)])
results$model <- factor(results$model, levels = c("Cholesky-VL", "Iterative-VL", "GPVecchia"), ordered=TRUE)

##Log-likelihood relative to Cholesky-VL########################################
results$relative_L[results$model=="Cholesky-VL"] <- NA
results$relative_L[results$model=="Iterative-VL"] <- 
  results$negLL_at_true_covpars[results$model=="Iterative-VL"] / results$negLL_at_true_covpars[results$model=="Cholesky-VL"]
results$relative_L[results$model=="GPVecchia"] <- 
  results$negLL_at_true_covpars[results$model=="GPVecchia"] / results$negLL_at_true_covpars[results$model=="Cholesky-VL"]

ggplot(results[results$model!="Cholesky-VL", ], aes(x=model, y=relative_L)) + 
  geom_hline(yintercept=1, linetype="dashed") + geom_boxplot() + 
  ylab("relative log-likelihood") + xlab("") + theme_bw() + 
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE)

##Parameter estimates###########################################################
p_var <- ggplot(results, aes(x=model, y=sigma2)) + 
  geom_hline(yintercept=sigma2_true, linetype="dashed") + geom_boxplot() +
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE) + 
  ylab(expression(sigma[1]^2))+ xlab("") + theme_bw()

p_range <- ggplot(results, aes(x=model, y=rho)) + 
  geom_hline(yintercept=rho_true, linetype="dashed") + geom_boxplot() +
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE) +
  xlab("") + ylab(expression(rho)) + theme_bw()

grid.newpage()
grid.draw(cbind(ggplotGrob(p_var), ggplotGrob(p_range), size = "last"))

##Prediction accuracy###########################################################
#RMSE
ggplot(results, aes(x=model, y=RMSE_latent_mu)) + geom_boxplot() +
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE) + 
  ylab("RMSE") + theme_bw() + xlab("")

#log-score
ggplot(results[results$model!="GPVecchia", ], aes(x=model, y=log_score)) + geom_boxplot() +
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE) + 
  xlab("") + ylab("LS") + theme_bw()

################################################################################
# Results for tables
################################################################################

##Parameter estimates###########################################################

##RMSE sigma^2
cat("GPVecchia: ", round(sqrt(mean((results_GPVecchia$sigma2 - sigma2_true)^2)), digits = 4))
cat("Cholesky-VL: ", round(sqrt(mean((results_VL$sigma2 - sigma2_true)^2)), digits = 4))
cat("Iterative-VL: ", round(sqrt(mean((results_iterative_VL$sigma2 - sigma2_true)^2)), digits = 4))

##RMSE rho
cat("GPVecchia: ", round(sqrt(mean((results_GPVecchia$rho - rho_true)^2)), digits = 4))
cat("Cholesky-VL: ", round(sqrt(mean((results_VL$rho - rho_true)^2)), digits = 4))
cat("Iterative-VL: ", round(sqrt(mean((results_iterative_VL$rho - rho_true)^2)), digits = 4))

##Bias sigma^2
cat("GPVecchia: ", round(mean(results_GPVecchia$sigma2) - sigma2_true, digits = 4))
cat("Cholesky-VL: ", round(mean(results_VL$sigma2) - sigma2_true, digits = 4))
cat("Iterative-VL: ", round(mean(results_iterative_VL$sigma2) - sigma2_true, digits = 4))

##Bias rho
cat("GPVecchia: ", round(mean(results_GPVecchia$rho) - rho_true, digits = 4))
cat("Cholesky-VL: ", round(mean(results_VL$rho) - rho_true, digits = 4))
cat("Iterative-VL: ", round(mean(results_iterative_VL$rho) - rho_true, digits = 4))

##Prediction accuracy###########################################################

##Mean RMSE
cat("GPVecchia: ", round(mean(results_GPVecchia$RMSE_latent_mu), digits = 4))
cat("Cholesky-VL: ", round(mean(results_VL$RMSE_latent_mu), digits = 4))
cat("Iterative-VL: ", round(mean(results_iterative_VL$RMSE_latent_mu), digits = 4))

##Mean log-score
#cat("GPVecchia: ", round(mean(results_GPVecchia$log_score, na.rm=TRUE), digits = 4))
cat("Cholesky-VL: ", round(mean(results_VL$log_score), digits = 4))
cat("Iterative-VL: ", round(mean(results_iterative_VL$log_score), digits = 4))

##SD mean RMSE
cat("GPVeccia: ", round(sqrt(mean((results_GPVecchia$RMSE_latent_mu - mean(results_GPVecchia$RMSE_latent_mu))^2))/sqrt(length(results_GPVecchia$RMSE_latent_mu)), digits=4))
cat("Cholesky-VL: ", round(sqrt(mean((results_VL$RMSE_latent_mu - mean(results_VL$RMSE_latent_mu))^2))/sqrt(length(results_VL$RMSE_latent_mu)), digits=4))
cat("Iterative-VL: ", round(sqrt(mean((results_iterative_VL$RMSE_latent_mu - mean(results_iterative_VL$RMSE_latent_mu))^2))/sqrt(length(results_iterative_VL$RMSE_latent_mu)), digits=4))

##SD mean log-score
cat("Cholesky-VL: ", round(sqrt(mean((results_VL$log_score - mean(results_VL$log_score))^2))/sqrt(length(results_VL$log_score)), digits=4))
cat("Iterative-VL: ", round(sqrt(mean((results_iterative_VL$log_score - mean(results_iterative_VL$log_score))^2))/sqrt(length(results_iterative_VL$log_score)), digits=4))

##Runtime estimation############################################################

cat("GPVecchia: ", round(mean(results_GPVecchia$time_estimation), digits = 1))
cat("Cholesky-VL: ", round(mean(results_VL$time_estimation), digits = 1))
cat("Iterative-VL: ", round(mean(results_iterative_VL$time_estimation), digits = 1))
