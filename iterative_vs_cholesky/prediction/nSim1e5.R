################################################################################
# Simulation-based prediction
################################################################################
library(gpboost)

################################################################################
#Generate data

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../../data_sets/simulated/make_data.R")

sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
n <- 1e5
n_rep <- 10
NN <- 20

set.seed(1)
mydata <- make_data(n=2*n,
                    likelihood = "bernoulli_logit",
                    sigma2=sigma2_true,
                    rho=rho_true,
                    cov_function="matern",
                    matern_nu=1.5,
                    with_fixed_effects=FALSE)

################################################################################
#Subsample train and test set

set.seed(123)
i_train <- sample((2*n), n, replace = FALSE)
i_test <- 1:(2*n)
i_test <- i_test[!(i_test %in% i_train)]

coords_train <- mydata$coords[i_train,]
y_train <- mydata$y[i_train]
b_train <- mydata$b[i_train]

coords_test <- mydata$coords[i_test,]
y_test <- mydata$y[i_test]
b_test <- mydata$b[i_test]

################################################################################
#Import Cholesky-based prediction

chol_results <- readRDS("./data/cholesky1e5.rds")

################################################################################
#Simulation-based prediction
nSim <- c(50, 100, 200, 500, 1000, 2000, 3000, 4000)

results <- data.frame(matrix(nrow=2*n_rep*length(nSim), ncol = 5))
colnames(results) <- c("nSim", "log_score", "RMSE", "time", "preconditioner")

pred_var <- data.frame(matrix(nrow=2*n_rep*length(nSim), ncol = n + 2))
colnames(pred_var) <- c("preconditioner", "nSim", paste0("var_", 1:n))

i <- 1
for(s in 1:length(nSim)){
  for(r in 1:n_rep){
    
    #Sigma_inv_plus_BtWB
    model_Sigma_inv_plus_BtWB <- GPModel(gp_coords = coords_train,
                                         cov_function = "matern",
                                         cov_fct_shape=1.5,
                                         likelihood="bernoulli_logit",
                                         gp_approx = "vecchia",
                                         num_neighbors=NN,
                                         vecchia_ordering="random",
                                         matrix_inversion_method = "iterative")
  
    model_Sigma_inv_plus_BtWB$set_optim_params(params = list(cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                                             seed_rand_vec_trace = i))
    
    model_Sigma_inv_plus_BtWB$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                                                  num_neighbors_pred=NN,
                                                  nsim_var_pred = nSim[s])
  
    time <- system.time(pred_Sigma_inv_plus_BtWB <- predict(model_Sigma_inv_plus_BtWB, 
                                                            y=y_train,
                                                            gp_coords_pred=coords_test,
                                                            cov_pars = true_covpars,
                                                            predict_var = TRUE, 
                                                            predict_cov_mat = FALSE,
                                                            predict_response = FALSE))[3]
    
    results$preconditioner[i] <- "P[1]"
    results$time[i] <- time
    results$nSim[i] <- nSim[s]
    results$RMSE[i] <- sqrt(mean((pred_Sigma_inv_plus_BtWB$var-chol_results$pred_var)^2))
    results$log_score[i] <- -sum(dnorm(b_test, pred_Sigma_inv_plus_BtWB$mu, sqrt(pred_Sigma_inv_plus_BtWB$var), log = TRUE))
    
    pred_var$preconditioner[i] <- "P[1]"
    pred_var$nSim[i] <- nSim[s]
    pred_var[i,3:(n+2)] <- pred_Sigma_inv_plus_BtWB$var
    i <- i + 1
    
    #piv_chol_on_Sigma
    model_piv_chol_on_Sigma <- GPModel(gp_coords = coords_train,
                                       cov_function = "matern",
                                       cov_fct_shape=1.5,
                                       likelihood="bernoulli_logit",
                                       gp_approx = "vecchia",
                                       num_neighbors=NN,
                                       vecchia_ordering="random",
                                       matrix_inversion_method = "iterative")
    
    
    model_piv_chol_on_Sigma$set_optim_params(params = list(cg_preconditioner_type = "piv_chol_on_Sigma",
                                                           seed_rand_vec_trace = i))
    
    model_piv_chol_on_Sigma$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                                                num_neighbors_pred=NN,
                                                nsim_var_pred = nSim[s])
  
    time <- system.time(pred_piv_chol_on_Sigma <- predict(model_piv_chol_on_Sigma, 
                                                          y=y_train,
                                                          gp_coords_pred=coords_test,
                                                          cov_pars = true_covpars,
                                                          predict_var = TRUE, 
                                                          predict_cov_mat = FALSE,
                                                          predict_response = FALSE))[3]
    
    results$preconditioner[i] <- "P[2]"
    results$time[i] <- time
    results$nSim[i] <- nSim[s]
    results$RMSE[i] <- sqrt(mean((pred_piv_chol_on_Sigma$var-chol_results$pred_var)^2))
    results$log_score[i] <- -sum(dnorm(b_test, pred_piv_chol_on_Sigma$mu, sqrt(pred_piv_chol_on_Sigma$var), log = TRUE))
    
    pred_var$preconditioner[i] <- "P[2]"
    pred_var$nSim[i] <- nSim[s]
    pred_var[i,3:(n+2)] <- pred_piv_chol_on_Sigma$var
    
    i <- i + 1
    
    print(i)
    saveRDS(list(results=results,
                 pred_var=pred_var), "./data/nSim1e5.rds")
    gc()
  }
}
