library(fields)
library(gpboost)

################################################################################
# Iterative-VL: Estimation on training data set with fixed shape 
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
res_n5e4 <- readRDS("./data/IterativeVL.rds")
L_MYD05 <- readRDS("../../data_sets/MYD05/L_MYD05.rds")

#glm - for initial fixed effects
glm_train <- glm(L_MYD05$L_vapor_train ~ -1 + L_MYD05$L_X_train, family = Gamma(link = "log"))
betas_init <- glm_train$coefficients

(init_shape <- res_n5e4$estimates$alpha) #optimal from estimation on subsample

##Estimation
gp_model <- GPModel(gp_coords = L_MYD05$L_locations_train,
                    cov_function = "matern",
                    likelihood="gamma",
                    cov_fct_shape=1.5,
                    matrix_inversion_method = "iterative",
                    gp_approx = "vecchia",
                    vecchia_ordering = "random",
                    num_neighbors=20)

gp_model$set_optim_params(params = list(trace=TRUE,
                                        optimizer_cov="lbfgs",
                                        estimate_aux_pars = FALSE,
                                        init_coef=betas_init,
                                        init_aux_pars=init_shape,
                                        cg_preconditioner_type = "Sigma_inv_plus_BtWB"))

gp_model$fit(y=L_MYD05$L_vapor_train, X = L_MYD05$L_X_train)

#Estimates
estimates <- list()
(estimates$betas  <- gp_model$get_coef())
(estimates$cov_pars <- gp_model$get_cov_pars())
(estimates$alpha <- init_shape)

################################################################################
# Prediction
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
L_MYD05 <- readRDS("../../data_sets/MYD05/L_MYD05.rds")

L_X_train_betas <- c(L_MYD05$L_X_train %*% estimates$betas)
L_X_test_betas <- c(L_MYD05$L_X_test %*% estimates$betas)

gp_model <- GPModel(gp_coords = L_MYD05$L_locations_train,
                    cov_function = "matern",
                    likelihood="gamma",
                    cov_fct_shape=1.5,
                    matrix_inversion_method = "iterative",
                    gp_approx = "vecchia",
                    vecchia_ordering = "random",
                    num_neighbors = 20)

gp_model$set_optim_params(params = list(trace=TRUE,
                                        optimizer_cov="lbfgs",
                                        init_coef = estimates$betas,
                                        init_aux_pars = estimates$alpha,
                                        cg_preconditioner_type = "Sigma_inv_plus_BtWB"))

gp_model$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                             num_neighbors_pred = 20,
                             nsim_var_pred = 1000)

pred_latent_gp_model <- predict(gp_model,
                                y = L_MYD05$L_vapor_train,
                                offset = L_X_train_betas,
                                offset_pred = L_X_test_betas,
                                gp_coords_pred = L_MYD05$L_locations_test,
                                cov_pars = estimates$cov_pars,
                                predict_var = TRUE,
                                predict_cov_mat = FALSE,
                                predict_response = FALSE)

#Response mean
pred_response_mu <- exp(pred_latent_gp_model$mu+1/2*pred_latent_gp_model$var)

#Response variance
pred_response_var <- exp(2*pred_latent_gp_model$mu+2*pred_latent_gp_model$var)/estimates$alpha +
  (exp(pred_latent_gp_model$var)-1)*exp(2*pred_latent_gp_model$mu+pred_latent_gp_model$var)

#RMSE: Response Prediction
RMSE <- sqrt(mean((pred_response_mu-L_MYD05$L_vapor_test)^2))

#CRPS: Latent Prediction
pred_latent_mu <- pred_latent_gp_model$mu #fixed effects included
pred_latent_sd <- sqrt(pred_latent_gp_model$var)

n_samples <- 100 #number of samples
sample_mat <- matrix(0, ncol=n_samples, nrow = length(L_MYD05$L_vapor_test))

set.seed(1)
for(s in 1:n_samples){
  sample_latent_mu = rnorm(n = length(L_MYD05$L_vapor_test), mean = pred_latent_mu, sd = pred_latent_sd)
  sample_response_mean = exp(sample_latent_mu)
  sample_mat[,s]= rgamma(length(L_MYD05$L_vapor_test), shape = estimates$alpha, rate = estimates$alpha/sample_response_mean)
}

CRPS <- mean(scoringRules::crps_sample(y = L_MYD05$L_vapor_test, dat = sample_mat))
  
saveRDS(list(estimates=estimates,
             pred_response_mu=pred_response_mu,
             pred_response_var=pred_response_var,
             RMSE=RMSE,
             CRPS=CRPS), "./data/IterativeVL_25e4.rds")
