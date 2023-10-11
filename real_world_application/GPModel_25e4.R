require(fields)
require(gpboost)

################################################################################
# Estimation on large data set
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
L_MYD05 <- readRDS("../data_sets/MYD05/L_MYD05.rds")

#glm - for initial fixed effects
glm_train <- glm(L_MYD05$L_vapor_train ~ -1 + L_MYD05$L_X_train, family = Gamma(link = "log"))
betas_init = glm_train$coefficients

init_shape = 27.98258 #optimal from GPModel estimation on small data set

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
                                        estimate_aux_pars = FALSE,
                                        init_coef=betas_init,
                                        init_aux_pars=init_shape,
                                        cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                        cg_max_num_it=100000,
                                        cg_max_num_it_tridiag=100000))

system.time(gp_model$fit(y=L_MYD05$L_vapor_train, X = L_MYD05$L_X_train))[3]

#Estimates
estimates <- list()
gp_model$get_coef()
# Intercept     Longitude      Latitude 
# -2.550756e-02 -4.657731e-05 -8.489642e-04

gp_model$get_cov_pars()
# GP_var   GP_range 
# 0.2871984 11.9766092 

gp_model$get_current_neg_log_likelihood()
#-242545.3

################################################################################
# Prediction
################################################################################
require(fields)
require(gpboost)

estimates <- list()
estimates$betas <- c(-2.550756e-02, -4.657731e-05, -8.489642e-04)
estimates$cov_pars <- c(0.2871984, 11.9766092)
estimates$alpha <- c(27.98258)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
L_MYD05 <- readRDS("../data_sets/MYD05/L_MYD05.rds")

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
                                        init_coef = estimates$betas,
                                        init_aux_pars = estimates$alpha,
                                        cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                        cg_max_num_it=100000,
                                        cg_max_num_it_tridiag=100000))

gp_model$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                             num_neighbors_pred = 20,
                             nsim_var_pred = 2000)

#RMSE: Response Prediction
system.time(pred_response_gp_model <- predict(gp_model,
                                              y = L_MYD05$L_vapor_train,
                                              fixed_effects = L_X_train_betas,
                                              fixed_effects_pred = L_X_test_betas,
                                              gp_coords_pred = L_MYD05$L_locations_test,
                                              cov_pars = estimates$cov_pars,
                                              predict_var = FALSE,
                                              predict_cov_mat = FALSE,
                                              predict_response = TRUE))[3]

sqrt(mean((pred_response_gp_model$mu-L_MYD05$L_vapor_test)^2)) 
#RMSE: 0.07910425

#CRPS: Probabilistic Prediction
system.time(pred_latent_gp_model <- predict(gp_model,
                                            y = L_MYD05$L_vapor_train,
                                            fixed_effects = L_X_train_betas,
                                            fixed_effects_pred = L_X_test_betas,
                                            gp_coords_pred = L_MYD05$L_locations_test,
                                            cov_pars = estimates$cov_pars,
                                            predict_var = TRUE,
                                            predict_cov_mat = FALSE,
                                            predict_response = FALSE))[3]

pred_latent_mu <- pred_latent_gp_model$mu
pred_latent_sd <- sqrt(pred_latent_gp_model$var)

n_samples <- 100 #number of samples
sample_mat = matrix(0, ncol=n_samples, nrow = length(L_MYD05$L_vapor_test))

set.seed(1)
for(s in 1:n_samples){
  sample_latent_mu = rnorm(n = length(L_MYD05$L_vapor_test), mean = pred_latent_mu, sd = pred_latent_sd)
  sample_response_mean = exp(sample_latent_mu)
  sample_mat[,s]= rgamma(length(L_MYD05$L_vapor_test), shape = estimates$alpha, rate = estimates$alpha/sample_response_mean)
}

mean(scoringRules::crps_sample(y = L_MYD05$L_vapor_test, dat = sample_mat))
#CRPS: 0.04415603

#Plotting
#quilt.plot(L_MYD05$L_locations_test, pred_response_gp_model$mu, nx=150, ny=300, main="Iterative-VL trained on large data set", add.legend = FALSE)

saveRDS(pred_response_gp_model$mu, "./data/pred_response_GPModel_25e4.rds")
