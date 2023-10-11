require(fields)
require(gpboost)

################################################################################
# Estimation
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
M_MYD05 <- readRDS("../data_sets/MYD05/M_MYD05.rds")

#glm - for initial fixed effects
glm_train <- glm(M_MYD05$M_vapor_train ~ -1 + M_MYD05$M_X_train, family = Gamma(link = "log"))
betas_init = glm_train$coefficients

##Estimation
gp_model <- GPModel(gp_coords = M_MYD05$M_locations_train,
                    cov_function = "matern",
                    likelihood="gamma",
                    cov_fct_shape=1.5,
                    matrix_inversion_method = "iterative",
                    gp_approx = "vecchia",
                    vecchia_ordering = "random",
                    num_neighbors=20)

gp_model$set_optim_params(params = list(maxit=1000,
                                        trace=TRUE,
                                        cg_max_num_it=10000,
                                        cg_max_num_it_tridiag=10000,
                                        estimate_aux_pars = TRUE,
                                        init_coef=betas_init,
                                        cg_preconditioner_type = "Sigma_inv_plus_BtWB"))

system.time(gp_model$fit(y=M_MYD05$M_vapor_train, X = M_MYD05$M_X_train))[3] #4607.5 s

#Estimates
estimates <- list()
gp_model$get_coef()
# Intercept     Longitude      Latitude 
# -2.091357e-02 -5.101541e-05 -8.490558e-04 

gp_model$get_cov_pars()
# GP_var   GP_range 
# 0.2967264 15.3465348 

gp_model$get_aux_pars()
# shape 
# 27.98258 

gp_model$get_current_neg_log_likelihood()
# -33510.74

################################################################################
# Prediction
################################################################################
require(fields)
require(gpboost)

estimates <- list()
estimates$betas <- c(-2.091357e-02, -5.101541e-05, -8.490558e-04 )
estimates$cov_pars <- c(0.2967264, 15.3465348)
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
                                        cg_max_num_it=100000,
                                        cg_max_num_it_tridiag=100000, 
                                        init_coef = estimates$betas,
                                        init_aux_pars = estimates$alpha,
                                        cg_preconditioner_type = "Sigma_inv_plus_BtWB"))

gp_model$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                             num_neighbors_pred = 20,
                             nsim_var_pred = 2000)

system.time(pred_latent_gp_model <- predict(gp_model,
                                            y = L_MYD05$L_vapor_train,
                                            fixed_effects = L_X_train_betas,
                                            fixed_effects_pred = L_X_test_betas,
                                            gp_coords_pred = L_MYD05$L_locations_test,
                                            cov_pars = estimates$cov_pars,
                                            predict_var = TRUE,
                                            predict_cov_mat = FALSE,
                                            predict_response = FALSE))[3]

#Response mean
pred_response_mu <- exp(pred_latent_gp_model$mu+1/2*pred_latent_gp_model$var)
quilt.plot(L_MYD05$L_locations_test, pred_response_mu, nx=100, ny=200, main="Iterative-VL", add.legend = FALSE)
saveRDS(pred_response_mu, "./data/pred_response_mean_GPModel.rds")

#Response variance
pred_response_var <- exp(2*pred_latent_gp_model$mu+2*pred_latent_gp_model$var)/estimates$alpha +
  (exp(pred_latent_gp_model$var)-1)*exp(2*pred_latent_gp_model$mu+pred_latent_gp_model$var)
quilt.plot(L_MYD05$L_locations_test, pred_response_var, nx=100, ny=200, main="Iterative-VL", add.legend = FALSE)
saveRDS(pred_response_var, "./data/pred_response_var_GPModel.rds")

#RMSE: Response Prediction
sqrt(mean((pred_response_mu-L_MYD05$L_vapor_test)^2)) 
#RMSE: 0.08121711

#CRPS: Probabilistic Prediction
pred_latent_mu <- pred_latent_gp_model$mu #fixed effects included.
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
#CRPS:0.04449582
