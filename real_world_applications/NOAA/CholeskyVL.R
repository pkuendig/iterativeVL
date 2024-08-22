library(fields)
library(gpboost)

################################################################################
# Cholesky-VL: Estimation
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
NOAA <- readRDS("../../data/NOAA/NOAA.rds")

#glm - for initial fixed effects
glm_train <- glm(NOAA$Y_train ~ -1 + NOAA$X_train, family = Gamma(link = "log"))
betas_init <- glm_train$coefficients

##Estimation
gp_model <- GPModel(gp_coords = NOAA$coords_train,
                    cov_function = "matern",
                    likelihood="gamma",
                    cov_fct_shape=1.5,
                    matrix_inversion_method = "cholesky",
                    gp_approx = "vecchia",
                    vecchia_ordering = "random",
                    num_neighbors=20)

gp_model$set_optim_params(params = list(trace=TRUE,
                                        optimizer_cov="lbfgs",
                                        estimate_aux_pars = TRUE,
                                        init_coef=betas_init))

gp_model$fit(y=NOAA$Y_train, X = NOAA$X_train)

#Store estimates
estimates <- list()
(estimates$betas  <- gp_model$get_coef())
(estimates$cov_pars <- gp_model$get_cov_pars())
(estimates$alpha <- gp_model$get_aux_pars())

################################################################################
# Prediction
################################################################################

gp_model$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                             num_neighbors_pred = 20)

pred_latent_gp_model <- predict(gp_model,
                                X_pred = NOAA$X_test,
                                gp_coords_pred = NOAA$coords_test,
                                predict_var = TRUE,
                                predict_cov_mat = FALSE,
                                predict_response = FALSE)

#Response mean
pred_response_mu <- exp(pred_latent_gp_model$mu+1/2*pred_latent_gp_model$var)

#Response variance
pred_response_var <- exp(2*pred_latent_gp_model$mu+2*pred_latent_gp_model$var)/estimates$alpha +
  (exp(pred_latent_gp_model$var)-1)*exp(2*pred_latent_gp_model$mu+pred_latent_gp_model$var)

#RMSE: Response Prediction
RMSE <- sqrt(mean((pred_response_mu - NOAA$Y_test)^2))

#CRPS: Latent Prediction
pred_latent_mu <- pred_latent_gp_model$mu #fixed effects included.
pred_latent_sd <- sqrt(pred_latent_gp_model$var)

n_samples <- 100 #number of samples
sample_mat <- matrix(0, ncol=n_samples, nrow = length(NOAA$Y_test))

set.seed(1)
for(s in 1:n_samples){
  sample_latent_mu = rnorm(n = length(NOAA$Y_test), mean = pred_latent_mu, sd = pred_latent_sd)
  sample_response_mean = exp(sample_latent_mu)
  sample_mat[,s]= rgamma(length(NOAA$Y_test), shape = estimates$alpha, rate = estimates$alpha/sample_response_mean)
}

CRPS <- mean(scoringRules::crps_sample(y = NOAA$Y_test, dat = sample_mat))

################################################################################
saveRDS(list(estimates=estimates,
             pred_response_mu=pred_response_mu,
             pred_response_var=pred_response_var,
             RMSE=RMSE,
             CRPS=CRPS), "./CholeskyVL.rds")