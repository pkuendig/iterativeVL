library(fields)
library(gpboost)

################################################################################
# LowRank (FITC): Estimation on subsample
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
M_MYD05 <- readRDS("../../data/MYD05/M_MYD05.rds")

#glm - for initial fixed effects
glm_train <- glm(M_MYD05$M_vapor_train ~ -1 + M_MYD05$M_X_train, family = Gamma(link = "log"))
betas_init <- glm_train$coefficients

##Estimation
gp_model <- GPModel(gp_coords = M_MYD05$M_locations_train,
                    cov_function = "matern",
                    likelihood="gamma",
                    cov_fct_shape=1.5,
                    matrix_inversion_method = "cholesky",
                    gp_approx = "fitc",
                    num_ind_points=500)

gp_model$set_optim_params(params = list(trace=TRUE,
                                        optimizer_cov="lbfgs",
                                        estimate_aux_pars = TRUE,
                                        init_coef=betas_init))

gp_model$fit(y=M_MYD05$M_vapor_train, X = M_MYD05$M_X_train)

#Estimates
estimates <- list()
(estimates$betas  <- gp_model$get_coef())
(estimates$cov_pars <- gp_model$get_cov_pars())
(estimates$alpha <- gp_model$get_aux_pars())

################################################################################
# Prediction
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
L_MYD05 <- readRDS("../../data/MYD05/L_MYD05.rds")

L_X_train_betas <- c(L_MYD05$L_X_train %*% estimates$betas)
L_X_test_betas <- c(L_MYD05$L_X_test %*% estimates$betas)

gp_model <- GPModel(gp_coords = L_MYD05$L_locations_train,
                    cov_function = "matern",
                    likelihood="gamma",
                    cov_fct_shape=1.5,
                    matrix_inversion_method = "cholesky",
                    gp_approx = "fitc",
                    num_ind_points=500)

gp_model$set_optim_params(params = list(trace=TRUE,
                                        optimizer_cov="lbfgs",
                                        init_coef = estimates$betas,
                                        init_aux_pars = estimates$alpha))

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
pred_latent_mu <- pred_latent_gp_model$mu #fixed effects included.
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

################################################################################
saveRDS(list(estimates=estimates,
             pred_response_mu=pred_response_mu,
             pred_response_var=pred_response_var,
             RMSE=RMSE,
             CRPS=CRPS), "./FITC.rds")
