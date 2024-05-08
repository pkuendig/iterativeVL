################################################################################
# Cholesky-based prediction
################################################################################

library(gpboost)

################################################################################
#Generate data

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../../../../data_sets/simulated/make_data.R")

sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
n <- 1e5
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
#Cholesky

model_chol <- GPModel(gp_coords = coords_train,
                      cov_function = "matern",
                      cov_fct_shape = 1.5,
                      likelihood = "bernoulli_logit",
                      gp_approx = "vecchia",
                      num_neighbors=NN,
                      vecchia_ordering="random",
                      matrix_inversion_method = "cholesky")

model_chol$set_prediction_data(vecchia_pred_type = "latent_order_obs_first_cond_obs_only",
                               num_neighbors_pred=NN)

time_chol <- system.time(pred_chol <- predict(model_chol, 
                                              y=y_train,
                                              gp_coords_pred=coords_test,
                                              cov_pars = true_covpars,
                                              predict_var = TRUE, 
                                              predict_cov_mat = FALSE,
                                              predict_response = FALSE))[3]

log_score_chol <- -sum(dnorm(b_test, pred_chol$mu, sqrt(pred_chol$var), log = TRUE))

saveRDS(list(pred_mu=pred_chol$mu,
             pred_var=pred_chol$var,
             time_chol=time_chol,
             log_score_chol=log_score_chol), "./data/cholesky1e5.rds")
