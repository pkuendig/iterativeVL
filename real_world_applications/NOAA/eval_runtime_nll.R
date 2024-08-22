library(gpboost)
library(GPvecchia)

################################################################################
# Runtime for likelihood evaluation on training data set
################################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
NOAA <- readRDS("../../data/NOAA/NOAA.rds")

set.seed(2)

sigma2 <- 0.4
rho <- 0.86
alpha <- 5
NN <- 20 
cov_pars <- c(sigma2, rho) 

res_cols <- c("model", "time_negLL", "mll")
runtime_results <- data.frame(matrix(nrow=4, ncol = length(res_cols)))
colnames(runtime_results) <- res_cols

############################################################################  
# Sample Data
############################################################################  
coords_train <- NOAA$coords_train
y_train <- NOAA$Y_train
X_train <- NOAA$X_train

##############################################################################
# glm for fixed effects
##############################################################################
glm_train <- glm(y_train ~ -1 + X_train, family = Gamma(link = "log"))
betas <- glm_train$coefficients
fixed_effects <- as.numeric(X_train %*% betas)

############################################################################
# Calculate Likelihood
############################################################################
i <- 1

##Cholesky-VL############################################################
model_cholVL <- GPModel(gp_coords = coords_train,
                    cov_function = "matern",
                    likelihood="gamma",
                    cov_fct_shape=1.5,
                    matrix_inversion_method = "cholesky",
                    gp_approx = "vecchia",
                    vecchia_ordering = "random",
                    num_neighbors = 20)

#Likelihood
runtime_results$time_negLL[i] <- system.time(mll <- model_cholVL$neg_log_likelihood(cov_pars=cov_pars, y=y_train, aux_pars=c(alpha),
                                                                                    fixed_effects = fixed_effects))[3]

runtime_results$model[i] <- "CholeskyVL"
runtime_results$mll[i] <- mll
i <- i + 1

##Iterative-VL################################################################
model_itVL <- GPModel(gp_coords = coords_train,
                    cov_function = "matern",
                    likelihood="gamma",
                    cov_fct_shape=1.5,
                    matrix_inversion_method = "iterative",
                    gp_approx = "vecchia",
                    vecchia_ordering = "random",
                    num_neighbors=20)


model_itVL$set_optim_params(params = list(cg_preconditioner_type = "Sigma_inv_plus_BtWB"))

#Likelihood
runtime_results$time_negLL[i] <- system.time(mll <- model_itVL$neg_log_likelihood(cov_pars=cov_pars, y=y_train, aux_pars=c(alpha),
                                                                                  fixed_effects = fixed_effects))[3]

runtime_results$model[i] <- "IterativeVL"
runtime_results$mll[i] <- mll
i <- i + 1

##FITC########################################################################
model_fitc <- GPModel(gp_coords = coords_train,
                      cov_function = "matern",
                      likelihood="gamma",
                      cov_fct_shape=1.5,
                      matrix_inversion_method = "cholesky",
                      gp_approx = "fitc",
                      num_ind_points=500)


#Likelihood
runtime_results$time_negLL[i] <- system.time(mll <- model_fitc$neg_log_likelihood(cov_pars=cov_pars, y=y_train, aux_pars=c(alpha),
                                                                                  fixed_effects = fixed_effects))[3]

runtime_results$model[i] <- "FITC"
runtime_results$mll[i] <- mll
i <- i + 1

##GPVeccia##################################################################
vecchia.approx <- vecchia_specify(coords_train, m=NN, cond.yz = "zy") #for posterior
vecchia.approx.IW <- vecchia_specify(coords_train, m=NN) #for integrated likelihood
covparms <- c(cov_pars[1], cov_pars[2]/sqrt(2*1.5), 1.5) #sigma^2, range, smoothness -> range is divided by sqrt(2*smoothness) such that Matern covariance is identical with GPModel
default_lh_params = list("alpha"=alpha, "sigma"=sqrt(.1))

#Likelihood
runtime_results$time_negLL[i] <- system.time(mll <- vecchia_laplace_likelihood(y_train,
                                                                               vecchia.approx,
                                                                               likelihood_model="gamma",
                                                                               covmodel = "matern",
                                                                               covparms = covparms,
                                                                               return_all = FALSE,
                                                                               likparms = default_lh_params,
                                                                               vecchia.approx.IW=vecchia.approx.IW,
                                                                               prior_mean = fixed_effects))[3]

runtime_results$model[i] <- "GPVecchia"
runtime_results$mll[i] <- mll

#############################################################################
saveRDS(runtime_results, "./runtime.rds")