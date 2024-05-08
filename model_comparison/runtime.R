################################################################################
# Runtime for likelihood
# Compare models: Laplace | Cholesky-VL | Iterative-VL | GPVeccia
################################################################################

library(gpboost)
library(GPvecchia)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data_sets/simulated/make_data.R")

N <- c(1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000)
sigma2_true <- 1
rho_true <- 1/20
NN <- 20 
cov_pars <- c(sigma2_true, rho_true) #Evaluation at the optimum

res_cols <- c("n", "model", "time_negLL")
runtime_results <- data.frame(matrix(nrow=length(N)*4, ncol = length(res_cols)))
colnames(runtime_results) <- res_cols

i <- 1
for(q in 1:length(N)){
  
  ############################################################################  
  # Generate Data
  ############################################################################  
  
  n <- N[q]
  print(n)
  set.seed(i)
  mydata <- make_data(n=n,
                      likelihood = "bernoulli_logit",
                      sigma2=sigma2_true,
                      rho=rho_true,
                      cov_function="matern",
                      matern_nu=1.5,
                      with_fixed_effects=FALSE)
  
  coords_train <- mydata$coords
  y_train <- mydata$y
  b_train <- mydata$b

  ############################################################################
  # Calculate Likelihood + Gradient 
  ############################################################################
  
  ##Laplace###################################################################
  if(n <= 10000) {
    model_Laplace <- GPModel(gp_coords = coords_train,
                             cov_function = "matern",
                             cov_fct_shape=1.5,
                             likelihood="bernoulli_logit",
                             matrix_inversion_method = "cholesky",
                             gp_approx="none")

    model_Laplace$set_optim_params(params = list(maxit=1,
                                            trace=TRUE,
                                            lr_cov=1e-100, #small -> no LR-Update
                                            init_cov_pars = cov_pars))

    #Likelihood
    runtime_results$time_negLL[i] <- system.time(model_Laplace$neg_log_likelihood(cov_pars=cov_pars, y=y_train))[3]
  }
  
  runtime_results$n[i] <- n
  runtime_results$model[i] <- "Laplace"
  i <- i + 1
  
  ##Cholesky-VL############################################################
  if(n <= 50000) {

    model_VL <- GPModel(gp_coords = coords_train,
                        cov_function = "matern",
                        cov_fct_shape=1.5,
                        likelihood="bernoulli_logit",
                        matrix_inversion_method = "cholesky",
                        gp_approx = "vecchia",
                        vecchia_ordering = "random",
                        num_neighbors=NN)

    model_VL$set_optim_params(params = list(maxit=1,
                                           trace=TRUE,
                                           lr_cov=1e-100, #small -> no LR-Update
                                           init_cov_pars = cov_pars))

    #Likelihood
    runtime_results$time_negLL[i] <- system.time(model_VL$neg_log_likelihood(cov_pars=cov_pars, y=y_train))[3]
  }
  
  runtime_results$n[i] <- n
  runtime_results$model[i] <- "Cholesky-VL"
  i <- i + 1
  
  ##Iterative-VL################################################################
  # cg_delta_conv = 1e-2
  # num_rand_vec_trace = 50
  # cg_preconditioner_type = "Sigma_inv_plus_BtWB"
  
  model_VL_CG <- GPModel(gp_coords = coords_train,
                         cov_function = "matern",
                         cov_fct_shape=1.5,
                         likelihood="bernoulli_logit",
                         matrix_inversion_method = "iterative",
                         gp_approx = "vecchia",
                         vecchia_ordering = "random",
                         num_neighbors=NN)
  
  model_VL_CG$set_optim_params(params = list(maxit=1,
                                             trace=TRUE,
                                             lr_cov=1e-100, #small -> no LR-Update
                                             cg_max_num_it=100000,
                                             cg_max_num_it_tridiag=100000,
                                             cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                             seed_rand_vec_trace = i))
  
  #Likelihood
  runtime_results$time_negLL[i] <- system.time(model_VL_CG$neg_log_likelihood(cov_pars=cov_pars, y=y_train))[3]
  
  runtime_results$n[i] <- n
  runtime_results$model[i] <- "Iterative-VL"
  i <- i + 1
  
  ##GPVeccia##################################################################
  if(n <= 100000) {
    vecchia.approx <- vecchia_specify(coords_train, m=NN, cond.yz = "zy") #for posterior
    vecchia.approx.IW <- vecchia_specify(coords_train, m=NN) #for integrated likelihood
    covparms <- c(cov_pars[1], cov_pars[2]/sqrt(2*1.5), 1.5) #sigma^2, range, smoothness -> range is divided by sqrt(2*smoothness) such that Matern covariance is identical with GPModel

    #Likelihood
    runtime_results$time_negLL[i] <- system.time(mll <- vecchia_laplace_likelihood(y_train,
                                                                                   vecchia.approx,
                                                                                   likelihood_model="logistic",
                                                                                   covmodel = "matern",
                                                                                   covparms = covparms,
                                                                                   return_all = FALSE,
                                                                                   vecchia.approx.IW=vecchia.approx.IW))[3]
  }
  
  runtime_results$n[i] <- n
  runtime_results$model[i] <- "GPVecchia"
  i <- i + 1
  #############################################################################
  saveRDS(runtime_results, "./data/runtime_results.rds")
  gc()
}
