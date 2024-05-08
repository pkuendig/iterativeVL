library(gpboost)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../../data_sets/simulated/make_data.R")

N <- c(500,1000,2000,5000,10000,20000,50000)
n_rep <- 100
sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
NN <- 20 

res_cols <- c("n", "sigma2", "rho", "num_optim_iter")
results_Sigma_inv_plus_BtWB <- data.frame(matrix(nrow=length(N)*n_rep, ncol = length(res_cols)))
colnames(results_Sigma_inv_plus_BtWB) <- res_cols

i <- 1
for(n in 1:length(N)){
  for(r in 1:n_rep){
    
    ############################################################################  
    # Generate Data
    ############################################################################  
    
    set.seed(i)
    mydata <- make_data(n=N[n],
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
    # Model Estimation + Prediction on Test Set
    ############################################################################
    
    ##Estimation
    model_Sigma_inv_plus_BtWB <- GPModel(gp_coords = coords_train,
                                         cov_function = "matern",
                                         cov_fct_shape=1.5,
                                         likelihood="bernoulli_logit",
                                         matrix_inversion_method = "iterative",
                                         gp_approx = "vecchia",
                                         vecchia_ordering = "random",
                                         num_neighbors=NN)
    
    model_Sigma_inv_plus_BtWB$set_optim_params(params = list(maxit=1000,
                                                            trace=TRUE,
                                                            cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                                            seed_rand_vec_trace = i))
    
    model_Sigma_inv_plus_BtWB$fit(y=y_train)
    
    results_Sigma_inv_plus_BtWB$n[i] <- N[n]
    results_Sigma_inv_plus_BtWB$num_optim_iter[i] <- model_Sigma_inv_plus_BtWB$get_num_optim_iter()
    results_Sigma_inv_plus_BtWB$sigma2[i] <- model_Sigma_inv_plus_BtWB$get_cov_pars()[1]
    results_Sigma_inv_plus_BtWB$rho[i] <- model_Sigma_inv_plus_BtWB$get_cov_pars()[2]
    
    i <- i + 1
    print(i)
    saveRDS(results_Sigma_inv_plus_BtWB, "./data/bias_analysis.rds")
    gc()
  }
}
