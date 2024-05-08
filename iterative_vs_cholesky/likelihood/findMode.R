################################################################################
# Preconditioner comparison based on likelihood
################################################################################
library(gpboost)

################################################################################
#Generate data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../../data_sets/simulated/make_data.R")

sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
n <- 100000
NN <- 20

set.seed(1)
mydata <- make_data(n=n,
                    likelihood = "bernoulli_logit",
                    sigma2=sigma2_true,
                    rho=rho_true,
                    cov_function="matern",
                    matern_nu=1.5,
                    with_fixed_effects=FALSE)

################################################################################
NUM_RAND_VEC_TRACE <- c(10, 20, 50, 100)
PRECONDITIONER <- c("piv_chol_on_Sigma", "Sigma_inv_plus_BtWB")
n_rep <- 100

VLresult <- NA
VLtime <- NA

Itresults <- data.frame(matrix(nrow=length(NUM_RAND_VEC_TRACE)*length(PRECONDITIONER)*n_rep,ncol = 4))
colnames(Itresults) <- c("preconditioner", "t", "negLL", "time")

################################################################################
#Iterative-VL

i = 1
for(p in 1:length(PRECONDITIONER)){
  for(t in 1:length(NUM_RAND_VEC_TRACE)){
    for(r in 1:n_rep){
      
      Itmodel <- GPModel(gp_coords = mydata$coords,
                         cov_function = "matern",
                         cov_fct_shape = 1.5,
                         likelihood="bernoulli_logit",
                         matrix_inversion_method = "iterative",
                         gp_approx = "vecchia",
                         vecchia_ordering = "random",
                         num_neighbors=NN)

      Itmodel$set_optim_params(params = list(maxit=1,
                                            trace=TRUE,
                                            num_rand_vec_trace=NUM_RAND_VEC_TRACE[t],
                                            cg_preconditioner_type = PRECONDITIONER[p],
                                            seed_rand_vec_trace = i))
      
      Itresults$preconditioner[i] <- PRECONDITIONER[p]
      Itresults$t[i] <- NUM_RAND_VEC_TRACE[t]
      Itresults$time[i] <- system.time(Itresults$negLL[i] <- Itmodel$neg_log_likelihood(cov_pars=true_covpars, y=mydata$y))[3]
      
      i = i+1
      saveRDS(list(VLresult=VLresult, VLtime=VLtime, Itresults=Itresults), "./data/findMode.rds")
      gc()
    }
  }
}

################################################################################
#Cholesky-VL

VLmodel <- GPModel(gp_coords = mydata$coords,
                   cov_function = "matern",
                   cov_fct_shape=1.5,
                   likelihood="bernoulli_logit",
                   matrix_inversion_method = "cholesky",
                   gp_approx = "vecchia",
                   vecchia_ordering = "random",
                   num_neighbors=NN)

VLmodel$set_optim_params(params = list(maxit=1,
                                       trace=TRUE))

VLtime <- system.time(VLresult <- VLmodel$neg_log_likelihood(cov_pars=true_covpars, y=mydata$y))[3]

saveRDS(list(VLresult=VLresult, VLtime=VLtime, Itresults=Itresults), "./data/findMode.rds")
