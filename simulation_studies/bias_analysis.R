################################################################################
# Parameter estimation for different sample sizes with iterative Vecchia-Laplace
################################################################################

library(gpboost)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/make_data.R")

N <- c(500,1000,2000,5000,10000,20000,50000)
n_rep <- 100
sigma2_true <- 1
rho_true <- 1/20
true_covpars <- c(sigma2_true, rho_true)
NN <- 20 

res_cols <- c("n", "sigma2", "rho", "num_optim_iter")
results <- data.frame(matrix(nrow=length(N)*n_rep, ncol = length(res_cols)))
colnames(results) <- res_cols

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
    model_iterativeVL <- GPModel(gp_coords = coords_train,
                                         cov_function = "matern",
                                         cov_fct_shape=1.5,
                                         likelihood="bernoulli_logit",
                                         matrix_inversion_method = "iterative",
                                         gp_approx = "vecchia",
                                         vecchia_ordering = "random",
                                         num_neighbors=NN)
    
    model_iterativeVL$set_optim_params(params = list(maxit=1000,
                                                            trace=TRUE,
                                                            cg_preconditioner_type = "Sigma_inv_plus_BtWB",
                                                            seed_rand_vec_trace = i))
    
    model_iterativeVL$fit(y=y_train)
    
    results$n[i] <- N[n]
    results$num_optim_iter[i] <- model_iterativeVL$get_num_optim_iter()
    results$sigma2[i] <- model_iterativeVL$get_cov_pars()[1]
    results$rho[i] <- model_iterativeVL$get_cov_pars()[2]
    
    i <- i + 1
    print(i)
    gc()
  }
}

############################################################################
# Plotting
############################################################################
require(ggplot2)

options(scipen = 999)
results$n <- as.factor(results$n)

#sigma^2 with Errorbar
ggplot(results, aes(x=n, y=sigma2, group=n)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               linewidth=0.5,
               width = 0.8,
               fun.max = function(x) mean(x) + 2*sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - 2*sd(x) / sqrt(length(x))) + 
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE)  +
  theme_bw() + ylab(expression(sigma[1]^2)) + geom_hline(yintercept=sigma2_true, linetype="dashed")
