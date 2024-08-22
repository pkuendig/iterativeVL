################################################################################
# Marginal log-likelihood for different ranks for the LRAC preconditioner
################################################################################
library(gpboost)

################################################################################
#Generate data

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/make_data.R")

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
RANK <- c(10,20,50,100,200,500,1000)
n_rep <- 100

VLresult <- NA
VLtime <- NA

Itresults <- data.frame(matrix(nrow=length(RANK)*n_rep,ncol = 4))
colnames(Itresults) <- c("preconditioner", "rank", "negLL", "time")

################################################################################
#Iterative-VL

i = 1
for(rank in RANK){
  for(r in 1:n_rep){
    
    Itmodel <- GPModel(gp_coords = mydata$coords,
                       cov_function = "matern",
                       cov_fct_shape=1.5,
                       likelihood="bernoulli_logit",
                       matrix_inversion_method = "iterative",
                       gp_approx = "vecchia",
                       vecchia_ordering = "random",
                       num_neighbors=NN)

    Itmodel$set_optim_params(params = list(maxit=1,
                                          trace=TRUE,
                                          num_rand_vec_trace=50,
                                          piv_chol_rank=rank,
                                          cg_preconditioner_type = "piv_chol_on_Sigma",
                                          seed_rand_vec_trace = i))
    
    Itresults$preconditioner[i] <- "P[LRAC]"
    Itresults$rank[i] <- rank
    Itresults$time[i] <- system.time(Itresults$negLL[i] <- Itmodel$neg_log_likelihood(cov_pars=true_covpars, y=mydata$y))[3]
    
    i = i+1
    gc()
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

################################################################################
# Plotting
################################################################################

library(ggplot2)
library(grid)

Itresults$rank <- as.factor(Itresults$rank)

p1 <- ggplot(Itresults, aes(x=rank, y=negLL)) + 
  geom_hline(yintercept=VLresult, linetype = "dashed") +  
  geom_boxplot(fill="#377EB8") + 
  scale_fill_brewer(type = "qual", palette=6, labels = scales::parse_format()) + 
  labs(fill  = "") + ylab("log-likelihood") + theme_bw() + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "top", 
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) + 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

p2 <- ggplot(Itresults, aes(x=rank, y=time, group=1)) +
  stat_summary(fun = mean, geom = 'line', size=1, color="#377EB8", alpha=0.9) + 
  stat_summary(fun = mean, geom = 'point', size=2, color="#377EB8", shape=1) +
  scale_y_log10() + ylab("Time (s)") + theme_bw() + xlab("k") +
  theme(legend.position = "none", 
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

grid.newpage()
grid.draw(rbind(ggplotGrob(p1), 
                ggplotGrob(p2), size = "first"))
