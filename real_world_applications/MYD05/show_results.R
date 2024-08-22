library(fields)
library(xtable)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
L_MYD05 <- readRDS("../../data/MYD05/L_MYD05.rds")

models <- c("CholeskyVL", "IterativeVL", "GPVecchia", "FITC")

results <- vector(mode = "list", length = length(models))
names(results) <- models
for(model in models){
  results[[model]] <- readRDS(paste0("./", model, ".rds"))
}

#Runtime nll
runtime <-  readRDS("./runtime.rds")
runtime <- runtime[runtime$n==5e4,] #for n=50'000

#GPVecchia
results[["GPVecchia"]]$estimates$cov_pars <- results[["GPVecchia"]]$estimates$cov_pars[1:2] #remove smoothness param
results[["GPVecchia"]]$estimates$cov_pars[2] <- results[["GPVecchia"]]$estimates$cov_pars[2]*sqrt(2*1.5) #adjust range

##Table performance#############################################################
p_cols <- c("RMSE", "CRPS", "TimeLL")
performance <- matrix(nrow = length(models), ncol = length(p_cols))
rownames(performance) <- models
colnames(performance) <- p_cols

for(model in models){
  performance[rownames(performance) == model,c("RMSE", "CRPS")] <- unlist(results[[model]][c("RMSE", "CRPS")])
  performance[rownames(performance) == model,c("TimeLL")] <- runtime$time_negLL[runtime$model==model]
}
xtable(performance)

##Table estimates###############################################################
e_cols <- c("alpha", "variance", "range", "beta_0", "beta_1", "beta_2")
estimates <- matrix(nrow = length(models), ncol = length(e_cols))
rownames(estimates) <- models
colnames(estimates) <- e_cols
for(model in models){
  estimates[rownames(estimates) == model,] <- unlist(results[[model]][["estimates"]][c("alpha", "cov_pars", "betas")])
}
xtable(estimates)

##Plot response#################################################################
coltab <- tim.colors(64)
zr<- range(c(L_MYD05$L_vapor_test, 
             results$CholeskyVL$pred_response_mu,
             results$IterativeVL$pred_response_mu,
             results$GPVecchia$pred_response_mu,
             results$FITC$pred_response_mu), na.rm=TRUE)

#Plot first line
par(mar = c(2, 0.5, 2, 1))
par(oma=c(0,1.5,0,1))
set.panel(1,3)
quilt.plot(L_MYD05$L_locations_test, L_MYD05$L_vapor_test, col=coltab, zlim=zr, nx=100, ny=200, main="Test Data", add.legend = FALSE)
quilt.plot(L_MYD05$L_locations_test, results$CholeskyVL$pred_response_mu, zlim=zr, col=coltab, nx=100, ny=200, main="Cholesky-VL", add.legend = FALSE, yaxt='n')
quilt.plot(L_MYD05$L_locations_test, results$IterativeVL$pred_response_mu, zlim=zr, col=coltab, nx=100, ny=200, main="Iterative-VL", add.legend = FALSE, yaxt='n')

#Plot second line
set.panel(1,2)
par(mar = c(2, 0.5, 2, 1))
par(oma=c(0,1.5,0,3))
quilt.plot(L_MYD05$L_locations_test, results$GPVecchia$pred_response_mu, zlim=zr, col=coltab, nx=100, ny=200, main="GPVecchia", add.legend = FALSE, 
           cex.main=0.85, cex.axis=0.7, mgp = c(0.5, .75, 0))
quilt.plot(L_MYD05$L_locations_test, results$FITC$pred_response_mu, zlim=zr, col=coltab, nx=100, ny=200, main="LowRank", add.legend = FALSE, 
           cex.main=0.85, cex.axis=0.7, mgp = c(0.5, .75, 0), yaxt='n')

par(oma=c(0,0,0,1))
image.plot(zlim=zr,legend.only=TRUE, axis.args = list(cex.axis = .7))

##Plot variance#################################################################
coltab <- tim.colors(64)
par(mar=c(1.75,1.4,1.5,1.5))
par(oma=c(0,0.5,0,0.5))
set.panel(2,2)
quilt.plot(L_MYD05$L_locations_test, results$CholeskyVL$pred_response_var, col=coltab, nx=100, ny=200, main="Cholesky-VL",
           cex.main=0.95, cex.axis=0.8, mgp = c(0.5, .75, 0), axis.args = list(cex.axis = .8))
quilt.plot(L_MYD05$L_locations_test, results$IterativeVL$pred_response_var, col=coltab, nx=100, ny=200, main="Iterative-VL", yaxt='n',
           cex.main=0.95, cex.axis=0.8, mgp = c(0.5, .75, 0), axis.args = list(cex.axis = .8))
quilt.plot(L_MYD05$L_locations_test, results$GPVecchia$pred_response_var, col=coltab, nx=100, ny=200, main="GPVecchia",
           cex.main=0.95, cex.axis=0.8, mgp = c(0.5, .75, 0), axis.args = list(cex.axis = .8))
quilt.plot(L_MYD05$L_locations_test, results$FITC$pred_response_var, col=coltab, nx=100, ny=200, main="LowRank", yaxt='n',
           cex.main=0.95, cex.axis=0.8, mgp = c(0.5, .75, 0), axis.args = list(cex.axis = .8))
