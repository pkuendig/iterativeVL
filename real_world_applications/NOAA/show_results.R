library(fields)
library(xtable)
library(RColorBrewer)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
NOAA <- readRDS("../../data/NOAA/NOAA.rds")

models <- c("CholeskyVL", "IterativeVL", "GPVecchia", "FITC")

results <- vector(mode = "list", length = length(models))
names(results) <- models
for(model in models){
  results[[model]] <- readRDS(paste0("./", model, ".rds"))
}

#Runtime nll
runtime <-  readRDS("./runtime.rds")

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
par(mar = c(2, 0.5, 2, 1))
par(oma=c(0,1.5,0,0))
set.panel(3,2)
zr<- range(c(NOAA$Y_test, 
             results$IterativeVL$pred_response_mu,
             results$CholeskyVL$pred_response_mu, 
             results$GPVecchia$pred_response_mu,
             results$FITC$pred_response_mu), na.rm=TRUE)
nx <- 50
ny <- 24

quilt.plot(NOAA$coords_test, NOAA$Y_test, col=coltab, zlim=zr, nx=nx, ny=ny, main="Test Data", add.legend = FALSE)
quilt.plot(NOAA$coords_test, results$CholeskyVL$pred_response_mu, zlim=zr, col=coltab, nx=nx, ny=ny, main="Cholesky-VL", add.legend = FALSE, yaxt='n')
quilt.plot(NOAA$coords_test, results$IterativeVL$pred_response_mu, zlim=zr, col=coltab, nx=nx, ny=ny, main="Iterative-VL", add.legend = FALSE)
quilt.plot(NOAA$coords_test, results$GPVecchia$pred_response_mu, zlim=zr, col=coltab, nx=nx, ny=ny, main="GPVecchia", add.legend = FALSE, yaxt='n')
quilt.plot(NOAA$coords_test, results$FITC$pred_response_mu, zlim=zr, col=coltab, nx=nx, ny=ny, main="FITC", add.legend = FALSE)

par(oma=c(0,0,0,31.5))
image.plot(zlim=zr,legend.only=TRUE)

##Plot variance#################################################################
coltab <- tim.colors(64)
par(mar=c(1.75,1.4,1.5,1.5))
par(oma=c(0,0.1,0,0.5))
set.panel(2,2)

nx <- 50
ny <- 25

quilt.plot(NOAA$coords_test, results$CholeskyVL$pred_response_var, col=coltab, nx=nx, ny=ny, main="Cholesky-VL",
           cex.main=0.95, cex.axis=0.8, mgp = c(0.5, .75, 0), axis.args = list(cex.axis = .8))
quilt.plot(NOAA$coords_test, results$IterativeVL$pred_response_var, col=coltab, nx=nx, ny=ny, main="Iterative-VL", yaxt='n',
           cex.main=0.95, cex.axis=0.8, mgp = c(0.5, .75, 0), axis.args = list(cex.axis = .8))
quilt.plot(NOAA$coords_test, results$GPVecchia$pred_response_var, col=coltab, nx=nx, ny=ny, main="GPVecchia",
           cex.main=0.95, cex.axis=0.8, mgp = c(0.5, .75, 0), axis.args = list(cex.axis = .8))
quilt.plot(NOAA$coords_test, results$FITC$pred_response_var, col=coltab, nx=nx, ny=ny, main="FITC", yaxt='n',
           cex.main=0.95, cex.axis=0.8, mgp = c(0.5, .75, 0), axis.args = list(cex.axis = .8))
