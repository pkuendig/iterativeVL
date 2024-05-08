require(fields)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
MYD05 <- readRDS("./MYD05.rds")

################################################################################
#Large train+test
n <- length(MYD05$vapor)

n_train <- 25e4
n_test <- 25e4

set.seed(123)
i_train <- sample(n, n_train, replace = FALSE)
i_remain <- 1:n
i_remain <- i_remain[!(i_remain %in% i_train)]
i_test <- sample(i_remain, n_test, replace=FALSE)

L_locations_train <- as.matrix(MYD05$locations[i_train,])
L_vapor_train <- MYD05$vapor[i_train]
L_X_train <- cbind(rep(1, n_train), L_locations_train)
colnames(L_X_train) <- c("Intercept", "Longitude", "Latitude")

L_locations_test <- as.matrix(MYD05$locations[i_test,])
L_vapor_test <- MYD05$vapor[i_test]
L_X_test <- cbind(rep(1, n_test), L_locations_test)
colnames(L_X_test) <- c("Intercept", "Longitude", "Latitude")

#quilt.plot(L_locations_test, L_vapor_test, nx=150, ny=300, main="Test Data")

saveRDS(list(L_vapor_train=L_vapor_train,
             L_locations_train=L_locations_train,
             L_X_train=L_X_train,
             L_vapor_test=L_vapor_test,
             L_locations_test=L_locations_test,
             L_X_test=L_X_test), "L_MYD05.rds")

################################################################################
#Small train+test (Subset from large train+test)
n <- length(MYD05$vapor)

n_train <- 1e4
n_test <- 1e4

set.seed(123)
i_train <- sample(25e4, n_train, replace = FALSE)
i_test <- sample(25e4, n_test, replace=FALSE)

S_locations_train <- L_locations_train[i_train,]
S_vapor_train <- L_vapor_train[i_train]
S_X_train <- L_X_train[i_train,]

S_locations_test <- L_locations_test[i_test,]
S_vapor_test <- L_vapor_test[i_test]
S_X_test <- L_X_test[i_test,]

#quilt.plot(S_locations_test, S_vapor_test, nx=100, ny=200, main="Test Data")

saveRDS(list(S_vapor_train=S_vapor_train,
             S_locations_train=S_locations_train,
             S_X_train=S_X_train,
             S_vapor_test=S_vapor_test,
             S_locations_test=S_locations_test,
             S_X_test=S_X_test), "S_MYD05.rds")

################################################################################
#Medium train+test (Subset from large train+test)
n <- length(MYD05$vapor)

n_train <- 5e4
n_test <- 5e4

set.seed(123)
i_train <- sample(25e4, n_train, replace = FALSE)
i_test <- sample(25e4, n_test, replace=FALSE)

M_locations_train <- L_locations_train[i_train,]
M_vapor_train <- L_vapor_train[i_train]
M_X_train <- L_X_train[i_train,]

M_locations_test <- L_locations_test[i_test,]
M_vapor_test <- L_vapor_test[i_test]
M_X_test <- L_X_test[i_test,]

#quilt.plot(M_locations_test, M_vapor_test, nx=100, ny=200, main="Test Data")

saveRDS(list(M_vapor_train=M_vapor_train,
             M_locations_train=M_locations_train,
             M_X_train=M_X_train,
             M_vapor_test=M_vapor_test,
             M_locations_test=M_locations_test,
             M_X_test=M_X_test), "M_MYD05.rds")
