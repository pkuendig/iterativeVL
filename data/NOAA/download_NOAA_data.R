library(lubridate)
library(fields)

# library(remotes)
# install_version("rnoaa", "1.4.0")
library(rnoaa)

################################################################################
month_1day <- as.Date("2019-11-01")
month_days <- seq(floor_date(month_1day, 'month'), 
                  ceiling_date(month_1day, 'month') - 1, by = 1)

month_data <- c()
for(d in 1:length(month_days)){
  day_data <- cpc_prcp(date = as.character(month_days[d]), us=TRUE, drop_undefined = TRUE)
  print(dim(day_data))
  if(d == 1){
    month_data <- day_data
  }else{
    month_data <- merge(x=month_data, y=day_data, by=c("lon", "lat"), all=TRUE)
  }
  colnames(month_data)[d+2] <- paste0("precip", d)
}

if(any(is.na(month_data))) stop("NA's found")

month_data$total <- rowSums(month_data[,3:ncol(month_data)])

Y <- month_data$total
coords <- as.matrix(month_data[,c("lon", "lat")])
X <- as.matrix(cbind(1, coords))
colnames(X) <- c("Intercept", "lon", "lat")

##Split 20% Test/ 80% Training##################################################

n <- length(Y)
n_test <- round(0.2*n)

set.seed(1)
i_test <- sample(1:n, n_test, replace=FALSE)
i_train <- 1:n
i_train <- i_train[!(i_train %in% i_test)]

Y_test <- Y[i_test]
coords_test <- coords[i_test,]
X_test <- X[i_test,]

Y_train <- Y[i_train]
coords_train <- coords[i_train,]
X_train <- X[i_train,]

#Train
quilt.plot(cbind(coords_train, Y_train),
           ny=length(unique(month_data$lat)), nx=length(unique(month_data$lon)))
#Test
quilt.plot(cbind(coords_test, Y_test),
           ny=length(unique(month_data$lat)), nx=length(unique(month_data$lon)))

saveRDS(list(Y_train=Y_train,
             coords_train=coords_train,
             X_train=X_train,
             Y_test=Y_test,
             coords_test=coords_test,
             X_test=X_test), "NOAA.rds")
