library(ncdf4)
library(fields)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
nc <- nc_open("MYD05_L2.A2019087.1345.061.2019088165734.hdf")

#Read vapor data (1km granularity)
v1 <- nc$var[["Water_Vapor_Near_Infrared"]]
vapor_original_mat <- ncvar_get(nc, v1)

#Read Longitude + Latitude (only 5km granularity) https://atmosphere-imager.gsfc.nasa.gov/sites/default/files/ModAtmo/MYD05_L2.C6.CDL.fs
v2 <- nc$var[["Latitude"]]
lat_original_mat <- ncvar_get(nc, v2)
max(lat_original_mat) #top: 67.4436
min(lat_original_mat) #button: 45.25951
v3 <- nc$var[["Longitude"]]
long_original_mat <- ncvar_get(nc, v3)
max(long_original_mat) #right: 4.26105
min(long_original_mat) #left: -42.24899

#image.plot(vapor_original_mat[nrow(vapor_original_mat):1,], x=1:nrow(vapor_original_mat), y=1:ncol(vapor_original_mat)) #revert row-order -> most accurate way

#Create linear artifical coordinates of 1 unit via raster:
library(raster)
vapor_raster <- raster(t(vapor_original_mat[nrow(vapor_original_mat):1,ncol(vapor_original_mat):1]), xmn=0, xmx=1354, ymn=0, ymx=2030)
vapor_df_na = as.data.frame(vapor_raster, xy=TRUE) #create artifical coordinates
vapor_df = vapor_df_na[!is.na(vapor_df_na$layer),]
vapor <- vapor_df$layer
locations <- vapor_df[,c("x","y")]

#quilt.plot(locations, vapor, nx=ncol(vapor_original_mat)/2, ny=nrow(vapor_original_mat)/2)

MYD05 = list(vapor=vapor, locations=locations)
saveRDS(MYD05, file="MYD05.rds")
