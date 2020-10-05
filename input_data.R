rand <- NLMR::nlm_random(5,5,1,TRUE)
landscapetools::show_landscape(rand)
rand <- rand*0 + 0.2

#pull precip data from noaa
#ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/precip/total/daily/p.full.1stday_month_yearmtdy.tif
#ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/temp/total/daily/t.full.1stday_month_yearmtdy.tif
download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/temp/total/daily/t.full.1stday_month_20130611.tif", "temp.tif", "auto", mode = "wb")
download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/precip/total/daily/p.full.1stday_month_20130611.tif", "precip.tif", "auto", mode = "wb")

#import temp and precip data
tif_name <- 'temp.tif'
temp <- raster(tif_name, layer=0)
landscapetools::show_landscape(temp)

tif_name <- 'precip.tif'
precip <- raster(tif_name, layer=0)

#import clip extent
nc_cnty <- rgdal::readOGR("NC_wake_county.shp", "NC_wake_county")
plot(nc_cnty)

str(temp)
str(nc_cnty)
crs(temp)
crs(nc_cnty)
extent(temp)
extent(nc_cnty)

#match projections
nc_cnty_reproject <- spTransform(nc_cnty, crs(temp))
crs(nc_cnty_reproject)
extent(nc_cnty_reproject)

#crop to same area
crop_temp <- crop(temp, nc_cnty_reproject)
plot(crop_temp)

crop_precip <- crop(precip, nc_cnty_reproject)
plot(crop_precip)
extent(crop_temp)
extent(crop_precip)

#give input maps same extent and projection as temp & precip
rand <- NLMR::nlm_random(10,10,1,TRUE)
crs(rand) <- crs(temp)
extent(rand) <- extent(nc_cnty_reproject)

#create a raster of desired extent
ras <- raster(nrow = 3, ncol = 3)
values(ras) <- 3
crs(ras) <- crs(temp)
extent(ras) <- extent(nc_cnty_reproject)
plot(ras)
# convert to polygon to crop
poly <- rasterToPolygons(ras, fun = function(x){x==3})
plot(poly)
precip_second_crop <- crop(precip, poly)
plot(precip_second_crop)
plot(precip)

# testing new library option
library(sensobol)
# Define settings:
n <- 1000; k <- 8; R <- 100
# Design the sample matrix:
# outputs a matrix
A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
# Compute the model output:
# outputs a numeric vector with model output
Y <- sobol_Fun(A)
# Compute the Sobol' indices:
# outputs data.table with results 
sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)),
                      R = R, n = n, parallel = "no", ncpus = 1, second = TRUE, third = TRUE)
