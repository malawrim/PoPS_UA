library(PoPS)
library(raster)
# on each loop get new temp/host/precip/etc map

# create needed inputs (Rasters)
# all constants #
# pull precip & temp data from noaa ** randomize day **
download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/temp/total/daily/t.full.1stday_month_20130611.tif", "temp.tif", "auto", mode = "wb")
download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/precip/total/daily/p.full.1stday_month_20130611.tif", "precip.tif", "auto", mode = "wb")

# import temp and precip data
tif_name <- 'temp.tif'
temp <- raster(tif_name, layer=0)

download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/temp/total/daily/t.full.1stday_month_20130612.tif", "temp2.tif", "auto", mode = "wb")
tif_name <- 'temp2.tif'
temp2 <- raster(tif_name, layer=0)

temp_merge <- merge(temp, temp2)
plot(temp)

tif_name <- 'precip.tif'
precip <- raster(tif_name, layer=0)

# SOD projection: +proj=utm +zone=10 +datum=WGS84 +units=m +no_defs 

# call sobol_matrices or create own matrix

# call pops_multirun

# record results in numeric vector

# call sobol_indices with results from pops_multirun