library(PoPS)
library(raster)
# call sobol_indices with results from pops_multirun
library(sensobol)

# create matrix of potential sd inputs for infected
# inputs duplicated below
# sample size
pops_n <- 4
# number of inputs
pops_k <- 5
# matrix is of size n * 2k
# matrix should be created in respect to pdf of param
pops_matrices <- sobol_matrices(n = pops_n, k = pops_k, second = TRUE, third = TRUE)
l <- pops_n * (2 + pops_k)
count <- c(1:108)

data_list <- list(list())
# access element in 2D list data_list[[1]][1]
# access whole list in 2D list data_list[[1]]
#infected_old <- infected
infected <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                   crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                   vals=as.integer(stats::rnorm(16200, mean=0, sd=2)))
infected[infected < 0] <- 0
values(infected) <- round(values(infected), 0)
# plot(infected)
writeRaster(infected, "infected_file.tif", overwrite=TRUE)
# plot(infected_file)
infected_file <- "infected_file.tif"

host <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
               crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
               vals=as.integer(stats::rnorm(16200, mean=4, sd=1)))
# plot(host)
host_file <- writeRaster(host, "host_file.tif", overwrite=TRUE)
# plot(host_file)
host_file <- "host_file.tif"

total_populations <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                            crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                            vals=9)
writeRaster(total_populations, "total_populations_file.tif", overwrite=TRUE)
total_populations_file <- "total_populations_file.tif"
plot(total_populations)

parameter_means <- c(2.55, 1.45, 1, 0, 0, 0)
parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
parameter_cov_matrix[1,1] <- 2.921333
parameter_cov_matrix[2,1] <- 8.265058
parameter_cov_matrix[1,2] <- 8.265058
parameter_cov_matrix[2,2] <- 23.469279
# constant #
model_type <- "SI"
latency_period <- 0
time_step <- "month"
start_date <- '2003-01-01'
end_date <- '2003-12-31'
season_month_start <- 1
season_month_end <- 12
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
anthropogenic_dir <- "NONE"
number_of_iterations <- 2
number_of_cores <- NA

random_seed <- NULL
output_frequency <- "year"
output_frequency_n <- 1
movements_file <- ""
use_movements <- FALSE
start_exposed <- FALSE
generate_stochasticity <- FALSE
establishment_stochasticity <- FALSE
movement_stochasticity <- FALSE
deterministic <- TRUE
establishment_probability <- 0.5
dispersal_percentage <- 0.99
quarantine_areas_file <- ""
use_quarantine <- FALSE
use_spreadrates <- FALSE

for ( i in count ) {
  # download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/temp/total/daily/t.full.1stday_month_20130611.tif", "temp.tif", "auto", mode = "wb")
  # download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/precip/total/daily/p.full.1stday_month_20130611.tif", "precip.tif", "auto", mode = "wb")
  # 
  # # import temp and precip data
  # tif_name <- 'temp.tif'
  # temp <- raster(tif_name, values=TRUE)
  temp <- TRUE
  temp_raster <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                    crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                    vals=pops_matrices[i,1])
  temp_stack <- stack(temp_raster, temp_raster, temp_raster, temp_raster, temp_raster, temp_raster, temp_raster, temp_raster, temp_raster, temp_raster, temp_raster, temp_raster)
  temp_file <- writeRaster(temp_stack, "temp_file.tif", type= "GTIFF", overwrite=TRUE)
  temperature_coefficient_file <- 'temp_file.tif'
  # values(temp2)
  # 
  # download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/temp/total/daily/t.full.1stday_month_20130612.tif", "temp2.tif", "auto", mode = "wb")
  # tif_name <- 'temp2.tif'
  # temp2 <- raster(tif_name, layer=0)
  # 
  # temp_merge <- merge(temp, temp2)
  # plot(temp2)
  # 
  # tif_name <- 'precip.tif'
  # precip <- raster(tif_name, layer=0)
  precip <- TRUE
  precip_raster <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                 crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                 vals=pops_matrices[i,2])
  precip_stack <- stack(precip_raster, precip_raster, precip_raster, precip_raster, precip_raster, precip_raster, precip_raster, precip_raster, precip_raster, precip_raster, precip_raster, precip_raster)
  precip_file <- writeRaster(precip_raster, "precip_file.tif", type= "GTIFF", overwrite=TRUE)
  precipitation_coefficient_file <- 'precip_file.tif'
  
  use_lethal_temperature <- FALSE
  # min temp map
  temperature_file <- ""
  lethal_temperature <- -12.87
  lethal_temperature_month <- 1
  mortality_on <- TRUE
  mortality_rate <- pops_matrices[i,3]
  # multiply time lag by number of years simulation - needs to be greater than 1
  num_years <- 10
  mortality_time_lag <- pops_matrices[i,4] * num_years

  management <- TRUE
  treatment_dates <- c("2003-01-01")
  # one layer per timestep (0 or 1)
  treatment <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                   crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                   vals=1)
  
  treatments_raster <- writeRaster(treatment, "treatments_file.tif", overwrite=TRUE)
  treatments_file <- c("treatments_file.tif")
  treatment_method <- "ratio"
  # may just keep constant to match timestep of treatments
  pesticide_duration <- c(1)
  pesticide_efficacy <- pops_matrices[i,5]
  
  # call sobol_matrices or create own matrix
  # since the only thing we are altering for initial condition is std-dev of
  # infected then matrix shouldn't be necessary?
  
  data <- PoPS::pops_multirun(infected_file, 
                              host_file, 
                              total_populations_file,
                              parameter_means,
                              parameter_cov_matrix,
                              temp, 
                              temperature_coefficient_file, 
                              precip, 
                              precipitation_coefficient_file,
                              model_type,
                              latency_period,
                              time_step,
                              season_month_start,
                              season_month_end,
                              start_date,
                              end_date,
                              use_lethal_temperature,
                              temperature_file,
                              lethal_temperature, 
                              lethal_temperature_month,
                              mortality_on, 
                              mortality_rate, 
                              mortality_time_lag, 
                              management,
                              treatment_dates, 
                              treatments_file,
                              treatment_method,
                              natural_kernel_type, 
                              anthropogenic_kernel_type,
                              natural_dir, 
                              anthropogenic_dir,
                              number_of_iterations, 
                              number_of_cores,
                              pesticide_duration,
                              pesticide_efficacy,
                              random_seed,
                              output_frequency,
                              output_frequency_n,
                              movements_file, 
                              use_movements,
                              start_exposed,
                              generate_stochasticity,
                              establishment_stochasticity,
                              movement_stochasticity,
                              deterministic,
                              establishment_probability,
                              dispersal_percentage)
  
  # cellStats(data$simulation_mean, 'mean')
  # cellStats(data$simulation_sd, 'mean')
  
  # maybe just use these?
  # both have mean and standard deviation
  # order: number of infected mean, number of infected sd, area infected mean, area infected sd
  data_list[[i]] <- c(data$number_infecteds[1], data$number_infecteds[2], data$infected_areas[1], data$infected_areas[2])
  
  # save(data_1, file="data_1.RData")
  # save(data_sd1, file="data_sd1.RData")
  # save(data_sd0, file='data_sd0.RData')
  # save(data_sd0.5, file="data_sd0.5.RData")
}
matrix_data_list <- matrix(unlist(data_list), nrow=length(data_list), byrow=TRUE)

pops_output <- matrix_data_list[,1]
pops_params <- c("temp", "precip", "mortality_rate", "mortality_time_lag", "pesticide_efficacy")

# number of bootstrap replicas
pops_R <- 5000

plot_uncertainty(pops_output, pops_n)
# Compute the Sobol' indices:
# will have to separate indices for each of the four results in the output list (num infect mean vs sd)
pops_sens <- sobol_indices(Y = pops_output, params = pops_params, type= "saltelli",
                           R = pops_R, n = pops_n, parallel = "no", ncpus = 1, second = FALSE, third = FALSE)

# pops_replicas <- sobol_replicas(pops_sens, pops_k, second=FALSE, third=FALSE)

pops_dummy <- sobol_dummy(pops_output, pops_params, pops_R, pops_n)
pops_dummy_ci <- sobol_ci_dummy(pops_dummy, type= "norm", conf = 0.95)
# compute confidence intervals
# only works with 2+ params
pops_ci <- sobol_ci(pops_sens, params = pops_params, type = "norm", conf = 0.95, second = FALSE, third = FALSE)

plot_scatter(pops_matrices, pops_output, pops_n, pops_params)

plot_sobol(pops_ci, dummy = pops_dummy_ci, type = 1)
