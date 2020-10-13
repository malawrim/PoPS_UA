library(PoPS)
library(raster)
#on each loop get new host/precip/etc map

# create needed inputs (Rasters)

#give input maps same extent and projection as temp & precip
#infected <- NLMR::nlm_random(30,30,1,TRUE)
#values(infected) <- values(infected) * 0
#crs(infected) <- "+proj=longlat +datum=WGS84 +no_defs"

# infected_file <-  system.file("extdata", "SODexample", "initial_infections.tif", package = "PoPS")
# host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
# #total_populations_file <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
# total_populations_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")

# infected_SOD <- raster(infected_file)
# plot(infected_SOD)
# host_SOD <- raster(host_file) 
# plot(host_SOD)
# total_pop_SOD <- raster(total_populations_file)
# plot(total_pop_SOD)
# call pops_multirun

infected <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                   crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                   vals=as.integer(stats::rnorm(16200, mean=0, sd=2)))
infected[infected < 0] <- 0
values(infected) <- round(values(infected), 0)
plot(infected)
# writeRaster(infected, "infected_file.tif", overwrite=TRUE)
# probably don't want sd randomized #
infected_sd <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                   crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                   vals=runif(1, min=0, max=1))
plot(infected_sd)
infected_stack <- stack(infected, infected_sd)
infected_brick <- brick(infected_stack)
plot(infected_brick)
infected_file <- writeRaster(infected_brick, "infected_file.tif", format="GTiff", overwrite=TRUE)
plot(infected_file)
infected_file <- "infected_file.tif"

host <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                   crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                   vals=as.integer(stats::rnorm(16200, mean=4, sd=1)))
writeRaster(host, "host_file.tif", overwrite=TRUE)
host_file <- "host_file.tif"
plot(host)

total_populations <- raster(nrows=100, ncols=100, xmn=0, ymn=0,
                            crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                            vals=9)
writeRaster(total_populations, "total_populations_file.tif", overwrite=TRUE)
total_populations_file <- "total_populations_file.tif"
plot(total_populations)

#host <- NLMR::nlm_random(20,20,100,TRUE)
#values(host) <- values(host) * 5
#crs(host) <- crs(infected)
#extent(host) <- extent(infected)

#host_sd <- NLMR::nlm_random(20,20,100,TRUE)
#crs(host_sd) <- crs(host)
#extent(host_sd) <- extent(host)

#host_stack <- stack(host, host_sd)
#plot(host_stack)
#crs(host_stack) <- crs(host)
#extent(host_stack) <- extent(host)

#ex <- landscapemetrics::landscape
#plot(ex)
#crs(ex) <- crs(infected)
#extent(ex) <- extent(infected)
#res(ex) <- 1
#writeRaster(ex, "host_file.tif", overwrite=TRUE)
#host_file <- "host_file.tif"

#total_populations <- ex
#writeRaster(total_populations, "total_populations_file.tif", overwrite=TRUE)
#total_populations_file <- "total_populations_file.tif"

parameter_means <- c(2.55, 1.45, 1, 0, 0, 0)
parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
parameter_cov_matrix[1,1] <- 2.921333
parameter_cov_matrix[2,1] <- 8.265058
parameter_cov_matrix[1,2] <- 8.265058
parameter_cov_matrix[2,2] <- 23.469279
# parameter_means <- read.csv("means_control.csv")
# parameter_means <- c(unlist(parameter_means, use.names=FALSE))
# 
# parameter_cov_matrix <- read.csv("cov_control.csv")
# parameter_cov_matrix <- matrix(unlist(parameter_cov_matrix), ncol = 6, byrow = TRUE)
# parameter_cov_matrix <- t(parameter_cov_matrix)
# parameter_cov_matrix <- parameter_cov_matrix[,-1]
temp <- FALSE
temperature_coefficient_file <- ""
precip <- FALSE
precipitation_coefficient_file <- ""
model_type <- "SI"
latency_period <- 0

time_step <- "month"
start_date <- '2003-01-01'
end_date <- '2003-12-31'

season_month_start <- 1
season_month_end <- 12
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -12.87
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c(0)
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
anthropogenic_dir <- "NONE"
number_of_iterations <- 10
number_of_cores <- NA
pesticide_duration <- 0
pesticide_efficacy <- 1
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
# call sobol_matrices or create own matrix
# since the only thing we are altering for initial condition is std-dev of
# infected then matrix shouldn't be necessary?

data_1_new <- PoPS::pops_multirun(infected_file, 
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

# save(data_1, file="data_1.RData")
# save(data_sd1, file="data_sd1.RData")
# save(data_sd0, file='data_sd0.RData')
# save(data_sd0.5, file="data_sd0.5.RData")

# temp_coefficient_file <- system.file("extdata", "simple2x2", "temperature_coefficient.tif", package = "PoPS")
# temp <- raster(temp_coefficient_file)
# plot(temp$temperature_coefficient)
# 
# temp_coefficient_days <- system.file("extdata", "simple2x2", "temperature_coefficient_days.tif", package = "PoPS")
# temp_days <- raster(temp_coefficient_days)
# plot(temp_days)
# 
# temp_coefficient_weeks <- system.file("extdata", "simple2x2", "temperature_coefficient_weeks.tif", package = "PoPS")
# temp_weeks <- raster(temp_coefficient_weeks)
# plot(temp_weeks)
# record results in numeric vector

pops_output <- matrix(c(4.5, 6.7, 9.0, 10.0, 12.3, 3.2, 6.7, 9.1), ncol=2, byrow=FALSE)
pops_params <- c("V1", "V2")
pops_n <- 2
pops_R <- 100

library(sensobol)
# # Define settings:
# n <- 1000 #sample size of the sample matrix.
# # k <- 8
# R <- 100 #number of bootstrap replicas.
# # Design the sample matrix:
# A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
# # Compute the model output:
# Y <- sobol_Fun(A)
# # Compute the Sobol' indices:
# sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)),
#                       R = R, n = n, parallel = "no", ncpus = 1, second = TRUE, third = TRUE)
sens <- sobol_indices(Y = pops_output, params = pops_params,
                      R = pops_R, n = pops_n, parallel = "no", ncpus = 1, second = TRUE, third = TRUE)

# call sobol_indices with results from pops_multirun

