# install.packages("devtools")
library(devtools)
# devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "v1.0.2")
library(PoPS)
# install.packages("raster")
library(raster)
# call sobol_indices with results from pops_multirun
# install.packages("sensobol")
library(sensobol)
library(foreach)
library(doParallel)

# create matrix of potential sd inputs for infected
# inputs duplicated below
# sample size
pops_n <- 1000
# number of inputs
pops_k <- 10
# matrix is of size n * 2k
# matrix should be created in respect to pdf of param
pops_matrices <- sobol_matrices(n = pops_n, k = pops_k, second = TRUE, third = TRUE)
count <- nrow(pops_matrices)

list_output <- list(list())
# access element in 2D list data_list[[1]][1]
# access whole list in 2D list data_list[[1]]
#infected_old <- infected
infected <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                   crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                   vals=as.integer(stats::rnorm(16200, mean=0, sd=2)))
infected[infected < 0] <- 0
values(infected) <- round(values(infected), 0)
# plot(infected)
raster::writeRaster(infected, "infected_file.tif", overwrite=TRUE)
# plot(infected_file)
infected_file <- "infected_file.tif"

total_populations <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                            crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                            vals=9)
raster::writeRaster(total_populations, "total_populations_file.tif", overwrite=TRUE)
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
number_of_iterations <- 10
number_of_cores <- 1

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

# param_on_off <- matrix(c(FALSE, FALSE, FALSE, FALSE,
#                          TRUE, FALSE, FALSE, FALSE,
#                          FALSE, TRUE, FALSE, FALSE,
#                          FALSE, FALSE, TRUE, FALSE,
#                          FALSE, FALSE, FALSE, TRUE,
#                          TRUE, TRUE, FALSE, FALSE,
#                          TRUE, FALSE, TRUE, FALSE,
#                          TRUE, FALSE, FALSE, TRUE,
#                          FALSE, TRUE, TRUE, FALSE,
#                          FALSE, TRUE, FALSE, TRUE,
#                          FALSE, FALSE, TRUE, TRUE,
#                          TRUE, TRUE, TRUE, FALSE,
#                          TRUE, TRUE, FALSE, TRUE,
#                          TRUE, FALSE, TRUE, TRUE,
#                          FALSE, TRUE, TRUE, TRUE,
#                          TRUE, TRUE, TRUE, TRUE), nrow=16, ncol=4, byrow=TRUE)
numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

data_list <- foreach (i=1:count) %dopar% {
  
  host <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                 crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                 vals=as.integer(stats::rnorm(16200, mean=4, sd=1)))
  host_sd <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                    crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                    vals=pops_matrices[i,1])
  # plot(host_sd)
  host_stack <- raster::stack(host, host_sd)
  host_brick <- raster::brick(host_stack)
  host_file <- paste("driver_host_file_", i, ".tif", sep="")
  # plot(host_brick)
  raster::writeRaster(host_brick, host_file, format="GTiff", overwrite=TRUE)
  # plot(host_file)
  rm(host)
  rm(host_sd)
  rm(host_stack)
  rm(host_brick)
  # download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/temp/total/daily/t.full.1stday_month_20130611.tif", "temp.tif", "auto", mode = "wb")
  # download.file("ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/precip/total/daily/p.full.1stday_month_20130611.tif", "precip.tif", "auto", mode = "wb")
  # 
  # # import temp and precip data
  # tif_name <- 'temp.tif'
  # temp <- raster::raster(tif_name, values=TRUE)
  temp <-  as.logical(round(pops_matrices[i,2]))
  temp_raster <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                    crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                    vals=runif(16200, 0, 1))
  temp_stack <- raster::stack(replicate(12, temp_raster))
  temperature_coefficient_file <- paste("driver_temp_file_", i, ".tif", sep="")
  raster::writeRaster(temp_stack, temperature_coefficient_file, type= "GTIFF", overwrite=TRUE)
  rm(temp_raster)
  rm(temp_stack)
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
  precip <-  as.logical(round(pops_matrices[i,3]))
  precip_raster <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                 crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                 vals=runif(16200, 0, 1))
  precip_stack <- raster::stack(replicate(12, precip_raster))
  precipitation_coefficient_file <- paste("driver_precip_file_", i, ".tif", sep="")
  raster::writeRaster(precip_stack, precipitation_coefficient_file, type= "GTIFF", overwrite=TRUE)
  rm(precip_raster)
  rm(precip_stack)
  
  use_lethal_temperature <-  as.logical(round(pops_matrices[i,4]))
  # min temp map
  # max minimum temperature
  max_min_temp <- 0
  # min minimum temperature
  min_min_temp <- -20
  # denormalize the time lag value from pops_matrices
  min_temp_raster <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                          crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                          vals=(pops_matrices[i,5] * (max_min_temp - min_min_temp) + min_min_temp))
  min_temp_stack <- raster::stack(replicate(12, min_temp_raster))
  temperature_file <- paste("driver_min_temperature_file_", i, ".tif", sep="")
  raster::writeRaster(min_temp_stack, temperature_file, type= "GTIFF", overwrite=TRUE)
  rm(min_temp_raster)
  rm(min_temp_stack)
  
  lethal_temperature <- -10
  lethal_temperature_month <- 1
  mortality_on <- as.logical(round(pops_matrices[i,6]))
  mortality_rate <- pops_matrices[i,7]
  # max lag = max years of simulation run
  max_years <- 2
  # min is always 1
  min_years <- 1 
  # denormalize the time lag value from pops_matrices
  mortality_time_lag <- pops_matrices[i,8] * (max_years - min_years) + min_years
  if(mortality_time_lag < 1) {
    mortality_time_lag <- mortality_time_lag * 10
  }

  management <-  as.logical(round(pops_matrices[i,9]))
  treatment_dates <- c("2003-01-01")
  # one layer per timestep (0 or 1)
  treatment <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                   crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                   vals=1)
  treatments_file <- c(paste("driver_treatments_file_", i, ".tif", sep=""))
  raster::writeRaster(treatment, treatments_file[[1]], overwrite=TRUE)
  rm(treatment)
  treatment_method <- "ratio"
  # may just keep constant to match timestep of treatments
  pesticide_duration <- c(1)
  pesticide_efficacy <- pops_matrices[i,10]
  
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
  file.remove(host_file)
  file.remove(temperature_coefficient_file)
  file.remove(precipitation_coefficient_file)
  file.remove(temperature_file)
  file.remove(treatments_file)
  # cellStats(data$simulation_mean, 'mean')
  # cellStats(data$simulation_sd, 'mean')
  
  # maybe just use these?
  # both have mean and standard deviation
  # order: number of infected mean, number of infected sd, area infected mean, area infected sd
  list_output[[i]] <- c(data$number_infecteds[1], data$number_infecteds[2], data$infected_areas[1], data$infected_areas[2])
}

stopCluster(cl)

matrix_data_list <- matrix(unlist(data_list), nrow=length(data_list), byrow=TRUE)
# colnames(matrix_data_list) <- c("host", "temp", "precip", "lethal_temp", "min_temp", "mortality", "mortality_rate", "mortality_time_lag", "management", "pesticide_efficacy")
# write.table(matrix_data_list, file = "matrix_data_list.csv")

# pops_params <- c("host", "temp", "precip", "lethal_temp", "min_temp", "mortality", "mortality_rate", "mortality_time_lag", "management", "pesticide_efficacy")

pops_params <-c("host", "temp", "precip", "lethal_temp", "min_temp", "mortality", "mortality_rate", "mortality_time_lag", "management", "pesticide_efficacy")
params <- c(1:4)
indices <- list(data.frame())
pops_dummy <- list(data.frame())
pops_dummy_ci <- list(data.frame())
pops_ci <- list(data.frame())
# number of bootstrap replicas
pops_R <- 5000
for ( i in params ) {
  pops_output <- matrix_data_list[,i]

#  plot_name <- paste("driver_plot_uncertainty_", i,".pdf", sep="")
#  pdf(file = plot_name)
#  plot_uncertainty(pops_output, pops_n)
#  dev.off()
  
  # Compute the Sobol' indices:
  # will have to separate indices for each of the four results in the output list (num infect mean vs sd)
  indices[[i]] <- sobol_indices(Y = pops_output, params = pops_params, type= "saltelli",
                           R = pops_R, n = pops_n, parallel = "no", ncpus = 1, second = TRUE, third = TRUE)

  # pops_replicas <- sobol_replicas(pops_sens, pops_k, second=FALSE, third=FALSE)
  pops_dummy[[i]] <- sobol_dummy(pops_output, pops_params, pops_R, pops_n)
  pops_dummy_ci[[i]] <- sobol_ci_dummy(pops_dummy[[i]], type= "norm", conf = 0.95)
  # compute confidence intervals
  # only works with 2+ params
  pops_ci[[i]] <- sobol_ci(indices[[i]], params = pops_params, type = "norm", conf = 0.95, second = TRUE, third = TRUE)
#  plot_name_1 <- paste("driver_plot_scatter_", i,".pdf", sep="")
#  pdf(file = plot_name_1)
#  plot_scatter(pops_matrices, pops_output, pops_n, pops_params)
#  dev.off()

#  plot_name_2 <- paste("driver_plot_sobol_", i,".pdf", sep="")
#  pdf(file = plot_name_2)
#  plot_sobol(pops_ci[[i]], dummy = pops_dummy_ci[[i]], type = 1)
#  dev.off()
}

save(indices, file="driver_indices.Rdata")
save(pops_dummy, file="driver_pops_dummy.Rdata")
save(pops_dummy_ci, file="driver_pops_dummy_ci.Rdata")
save(pops_ci, file="driver_pops_ci.Rdata")
