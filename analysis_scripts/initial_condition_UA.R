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
# library(parallel) - may not need

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

# create matrix of potential sd inputs for infected
# inputs duplicated below
# sample size
pops_n <- 500
# number of inputs
pops_k <- 2
# # matrix is of size n * 2k
# matrix should be created in respect to pdf of param
pops_matrices <- sobol_matrices(n = pops_n, k = pops_k, second = FALSE, third = FALSE)
count <- nrow(pops_matrices)

list_output <- list(list())
# access element in 2D list data_list[[1]][1]
# access whole list in 2D list data_list[[1]]
numCores <- detectCores()
cl <- makeCluster(numCores - 1)
registerDoParallel(cl)
data_list <- foreach (i=1:count) %dopar% {

  # on each loop get new host/precip/etc map
  
  # create needed inputs (Rasters)
  
  # give input maps same extent and projection as temp & precip
  # infected <- NLMR::nlm_random(30,30,1,TRUE)
  # values(infected) <- values(infected) * 0
  # crs(infected) <- "+proj=longlat +datum=WGS84 +no_defs"
  
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
  
  infected <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                     crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                     vals=as.integer(stats::rnorm(16200, mean=0, sd=2)))
  infected[infected < 0] <- 0
  raster::values(infected) <- round(raster::values(infected), 0)
  # plot(infected)
  # writeRaster(infected, "infected_file.tif", overwrite=TRUE)
  # probably don't want sd randomized #
  infected_sd <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                        crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                        vals=pops_matrices[i,1])
  # plot(infected_sd)
  infected_stack <- raster::stack(infected, infected_sd)
  infected_brick <- raster::brick(infected_stack)
  infected_file <- paste("ic_infected_file_", i, ".tif", sep="")
  # plot(infected_brick)
  raster::writeRaster(infected_brick, infected_file, format="GTiff", overwrite=TRUE)
  # plot(infected_file)
  
  rm(infected)
  rm(infected_sd)
  rm(infected_stack)
  rm(infected_brick)
  
  host <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                 crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                 vals=as.integer(stats::rnorm(16200, mean=4, sd=1)))
  host_sd <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                        crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                        vals=0)
  # plot(host_sd)
  host_stack <- raster::stack(host, host_sd)
  host_brick <- raster::brick(host_stack)
  # plot(host_brick)
  host_file <- paste("ic_host_file_", i, ".tif", sep="")
  raster::writeRaster(host_brick, host_file, format="GTiff", overwrite=TRUE)
  # plot(host_file)
  
  rm(host)
  rm(host_sd)
  rm(host_stack)
  rm(host_brick)
  total_populations <- raster::raster(nrows=100, ncols=100, xmn=0, ymn=0,
                              crs="+proj=longlat +datum=WGS84 +no_defs", resolution=1,
                              vals=9)
  total_populations_file <- paste("ic_total_populations_file_", i, ".tif", sep="")
  raster::writeRaster(total_populations, total_populations_file, overwrite=TRUE)
  rm(total_populations)
  # plot(total_populations)
  
  # host <- NLMR::nlm_random(20,20,100,TRUE)
  # values(host) <- values(host) * 5
  # crs(host) <- crs(infected)
  # extent(host) <- extent(infected)
  # 
  # host_sd <- NLMR::nlm_random(20,20,100,TRUE)
  # crs(host_sd) <- crs(host)
  # extent(host_sd) <- extent(host)
  # 
  # host_stack <- stack(host, host_sd)
  # plot(host_stack)
  # crs(host_stack) <- crs(host)
  # extent(host_stack) <- extent(host)
  # 
  # ex <- landscapemetrics::landscape
  # plot(ex)
  # crs(ex) <- crs(infected)
  # extent(ex) <- extent(infected)
  # res(ex) <- 1
  # writeRaster(ex, "host_file.tif", overwrite=TRUE)
  # host_file <- "host_file.tif"
  # 
  # total_populations <- ex
  # writeRaster(total_populations, "total_populations_file.tif", overwrite=TRUE)
  # total_populations_file <- "total_populations_file.tif"
  
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
  
  file.remove(infected_file)
  file.remove(host_file)
  file.remove(total_populations_file)
  
  # cellStats(data$simulation_mean, 'mean')
  # cellStats(data$simulation_sd, 'mean')
  
  # maybe just use these?
  # both have mean and standard deviation
  # order: number of infected mean, number of infected sd, area infected mean, area infected sd
  list_output[[i]] <- c(data$number_infecteds[1], data$number_infecteds[2], data$infected_areas[1], data$infected_areas[2])
}
stopCluster(cl)
matrix_data_list <- matrix(unlist(data_list), nrow=length(data_list), byrow=TRUE)

pops_params <- c("infected_mean", "infected_sd", "area_infected_mean", "area_infected_sd")

params <- c(1:4)
indices <- list(data.frame())
pops_dummy <- list(data.frame())
pops_dummy_ci <- list(data.frame())
pops_ci <- list(data.frame())
# number of bootstrap replicas
pops_R <- 5000
for ( i in params ) {
  pops_output <- matrix_data_list[,i]
  
  # plot_name <- paste("ic_plot_uncertainty_", i,".pdf", sep="")
  # pdf(file = plot_name)
  # plot_uncertainty(pops_output, pops_n)
  # dev.off()
  
  # Compute the Sobol' indices:
  # will have to separate indices for each of the four results in the output list (num infect mean vs sd)
  indices[[i]] <- sobol_indices(Y = pops_output, params = pops_params, type= "saltelli",
                                R = pops_R, n = pops_n, parallel = "no", ncpus = 1, second = FALSE, third = FALSE)
  
  # pops_replicas <- sobol_replicas(pops_sens, pops_k, second=FALSE, third=FALSE)
  pops_dummy[[i]] <- sobol_dummy(pops_output, pops_params, pops_R, pops_n)
  
  pops_dummy_ci[[i]] <- sobol_ci_dummy(pops_dummy[[i]], type= "norm", conf = 0.95)
  # compute confidence intervals
  # only works with 2+ params
  pops_ci[[i]] <- sobol_ci(indices[[i]], params = pops_params, type = "norm", conf = 0.95, second = FALSE, third = FALSE)
  
  # plot_name_1 <- paste("ic_plot_scatter_", i,".pdf", sep="")
  # pdf(file = plot_name_1)
  # plot_scatter(pops_matrices, pops_output, pops_n, pops_params)
  # dev.off()
  # 
  # plot_name_2 <- paste("ic_plot_sobol_", i,".pdf", sep="")
  # pdf(file = plot_name_2)
  # plot_sobol(pops_ci[[i]], dummy = pops_dummy_ci[[i]], type = 1)
  # dev.off()
}

save(indices, file="ic_indices.Rdata")
save(pops_dummy, file="ic_pops_dummy.Rdata")
save(pops_dummy_ci, file="ic_pops_dummy_ci.Rdata")
save(pops_ci, file="ic_pops_ci.Rdata")
