library(PoPS)
library(raster)
#on each loop get new temp/host/precip/etc map

# create needed inputs (Rasters)
# all constants #

infected_file <- system.file

#give input maps same extent and projection as temp & precip
infected <- NLMR::nlm_random(20,20,100,TRUE)
values(infected) <- values(infected) * 0
crs(infected) <- "+proj=longlat +datum=WGS84 +no_defs"
writeRaster(infected, "infected_file.tif", overwrite=TRUE)
infected_file <- "infected_file.tif"

host <- NLMR::nlm_random(20,20,100,TRUE)
values(host) <- values(host) * 5
crs(host) <- crs(infected)
extent(host) <- extent(infected)


host_sd <- NLMR::nlm_random(20,20,100,TRUE)
crs(host_sd) <- crs(host)
extent(host_sd) <- extent(host)

host_stack <- stack(host, host_sd)

plot(host_stack)
crs(host_stack) <- crs(host)
extent(host_stack) <- extent(host)
writeRaster(host_stack, "host_file.tif", overwrite=TRUE)
host_file <- "host_file.tif"

## TODO should total_populations be created from infected and/or host? ##
total_populations <- host
writeRaster(total_populations, "total_populations_file.tif", overwrite=TRUE)
total_populations_file <- "total_populations_file.tif"

parameter_means <- read.csv("2018_2019_means.csv")$X
parameter_cov_matrix  <- read.csv("2018_2019_cov_matrix.csv")

temp <- FALSE
temperature_coefficient_file <- ""
precip <- FALSE
precipitation_coefficient_file <- ""
model_type <- "SI"
latency_period <- 0

## TODO ##
time_step <- "month"
start_date <- '2003-01-01'
end_date <- '2003-12-31'
## ##
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
use_movement <- FALSE
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
infected_file <-  system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
host_file <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
total_populations_file <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
total_populations_file <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")

infected <- raster(infected_file)
plot(infected)
host <- raster(host_file)
plot(host)
total_pop <- raster(total_populations_file)
plot(total_pop)
# call pops_multirun
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
                      use_movement,
                      start_exposed,
                      generate_stochasticity,
                      establishment_stochasticity,
                      movement_stochasticity,
                      deterministic,
                      establishment_probability,
                      dispersal_percentage)

data <- PoPS::pops(
  infected_file, 
  host_file, 
  total_populations_file,
  parameter_means,
  parameter_cov_matrix
)

# record results in numeric vector

# call sobol_indices with results from pops_multirun

