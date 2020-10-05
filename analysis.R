##
#calibrate
##

devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "master", force = TRUE)
library(PoPS)
library(raster)
infected_years_file <- system.file("extdata", "simple20x20", "infected_single.tif", package = "PoPS")
infected_file <-  system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
host_file <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
total_plants_file <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
#temperature_coefficient_file <- system.file("extdata", "SODexample", "weather.tif", package = "PoPS")
#treatments_file <- system.file("extdata", "SODexample", "management.tif", package = "PoPS")
temperature_file <- ""
temperature_coefficient_file <- ""
precipitation_coefficient_file <-""
use_lethal_temperature <- FALSE
temp <- FALSE
precip <- FALSE
season_month_start <- 1
season_month_end <- 12
time_step <- "month"
start_date <- '2003-01-01'
end_date <- '2003-12-31'
lethal_temperature <- -35
lethal_temperature_month <- 1
random_seed <- 42
treatments_file <- ""
treatment_dates <- c('2003-01-24')
treatment_method <- "ratio"
management <- FALSE
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
natural_kernel_type <- "cauchy"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
anthropogenic_dir <- "NONE"
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
random_seed = NULL
output_frequency = "year"
movements_file <- ""
use_movements <- FALSE
number_of_iterations = 2
number_of_cores <- 2
model_type <- "SI"
latency_period <- 0
number_of_observations <- 68
prior_number_of_observations <- 0
prior_means <- c(0, 0, 0, 0, 0, 0)    ### leave as 0 for now, means that you are giving them no weight
prior_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
params_to_estimate <- c(T, T, T, T, F, F)  ### 1st: reproductive rates, 2nd: natural distance, 3rd: percent natural, 4tH: anthropogenic distance, 5th Natural Kappa, 6th anthropogenic kappa
number_of_generations <- 6
generation_size <- 10
checks = c(1200, 90000, 900, 1000)
success_metric <- "number of locations and total distance"
natural_kappa <- 0
anthropogenic_kappa <- 0
mask <- NULL
percent_natural_dispersal <- 1.0
anthropogenic_distance_scale <- 0.0

data_calibrate <- abc_calibration(infected_years_file, 
                        number_of_observations, 
                        prior_number_of_observations,
                        prior_means, prior_cov_matrix, 
                        params_to_estimate,
                        number_of_generations,
                        generation_size,
                        checks,
                        infected_file, 
                        host_file, 
                        total_plants_file, 
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
                        natural_kappa, 
                        anthropogenic_dir, 
                        anthropogenic_kappa,
                        pesticide_duration, 
                        pesticide_efficacy,
                        mask, 
                        success_metric, 
                        output_frequency,
                        movements_file, 
                        use_movements)

##
#validate
##

##################
short_distance_scale <- getmode(data_calibrate$short_distance_scale)
reproductive_rate <- getmode(data_calibrate$reproductive_rate)

infected_years_file <- system.file("extdata", "simple20x20", "infected_years.tif", package = "PoPS")
parameter_means <- c(1.8, 16.4, 0.973, 7803, 0, 0)
parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
#TODO should be these (below) 2018-2019 so need dif infected etc. files
#parameter_means <- read.csv("2018_2019_means.csv")$X
#parameter_cov_matrix  <- read.csv("2018_2019_cov_matrix.csv")
infected_file <- system.file("extdata", "simple20x20", "initial_infection.tif", package = "PoPS")
host_file <- system.file("extdata", "simple20x20", "host.tif", package = "PoPS")
total_plants_file <- system.file("extdata", "simple20x20", "all_plants.tif", package = "PoPS")
temp <- FALSE
temperature_coefficient_file <- ""
precip <- FALSE
precipitation_coefficient_file <- ""
model_type = "SEI"
latency_period = 14
time_step <- "day"
season_month_start <- 1
season_month_end <- 12
start_date <- '2003-01-01'
end_date <- '2003-02-11'
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2003-01-24')
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "exponential"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
natural_kappa <- 0
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
success_metric <- "quantity and configuration"
output_frequency <- "week"
movements_file = ""
use_movements = FALSE
percent_natural_dispersal <- 1.0
anthropogenic_distance_scale <- 0.0
number_of_iterations = 10
number_of_cores = 2

data <- abc_validate(infected_years_file, 
                     number_of_iterations, 
                     number_of_cores,
                     parameter_means,
                     parameter_cov_matrix,
                     infected_file, 
                     host_file, 
                     total_plants_file, 
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
                     pesticide_duration, 
                     pesticide_efficacy,
                     mask, 
                     success_metric, 
                     output_frequency,
                     movements_file, 
                     use_movements)

############
## PoPS multi-run
###########

infected_file <-  system.file("extdata", "SODexample", "initial_infection2001.tif", package = "PoPS")
host_file <- system.file("extdata", "SODexample", "host.tif", package = "PoPS")
total_plants_file <- system.file("extdata", "SODexample", "all_plants.tif", package = "PoPS")
parameter_means <- read.csv("2018_2019_means.csv")$X
parameter_cov_matrix  <- read.csv("2018_2019_cov_matrix.csv")

data <- PoPS::pops(
infected_file, 
host_file, 
total_plants_file,
parameter_means,
parameter_cov_matrix
)
nrow(parameter_cov_matrix)

data <- PoPS::pops_multirun(infected_file, 
                      host_file, 
                      total_plants_file,
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
                      FALSE,
                      FALSE,
                      FALSE,
                      FALSE,
                      TRUE,
                      0.5,
                      0.99)


