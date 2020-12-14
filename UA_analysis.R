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
library(ggplot2)

load("ic_indices.Rdata")
ic_indicies <- indices
load("ic_pops_ci.Rdata")
ic_pops_ci <- pops_ci
load("ic_pops_dummy.Rdata")
ic_pops_dummy <- pops_dummy
load("ic_pops_dummy_ci.Rdata")
ic_pops_dummy_ci <- pops_dummy_ci

load("driver_indices.Rdata")
driver_indices <- indices
load("driver_pops_ci.Rdata")
driver_pops_ci <- pops_ci
load("driver_pops_dummy.Rdata")
driver_pops_dummy <- pops_dummy
load("driver_pops_dummy_ci.Rdata")
driver_pops_dummy_ci <- pops_dummy_ci

load("param_indices.Rdata")
param_indices <- indices
load("param_pops_ci.Rdata")
param_pops_ci <- pops_ci
load("param_pops_dummy.Rdata")
param_pops_dummy <- pops_dummy
load("param_pops_dummy_ci.Rdata")
param_pops_dummy_ci <- pops_dummy_ci

rm(indices)
rm(pops_ci)
rm(pops_dummy)
rm(pops_dummy_ci)

ic_indicies[[1]]$V1[[1]]$t0
ic_indicies[[1]]$V1[[2]]$t0

ic_indicies[[2]]$V1[[1]]$t0
ic_indicies[[2]]$V1[[2]]$t0

ic_indicies[[3]]$V1[[1]]$t0
ic_indicies[[3]]$V1[[2]]$t0

ic_indicies[[4]]$V1[[1]]$t0
ic_indicies[[4]]$V1[[2]]$t0


driver_indices[[1]]$V1[[1]]$t0
driver_indices[[1]]$V1[[2]]$t0
driver_indices[[1]]$V1[[3]]$t0
driver_indices[[1]]$V1[[4]]$t0
driver_indices[[1]]$V1[[5]]$t0
driver_indices[[1]]$V1[[6]]$t0
driver_indices[[1]]$V1[[7]]$t0
driver_indices[[1]]$V1[[8]]$t0
driver_indices[[1]]$V1[[9]]$t0
driver_indices[[1]]$V1[[10]]$t0

driver_indices[[2]]$V1[[1]]$t0
driver_indices[[2]]$V1[[2]]$t0
driver_indices[[2]]$V1[[3]]$t0
driver_indices[[2]]$V1[[4]]$t0
driver_indices[[2]]$V1[[5]]$t0
driver_indices[[2]]$V1[[6]]$t0
driver_indices[[2]]$V1[[7]]$t0
driver_indices[[2]]$V1[[8]]$t0
driver_indices[[2]]$V1[[9]]$t0
driver_indices[[2]]$V1[[10]]$t0

driver_indices[[3]]$V1[[1]]$t0
driver_indices[[3]]$V1[[2]]$t0
driver_indices[[3]]$V1[[3]]$t0
driver_indices[[3]]$V1[[4]]$t0
driver_indices[[3]]$V1[[5]]$t0
driver_indices[[3]]$V1[[6]]$t0
driver_indices[[3]]$V1[[7]]$t0
driver_indices[[3]]$V1[[8]]$t0
driver_indices[[3]]$V1[[9]]$t0
driver_indices[[3]]$V1[[10]]$t0

driver_indices[[4]]$V1[[1]]$t0
driver_indices[[4]]$V1[[2]]$t0
driver_indices[[4]]$V1[[3]]$t0
driver_indices[[4]]$V1[[4]]$t0
driver_indices[[4]]$V1[[5]]$t0
driver_indices[[4]]$V1[[6]]$t0
driver_indices[[4]]$V1[[7]]$t0
driver_indices[[4]]$V1[[8]]$t0
driver_indices[[4]]$V1[[9]]$t0
driver_indices[[4]]$V1[[10]]$t0



param_indices[[1]]$V1[[1]]$t0
param_indices[[1]]$V1[[2]]$t0
param_indices[[1]]$V1[[3]]$t0
param_indices[[1]]$V1[[4]]$t0

param_indices[[2]]$V1[[1]]$t0
param_indices[[2]]$V1[[2]]$t0
param_indices[[2]]$V1[[3]]$t0
param_indices[[2]]$V1[[4]]$t0

param_indices[[3]]$V1[[1]]$t0
param_indices[[3]]$V1[[2]]$t0
param_indices[[3]]$V1[[3]]$t0
param_indices[[3]]$V1[[4]]$t0

param_indices[[4]]$V1[[1]]$t0
param_indices[[4]]$V1[[2]]$t0
param_indices[[4]]$V1[[3]]$t0
param_indices[[4]]$V1[[4]]$t0


# plot_name <- paste("ic_plot_uncertainty_", 1,".pdf", sep="")
# pdf(file = plot_name)
# plot_uncertainty(pops_output, pops_n)
# dev.off()
# 
# plot_name_1 <- paste("ic_plot_scatter_", 1,".pdf", sep="")
# pdf(file = plot_name_1)
# plot_scatter(pops_matrices, pops_output, pops_n, pops_params)
# dev.off()

# plot_name_2 <- paste("ic_plot_sobol_", 1,".pdf", sep="")
# pdf(file = plot_name_2)
# type determines what is plotted
# type = 1 is First and total-order Sobol' indies.
plot_sobol(ic_pops_ci[[3]], dummy = ic_pops_dummy_ci[[3]], type = 1)
# dev.off()

params <- 4
param_names <- c("reproductive rate", "natural dispersal distance", "percent natural dispersal", "anthropogenic dispersal distance")
for ( i in params ) {
  
  ic_sobol <- plot_sobol(ic_pops_ci[[i]], dummy = ic_pops_dummy_ci[[i]], type = 1)
  graph_name <- paste("Initial Condition UA -", param_names[[i]], sep = " ")
  print(ic_sobol + ggtitle(graph_name) + theme(legend.position="right")) + theme(plot.title = element_text(hjust = 0.5))
   
  param_sobol <- plot_sobol(param_pops_ci[[i]], dummy = param_pops_dummy_ci[[i]], type = 1)
  graph_name <- paste("Parameter UA -", param_names[[i]], sep = " ")
  print(param_sobol + ggtitle(graph_name) + theme(legend.position="right")) + theme(plot.title = element_text(hjust = 0.5))
  

  # plot_name <- paste("driver_plot_uncertainty_", i,".pdf", sep="")
  # pdf(file = plot_name)
  # plot_uncertainty(pops_output, pops_n)
  # dev.off()
  # 
  # plot_name_1 <- paste("driver_plot_scatter_", i,".pdf", sep="")
  # pdf(file = plot_name_1)
  # plot_scatter(pops_matrices, pops_output, pops_n, pops_params)
  # dev.off()

  # plot_name_2 <- paste("driver_plot_sobol_", i,".pdf", sep="")
  # png(file = plot_name_2)
  # type determines what is plotted
  # type = 1 is First and total-order Sobol' indices.
  driver_sobol <- plot_sobol(driver_pops_ci[[i]], dummy = driver_pops_dummy_ci[[i]], type = 1)
  graph_name <- paste("Driver UA -", param_names[[i]], sep = " ")
  print(driver_sobol + ggtitle(graph_name) + theme(legend.position="right")) + theme(plot.title = element_text(hjust = 0.5))
  # type = 2 second order Sobol' indices.
  # plot_sobol(driver_pops_ci[[i]], dummy = driver_pops_dummy_ci[[i]], type = 2)
  # type = 3 third order Sobol' indices.
  # plot_sobol(driver_pops_ci[[i]], dummy = driver_pops_dummy_ci[[i]], type = 3)
  # dev.off()
  
  # ggplot2::ggplot(b.rep, aes(value)) +
  #   geom_histogram() +
  #   labs(x = "Y",
  #        y = "Count") +
  #   facet_wrap(parameters~variable, 
  #              scales = "free") +
  #   labs(x = "Variance", 
  #        y = "Count") +
  #   theme_bw() +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         legend.background = element_rect(fill = "transparent",
  #                                          color = NA),
  #         legend.key = element_rect(fill = "transparent",
  #                                   color = NA))
}

