rm(list = ls())
setwd("")
sim_id <- 1                         # 1 or 2 corresponding to simulation 1 and simulation 2

cmdArgs <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(cmdArgs[1])

library(cmdstanr)

dat_all <- readRDS(paste0("../Data/sim", sim_id, "_rep50.rds"))
dat <- dat_all[[id]]
dat$scale_icept <- 10
dat$L <- 20  

file <- file.path("AFT_GEV_SIM_GP_GLT_scaleGEVtrunc_new.stan")
mod <- cmdstan_model(file)
fit <- mod$sample(
  data = dat, 
  iter_sampling = 5000, 
  iter_warmup = 10000, 
  seed = 100*id, 
  parallel_chains = getOption("mc.cores", 2L),
  chains = 2,
  refresh = 100)

fit$save_object(file = paste0("../Result/simulation/sim", sim_id, "/res_id",id,".RDS"))


