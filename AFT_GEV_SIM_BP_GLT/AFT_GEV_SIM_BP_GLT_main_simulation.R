rm(list = ls())
setwd("")
sim_id <- 1                         # 1 or 2 corresponding to simulation 1 and simulation 2

cmdArgs <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(cmdArgs[1])

library(cmdstanr)

dat_all <- readRDS(paste0("../Data/sim", sim_id, "_rep50.rds"))
dat <- dat_all[[id]]
dat$scale_icept <- 10
dat$M <- 20
dat$A <- matrix(1, dat$M+1, dat$M+1)
dat$A[upper.tri(dat$A)] <- 0 

file <- file.path("AFT_GEV_SIM_BP_GLT_scaleGEVtrunc_new.stan")
mod <- cmdstan_model(file)
fit <- mod$sample(
  data = dat, 
  iter_sampling = 5000, 
  iter_warmup = 10000, 
  seed = 100*id, 
  chains = 1,
  refresh = 100)

fit$save_object(file = paste0("../Result/simulation/sim", sim_id, "/res_id",id,"_GEVBP.RDS"))


