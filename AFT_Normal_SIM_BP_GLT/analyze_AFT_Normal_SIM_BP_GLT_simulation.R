rm(list = ls())
setwd("")
sim_id <- 1                         # 1 or 2 corresponding to simulation 1 and simulation 2

library(loo)
library(ggplot2)
library(abind)
source("../Sequential-2-Means/S2M.R")

## local robust minima
local_min <- function(x, robust_l = 4){
  x_rank <- rank(x)
  local <- 1
  if(length(x_rank) < (robust_l+1)) return(NULL)
  for(tt in 2:(length(x_rank)-robust_l+1)){
    if((x_rank[tt] < x_rank[tt-1])&(x_rank[tt]==x_rank[tt+robust_l-1])) {
      local <- tt
      break
    }
  }
  local
}

res <- list()

for(id in 1:50){
  fit <- readRDS(paste0("../Result/simulation/sim", sim_id, "/res_id",id,"_NormalBP.RDS"))
  chain_id <- 1                 #which chain is used for the posterior summaries
  
  sc_fix <- apply(fit$draws("sc_fix")[,chain_id,],3,cbind)
 
  ################################################################################
  # Variable selection via sequential 2-means, Li and Debdeep (2017)             #
  ################################################################################
  sc_fix_max <- matrix(0, nrow=nrow(sc_fix), ncol=ncol(sc_fix)/4)
  sd_max <- rep(0, nrow(sc_fix))
  for(ii in 1:ncol(sc_fix_max)){
    gg <- 1:4
    ii_total <- (gg-1)*ncol(sc_fix)/4 + ii
    sc_fix_max[,ii] <- apply(abs(sc_fix[,ii_total]), 1, max)
    sd_max <- mapply(max, sd_max, apply(sc_fix[,ii_total],1,sd))
    
  }
  sd_sc_fix <- max(sd_max)
  S2M.ob <- S2M(sc_fix_max, lower=0.01, upper=1, l=50)

  H_hat <- local_min(S2M.ob$H.b.i, robust_l=10)
  ind_select <- S2M.vs(S2M.ob, H=S2M.ob$H.b.i[H_hat])
  res[[id]] <- ind_select
}

saveRDS(res, paste0("../Result/simulation/sim", sim_id, "/ind_selected_all_NormalBP.rds"))

