rm(list = ls())

setwd("")
sim_id <- 1                         # 1 or 2 corresponding to simulation 1 and simulation 2

library(ggplot2)
library(abind)

res <- list()

for(id in 1:50){
  
  fit <- readRDS(paste0("../Result/simulation/sim", sim_id, "/res_id",id,"_SelectedVariables.RDS"))
  chain_id <- c(1,2)

  
  intercept <- apply(fit$draws("intercept")[,chain_id,],3,cbind)
  shape_GEV <- apply(fit$draws("shape_GEV")[,chain_id,],3,cbind)
  scale_GEV <- apply(fit$draws("scale_GEV")[,chain_id,],3,cbind)
  sc_fix <- apply(fit$draws("sc_fix")[,chain_id,],3,cbind)
  alpha <- apply(fit$draws("alpha")[,chain_id,],3,cbind)

  res[[id]] <- list(icept = intercept, sh_GEV = shape_GEV, sc_GEV = scale_GEV,
                    sc = sc_fix, alpha = alpha)
}

## check convergence of the single index function at test points
#plot the individual trajectory estimated from the method
# the base phi function
# eta: x*beta
phi <- function(eta) {
  if (eta < -1.0) {
    output = 0.0
  }
  else if (eta < 0.0) {
    output = 0.5*(1+eta)^2;
  }
  else if (eta < 1.0) {
    output = 0.5+0.5*(2-eta)*eta
  }
  else {
    output = 1.0
  }
  return(output)
}
# phi function for l = 1
# eta: x*beta
# u: vector of specified knots from -1 to 1
# N: sample size
phi0 <- function(eta, u, N) {
  output <- rep(0, N)
  for (i in 1:N) {
    if (eta[i] < u[2]) {
      output[i] = (u[2]*(eta[i]+1.0)+0.5-0.5*eta[i]^2.0)/(u[2]+1)
    }
    else {
      output[i] = 0.5 + 0.5*u[2]
    }
  }
  return(output)
}
# phi function for 2 <= l <= (L)
# eta: x*beta
# u: vector of specified knots from -1 to 1
# N: sample size
# l : index ranging from 1 to #knots
phil <- function(eta, u, N, l){
  output <- rep(0, N)
  for (i in 1:N) {
    if (eta[i] < u[l]) {
      output[i] = (u[l]-u[l-1])*phi((eta[i]-u[l])/(u[l]-u[l-1]))
    }
    else {
      output[i] = 0.5*(u[l]-u[l-1])+(u[l+1]-u[l])*
        (phi((eta[i]-u[l])/(u[l+1]-u[l]))-0.5)
    }
  }
  return(output)
}
# phi function for l = L + 1
# eta: x*beta
# u: vector of specified knots from -1 to 1
# N: sample size
# l : index ranging from 1 to #knots
phiL <- function(eta, u, N, l){
  output <- rep(0, N)
  for (i in 1:N) {
    if (eta[i] < u[l]) {
      output[i] = (u[l]-u[l-1])*phi((eta[i]-u[l])/(u[l]-u[l-1]))
    }
    else {
      output[i] = 0.5*(u[l]-u[l-1])+(u[l]-u[l-1])*
        (phi((eta[i]-u[l])/(u[l]-u[l-1]))-0.5)
    }
  }
  return(output)
}
# return matrix whose columns evaluate each basis (GP of L+1 basis) at the input vector eta
phi_mat <- function (eta, u, N, L){
  phiX <- matrix(0,nrow=N,ncol=L+1)  
  phiX[,1] = phi0(eta,u,N)
  for (j in 2:L) {
    phiX[,j] = phil(eta,u,N,j)
  }
  phiX[,L+1] = phiL(eta,u,N,L+1)
  return(phiX)
}


### plotting of estimated single index function and credible bands
a_si <- 5
my_sin <- function(x, a_si){
  if(max(abs(x)) > 1) return(0)
  return((3*a_si/8)*(sin(x*pi/2)+1)+(a_si/8)*(x+1))
}
mse_fun <- rep(0, 50)
x_test <- seq(-0.999,0.999,length.out = 50)
l_test <- length(x_test)
fun_est <- matrix(0, l_test, ncol=50)
L <- 20                           # the number of knots (minus 1)
u <- seq(-1,1,length=L+1)        # equal distanced grid


chain_id_final <- rep(0, 50)
for(id in (1:50)){
  alpha <- res[[id]]$alpha
  mse_twochain <- rep(0,2)
  for(jj in 1:2){
    temp <- phi_mat(x_test, u, N = l_test, L)%*%apply(abs(t(alpha[((jj-1)*(nrow(alpha)/2)+1):(jj*nrow(alpha)/2),])), 1, median)   #for better/robust results, posterior median is used
    mse_twochain[jj] <- sum((temp - my_sin(x_test, a_si))^2*(x_test[2]-x_test[1]))
  }
  chain_id_final[id] <- which(mse_twochain == min(mse_twochain))
  fun_est[,id] <- phi_mat(x_test, u, N = l_test, L)%*%apply(abs(t(alpha[((chain_id_final[id]-1)*(nrow(alpha)/2)+1):(chain_id_final[id]*nrow(alpha)/2),])), 1, median)
  mse_fun[id] <- mse_twochain[chain_id_final[id]] 
}

fun_summary <- apply(fun_est, 1, function(x) c(mean(x),sd(x),quantile(x,probs=c(.025,.975))))

df_fun <- data.frame(x = x_test, fun_est = fun_summary[1,], 
                     fun_2half5 = fun_summary[3,], fun_97half5 = fun_summary[4,],
                     true_fun = my_sin(x_test, a_si))
pdf(paste0("../Result/figure/simulation/sim", sim_id, "/function_estimate_GEV_GP_SelectedVariables.pdf"))
par(mfrow=c(1,1),mar=c(2,2,2,1))
p1 <- ggplot(df_fun, aes(x = x, y = fun_est))+
  geom_line()+
  geom_ribbon(data=df_fun, aes(ymin=fun_2half5, ymax=fun_97half5), alpha=0.3) +
  geom_line(data=df_fun, aes(x=x, y=true_fun), color="red")
p1 + labs(y=bquote(f(u)), x="u", title="GEV-GP") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### mse for beta estimates
mse_beta <- rep(0, 50)
beta_true <- matrix(c(-0.4,-0.4,-0.5,-0.5,-0.4,0.4,-0.5,-0.5,-0.4,-0.4,-0.2,0.2), 
                    nrow=4, ncol=3)             #coefficient matrix, column is the nonzero coefficient, row the group
group_norm <- apply(beta_true, 1, function(x) sqrt(sum(x^2)))
beta_true <- beta_true/group_norm
ind_all <- readRDS(paste0("../Result/simulation/sim", sim_id, "/ind_selected_all.rds"))

for(id in 1:50){
  p_onesim <- (ncol(res[[id]]$sc)/4)
  beta_true_vec <- NULL
  if(max(ind_all[[id]]) < 4){
    for(gg in 1:4){
      beta_true_vec <- c(beta_true_vec, beta_true[gg,sort(ind_all[[id]])])
    }
  }
  else{
    for(gg in 1:4){
      temp <- rep(0, p_onesim)
      ind_nonzero <- ind_all[[id]]<4
      temp[1:sum(ind_nonzero)] <- beta_true[gg,sort(ind_all[[id]][ind_nonzero])]
      beta_true_vec <- c(beta_true_vec, temp)
    }
  }
  temp <- apply(res[[id]]$sc[((chain_id_final[id]-1)*(nrow(alpha)/2)+1):(chain_id_final[id]*nrow(alpha)/2),], 2, mean)
  mse_beta[id] <- sum((temp-beta_true_vec)^2)  
}
