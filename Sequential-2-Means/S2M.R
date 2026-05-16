###################
#S2M main function#
###################


#This function creates a sequence $\{b_i\}$ of $l$ values of the tuning parameter and implements Sequential-2-Means on the posterior samples to obtain variable selection results corresponding to each value in the sequence $\{b_i\}$. }
#\arguments{
# \item{Beta}{
#   it is an $N$ by $p$ matrix consisting of $N$ posterior samples of $p$ variables
# }
# \item{lower}{
#   the lower bound of the chosen values of the tuning parameter
# }
# \item{upper}{
#   the upper bound of the chosen values of the tuning parameter
# }
# \item{l}{
#   the number of chosen values of the tuning parameter
# }


# \value{
#   \item{N }{the posterior sample size}
#   \item{p }{the total number of variables}
#   
#   \item{H.b.i }{the estimated number of signals corresponding to each $b_i$}
#   \item{b.i }{the values of the tuning parameter}
#   
#   \item{abs.post.median }{medians of the absolute values of the posterior samples of each variable}
#   
S2M=function(Beta,lower,upper,l)
{
  
  N=dim(Beta)[1]
  p=dim(Beta)[2]
  b.i=seq(lower,upper,length=l)
  H.b.i=NULL
  for(r in 1:l)
  {
    KK=NULL
    for (i in 1:N)
    {
      fit=kmeans(abs(Beta[i,]),2)
      cen1=min(fit$centers)
      cen2=max(fit$centers)
      temp1=Beta[i,]

      
      fit1=fit
      
      while(cen2-cen1>b.i[r])
      {
        fit1=fit
        temp=which.min(fit$centers)
        temp1=temp1[which(fit$cluster==temp)]
        if(length(temp1)<=2){break;}
        fit=kmeans(abs(temp1),2)
        cen2=max(fit$centers)
        cen1=min(fit$centers)
      }
      temp=which.min(fit1$centers)
      KK[i]=p-length(which(fit1$cluster==temp))
    }
    H.b.i[r]=as.numeric(names(sort(-table(KK)))[1])
  }
  abs.post.median=NULL
  for(i in 1:p)
  {
    abs.post.median[i]=median(abs(Beta[,i]))
  }
  return(list(N=N,
              p=p,
              H.b.i=H.b.i,
              b.i=b.i,
              abs.post.median=abs.post.median
    ))
}



############################
#make a H_b_i v.s. b_i plot#
############################
#This function is to make the $H_{b_{i}}$-v.s.-$b_i$ plot. The estimated number of signals, or equivalently, the optimal value of the tuning parameter can be decided by looking at this plot.
# \arguments{
#   \item{S2M}{a list obtained from the function \code{S2M}.}
# }
Hbi.vs.bi=function(S2M)
{
  plot(S2M$b.i,S2M$H.b.i)
}




###############################
#variable selection results   #
#given a particular value of H#
###############################
#This function is to print out the variable selection results given the estimated number of signals.}
# \arguments{
#   \item{S2M}{a list obtained from the function \code{S2M}.}
#   \item{H}{the estimated number of signals}
# }
# \value{
#   \item{var.selected}{The indices of selected variables}
#  }

S2M.vs=function(S2M,H)
{

  p=S2M$p
  abs.post.median=S2M$abs.post.median
  var.selected=order(abs.post.median)[p:(p-H+1)]
  return(var.selected)
}


