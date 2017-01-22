library(flare)
library(igraph)
library(lpSolve)
library(mnormt)
library(matrixcalc)
library(pracma)
library(MASS)
library(doParallel)
library(foreach)
library(parallel)
library(abind)
library(compiler)
enableJIT(3)


weight_epanechnikovker <- function(t, m, n, b)
{
  v = (t-m/n)/b
  if(abs(v) > 1)
  {
    return(0)
  }
  else
  {
    #Epanechnikov
    numerator = 0.75 * (1 - v^2)
    vi=sapply(1:n,function(x){(t-x/n)/b})
    indx=abs(vi)<=1
    denumerator=sum(0.75*(1-vi[indx]^2))
    return(numerator/denumerator)
  }
}
vc_weight_epanechnikovker<-Vectorize(weight_epanechnikovker,vectorize.args='m')#vectorization of a scalar function
vc_weight_epanechnikovker=cmpfun(vc_weight_epanechnikovker)


#' fit a TV_VAR model
#' @param it Specify the time point we are researching
#' @param X T*d array which stores the data of the time series. T is the length and d is the dimension.
#' @param tau penalty parameter
#' @param b bandwidth
#' @return Estimation of transition matrix at time t
tvvar_fit<-function(it,X,tau,b){
  T=dim(X)[1]
  d=dim(X)[2]
  #compute the kernel weights for t_i
  ker_weight=vc_weight_epanechnikovker(it/T,1:T,T,b)
  #compute the kernel weights for t_{i-1}
  ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
  #sigma_{i-1,0}
  s10=crossprod(X*drop(ker_weight_1),X)
  #sigma_{i,0}
  s00=crossprod(X*drop(ker_weight),X)
  #sigma_{i-1,1}
  s11=crossprod(X[1:(T-1),]*drop(ker_weight_1[1:(T-1)]),X[2:T,])
  #sigma_{i,1}
  s01=crossprod(X[1:(T-1),]*drop(ker_weight[1:(T-1)]),X[2:T,])
  #sigma_{i-1,-1}
  s1_1=crossprod(X[2:T,]*drop(ker_weight_1[2:T]),X[1:(T-1),])
  #sigma_{i,-1}
  s0_1=crossprod(X[2:T,]*drop(ker_weight[2:T]),X[1:(T-1),])
  #solve linear programming
  f.obj=rep(1,2*d)
  # f.con1=cbind(-s10,s10)
  # f.con2=cbind(s00,-s10)
  # f.con=rbind(f.con1,f.con2)
  f.con=matrix(rep(0,4*d^2),nrow=2*d)
  f.con[1:d,1:d]=-s10
  f.con[1:d,(d+1):(2*d)]=s10
  f.con[(d+1):(2*d),1:d]=s00
  f.con[(d+1):(2*d),(d+1):(2*d)]=-s10
  f.dir=rep("<=",2*d)
  f.rhs=matrix(0,2*d,d)
  A_hat_tmp=sapply(1:d,FUN=function(x){
    f.rhs[,x]=matrix(c(tau-pmax(s11[,x],s1_1[x,]),
                       tau+pmin(s01[,x],s0_1[x,])),2*d,1)
    lp.out=lp('min',f.obj,f.con,f.dir,f.rhs[,x])$solution
    return(matrix(lp.out[1:d]-lp.out[(1+d):(d+d)],d,1))
  })
  A_hat=t(A_hat_tmp)
  # A_hat=(A_hat_tmp)
  return(A_hat)
}


#' fit a Ridge_VAR model based on LS with l2 penalty
#' @param it Specify the time point we are researching
#' @param X T*d array which stores the data of the time series. T is the length and d is the dimension.
#' @param lam penalty parameter
#' @param b bandwidth
#' @return Estimation of transition matrix at time t
lsvar_fit<-function(it,X,lam,b){
  T=dim(X)[1]
  d=dim(X)[2]
  #compute the kernel weights for t_i
  ker_weight=vc_weight_epanechnikovker(it/T,1:T,T,b)
  # #compute the kernel weights for t_{i-1}
  # ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
  W1=crossprod(X[2:T,]*drop(ker_weight[2:T]),X[1:(T-1),])
  W2=crossprod(X[1:(T-1),]*drop(ker_weight[2:T]),X[1:(T-1),])
  # A_hat_LS=solve(W2+lam*diag(d))%*%W1
  # A_hat_LS=crossprod(chol2inv(chol(W2+lam*diag(d))),W1)#W2 is symmetric
  A_hat_LS=tcrossprod(W1,chol2inv(chol(W2+lam*diag(d))))#W2 is symmetric
  return(A_hat_LS)
}


#' fit a MLE_VAR model based on MLE
#' @param it Specify the time point we are researching
#' @param X T*d array which stores the data of the time series. T is the length and d is the dimension.
#' @param b bandwidth
#' @return Estimation of transition matrix at time t
mlevar_fit<-function(it,X,b){
  T=dim(X)[1]
  d=dim(X)[2]
  #compute the kernel weights for t_i
  ker_weight=vc_weight_epanechnikovker(it/T,1:T,T,b)
  # #compute the kernel weights for t_{i-1}
  # ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
  W1=crossprod(X[2:T,]*drop(ker_weight[2:T]),X[1:(T-1),])
  W2=crossprod(X[1:(T-1),]*drop(ker_weight[2:T]),X[1:(T-1),])
  A_hat_MLE=tcrossprod(W1,chol2inv(chol(W2)))
  # A_hat_MLE=tcrossprod(W1,solve(W2))
  # A_hat_MLE=crossprod(chol2inv(chol(W2)),W1)
  # A_hat_MLE=W1%*%chol2inv(chol(W2))
  return(A_hat_MLE)
}



abso<-function(x){
  x[x<0]=0
  return(x)
}
#' fit a Lasso_VAR model based on LS with l1 penalty
#' @param it Specify the time point we are researching
#' @param X T*d array which stores the data of the time series. T is the length and d is the dimension.
#' @param lam shrinkage parameter
#' @param b bandwidth
#' @param alpha stepsize
#' @param A_old initial matrix 1
#' @param A_oldold initial matrix 2
#' @return Estimation of transition matrix at time t
lavar_fit<-function(it,X,lam,b,alpha,A_old,A_oldold){
  # alpha=1/max(svd(crossprod(X,X))$d)
  T=dim(X)[1]
  d=dim(X)[2]
  #compute the kernel weights for t_i
  ker_weight=vc_weight_epanechnikovker(it/T,1:T,T,b)
  # #compute the kernel weights for t_{i-1}
  # ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
  # W1=crossprod(X[2:T,]*drop(ker_weight[1:(T-1)]),X[1:(T-1),])
  # W2=crossprod(X[1:(T-1),]*drop(ker_weight[1:(T-1)]),X[1:(T-1),])
  W1=crossprod(X[2:T,]*drop(ker_weight[2:T]),X[1:(T-1),])
  W2=crossprod(X[1:(T-1),]*drop(ker_weight[2:T]),X[1:(T-1),])
  dis=1
  step=2
  alpha=1/max(svd(2*W2)$d)
  # alpha=1/max(svd(crossprod(X,X))$d)
  while (dis>10^(-3)){
    print(c(step,dis))
    y=A_old+(step-1)/(step+2)*(A_old-A_oldold)
    tmp=y+2*alpha*(W1-tcrossprod(y,W2))
    # tmp=y+2*alpha*(W1-crossprod(W2,y))
    A_new=abso(abs(tmp)-alpha*lam)*sign(tmp)
    # dis=sum((A_new-A_old)^2)
    dis=norm(A_new-A_old,type="F")
    # dis=norm(A_new-y,type="F")
    # isallzero=(sum(A_new-A_old)==d*d)
    A_oldold=A_old
    A_old=A_new
    step=step+1
  }
  return(A_new)
}
