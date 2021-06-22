#' @title DP EM M-step Function
#'
#' @description  M-step to calculate the coefficients and covariance matrix in EM algorithm for DP
#'
#' @param Ymatrix matrix of dependent variables, multiple dimensions eligible
#' @param Xmatrix matrix of independent variables
#' @param t length of trajectory
#' @param N number of trajectories
#' @param pden.post probablity density from E-step
#' @param VarMatrix covariance matrix
#' @param Beta coefficient matrix
#' @param K number of clusters
#' @param Pcluster weights of clusters
#' @param VP sequence of i.i.d. random variables distributed as Beta(1,alpha), suppose to be infinite for DP, otherwise truncated by K
#' @param alpha hyperparameter with a constant value
#@param split index of columns with no penalty
#@param Lambda penalty coefficient
#' @return list of outputs used for M-step
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr slice n
#' @importFrom matrixcalc is.singular.matrix
# @examples


DP_EM_M_step <- function(Ymatrix,Xmatrix,#split,p,np,Lambda=0,
                   t,N,pden.post,VarMatrix,Beta,K,Pcluster,VP,alpha){
  # require(dplyr)
  # require(matrixcalc)
  Pdenmatrix.post <- c()
  varMatrix <- list()
  VarMatrix.post <- list()
  Pcluster.post <- c()
  VP.post <- c()
  Beta.post <- list()
  # update VP
  for(k in 1:K){

    if(k<K){
      VP_temp <- 0
      for(i in (k+1):K){
        VP_temp = VP_temp + sum(pden.post[[i]])
      }
      if(sum(pden.post[[k]])<=1e-10){
        VP.post[k] = 0
      }
      else {
        VP.post[k] <- sum(pden.post[[k]]) / (sum(pden.post[[k]]) + alpha - 1 + VP_temp)
      }
    }
    if(k==K){
      VP.post[k] = 1
    }
  }
  # use VP to update the weights of clusters
  for(k in 1:K){
    if(k==1){
      Pcluster.post[k] <- VP.post[k]
    }
    if(k>1){
      Pcluster.post[k] <- VP.post[k] * prod((1-VP.post[1:(k-1)]))
    }
    cat("pi_k is ",Pcluster.post[k],"\n")
  }

  # update coefficient matrix and corariance matrix
  for(k in 1:K){

    Pdenmatrix.post <- as.data.frame(pden.post[[k]])
    Pdenmatrix.post <- Pdenmatrix.post %>% slice(rep(1:n(), each = length(t)))
    Pdenmatrix.post <- as.numeric(as.matrix(Pdenmatrix.post))
    Pdenmatrix.post <- diag(Pdenmatrix.post)
    varMatrix[[k]] <- VarMatrix[[k]]
    cat("var_k is ",varMatrix[[k]],"\n")
    Beta.post[[k]] <- Beta[[k]]
    #if(Lambda == 0){
    #if(length(split)==dim(Xmatrix)[2]){
    if(Pcluster.post[k]<=1e-5){
      Beta.post[[k]] <- Beta[[k]]
      VarMatrix.post[[k]] <- varMatrix[[k]]
    }else{
      AYmatrix <- Ymatrix
      AXmatrix <- Xmatrix#[,split]
      # check if feature matrix is singular, if singular, not updating
      TF_singularA <- matrixcalc::is.singular.matrix(t(AXmatrix)%*%Pdenmatrix.post%*%AXmatrix)
      if(TF_singularA==TRUE){
        Beta.post[[k]] <- Beta[[k]]
        VarMatrix.post[[k]] <- varMatrix[[k]]
      }else{
        tempA <- solve(t(AXmatrix)%*%Pdenmatrix.post%*%AXmatrix)%*%t(AXmatrix)%*%Pdenmatrix.post%*%AYmatrix
        Beta.post[[k]] <- tempA
        VarMatrix.post[[k]] <-(1/sum(pden.post[[k]]))*t(Ymatrix - Xmatrix%*%Beta.post[[k]])%*%Pdenmatrix.post%*%(Ymatrix - Xmatrix%*%Beta.post[[k]])
        VarMatrix.post[[k]] <- VarMatrix.post[[k]] / length(t)
      }
    }
    #}

    #}
  }
  return(list(
    "Pcluster.post" = Pcluster.post,
    "VarMatrix.post" = VarMatrix.post,
    "Beta.post" = Beta.post,
    "pden.post" = pden.post,
    "VP.post" = VP.post))

}
