#' @title EM M-step Function
#'
#' @description  M-step to calculate the coefficients and covariance matrix in EM algorithm
#'
#' @param Ymatrix matrix of dependent variables, multiple dimensions eligible
#' @param Xmatrix matrix of independent variables
#' @param split index of columns with no penalty
#' @param t length of trajectory
#' @param N number of trajectories
#' @param pden.post probablity density from E-step
#' @param VarMatrix covariance matrix
#' @param Beta coefficient matrix
#' @param K number of clusters
#' @param Lambda penalty coefficient
#' @return list of outputs used for M-step
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr slice n
#' @importFrom matrixcalc is.singular.matrix
# @examples
# \donttest{
# EM_M_step(Ymatrix, Xmatrix, split=c(1,2,3,4), t=c(1:8), N=N, pden.post = EM_E$pden.post,
# VarMatrix=EM_E$VarMatrix.post,
# Beta=EM_E$Beta.post, K=2, Lambda=0)
# }




EM_M_step <- function(Ymatrix, Xmatrix, split, t, N, pden.post, VarMatrix, Beta, K, Lambda=0){
  # require(dplyr)
  # require(matrixcalc)
  Pdenmatrix.post <- c()
  varMatrix <- list()
  VarMatrix.post <- list()
  Pcluster.post <- c()
  Beta.post <- list()

  for(k in 1:K){

    Pdenmatrix.post <- as.data.frame(pden.post[[k]])
    Pdenmatrix.post <- Pdenmatrix.post %>% slice(rep(1:n(), each=length(t)))
    Pdenmatrix.post <- as.numeric(as.matrix(Pdenmatrix.post))
    Pdenmatrix.post <- diag(Pdenmatrix.post)

    varMatrix[[k]] <- VarMatrix[[k]]
    cat("var_k is ",varMatrix[[k]],"\n")

    #cat("pden.post is ",pden.post[[k]],"\n")
    Pcluster.post[k] <- sum(pden.post[[k]])/N
    cat("pi_k is ",Pcluster.post[k],"\n")

    Beta.post[[k]] <- Beta[[k]]

    if(Lambda == 0){
      if(length(split)==dim(Xmatrix)[2]){
        if(Pcluster.post[k]<=1e-3){

          Beta.post[[k]] <- Beta[[k]]
          #Beta.post[[k]] <- rep(0,np*p+4)

          VarMatrix.post[[k]] <- varMatrix[[k]]

        }else{
          AYmatrix <- Ymatrix
          AXmatrix <- Xmatrix[,split]
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
      }

    }
  }

  return(list(
    "Pcluster.post" = Pcluster.post,
    "VarMatrix.post" = VarMatrix.post,
    "Beta.post" = Beta.post,
    "pden.post" = pden.post))
}
