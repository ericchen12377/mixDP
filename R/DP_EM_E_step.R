#' @title DP EM E-step Function
#'
#' @description  E-step to calculate the expectation of log-likelihood in EM algorithm for DP
#'
#' @param Ymatrix matrix of dependent variables, multiple dimensions eligible
#' @param Xmatrix matrix of independent variables
#' @param t length of trajectory
#' @param N number of trajectories
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
#' @importFrom matrixStats logSumExp
#' @importFrom mixtools logdmvnorm
#' @importFrom stats na.omit rnorm runif sd


DP_EM_E_step <- function(Ymatrix, Xmatrix, #split,Lambda=0,
                   t, N, VarMatrix, Beta, K, Pcluster, VP, alpha){
  # require(matrixStats)
  # require(mixtools)
  # create intermediate variables
  # mu is the vector of mean value for each row
  mu <- list()
  log_den <- list()
  log_pden <- list()
  varMatrix <- list()
  # for each cluster, run though the loop
  for(k in 1:K){
    # input parameters for the kth cluster
    beta <- Beta[[k]]

    pcluster <- Pcluster[k]

    varMatrix[[k]] <- VarMatrix[[k]]

    # if any of the cluster has a invalid covariance matrix, put the log probability density to -Inf
    if(varMatrix[[k]][1,1]<=1e-3|varMatrix[[k]][2,2]<=1e-3){
      log_pden[[k]] <- rep(-Inf,N)
      next
    }
    # if the weight of cluster is valid, continue updating the parameters
    if(pcluster>=1e-5){
      mu[[k]] <- Xmatrix%*%as.matrix(beta)
      log_den[[k]] <- sapply(1:dim(Ymatrix)[1], function(row) mixtools::logdmvnorm(Ymatrix[row,], mu[[k]][row,],sigma = varMatrix[[k]]))

      dim(log_den[[k]]) <- c(length(t),N)

      log_den[[k]] <- t(log_den[[k]])
      log_pden[[k]]<- log(pcluster)+apply(log_den[[k]],1,sum)
    }else{
      # if the kth cluster has trival proportion, the probability density is trival
      log_pden[[k]] <- rep(-Inf,N)
    }

  }
  # vectorize the log probability density and use logSumExp technique to get the log sum
  vector_log_pden <- sapply(log_pden,as.vector)

  logsum <- apply(vector_log_pden,1,matrixStats::logSumExp)

  pden.post <- list()
  # btemp <- c()
  for(k in 1:K){
    pden.post[[k]] <- exp(log_pden[[k]] - logsum)
    # btemp[k]<- (-Lambda)*sum(Beta[[k]][-split]^2)
  }

  loglikelyhood <- sum(logsum) #+ sum(na.omit(btemp))
  return(list("pden.post" = pden.post,
              "VarMatrix.post" = VarMatrix,
              "Beta.post" = Beta,
              "loglikelyhood" = loglikelyhood))
}
