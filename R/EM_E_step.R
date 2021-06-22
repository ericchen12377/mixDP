#' @title EM E-step Function
#'
#' @description  E-step to calculate the expectation of log-likelihood in EM algorithm
#'
#' @param Ymatrix matrix of dependent variables, multiple dimensions eligible
#' @param Xmatrix matrix of independent variables
#' @param split index of columns with no penalty
#' @param t length of trajectory
#' @param N number of trajectories
#' @param VarMatrix covariance matrix
#' @param Beta coefficient matrix
#' @param K number of clusters
#' @param Pcluster weights of subgroups
#' @param Lambda penalty coefficient
#' @return list of outputs used for M-step
#' @export
#' @importFrom matrixStats logSumExp
#' @importFrom mixtools logdmvnorm
#' @importFrom stats na.omit rnorm runif sd
# @examples
# \donttest{
# initial <- list("K" = 2,
#                 "Pcluster" = c(0.5,0.5),
#                 "VarMatrix" = list(matrix(c(1,0.5,0.5,1),2,2),
#                                    matrix(c(1,0.5,0.5,1),2,2)),
#                 "Beta" = list(matrix(c(100,1,1,1,50,1,1,1),4,2),
#                               matrix(c(100,1,1,1,30,1,1,1),4,2)))
#
#
# EM_E = EM_E_step(Ymatrix, Xmatrix, split=c(1,2,3,4), t=c(1:8),
# N=N, VarMatrix=initial[['VarMatrix']],
# Beta=initial[['Beta']], K=2,
# Pcluster=initial[['Pcluster']], Lambda = 0)
# }

EM_E_step <- function(Ymatrix, Xmatrix, split, t, N, VarMatrix, Beta, K, Pcluster, Lambda = 0){
  # require(matrixStats)
  # require(mixtools)
  #create intermediate variables
  # mu is the vector of mean value for each row
  mu <- list()
  log_den <- list()
  log_pden <- list()
  varMatrix <- list()
  #for each cluster, run though the loop
  for(k in 1:K){
    # input parameters for the kth cluster
    beta <- Beta[[k]]
    pcluster <- Pcluster[k]
    varMatrix[[k]] <- VarMatrix[[k]]

    #check if covariance matrix has full rank. If not, the probability density is trival.
    if(any(diag(varMatrix[[k]])<= 1e-3)){
      log_pden[[k]] <- rep(-Inf, N)
      next
    }
    #check if the kth cluster has non-trival proportion of data
    if(pcluster >= 1e-3){
      mu[[k]] <- Xmatrix %*% as.matrix(beta)
      log_den[[k]] <- sapply(1:dim(Ymatrix)[1], function(row) mixtools::logdmvnorm(Ymatrix[row, ], mu[[k]][row, ],
                                                                         sigma = varMatrix[[k]]
                                                                         )
                             )
      #convert into matrix
      dim(log_den[[k]]) <- c(length(t), N)
      log_den[[k]] <- t(log_den[[k]])
      log_pden[[k]]<- log(pcluster) + apply(log_den[[k]], 1, sum)

    }
    else{
      #if the kth cluster has trival proportion, the probability density is trival
      log_pden[[k]] <- rep(-Inf, N)
    }
  }
  #vectorize the log probability density and use logSumExp technique to get the log sum
  vector_log_pden <- sapply(log_pden, as.vector)
  logsum <- apply(vector_log_pden, 1, matrixStats::logSumExp)


  #output probability density for next iteration
  pden.post <- list()
  btemp <- c()
  for(k in 1:K){
    pden.post[[k]] <- exp(log_pden[[k]] - logsum)
    btemp[k] <- (-Lambda) * sum(Beta[[k]][-split]^2)
  }

  #log likelyhood is used to check monotonous
  loglikelyhood <- sum(logsum) + sum(na.omit(btemp))
  return(list("pden.post" = pden.post,
              "VarMatrix.post" = VarMatrix,
              "Beta.post" = Beta,
              "loglikelyhood" = loglikelyhood))
}
