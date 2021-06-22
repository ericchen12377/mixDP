#' @title DP EM Algorithm Integrated Function
#'
#' @description  Integration of E-step and M-step to iteratively update the parameters for DP
#'
#' @param Ymatrix matrix of dependent variables, multiple dimensions eligible
#' @param Xmatrix matrix of independent variables
#' @param t length of trajectory
#' @param N number of trajectories
#' @param iter number of maximum iterations
#' @param initial list of initialization values
#' @return list of outputs used for M-step
#' @export
# @examples

DP_EM_algorithm <- function(Ymatrix,Xmatrix,
                         # split,,p,np,lambda=0
                         t,N, iter,initial){
  iter <- iter
  loglikelyhood <- c()
  Iteration <- c()
  Pcluster_matrix <- c()
  alpha = initial[["alpha"]]
  for( j in 1:iter){
    # if first step, initiate the parameters and go to e-step
    if(j==1){
      #time
      t = t
      #dim of samples
      N = N
      #dim of features
      #p = p
      #dim of clusters
      K = initial[["K"]]
      #dim of coefficients
      #sigma of clusters
      Beta <- initial[["Beta"]]
      VarMatrix <- initial[["VarMatrix"]]
      #
      #membership prob
      Pcluster <- c()
      VP <- initial[["VP"]]
      for(k in 1:K){
        if(k==1){
          Pcluster[k] <- VP[k]
        }
        if(k>1){
          Pcluster[k] <- VP[k] * prod((1-VP[1:(k-1)]))

        }
        cat("pi_k is ",Pcluster[k],"\n")
      }

      ##################################

      e.step <- DP_EM_E_step(Ymatrix,Xmatrix,
                       # split,Lambda=0,
                       t,N,VarMatrix,Beta,K,Pcluster,VP,alpha)
      {cat("\nThe iteration step is\n",j)
        cat("\nThe current value of loglikelyhood is\n", e.step[["loglikelyhood"]],"\n")
      }
      cur.loglikelyhood <- e.step[["loglikelyhood"]]
      loglikelyhood.vector <- e.step[["loglikelyhood"]]
      cur.VarMatrix <- VarMatrix
      cur.Beta <- Beta
      cur.VP <- VP
      Pcluster_matrix <- rbind(Pcluster_matrix,Pcluster)
    }
    # if not first step, use e-step outputs to update m-step
    else{
      #cat(j,'\n')
      alpha = initial[["alpha"]]
      m.step <- DP_EM_M_step(Ymatrix,Xmatrix,
                       # split,p,np,Lambda=0,
                       t,N,e.step[["pden.post"]],e.step[["VarMatrix.post"]],e.step[["Beta.post"]],K,VP,alpha=alpha)

      # Differences of coefficients and covariances as stopping criteria, stop if converges
      Enorm.diff <-
        sum(abs(sapply(cur.VarMatrix,as.vector) - sapply(m.step[["VarMatrix.post"]],as.vector))^2)+
        sum(abs(sapply(cur.Beta,as.vector) - sapply(m.step[["Beta.post"]],as.vector))^2)
      +sum(abs(cur.VP -  m.step[["VP.post"]])^2)
      cat("\nEnorm.diff is\n",Enorm.diff,'\n')
      # threshold need to be set and tuned
      if(Enorm.diff<0.1) {
        Iteration <- j
        cat("\nThe final iteration step is\n",Iteration)
        break
      }else{
        cur.VarMatrix <- m.step[["VarMatrix.post"]]
        cur.Beta <- m.step[["Beta.post"]]
        cur.Pcluster <- m.step[["Pcluster.post"]]
        Pcluster_matrix <- rbind(Pcluster_matrix,cur.Pcluster)
        cur.VP <- m.step[["VP.post"]]
      }

      # use m-step outputs to update e-step
      e.step <- DP_EM_E_step(Ymatrix,Xmatrix,
                       # split,Lambda=0,
                       t,N,m.step[["VarMatrix.post"]],m.step[["Beta.post"]],K,m.step[["Pcluster.post"]],m.step[["VP.post"]],alpha)

      # Difference of loglikelihood as the stopping criteria, stop if converges
      logLL.diff <- abs(cur.loglikelyhood - e.step[["loglikelyhood"]])

      # threshold need to be set and tuned
      if(logLL.diff<1) {
        Iteration <- j
        cat("\nThe final iteration step is\n",Iteration)
        break
      }else{
        cur.VarMatrix <- m.step[["VarMatrix.post"]]
        cur.Beta <- m.step[["Beta.post"]]
        cur.Pcluster <- m.step[["Pcluster.post"]]
        Pcluster_matrix <- rbind(Pcluster_matrix,cur.Pcluster)
        cur.VP <- m.step[["VP.post"]]

      }
      if(cur.loglikelyhood > e.step[["loglikelyhood"]]){
        cat("loglikelyhood over","\n")
        cur.loglikelyhood <- e.step[["loglikelyhood"]]
        loglikelyhood.vector <- c(loglikelyhood.vector,e.step[["loglikelyhood"]])

        #break
      }else{
        cur.loglikelyhood <- e.step[["loglikelyhood"]]
        loglikelyhood.vector <- c(loglikelyhood.vector,e.step[["loglikelyhood"]])
      }

      if(j%%1==0){
        cat("\nThe iteration step is\n",j)
        cat("\nThe current value of loglikelyhood is\n", e.step[["loglikelyhood"]],"\n")
      }
    }

  }

  num_cluster <- sum(m.step[["Pcluster.post"]]>1e-3)
  num_para <- (num_cluster*(dim(m.step[["Beta.post"]][[1]])[1]*dim(m.step[["Beta.post"]][[1]])[2]+ #np*p+
                              1+dim(m.step[["VarMatrix.post"]][[1]])[1]*dim(m.step[["VarMatrix.post"]][[1]])[2]))
  BICresult <- -2*loglikelyhood.vector[length(loglikelyhood.vector)] + num_para*log(N*length(t)*dim(m.step[["Beta.post"]][[1]])[2])
  AICresult <- -2*loglikelyhood.vector[length(loglikelyhood.vector)] + 2*num_para
  return(list(
    "loglikelyhood" = loglikelyhood.vector,
    "Pden.post" = e.step[["pden.post"]],
    "Beta.post" = m.step[["Beta.post"]],
    "Pcluster.post" = m.step[["Pcluster.post"]],
    "Pcluster_matrix" = Pcluster_matrix,
    "VP.post" = m.step[["VP.post"]],
    "VarMatrix.post" = m.step[["VarMatrix.post"]],
    "BIC" =   BICresult,
    "AIC" =   AICresult,
    "Iteration" = Iteration))
}


