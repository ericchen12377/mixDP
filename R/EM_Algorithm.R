#' @title EM Algorithm Integrated Function
#'
#' @description  Integration of E-step and M-step to iteratively update the parameters
#'
#' @param Ymatrix matrix of dependent variables, multiple dimensions eligible
#' @param Xmatrix matrix of independent variables
#' @param split index of columns with no penalty
#' @param t length of trajectory
#' @param N number of trajectories
#' @param p number of parameters
#' @param np number of groups of parameters
#' @param iter number of maximum iterations
#' @param initial list of initialization values
#' @param lambda penalty coefficient
#' @return list of outputs used for M-step
#' @export
# @examples
# \donttest{
# initial <- list("K" = 2,
# "Pcluster" = c(0.5,0.5),
# "VarMatrix" = list(matrix(c(1,0.5,0.5,1),2,2),matrix(c(1,0.5,0.5,1),2,2)),
# "Beta" = list(matrix(c(100,1,1,1,50,1,1,1),4,2), matrix(c(100,1,1,1,30,1,1,1),4,2)))
# EM_algorithm(Ymatrix=Ymatrix,Xmatrix=Xmatrix, split = c(1,2,3,4),t=c(1:8),N=N,
# p=0, np=0,iter=1000,initial=initial,lambda=0)
# }




EM_algorithm <- function(Ymatrix,Xmatrix,split,t,N,p,np,iter,initial,lambda=0){
  iter <- iter
  #Qloss.vector<-c()
  #Qloss.diff <- c()
  loglikelyhood <- c()
  Iteration <- c()
  for( j in 1:iter){
    if(j==1){
      #time
      t = t
      #dim of samples
      N = N
      #dim of features
      p = p
      #dim of clusters
      K = initial[["K"]]
      #dim of coefficients
      #sigma of clusters
      Beta <- initial[["Beta"]]
      VarMatrix <- initial[["VarMatrix"]]
      #
      #Qobj <- 0
      # Beta <- list()
      # Var <- c()
      # set.seed(initial[["seed"]])
      # for(i in 1:K){
      #   set.seed(initial[["seed"]]+i)
      #   Beta[[i]] = rnorm(np*p+4, mean = initial[["Beta.mean"]], sd = initial[["Beta.sd"]])
      #   Var[i] = runif(1,initial[["Var.lower"]],initial[["Var.upper"]])
      # }
      #membership prob
      Pcluster <- initial[["Pcluster"]]
      e.step <- EM_E_step(Ymatrix,Xmatrix,split,t,N,VarMatrix,Beta,K,Pcluster,Lambda=0)
      #m.step <- M_step(Ymatrix,Xmatrix,split,t,N,p,np,e.step[["pden.post"]],VarMatrix,Beta,K,Lambda=0)
      #cur.Qloss <- m.step[["Qloss"]]
      #Qloss.vector <- m.step[["Qloss"]]

      {cat("\nThe iteration step is\n",j)
        #cat("\nThe current value of objective is\n", m.step[["Qloss"]],"\n")
        cat("\nThe current value of loglikelyhood is\n", e.step[["loglikelyhood"]],"\n")
      }
      cur.loglikelyhood <- e.step[["loglikelyhood"]]
      loglikelyhood.vector <- e.step[["loglikelyhood"]]
      cur.VarMatrix <- VarMatrix
      cur.Beta <- Beta
      cur.Pcluster <- Pcluster
    }
    else{
      #cat(j,'\n')
      m.step <- EM_M_step(Ymatrix,Xmatrix,split,t,N,e.step[["pden.post"]],e.step[["VarMatrix.post"]],e.step[["Beta.post"]],K,Lambda=0)
      # cat("cur.VarMatrix",cur.VarMatrix[[1]],'\n')
      # cat("cur.VarMatrix",cur.VarMatrix[[2]],'\n')
      # print(m.step[["VarMatrix.post"]])
      # cat("m.step_VarMatrix.post",m.step[["VarMatrix.post"]][[1]],'\n')
      # cat("m.step_VarMatrix.post",m.step[["VarMatrix.post"]][[2]],'\n')

      Enorm.diff <-
        sum(abs(sapply(cur.VarMatrix,as.vector) - sapply(m.step[["VarMatrix.post"]],as.vector))^2)+
        sum(abs(sapply(cur.Beta,as.vector) - sapply(m.step[["Beta.post"]],as.vector))^2)
      +sum(abs(cur.Pcluster -  m.step[["Pcluster.post"]])^2)
      cat("\nEnorm.diff is\n",Enorm.diff,'\n')
      if(Enorm.diff<0.1 | min(m.step[["Pcluster.post"]])<1e-3) {
        Iteration <- j
        cat("\nThe final iteration step is\n",Iteration)
        break
      }else{
        cur.VarMatrix <- m.step[["VarMatrix.post"]]
        cur.Beta <- m.step[["Beta.post"]]
        cur.Pcluster <- m.step[["Pcluster.post"]]
      }

      e.step <- EM_E_step(Ymatrix,Xmatrix,split,t,N,m.step[["VarMatrix.post"]],m.step[["Beta.post"]],K,m.step[["Pcluster.post"]],Lambda=0)
      logLL.diff <- abs(cur.loglikelyhood - e.step[["loglikelyhood"]])
      if(logLL.diff<1) {
        Iteration <- j
        cat("\nThe final iteration step is\n",Iteration)
        break
      }else{
        cur.VarMatrix <- m.step[["VarMatrix.post"]]
        cur.Beta <- m.step[["Beta.post"]]
        cur.Pcluster <- m.step[["Pcluster.post"]]
      }
      if(cur.loglikelyhood > e.step[["loglikelyhood"]]){
        cat("loglikelyhood over","\n")
        cur.loglikelyhood <- e.step[["loglikelyhood"]]
        loglikelyhood.vector <- c(loglikelyhood.vector,e.step[["loglikelyhood"]])

        break
      }else{
        cur.loglikelyhood <- e.step[["loglikelyhood"]]
        loglikelyhood.vector <- c(loglikelyhood.vector,e.step[["loglikelyhood"]])
      }

      if(j%%1==0){
        cat("\nThe iteration step is\n",j)
        #cat("\nThe current value of objective is\n", m.step[["Qloss"]],"\n")
        cat("\nThe current value of loglikelyhood is\n", e.step[["loglikelyhood"]],"\n")
      }
    }

  }

  num_para <- (K*(np*p+dim(m.step[["Beta.post"]][[1]])[1]*dim(m.step[["Beta.post"]][[1]])[2]+
                    1+dim(m.step[["VarMatrix.post"]][[1]])[1]*dim(m.step[["VarMatrix.post"]][[1]])[2]))
  BICresult <- -2*loglikelyhood.vector[length(loglikelyhood.vector)] + num_para*log(N*length(t)*dim(m.step[["Beta.post"]][[1]])[2])
  AICresult <- -2*loglikelyhood.vector[length(loglikelyhood.vector)] + 2*num_para
  return(list(#"Qloss.vector" = Qloss.vector,
    "loglikelyhood" = loglikelyhood.vector,
    "Pden.post" = e.step[["pden.post"]],
    "Beta.post" = m.step[["Beta.post"]],
    "Pcluster.post" = m.step[["Pcluster.post"]],
    "VarMatrix.post" = m.step[["VarMatrix.post"]],
    "BIC" =   BICresult,
    "AIC" =   AICresult,
    "Iteration" = Iteration))
}

