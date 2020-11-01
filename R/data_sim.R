#' @title Data Simulation Function
#'
#' @description  Define the function to simulate the datasets with multiple dimensions and multivariate normal distributions.
#'
#' @param percetage Percentage of each cluster, sum to 1
#' @param SDlist List of standard deviation matrices for each cluster
#' @param t Length of trajectory, t=1 will be single point;
#' t can be used to specify the length of trajectory for training/validation/testing
#' @param Tmatrix Feature matrix for the trajectory.
#' @param N Number of trajectories
#' @param Coeflist List of coefficient matrices for each cluster
#' @param seed Set the random seed for simulation
#' @param random Specify whether the simulation data has exact percentage
#' or approximate percentage for each cluster
#' @return Simulation dataset
#' @export
#' @importFrom mvtnorm rmvnorm
#' @examples
#' \donttest{
#' ######single dimension simulation######
#' t <- c(1:8)
#' percentage <- c(0.7, 0.3)
#' N <- 200
#'
#' # generate the feature matrix from polynomial formula
#' T <- c(1:8)
#' Tmatrix <- cbind(1, T, T^2, T^3)
#'
#' # generate list of standard deviation matrices
#' SDlist = c(1, 1)
#'
#' # generate list of coefficient matrices
#' # cluster 1
#' A1 <- c(100, 1, 1, 1)
#' # cluster 2
#' A2 <- c(80, 1, 1, 1)
#' Coeflist <- list(A1, A2)
#'
#' # simulate the dataset
#' data_sim(percentage, SDlist, t, Tmatrix, N, Coeflist, seed=10, random = TRUE)
#'
#' ######multiple dimensions simulation######
#' t <- c(1:8)
#' percentage <- c(0.7, 0.3)
#' N <- 200
#'
#' # generate the feature matrix from polynomial formula
#' T <- c(1:8)
#' Tmatrix <- cbind(1, T, T^2, T^3)
#'
#' # generate list of standard deviation matrices
#' error_Sigma1 <- matrix(c(1,0,0,1),2,2)
#' error_Sigma2 <- matrix(c(1,0,0,1),2,2)
#' SDlist <- list(error_Sigma1,error_Sigma2)
#'

#' # generate list of coefficient matrices
#' # cluster 1
#' A1 <- c(100,1,1,1)
#' B1 <- c(40,1,1,1)
#' # cluster 2
#' A2 <- c(80,1,1,1)
#' B2 <- c(30,1,1,1)
#' Coeflist <- list(cbind(A1, B1), cbind(A2, B2))
#'
#' # simulate the dataset
#' data_sim(percentage, SDlist, t, Tmatrix, N, Coeflist, seed=10, random = TRUE)
#' }









data_sim <- function(percentage, SDlist, t, Tmatrix, N, Coeflist, seed, random = TRUE){
  coefdim <- ncol(Tmatrix)
  if(random == TRUE){
    print("Getting approximate percentages of each cluster!")
    #define the function to get random data label
    get_label <- function(percentage, rand){
      cumsump <- cumsum(percentage)
      i <- 1
      while(rand <= 1){
        if(rand <= cumsump[i]){
          break
        }
        i <- i + 1
      }
      return(i)
    }
    tempX = matrix(rep(t(Tmatrix[t, ]), N) , ncol =  coefdim, byrow = TRUE)
    tempY <- c()
    label <- c()
    set.seed(seed)
    Ulist <- runif(N, 0, 1)
    if(is.list(SDlist)){
      Ydim <- ncol(SDlist[[1]])
      for(i in 1:N){
        label[i] <- get_label(percentage, Ulist[i])
        tempY <- rbind(tempY,t(apply(Tmatrix %*% Coeflist[[label[i]]], 1,
                                     function(m) mvtnorm::rmvnorm(1, mean = m, sigma = SDlist[[label[i]]])
                                     )
                               )
                       )
      }
      label <- rep(label, each = length(t))
      XYmatrix <- cbind(tempY, tempX, label)
      colnames(XYmatrix)[1:Ydim] <- paste("Y", c(1:Ydim), sep = "")
      colnames(XYmatrix)[coefdim + Ydim + 1] <- c("label")
    } else {
      for(i in 1:N){
        label[i] <- get_label(percentage, Ulist[i])
        tempY <- rbind(tempY, matrix(apply(Tmatrix %*% cbind(Coeflist[[label[i]]]), 1,
                                           function(m) rnorm(1, mean = m, sd = SDlist[label[i]])
        )
        )
        )
      }
      label <- rep(label, each = length(t))
      XYmatrix <- cbind(tempY, tempX, label)
      colnames(XYmatrix)[1] <- paste("Y", c(), sep = "")
      colnames(XYmatrix)[coefdim + 1 + 1] <- c("label")
    }
  }
  else if(random == FALSE){
    print("Getting exact percentages of each cluster!")
    XYmatrix <- c()
    if(is.list(SDlist)){
      Ydim <- ncol(SDlist[[1]])
      for(i in 1:length(percentage)){
        tempX <- matrix(rep(t(Tmatrix[t, ]), percentage[i] * N) , ncol =  coefdim, byrow = TRUE)
        tempY <- t(apply(tempX %*% Coeflist[[i]], 1,
                         function(m) mvtnorm::rmvnorm(1, mean = m, sigma = SDlist[[i]])))
        XYmatrix <- rbind(XYmatrix, cbind(tempY, tempX, i))
      }
      label <- XYmatrix[, coefdim + Ydim + 1]
      colnames(XYmatrix)[1:Ydim] <- paste("Y", c(1:Ydim), sep = "")
      colnames(XYmatrix)[coefdim + Ydim + 1] <- c("label")
    } else {
      for(i in 1:length(percentage)){
        tempX <- matrix(rep(t(Tmatrix[t, ]), percentage[i] * N) , ncol =  coefdim, byrow = TRUE)
        tempY <- apply(tempX %*% cbind(Coeflist[[i]]), 1, function(m) rnorm(1, mean = m, sd = sd[i]))
        XYmatrix <- rbind(XYmatrix, cbind(tempY, tempX, i))
      }
      label <- XYmatrix[, coefdim + 1 + 1]
      colnames(XYmatrix)[1] <- c("Y")
      colnames(XYmatrix)[coefdim + 1 + 1] <- c("label")
    }

    #randomly shuffle the data
    set.seed(seed)
    id <- rep(sample(c(1:N)), each = length(t))
    XYmatrix <- cbind(id, XYmatrix)
    XYmatrix <- XYmatrix[order(id), -1]
  }

  return(list("Data" = XYmatrix,
              "Ground Truth" = list("Sample size" = N,
                                    "Number of clusters" = length(percentage),
                                    "Percentage" = table(label) / (N * length(t)),
                                    "Clusters" = table(label))
             )
        )
}


