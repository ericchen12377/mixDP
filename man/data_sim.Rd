% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_sim.R
\name{data_sim}
\alias{data_sim}
\title{Data Simulation Function}
\usage{
data_sim(percentage, SDlist, t, Tmatrix, N, Coeflist, seed, random = TRUE)
}
\arguments{
\item{percentage}{Percentage of each cluster, sum to 1}

\item{SDlist}{List of standard deviation matrices for each cluster}

\item{t}{Length of trajectory, t=1 will be single point;
t can be used to specify the length of trajectory for training/validation/testing}

\item{Tmatrix}{Feature matrix for the trajectory.}

\item{N}{Number of trajectories}

\item{Coeflist}{List of coefficient matrices for each cluster}

\item{seed}{Set the random seed for simulation}

\item{random}{Specify whether the simulation data has exact percentage
or approximate percentage for each cluster}
}
\value{
Simulation dataset
}
\description{
Define the function to simulate the datasets with multiple dimensions and multivariate normal distributions.

percentage, SDlist, t, Tmatrix, N, Coeflist, seed=10, random = TRUE
}
