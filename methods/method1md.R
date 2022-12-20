#-------------------------------------------------------------------------------   
#          Method 1: Matrix decomposition - "md"
#-------------------------------------------------------------------------------
# File:    method1md.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - This method requires a grfmodel of modeltype "covar".
#  - The simulation is split up into the precomputations (cholesky factor of
#    covariance matrix) using the function Precompute_md() and the
#    evaluation part in the main function.
#  - For multiple simulation of the same GRF and for the same points, the
#    Precomputation part has to be done only once and can be provided as
#    argument.
#
# DEPENDENCIES
#  - TrueCovMat() from "covariance.R"
#  
# CONTENTS
#   .SimGRF_md()       Main simulation function
#   Precompute_md()   Precomputations (cholesky factor of cov. matrix)
#-------------------------------------------------------------------------------

# Sources
source("R/simGRFS/models/covariance.R")



.SimGRF_md <- function(dummyvar, pts, grfmodel, precomputations = NULL){
  #-----------------------------------------------------------------------------
  # Simulates values of a GRF on the sphere given the points and a grfmodel.
  #
  # Precomputation option: The cholesky factor of the covariance matrix can
  # be precomputated and provided as an argument.
  #
  # DEPENDENCIES
  #  - Precompute_md()
  #
  # INPUT
  #        dummyvar:  dummy variable for parallel computing
  #             pts:  (3xN) matrix of points in columns format c(x,y,z)
  #        grfmodel:  instance of class GRFModel with modeltype "covar"
  # precomputations:  precomputed cholesky factor of covariance matrix
  #
  # OUTPUT
  #   vector (N) with simulated values in the same order as input
  #-----------------------------------------------------------------------------
  N <- length(pts) / 3
  
  # compute cholesky decomposition of covariance matrix if not provided
  if (is.null(precomputations)) {
    cholesky_factor <- Precompute_md(pts = pts, grfmodel = grfmodel)
  } else cholesky_factor <- precomputations
  
  # main computation
  values <- cholesky_factor %*% rnorm(N)
  
  return(as.vector(values))
}



Precompute_md <- function(pts, grfmodel) {
  #-----------------------------------------------------------------------------
  # Precomputation for matrix decomposition method.
  # Computes the cholesky factor of the true covariance matrix of given
  # points using the provided grfmodel.
  #
  # DEPENDENCIES
  #  - TrueCovMat() from "covariance.R"
  #
  # INPUT
  #      pts:  (3xN) matrix of points in columns format c(x,y,z)
  # grfmodel:  instance of class GRFmodel with modeltype "covar"
  #
  # OUTPUT
  #   (NxN) matrix (cholesky factor) of covariance matrix
  #-----------------------------------------------------------------------------
  
  # compute covariance matrix
  covmat <- TrueCovMat(pts = pts, grfmodel = grfmodel, variance = 1)

  # return cholesky factor
 return(t(chol(covmat)))
}