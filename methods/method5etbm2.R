#-------------------------------------------------------------------------------   
#          Method 5: Emery TBM, alternative Implementation - "etbm2"
#------------------------------------------------------------------------------- 
# File:    method5etbm2.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - This method requires a grfmodel of modeltype "spectral".
#  - The A^(d) matrices can be precomputed and passed as an argument.
#  - This method has significant (co)variance deficits.
#
# DEPENDENCIES
#  - UniformSpherePoints() from "helpfunctions.R"
#  - Legendre() from "legendre.R"
#  - package "gsl" for fast performance
#
# CONTENT
#   .SimGRF_etbm2()         Main simulation function
#   Precompute_etbm()       Precomputations for ETBM2
#-------------------------------------------------------------------------------

# Sources
gsl.installed <- require(gsl)
if (!gsl.installed) {
  source("R/simGRFS/helpers/no_gsl/legendre_no_gsl.R")
} else  {
  source("R/simGRFS/helpers/legendre.R")
}
source("R/simGRFS/helpers/helpfunctions.R")



.SimGRF_etbm2 <- function(dummyvar,
                          pts,
                          grfmodel,
                          L = 1,
                          Q = 30,
                          precomputations = NULL,
                          disp_progress = FALSE,
                          components_out = FALSE) {
  #-----------------------------------------------------------------------------
  # Simulates values of a GRF on the sphere given the points and a grfmodel.
  #
  # The matrices A^(d) can be precomputed and passed as an argument.
  #
  # DEPENDENCIES
  # - UniformSpherePoints() from "helpfunctions.R"
  # - Legendre() from "legendre.R"
  #
  # INPUT
  #        dummyvar:  dummy variable for parallel computing
  #             pts:  (3xN) matrix of points to simulate in format c(x,y,z)
  #        grfmodel:  instance of class GRFmodel, specifying covariance model
  #               L:  number of independent realizations to use CLT
  #               Q:  number of frequencies to use 
  # precomputations:  matrices A^(d) for all (d) up to D
  # disp_progress:  set TRUE for progress reports
  #  components_out:  set TRUE to return matrix of KL-components instead
  #
  # OUTPUT
  #   vector (N) with simulated values in the same order as input, or
  #   matrix ((D+1)xN) with simulated KL-components for each point
  #-----------------------------------------------------------------------------
  N <- dim(pts)[2]
  
  # get model parameters
  D <- grfmodel$D
  if (!components_out) coefs <- grfmodel$coefs
  else coefs <- rep(1, D + 1)
  
  # prepare progress reports
  if (disp_progress) progresslist <- floor(seq(0, 100, length = D + 1))
  
  # prepare values matrix (columns = points, rows = d's)
  ysim <- matrix(numeric((D + 1) * N), D + 1, N) 
  
  #------------------ MAIN COMPUTATION: with precomputations -------------------
  if (!is.null(precomputations)) {
    
    # get evenly distributed directions on unit sphere
    directions <- UniformSpherePoints(Q)
    
    # loop over all KL-components
    for (d in 0:D) {
      
      # set lambda parameter
      if (d == 0) lambda <- 1 / 2
      else lambda <- d / 2
      
      for (l in 1:L) {
        
        # random frequencies in R^3 distributed as g (chi-distribution)
        norms <- sqrt(2 * lambda * rchisq(Q, df = d + 3))
        freqs <- directions * rep(norms, each = 3)

        # random phase
        phi <- runif(1, min = 0, max = 2 * pi)
        
        # compute vector B depending on (precomputed) matrix A
        U <- rnorm(Q)
        B <- precomputations[ , ,d + 1] %*% U
        
        # main computation
        ysim[d + 1, ] <- ysim[d + 1, ] + cos(t(pts) %*% freqs + phi) %*% B
      }
      
      # correct for 1/sqrt(L) and sqrt(a_n) spectral measure
      ysim[d + 1, ] <- sqrt(coefs[d + 1] / L) * ysim[d + 1, ] / Q
      
      # report progress
      if (disp_progress) cat("\r Simulation progress: ",
                             progresslist[d + 1], "%")
    }
  }
  
  #------------------ MAIN COMPUTATION: without precomputations ----------------
  else {
    
    # prepare frequency directions
    directions <- UniformSpherePoints(Q)
    cosangles <- t(directions) %*% directions
    diag(cosangles) <- rep(1, Q)
    
    # loop over all Legendre components 
    for (d in 0:D) {
      
      # parameter
      if (d == 0) lambda <- 1 / 2
      else lambda <- d / 2
      
      # compute ratio of densities
      ratio <- exp(lambda) * sqrt(8 / (pi * lambda**d)) * gamma((d + 3) / 2)
      
      # get matrix of legendre polynomials evaluations at angles
      legendrematrix <- Legendre(d, cosangles)
      
      # get eigenvectors and values
      eigen_both <- eigen(legendrematrix)
      V <- eigen_both$vectors
      ev <- eigen_both$values
      ev[ev < 1e-14] <- 0      # handle numerical errors
      
      # get matrix A
      A <-  (V * rep(ratio * sqrt(ev), each = Q))
      
      # loop over L (L in paper)
      for (l in 1:L) {
        
        # random frequencies in R^3 
        norms <- sqrt(2 * lambda * rchisq(Q, df = d + 3))
        freqs <- directions * rep(norms, each = 3)
        
        # random phase
        phi <- runif(1, min = 0, max = 2 * pi)
  
        # compute vector B depending on matrix A
        U <- rnorm(Q)
        B <- A %*% U
        
        # main computation
        ysim[d + 1, ] <- ysim[d + 1, ] + cos(t(pts) %*% freqs + phi) %*% B
      }

      # correct for 1/sqrt(L) and sqrt(a_n) spectral measure
      ysim[d + 1, ] <- sqrt(coefs[d + 1] / L) * ysim[d + 1, ] / Q
      
      # report progress
      if (disp_progress) cat("\r Simulation progress: ",
                             progresslist[d + 1], "%")
    }
  }
  if (disp_progress) cat("\n")
  
  # compute output
  if (!components_out) return(colSums(ysim))
  else return(ysim)
}



Precompute_etbm <- function(Q, D) {
  #-----------------------------------------------------------------------------
  # Precomputation for Emery-TBM
  # Computes the matrices A^(d) for all d up to D
  #
  # DEPENDENCIES
  #  - UniformSpherePoints() from "helpfunctions.R"
  #  - Legendre() from "legendre.R"
  #
  # INPUT
  #  Q:  number of frequencies used in TBM-algorithm
  #  D:  maximum degree of KL-component
  #
  # OUTPUT
  #   (QxQx(D+1)) array containing A^(d) matrices (QxQ)
  #-----------------------------------------------------------------------------
  
  # setup output array
  a_array <- array(numeric(Q * Q * (D + 1)), c(Q, Q, D + 1) )
  
  # Compute cosangles for evenly distributed points
  directions <- UniformSpherePoints(Q)
  cosangles <- t(directions) %*% directions
  diag(cosangles) <- rep(1, Q)
  
  # loop over KL-components
  for (d in 0:D) {
    
    # parameter (dependent on n)
    if (d == 0) lambda <- 1 / 2
    else lambda <- d / 2
    
    # compute ratio of densities (dependent on n)
    ratio <- exp(lambda) * sqrt(8 / (pi * lambda**d)) * gamma((d + 3) / 2)
    
    # get matrix of legendre polynomials evaluations at angles
    legendrematrix <- Legendre(d, cosangles)
    
    # get eigenvectors and values
    eigen_both <- eigen(legendrematrix)
    V <- eigen_both$vectors
    ev <- eigen_both$values
    ev[ev < 1e-14] <- 0      # handle numerical errors
    
    # get matrix A
    a_array[ , ,d + 1] <-  (V * rep(ratio * sqrt(ev), each = Q))
  }
  
  return(a_array)
}
