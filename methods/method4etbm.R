#-------------------------------------------------------------------------------
#          Method 4: Emery Turning Bands Method - "etbm"
#------------------------------------------------------------------------------- 
# File:    method4etbm.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - This method requires a grfmodel of modeltype "spectral".
#  - There are no precomputation options for this implementation.
#  - This method has significant (co)variance deficits.
#
# DEPENDENCIES
#  - RandomSpherePoints() from "helpfunctions.R"
#  - Legendre() from "legendre.R"
#
# CONTENT
#   .SimGRF_etbm()         Main simulation function
#-------------------------------------------------------------------------------

# Sources
gsl.installed <- require(gsl)
if (!gsl.installed) {
  source("R/simGRFS/helpers/no_gsl/legendre_no_gsl.R")
} else  {
  source("R/simGRFS/helpers/legendre.R")
}
source("R/simGRFS/helpers/helpfunctions.R")



.SimGRF_etbm <- function(dummyvar,
                         pts,
                         grfmodel,
                         L = 30,
                         Q = 30,
                         disp_progress = FALSE,
                         components_out = FALSE) {
  #-----------------------------------------------------------------------------
  # Simulates values of a GRF on the sphere given the points and a grfmodel.
  #
  # No precomputation options.
  #
  # DEPENDENCIES
  # - RandomSpherePoints() from "helpfunctions.R"
  # - Legendre() from "legendre.R"
  #
  # INPUT
  #        dummyvar:  dummy variable for parallel computing
  #             pts:  (3xN) matrix of points to simulate in format c(x,y,z)
  #        grfmodel:  instance of class GRFmodel, specifying covariance model
  #               L:  number of independent realizations to use CLT 
  #               Q:  number of frequencies to use 
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
  
  # prepare ysim matrix (columns = points, rows = d's)
  ysim <- matrix(numeric((D + 1) * N), D + 1, N) 
  
  # loop over all Legendre components Y^d_sim
  for (n in 0:D) {
    
    # parameter (dependent on n)
    if (n == 0) lambda <- 1 / 2
    else lambda <- n / 2
    
    # compute ratio of densities (dependent on n)
    ratio <- exp(lambda) * sqrt(8 / (pi * lambda**n)) * gamma((n + 3) / 2)
    
    # loop over independent realizations
    for (l in 1:L) {
      
      # random directions in R^3
      directions <- RandomSpherePoints(Q)
      
      # random frequencies in R^3 distributed as g (chi^2-norms)
      norms <- sqrt(2 * lambda * rchisq(Q, df = n + 3))
      freqs <- directions * rep(norms, each = 3)
      
      # random phase
      phi <- runif(1, min = 0, max = 2 * pi)
      
      # get matrix of legendre polynomials evaluations at angles
      cosangles <- t(directions) %*% directions
      diag(cosangles) <- rep(1, Q)
      legendrematrix <- Legendre(n, cosangles)
      
      # get eigenvectors and values
      eigendecomp <- eigen(legendrematrix)
      V <- eigendecomp$vectors
      ev <- eigendecomp$values
      ev[ev < 1e-14] <- 0      # handle numerical errors
      
      # computation of vector B by computing matrix A
      # (more efficient code for (ratio * (V %*% diag(sqrt(ev)))))
      A <- (V * rep(ratio * sqrt(ev), each = Q))
      U <- rnorm(Q)
      B <- A %*% U

      # loop over all points and main computation
      ysim[n + 1, ] <- ysim[n + 1, ] + cos(t(pts) %*% freqs + phi) %*% B
    }
    
    # correct for 1/sqrt(L) and sqrt(a_n) spectral measure
    ysim[n + 1, ] <- sqrt(coefs[n + 1] / L) * ysim[n + 1, ] / Q
    
    # report progress
    if (disp_progress) cat("\r Simulation progress: ", progresslist[n + 1], "%")
  }
  if (disp_progress) cat("\n")
  
  # return output
  if (!components_out) return(colSums(ysim))
  else return(ysim)
}

