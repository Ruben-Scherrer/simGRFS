#-------------------------------------------------------------------------------   
#           Functions for computing covariances for spectral models
#-------------------------------------------------------------------------------  
# File:    spectralcovariance.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# DEPENDENCIES
#  - Symmetrize() from "helpfunctions.R"
#  - Legendres() from "legendre.R"
#  - default parameters from "defaultparameters.R"
#
# CONTENT
# Main function: 
#   ApproxCovMat()        Approximates true covariance matrix
#   ApproxCov()           Approximates covariances of multiple distances
#   SchoenbergExpansion() Approximates covariance via Schoenberg expansion
#
# Schoenberg coefficient generators:
#   SpectCovExp()         for covariance C(alpha) = exp(- param * alpha) 
#   SpectCovGeometric()   for covariance C(alpha) = (1-param)/
#                                        sqrt(1-2*param*cos(alpha)+param**2)
#   SpectCovPoisson()     for covariance C(alpha) = exp(param*(cos(alpha)-1))*
#                                                      J_0(param*sin(alpha))
#   SpectCovBessel()      for covariance C(alpha) = exp(param*(cos(alpha)-1))
#   SpectCovLinear()      for covariance C(alpha) = 1-2*alpha/pi
#-------------------------------------------------------------------------------   

# sources
source("R/simGRFS/models/defaultparameters.R")
source("R/simGRFS/helpers/helpfunctions.R")
source("R/simGRFS/helpers/legendre.R")



ApproxCovMat <- function(pts, grfmodel) {
  #-----------------------------------------------------------------------------
  # Approximate the covariance matrix using Schoenberg coefficients.
  # 
  # For values very close to 0 or to pi, some numerical instability occurs.
  # The algorithm (dependent on grfmodel$D) should be STABLE for the ranges:
  #       alpha in (0.0142,pi-0.0142)         up to D=85
  #       alpha in (0.0000142, pi -0.0000142) up to D=50
  #
  # DEPENDENCIES
  # - Symmetrize() from "helpfunctions.R"
  # - ApproxCov()
  #
  # INPUT
  #       pts: points to compute covariance matrix of in format c(x,y,z)
  #  grfmodel: instance of class GRFModel, specifying covariance model
  #
  # OUTPUT
  #   approximate covariance matrix (NxN)
  #-----------------------------------------------------------------------------
  N <- dim(pts)[2]
  
  # get angles between points
  proj <- t(pts) %*% pts
  diag(proj) <- rep(1, length(pts) / 3)   # handle numerical errors
  alphas <- suppressWarnings(acos(proj))
  diag(alphas) <- rep(0, length(pts) / 3) # handle numerical errors
  
  # set up output matrix
  out <- diag(rep(1, N))
  
  # compute only triangle (since matrix is symmetric)
  out[lower.tri(out)] <- ApproxCov(alphas[lower.tri(out)], grfmodel = grfmodel)
  
  # symmetrize and return
  return(Symmetrize(t(out)))
}



ApproxCov <- function(alphas, grfmodel) {
  #-----------------------------------------------------------------------------
  # Approximate the covariances of distances alphas using Schoenberg coefs.
  # 
  # alpha theoretical range: 0 <= alpha <= pi
  #
  # For values very close to 0 or to pi, some numerical instability occurs.
  # The algorithm (dependent on grfmodel$D) should be STABLE for the ranges:
  #       alpha in (0.0142,pi-0.0142)         up to D=85
  #       alpha in (0.0000142, pi -0.0000142) up to D=50
  #
  # DEPENDENCIES
  #  - SchoenbergExpansion()
  #
  # INPUT
  #    alpha:  distance to compute covariance of
  #  grfmodel: instance of class GRFmodel, specifying covariance model
  #
  # OUTPUT
  #   covariance of alpha given the covariance model
  #-----------------------------------------------------------------------------
  N <- length(alphas)
  
  # handle wrong input models
  if (grfmodel$modeltype == "covar") stop("Error: wrong input model")

  # main computation
  out <- numeric(N)
  for (n in 1:N) {
    out[n] <- SchoenbergExpansion(alpha = alphas[n], grfmodel = grfmodel)
  }
  
  return(out)
}



SchoenbergExpansion <- function(alpha, grfmodel) {
  #-----------------------------------------------------------------------------
  # Approximate the covariance of distance alpha using Schoenberg coefficients.
  # 
  # alpha range: 0 <= alpha <= pi
  #
  # DEPENDENCIES
  # - Legendres() from "legendre.R"
  #
  # INPUT
  #    alpha:  distance to compute covariance of
  #  grfmodel: instance of class GRFModel, specifying covariance model
  #
  # OUTPUT
  #   covariance of alpha given the covariance model
  #-----------------------------------------------------------------------------
  
  # handle wrong input models
  if (grfmodel$modeltype == "covar") stop("Error: wrong input model")
  
  # get Schoenberg coefficients and model parameter
  coefs <- grfmodel$coefs
  D <- grfmodel$D
  
  # compute truncated expansion
  out <- sum(coefs * Legendres(D, cos(alpha)))
  
  return(out)
}



#-------------------------------------------------------------------------------   
#                         Schoenberg Coefficients
#
# Each of the following generating functions has the same structure:
#
# No Dependencies.
#
# INPUT
#           D:  maximal degree of Schoenberg coefficient to generate
#  spectparam:  hyperparameter
#
# OUTPUT
#   D+1 Schoenberg coefficients as vector (starting with degree 0)
#   (Entry j is Schoenberg coefficient of degree j-1)
#------------------------------------------------------------------------------- 

SpectCovExp <- function(D, spectparam = def_spect_exp) {
  #-----------------------------------------------------------------------------
  # Parameter range: 0 < spectparam
  #
  # Corresponds to covariance:
  # C(alpha) = exp(- spectparam * alpha) 
  #-----------------------------------------------------------------------------
  
  # initialize list and compute first entries
  coefs <- numeric(D + 1)
  coefs[1] <- 0.5 * (1 + exp(-spectparam * pi)) / (1 + spectparam**2)
  coefs[2] <- 1.5 * (1 - exp(-spectparam * pi)) / (4 + spectparam**2)
  
  # recursively compute entries 
  for (d in 3:(D + 1)) {
    n <- d - 1
    coefs[d] <- coefs[d - 2] * (2 * n + 1) * (spectparam**2 + (n - 2)**2) /
                ((2 * n - 3) * (spectparam**2 + (n + 1)**2))
  }
  
  return(coefs)
}

SpectCovGeometric <- function(D, spectparam = def_spect_geometric) {
  #-----------------------------------------------------------------------------
  # Parameter range: 0 < spectparam < 1
  #
  # Corresponds to covariance:
  # C(alpha) = (1-spectparam) / 
  #            sqrt(1-2*spectparam*cos(alpha)+spectparam**2)
  #-----------------------------------------------------------------------------
  
  # initialize list and compute first entry
  coefs <- numeric(D + 1)
  coefs[1] <- 1 - spectparam
  
  # recursively compute entries
  for (d in 2:(D + 1)) coefs[d] <- coefs[d - 1] * spectparam
  
  return(coefs)
}

SpectCovPoisson <- function(D, spectparam = def_spect_poisson) {
  #-----------------------------------------------------------------------------
  # Parameter range: 0 < spectparam
  #
  # Corresponds to covariance:
  # C(alpha) = exp(spectparam*(cos(alpha)-1))*J_0(spectparam*sin(alpha))
  #-----------------------------------------------------------------------------
  
  # initialize list and compute first entry
  coefs <- numeric(D + 1)
  coefs[1] <- exp(-spectparam)
  
  # recursively compute entries
  for (d in 2:(D + 1)) coefs[d] <- coefs[d - 1] * spectparam / (d - 1)
  
  return(coefs)
}

SpectCovBessel <- function(D, spectparam = def_spect_bessel) {
  #-----------------------------------------------------------------------------
  # Parameter range: 0 < spectparam
  #
  # Corresponds to covariance:
  # C(alpha) = exp(spectparam*(cos(alpha)-1))
  #-----------------------------------------------------------------------------
  
  # initialize list and compute constant
  constant <- sqrt(pi / (2 * spectparam)) * exp(-spectparam)
  coefs <- rep(constant, D + 1)
  
  # compute entries
  for (d in 1:(D + 1)) {
    coefs[d] <- coefs[d] * (2 * (d - 1) + 1) * besselI(spectparam, d - 0.5)
  }
  
  return(coefs)
}

SpectCovLinear <- function(D, spectparam = def_spect_linear) {
  #-----------------------------------------------------------------------------
  # (No parameter used.)
  # (Remark: Even coefficients are zero.)
  #
  # Corresponds to covariance:
  # C(alpha) = 1-2*alpha/pi
  #-----------------------------------------------------------------------------
  
  # initialize list and compute first entries
  coefs <- numeric(D + 1)
  coefs[1] <- 0
  coefs[2] <- 3 / 4
  
  # compute entries
  for (d in 3:(D + 1)) {
    n <- d - 1
    coefs[d] <- coefs[d - 2] * (2 * n + 1) * (n - 2)**2 /
                ((2 * n - 3) * (n + 1)**2)
  }
  
  return(coefs)
}
