#-------------------------------------------------------------------------------   
#          Method 3: Lantuéjoul Spectral Method  - "ls"
#------------------------------------------------------------------------------- 
# File:    method3ls.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - This method requires a grfmodel of modeltype "spectral".
#  - The precomputation function computes much more values than are necessary
#    and hence is not used for a single GRF.
#  - However, if multiple simulations of the same GRF and the same  points are
#    to be realized, the precomputation can be done once and then be used for
#    all realizations.
#
# DEPENDENCIES 
#  - functions from "sphericalharmonics.R"
#  - Lat2Colat() from "helpfunctions.R"
#  - Package "gsl" for fast performance
#
# CONTENT
#   .SimGRF_ls()      Main simulation function
#   Precompute_ls()   Precomputations (evaluates SH parts of input)
#-------------------------------------------------------------------------------

# Libraries
gsl.installed <- require(gsl)

# Sources
source("R/simGRFS/helpers/helpfunctions.R")
if (!gsl.installed) {
  source("R/simGRFS/helpers/no_gsl/sphericalharmonics_no_gsl.R")
} else  {
  source("R/simGRFS/helpers/sphericalharmonics.R")
}



.SimGRF_ls <- function(dummyvar,
                             pts,
                             grfmodel,
                             L = 30,
                             precomputations = NULL) {
  #-----------------------------------------------------------------------------
  # Simulates values of a GRF on the sphere given the points and a grfmodel.
  #
  # Precomputation option: The spherical harmonics can be precomputated and
  # provided as an argument.
  #
  # DEPENDENCIES
  #  - function .SphericalHarmonic_part() from "sphericalharmonics.R"
  #  - function Lat2Colat() from "helpfunctions.R"
  #
  # INPUT
  #        dummyvar:  dummy variable for parallel computing
  #             pts:  (3xN) matrix of points in format c(lat,long,r)
  #        grfmodel:  instance of class GRFmodel with modeltype "spectral"
  #               L:  number of independent realizations to use CLT with
  # precomputations:  (optional) ((D+1)x(2D+1)xN) array of evaluated
  #                   spherical harmonics at given points.
  #
  # OUTPUT
  #   vector (N) with simulated values in the same order as input
  #-----------------------------------------------------------------------------
  N <- dim(pts)[2]
  
  # get Schoenberg coefficients
  coefs <- grfmodel$coefs
  D <- grfmodel$D
  
  # normalize schoenberg coefficients for using them as probability measure
  spect_measure <- coefs / sum(coefs)
  
  # get random components
  degrees <- sample(0:D, L, replace = TRUE, prob = spect_measure)
  orders <- as.vector(unlist(lapply(degrees, function(x) sample(-x:x, 1))))
  # faster this way than by using mapply...
  psis <- runif(L, 0, 2 * pi)
  
  
  # transform latitudes to colatitudes
  pts <- Lat2Colat(pts)
  
  #-------------------------- MAIN COMPUTATION ---------------------------------
  
  # without precomputations
  if (is.null(precomputations)) {

    # sum independent realizations for main computation
    zsim <- sqrt(8 * pi / L) *
            rowSums(mapply(FUN = function(deg, ord, psi)
                                 .SphericalHarmonic_part(deg, ord, pts[1, ]) *
                                 cos(ord * pts[2, ] + psi),
                           deg = degrees, ord = orders, psi = psis))
  }
  
  # with precomputations
  else {
    
    # sum independent realizations for main computation
    zsim <- sqrt(8 * pi / L) *
            rowSums(mapply(FUN = function(deg, ord, psi)
                                 precomputations[deg + 1,D + 1 + ord, ] *
                                 cos(ord * pts[2, ] + psi),
                           deg = degrees, ord = orders, psi = psis))
  }
  
  return(zsim)
}



Precompute_ls <- function(pts, D) {
  #-----------------------------------------------------------------------------
  # Precomputation for LS method. Computes the parts of spherical
  # harmonics (SH), needed for LS method, for all points up to deg D.
  #
  # DEPENDENCIES
  #  - .SphericalHarmonics_parts() from "sphericalharmonics.R"
  #  - Lat2Colat() from "helpfunctions.R"
  #
  # INPUT
  #   pts:  (3xN) matrix of points to simulate in format c(lat,long,r)
  #     D:  maximal degree of spherical harmonics part to compute
  #
  # OUTPUT
  #   ((D+1)x(2D+1)xN) array of evaluated SH parts at given points
  #-----------------------------------------------------------------------------
  N <- dim(pts)[2]
  
  # transform latitudes to colatitudes
  pts <- Lat2Colat(pts)

  # set up output array dim(rows: degree of SH; column: order of SH; N)
  spher_harm <- array(numeric((D + 1) * (2 * D + 1) * N),
                      c(D + 1, 2 * D + 1, N))
  
  # get spherical harmonics for every point
  for (n in 1:N) {
    spher_harm[ , ,n] <- .SphericalHarmonics_parts(D, pts[1,n])
  }
  
  return(spher_harm)
}
