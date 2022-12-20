#------------------------------------------------------------------------------- 
#          Method 2: Karhunen Loève Expansion  - "kle"
#------------------------------------------------------------------------------- 
# File:    method2kle.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - This method requires a grfmodel of modeltype "spectral".
#  - The precomputations function is not used in the method itself since it
#    uses a lot of memory, which is not necessary if one can just loop over
#    the different components.
#  - However, if multiple simulations of the same GRF and the same points are
#    to be realized, the precomputation has to be done only once and can be
#    passed as an argument.
#
# DEPENDENCIES
#  - functions from "sphericalharmonics.R"
#  - Lat2Colat() from "helpfunctions.R"
#
# CONTENT
#   .SimGRF_kle()            Main simulation function
#   Precompute_kle()   Precomputations (evaluated SH of given points)
#-------------------------------------------------------------------------------

# Libraries
gsl.installed <- require(gsl)

# sources
if (!gsl.installed) {
  source("R/simGRFS/helpers/no_gsl/sphericalharmonics_no_gsl.R")
} else  {
  source("R/simGRFS/helpers/sphericalharmonics.R")
}
source("R/simGRFS/helpers/helpfunctions.R")



.SimGRF_kle <- function(dummyvar,
                        pts,
                        grfmodel,
                        precomputations = NULL,
                        components_out = FALSE) {
  #-----------------------------------------------------------------------------
  # Simulates values of a GRF on the sphere given the points and a grfmodel.
  #
  # Precomputation option: The spherical harmonics can be precomputated and
  # provided as an argument.
  #
  # DEPENDENCIES
  # - .SphericalHarmonics_nn() from "sphericalharmonics.R"
  # - Lat2Colat() from "helpfunctions.R"
  #
  # INPUT
  #        dummyvar:  dummy variable for parallel computing
  #             pts:  (3xN) matrix of points to simulate in format c(lat,long,r)
  #        grfmodel:  instance of class GRFmodel with modeltype "spectral"
  # precomputations:  (optional) array of precomputed spherical harmonics 
  #  components_out:  if TRUE, returns Z^d_sim-values instead of Z^sim
  #
  # OUTPUT
  #   vector (N) with simulated values in the same order as input, or
  #   matrix ((D+1)xN) with Z^d_sim-values
  #-----------------------------------------------------------------------------
  N <- dim(pts)[2]
  
  # get model parameters
  coefs <- grfmodel$coefs
  D <- grfmodel$D

  # setup cmatrix
  cmatrix <- diag(D + 1)
  cmatrix[ ,1] <- sqrt(coefs * 4 * pi / (2 * (0:D) + 1)) * rnorm(D + 1)
  for (d in 2:(D + 1)) {
    cmatrix[d,2:d] <- sqrt(coefs[d] * 8 * pi / (2 * d - 1)) *
                      complex(real = rnorm(d - 1), imaginary = rnorm(d - 1))
  }

  # transform latitudes to colatitudes
  pts <- Lat2Colat(pts)
  
  #----------------- MAIN COMPUTATION: option components_out -------------------
  if (components_out) {
    
    # no precomputations provided
    if (is.null(precomputations)) {
      zsim <- mapply(FUN = function(lat, long)
                           rowSums(Re(.SphericalHarmonics_nn(D, lat, long) *
                                     cmatrix)),
                    lat = pts[1, ], long = pts[2, ])
    }
    
    # precomputations provided
    else {
      zsim <- apply(Re(precomputations * rep(cmatrix, N)), 3, rowSums)
    }
  }
  
  #----------------- MAIN COMPUTATION: standard option -------------------------
  else {
    
    # no precomputations provided
    if (is.null(precomputations)) {
      zsim <- mapply(FUN = function(lat, long)
                          sum(Re(.SphericalHarmonics_nn(D, lat, long) *
                                 cmatrix)),
                    lat = pts[1, ], long = pts[2, ])
    }
    
    # precomputations provided
    else {
      zsim <- apply(Re(precomputations * rep(cmatrix, N)), 3, sum)
    }
  }

  return(zsim)
}



Precompute_kle <- function(pts, D) {
  #-----------------------------------------------------------------------------
  # Precomputation for KLE-method.
  # Computes the spherical harmonics (SH) for all points of nonnegative orders
  # up to degree D.
  #
  # DEPENDENCIES
  #  - .SphericalHarmonics_nn() from "sphericalharmonics.R"
  #  - Lat2Colat() from "helpfunctions.R"
  #
  # INPUT
  #  pts:  (3xN) matrix of points to simulate in format c(lat,long,r)
  #    D:  maximal degree of spherical harmonics to compute
  #
  # OUTPUT
  #   ((D+1)x(D+1)xN) array of evaluated SH at given points
  #-----------------------------------------------------------------------------
  N <- dim(pts)[2]
  
  # transform latitudes to colatitudes
  pts <- Lat2Colat(pts)
  
  # set up output array dim(rows: degree of SH; column: order of SH; N)
  sh_array <- array(numeric((D + 1)**2 * N), c(D + 1, D + 1, N))
  
  # get spherical harmonics for every point
  for (n in 1:N) {
    sh_array[ , ,n] <- .SphericalHarmonics_nn(D, pts[1,n], pts[2,n])
  }
  
  return(sh_array)
}

