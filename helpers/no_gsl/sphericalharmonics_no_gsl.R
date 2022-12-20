#-------------------------------------------------------------------------------   
#          Functions for evaluating spherical harmonics without GSL
#-------------------------------------------------------------------------------  
# File:    sphericalharmonics_no_gsl.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - Using these functions instead of the GSL versions results in
#             similar performance of SphericalHarmonic()
#             0.15x slower performance for SphericalHarmonics()
#  - The functions are generally stable up to degree 85 (afterwards the numbers
#    get too small very quickly).
#  - As a consequence of using the AssocLegendre functions, the 
#             functions are stable only for 
#             range theta in (0.0142,pi-0.0142)         up to degree 85
#             range theta in (0.0000142, pi -0.0000142) up to degree 50
#
# DEPENDENCIES
#  - functions from "legendre_no_gsl.R"
#
# CONTENT
#   SphericalHarmonic()         Evaluate spherical harmonic (SH)
#   SphericalHarmonics()        Evaluate all spherical harmonics up to deg
#   .SphericalHarmonics_nn()    Evalute all SH of nonnegative order up to deg
#   .SphericalHarmonic_part()   Evaluates part of SH (used for spect method)
#   .SphericalHarmonics_parts() Evaluates part of SHs of nn ord up do deg
#------------------------------------------------------------------------------- 

# Sources
source("R/simGRFS/helpers/no_gsl/legendre_no_gsl.R")



SphericalHarmonic <- function(n, k, theta, phi) {
  #-----------------------------------------------------------------------------
  # Direct computation of spherical harmonic with degree n, order k
  # at data point (theta,phi) on the sphere.
  #
  # theoretical range theta: [0,pi]
  # practical range: see stability issues in comments.
  # range phi: [0,2*pi)
  # 
  #  DEPENDENCIES
  #  - AssocLegendre() from "legendre_no_gsl.R"
  #
  # INPUT
  #     n:  degree of spherical harmonic to be computed
  #     k:  order of spherical harmonic to be computed
  # theta:  polar coordinates of point to evaluate
  #   phi:  polar coordinates of point to evaluate
  #
  # OUPUT
  #  value of spherical harmonic at point (theta,phi)
  #-----------------------------------------------------------------------------
  return(sqrt((2 * n + 1) * factorial(n - k) /
              (4 * pi * factorial(n + k))) *
         AssocLegendre(n, k, cos(theta), csphase = TRUE) *
         exp(complex(imaginary = k * phi)))
}



SphericalHarmonics <- function(D, theta, phi) {
  #-----------------------------------------------------------------------------
  # Computation of spherical harmonics up to degree D
  # at data point (theta,phi) on the sphere.
  #
  # theoretical range theta [0,pi]
  # practical range: see stability issues in comments.
  # range phi [0,2*pi)
  #
  #  DEPENDENCIES
  #  - AssocLegendres() from "legendre_no_gsl.R"
  #
  # INPUT
  #     D:   maximal degree of spherical harmonic to be computed
  # theta:   polar coordinates of point to evaluate
  #   phi:   polar coordinates of point to evaluate
  #
  # OUPUT 
  #  ((D+1)x(2D+1)) matrix with generic entry [i,j] the evaluated 
  #         spherical harmonic of degree i-1 and order j-D-1
  #-----------------------------------------------------------------------------
  
  # take ALP-matrix of cos(theta) as starting point and compute first entry
  out <- AssocLegendres(D, cos(theta))
  out[1,D + 1] <- sqrt(1 / (4 * pi)) * out[1,D + 1]
  
  # set up exp-terms
  expvec <- exp(complex(imaginary = (1:D) * phi))
  
  #  loop over all degrees
  for (n in 1:D) {
    
    # prepare factor
    fa1 <- sqrt((2 * n + 1) / (4 * pi))
    fa2 <- 1
    
    # update all entries with order 0 (exp-term is 1)
    out[n + 1,D + 1] <- fa1 * out[n + 1,D + 1]
    
    for (k in 1:n) {
      
      # update factor 2
      fa2 <- fa2 / sqrt((n - k + 1) * (n + k))
      
      # update values of positive and negative order
      out[n + 1,D + 1 + k] <- out[n + 1,D + 1 + k] * expvec[k] * fa1 * fa2
      out[n + 1,D + 1 - k] <- (-1)**k * Conj(out[n + 1,D + 1 + k])
    }
  }
  
  return(out)
}



.SphericalHarmonics_nn <- function(D, theta, phi) {
  #-----------------------------------------------------------------------------
  # Computation of spherical harmonics of positive order up to degree D
  # at data point (theta,phi) on the sphere.
  #
  # (not optimized, since AssocLegendres() computes also negative orders)
  #
  # theoretical range theta [0,pi]
  # practical range: see stability issues in comments.
  # range phi [0,2*pi)
  #
  #  DEPENDENCIES
  #  - AssocLegendres() from "legendre_no_gsl.R"
  #
  # INPUT
  #     D:   maximal degree of spherical harmonic to be computed
  # theta:   polar coordinates of point to evaluate
  #   phi:   polar coordinates of point to evaluate
  #
  # OUPUT 
  #  ((D+1)x(2D+1)) matrix with generic entry [i,j] the evaluated 
  #         spherical harmonic of degree i-1 and order j-n-1
  #-----------------------------------------------------------------------------
  
  # take (part of) ALP-matrix of cos(theta) as starting point
  out <- AssocLegendres(D, cos(theta))[ ,(D + 1):(2 * D + 1)]
  out[1,1] <- sqrt(1 / (4 * pi)) * out[1,1]
  
  # set up exp-terms
  expvec <- exp(complex(imaginary = (1:D) * phi))
  
  #  loop over all degrees
  for (n in 1:D) {
    
    # prepare factors
    fa1 <- sqrt((2 * n + 1) / (4 * pi))
    fa2 <- 1
    
    # update all entries with order 0 (exp-term is 1)
    out[n + 1,1] <- fa1 * out[n + 1,1]
    
    for (k in 1:n) {
      
      # update factor 2
      fa2 <- fa2 / sqrt((n - k + 1) * (n + k))
      
      # update values of positive and negative order
      out[n + 1,1 + k] <- out[n + 1,1 + k] * expvec[k] * fa1 * fa2
    }
  }
  
  return(out)
}



.SphericalHarmonic_part <- function(n, k, theta) {
  #-----------------------------------------------------------------------------
  # Computes part of spherical harmonic (without the added angle):
  # sqrt((2n+1)/4pi*(n-k)!/(n+k)!)P_n,k(cos(theta))
  #
  # Remark: This is NOT exactly the real part of the spherical harmonic!
  #
  # theoretical range theta: [0,pi]
  # 
  #  DEPENDENCIES
  #  - AssocLegendre() from "legendre_no_gsl.R"
  #
  # INPUT
  #     n:   degree of spherical harmonic to be computed
  #     k:   nonnegative order of spherical harmonic to be computed
  # theta:   polar coordinates of point to evaluate
  #
  # OUPUT
  #  value of part of spherical harmonic at point (theta,phi)
  #-----------------------------------------------------------------------------
  
  return(sqrt((2 * n + 1) * factorial(n - k) /
              (4 * pi * factorial(n + k))) *
         AssocLegendre(n, k, cos(theta), csphase = TRUE))
}



.SphericalHarmonics_parts <- function(D, theta) {
  #-----------------------------------------------------------------------------
  # Computation of part of spherical harmonics of positive orders up to degree D
  # at data point (theta,phi) on the sphere. (see .SphericalHarmonic_part())
  #
  # theoretical range theta [0,pi]
  #
  # DEPENDENCIES
  # - function AssocLegendres() from "legendre_no_gsl.R"
  #
  # INPUT
  #      D:   maximal degree of spherical harmonic to be computed
  #  theta:   polar coordinates of point to evaluate
  #
  # OUPUT 
  #  ((D+1)x(D+1)) matrix with generic entry [i,j] the evaluated 
  #        part of spherical harmonic of degree i-1 and order j-1
  #-----------------------------------------------------------------------------
  
  # take ALP-matrix of cos(theta) as starting point and compute first entry
  out <- AssocLegendres(D, cos(theta))
  out[1,D + 1] <- sqrt(1 / (4 * pi)) * out[1,D + 1]
  
  #  loop over all degrees
  for (n in 1:D) {
    
    # prepare factor
    fa1 <- sqrt((2 * n + 1) / (4 * pi))
    fa2 <- 1
    
    # update all entries with order 0 (exp-term is 1)
    out[n + 1,D + 1] <- fa1 * out[n + 1,D + 1]
    
    for (k in 1:n) {
      
      # update factor 2
      fa2 <- fa2 / sqrt((n - k + 1) * (n + k))
      
      # update values of positive and negative order
      out[n + 1,D + 1 + k] <- out[n + 1,D + 1 + k] * fa1 * fa2
      out[n + 1,D + 1 - k] <- (-1)**k * out[n + 1,D + 1 + k]
    }
  }
  
  return(out)
}


