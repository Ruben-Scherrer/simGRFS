#-------------------------------------------------------------------------------
#          Functions for evaluating spherical harmonics
#------------------------------------------------------------------------------- 
# File:    sphericalharmonics.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - Due to normalized evaluation the algorithms are quite stable. However, for
#    large degrees and points close to the poles the values can get too
#    big or too small.
#
# DEPENDENCIES
#  - functions from package "gsl"
#
# CONTENT
#   SphericalHarmonic()         Evaluate spherical harmonic (SH)
#   SphericalHarmonics()        Evaluate all spherical harmonics up to deg
#   .SphericalHarmonics_nn()    Evalute all SH of nonnegative order up to deg
#   .SphericalHarmonic_part()   Evaluates part of SH (used for spect method)
#   .SphericalHarmonics_parts() Evaluates all parts of SH up to degree deg
#-------------------------------------------------------------------------------



SphericalHarmonic <- function(n, k, theta, phi) {
  #-----------------------------------------------------------------------------
  # Direct computation of spherical harmonic with degree n, order k
  # at data point (theta,phi) on the sphere.
  #
  # theoretical range theta: [0,pi]
  # range phi: [0,2*pi)
  # 
  #  DEPENDENCIES
  #  - gsl::legendre_sphPlm()
  #
  # INPUT
  #     n:   degree of spherical harmonic to be computed
  #     k:   order of spherical harmonic to be computed
  # theta:   polar coordinates of point to evaluate
  #   phi:   polar coordinates of point to evaluate
  #
  # OUPUT
  #  (complex) value of spherical harmonic at point (theta,phi)
  #-----------------------------------------------------------------------------
  return(exp(complex(imaginary = k * phi)) * 
         gsl::legendre_sphPlm(n, abs(k), cos(theta)) * sign(k)**k)
}



SphericalHarmonics <- function(D, theta, phi) {
  #-----------------------------------------------------------------------------
  # Computation of spherical harmonics up to degree D
  # at data point (theta,phi) on the sphere.
  #
  # theoretical range theta [0,pi]
  # range phi [0,2*pi)
  #
  #  DEPENDENCIES
  #  - gsl::legendre_array()
  #
  # INPUT
  #      D:   maximal degree of spherical harmonic to be computed
  #  theta:   polar coordinates of points to evaluate
  #    phi:   polar coordinates of points to evaluate
  #
  # OUPUT 
  #  ((D+1)x(2D+1)) matrix with generic entry [i,j] the evaluated 
  #         spherical harmonic of degree i-1 and order j-D-1
  #-----------------------------------------------------------------------------
  
  # prepare output matrix
  out <- matrix(numeric((D + 1) * (2 * D + 1)),
                D + 1, 2 * D + 1)
  
  # get positive orders and multiply by complex factor
  out[ ,(D + 1):(2 * D + 1)] <- gsl::legendre_array(cos(theta), D, norm = 3) *
                                exp(complex(imaginary = 1))**
                                (t(replicate(D + 1, 0:D)) * phi)
  
  # get negative orders
  out[ ,D:1] <- (-1)**t(replicate(D + 1, 1:D)) *
                Conj(out[ ,(D + 2):(2 * D + 1)])
  
  return(out)
}



.SphericalHarmonics_nn <- function(D, theta, phi) {
  #-----------------------------------------------------------------------------
  # Computation of spherical harmonics of nonnegative orders up to degree D
  # at data point (theta,phi) on the sphere.
  #
  # theoretical range theta [0,pi]
  # range phi [0,2*pi)
  #
  #  DEPENDENCIES
  #  - gsl::legendre_array()
  #
  # INPUT
  #      D:   maximal degree of spherical harmonic to be computed
  #  theta:   polar coordinates of point to evaluate
  #    phi:   polar coordinates of point to evaluate
  #
  # OUPUT
  #  ((D+1)x(D+1)) matrix with generic entry [i,j] the evaluated 
  #         spherical harmonic of degree i-1 and order j-1
  #-----------------------------------------------------------------------------
  return(gsl::legendre_array(cos(theta), D, norm = 3) *
         exp(complex(imaginary = 1))**(t(replicate(D + 1, 0:D)) * phi))
}



.SphericalHarmonic_part <- function(n, k, theta) {
  #-----------------------------------------------------------------------------
  # Wrapper for GSL-function legendre_sphPlm(). 
  # Computes part of spherical harmonic (without the added angle):
  # sqrt((2n+1)/4pi*(n-k)!/(n+k)!)P_n,k(cos(theta))
  #
  # Remark: This is NOT exactly the real part of the spherical harmonic!
  #
  # theoretical range theta: [0,pi]
  # range phi: [0,2*pi)
  # 
  #  DEPENDENCIES
  #  - gsl::legendre_sphPlm()
  #
  # INPUT
  #     n:   degree of spherical harmonic part to be computed
  #     k:   nonnegative order of spherical harmonic part to be computed
  # theta:   polar coordinates of point to evaluate
  #
  # OUPUT
  #  value of part of spherical harmonic at point (theta,phi)
  #-----------------------------------------------------------------------------
  return(gsl::legendre_sphPlm(n, abs(k), cos(theta)) * sign(k)**k)
}



.SphericalHarmonics_parts <- function(D, theta, phi) {
  #-----------------------------------------------------------------------------
  # Computation of all parts of spherical harmonics up to degree D at data
  # point (theta,phi) on the sphere. (see SphericalHarmonic_part()) 
  #
  # theoretical range theta [0,pi]
  # range phi [0,2*pi)
  #
  #  DEPENDENCIES
  #  - gsl::legendre_array()
  #
  # INPUT
  #      D:   maximal degree of spherical harmonic part to be computed
  #  theta:   polar coordinates of point to evaluate
  #    phi:   polar coordinates of point to evaluate
  #
  # OUPUT 
  #  ((D+1)x(2D+1)) matrix with generic entry [i,j] the evaluated 
  #       part of spherical harmonic of degree i-1 and order j-D-1
  #-----------------------------------------------------------------------------
  
  # setup output matrix
  out <- matrix(numeric((D + 1) * (2 * D + 1)), D + 1, 2 * D + 1)
  
  # compute values using GSL-function legendre_array() and (-1)**ord
  out[ ,(D + 1):(2 * D + 1)] <- gsl::legendre_array(cos(theta), D, norm = 3)
  out[ ,D:1] <- (-1)**t(replicate(D + 1, 1:D)) * out[ ,(D + 2):(2 * D + 1)]
  
  return(out)
}
