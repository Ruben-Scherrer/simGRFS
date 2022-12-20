#-------------------------------------------------------------------------------   
#     Functions for computing covariances on the sphere
#------------------------------------------------------------------------------- 
# File:    covariance.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - In some models, there are unnecessary parameters. They can be ignored.
#  - The parameters param1, param2 correspond to rho_1, rho_2 in the report.
#
# DEPENDENCIES 
#  - Symmetrize() from "helpfunctions.R"
#  - default parameters from "defaultparameters.R"
#
# CONTENT
# Main functions:
#   TrueCov()      Compute covariances using grfmodel
#   TrueCovMat()   Computes true covariance matrix for given points
#   SampleCovMat() Computes sample covariance matrix
# 
# Covariance functions:
#   CovExp()        Exponential (special case of "powexp")
#   CovPowexp()     Powered exponential
#   CovSinepow()    Sine power
#   CovSpherical()  Spherical
#   CovAskey()      Askey
#   CovWendc2()     Wendland-C2
#   CovWendc4()     Wendland-C4
#   CovCmatern()    Chordal Matèrn
#   CovGeometric()  Geometric
#   CovPoisson()    Poisson
#   CovBessel()     Bessel
#   CovLinear()     Linear
#-------------------------------------------------------------------------------   

source("R/simGRFS/models/defaultparameters.R")
source("R/simGRFS/helpers/helpfunctions.R")



TrueCov <- function(alphas, grfmodel) {
  #-----------------------------------------------------------------------------
  # Compute true covariances for angles alphas using the covariance
  # model specified by grfmodel. 
  #
  # alphas ranges 0 <= alpha <= 2*pi
  #
  # DEPENDENCIES
  # - covariance functions
  #
  # INPUT
  #   alphas:  angles between two points on unit sphere
  # grfmodel:  instance of class GRFmodel, specifying covariance model
  # 
  # OUTPUT
  #   resulting covariance
  #-----------------------------------------------------------------------------
  
  # get model parameters
  if (grfmodel$modeltype == "covar") {
  covmodel <- grfmodel$covmodel
  param1 <- grfmodel$param1
  param2 <- grfmodel$param2
  } else {
   spectmodel <- grfmodel$spectmodel
   spectparam <- grfmodel$spectparam 
  }
  
  # compute covariance
  if (grfmodel$modeltype == "covar") {
  if (covmodel == "exp")            out <- CovExp(alphas, param1, param2)
  else if (covmodel == "powexp")    out <- CovPowexp(alphas, param1, param2)
  else if (covmodel == "spherical") out <- CovSpherical(alphas, param1, param2)
  else if (covmodel == "sinepow")   out <- CovSinepow(alphas, param1, param2)
  else if (covmodel == "askey")     out <- CovAskey(alphas, param1, param2)
  else if (covmodel == "wendc2")    out <- CovWendc2(alphas, param1, param2)
  else if (covmodel == "wendc4")    out <- CovWendc4(alphas, param1, param2)
  else if (covmodel == "cmatern")   out <- CovCmatern(alphas, param1, param2)
  else if (covmodel == "geometric") out <- CovGeometric(alphas, param1, param2)
  else if (covmodel == "poisson")   out <- CovPoisson(alphas, param1, param2)
  else if (covmodel == "bessel")    out <- CovBessel(alphas, param1, param2)
  else if (covmodel == "linear")    out <- CovLinear(alphas, param1, param2)
  } else {
    if (spectmodel == "exp")        out <- CovExp(alphas, param1 = spectparam)
    else if (spectmodel == "geometric") out <- CovGeometric(alphas, spectparam)
    else if (spectmodel == "poisson")   out <- CovPoisson(alphas, spectparam)
    else if (spectmodel == "bessel")    out <- CovBessel(alphas, spectparam)
    else if (spectmodel == "linear")    out <- CovLinear(alphas, spectparam)
  }
  return(out)
}



TrueCovMat <- function(pts, grfmodel, variance = 1) {
  #-----------------------------------------------------------------------------
  # Compute true covariance matrix for points (in cart. coord) using the
  # covariance model specified by grfmodel of type "covar"
  #
  # DEPENDENCIES
  # - covariance functions
  #
  # INPUT
  #      pts:  matrix of points in column format c(x,y,z)
  # grfmodel:  instance of class GRFmodel, specifying covariance model
  #      var:  variance factor of GRF 
  # 
  # OUTPUT
  #   (NxN) covariance matrix of given points
  #-----------------------------------------------------------------------------

  # get model parameters
  if (grfmodel$modeltype == "covar") {
    covmodel <- grfmodel$covmodel
    param1 <- grfmodel$param1
    param2 <- grfmodel$param2
  } else {
    spectmodel <- grfmodel$spectmodel
    spectparam <- grfmodel$spectparam 
  }
  
  # get angles between points
  proj <- t(pts) %*% pts
  diag(proj) <- rep(1, length(pts) / 3)   # handle numerical errors 
  alphas <- suppressWarnings(acos(proj))
  diag(alphas) <- rep(0, length(pts) / 3) # handle numerical errors
  
  # compute covariance
  if (grfmodel$modeltype == "covar") {
    if (covmodel == "exp")            out <- CovExp(alphas, param1, param2)
    else if (covmodel == "powexp")    out <- CovPowexp(alphas, param1, param2)
    else if (covmodel == "spherical") out <- CovSpherical(alphas, param1, param2)
    else if (covmodel == "sinepow")   out <- CovSinepow(alphas, param1, param2)
    else if (covmodel == "askey")     out <- CovAskey(alphas, param1, param2)
    else if (covmodel == "wendc2")    out <- CovWendc2(alphas, param1, param2)
    else if (covmodel == "wendc4")    out <- CovWendc4(alphas, param1, param2)
    else if (covmodel == "cmatern")   out <- CovCmatern(alphas, param1, param2)
    else if (covmodel == "geometric") out <- CovGeometric(alphas, param1, param2)
    else if (covmodel == "poisson")   out <- CovPoisson(alphas, param1, param2)
    else if (covmodel == "bessel")    out <- CovBessel(alphas, param1, param2)
    else if (covmodel == "linear")    out <- CovLinear(alphas, param1, param2)
  } else {
    if (spectmodel == "exp")        out <- CovExp(alphas, param1 = spectparam)
    else if (spectmodel == "geometric") out <- CovGeometric(alphas, spectparam)
    else if (spectmodel == "poisson")   out <- CovPoisson(alphas, spectparam)
    else if (spectmodel == "bessel")    out <- CovBessel(alphas, spectparam)
    else if (spectmodel == "linear")    out <- CovLinear(alphas, spectparam)
  }

  # if non-unit variance, multiply with variance
  if (variance != 1) out <- variance * out
  
  return(out)
}



SampleCovMat <- function(values) {
  #-----------------------------------------------------------------------------
  # Computes sample covariance matrix, either for one experiment of 
  # multiple realizations, or for multiple experiments of mul. real.
  #
  # DEPENDENCIES
  #  - Symmetrize() from "helpfunctions.R"
  #
  # INPUT
  #  values: either matrix (M x N) or array (M x N x E)
  # 
  # OUTPUT
  #   (NxN) sample covariance matrix
  #-----------------------------------------------------------------------------
  
  # get size parameters and unify array
  M <- dim(values)[1]
  N <- dim(values)[2]
  if (is.na(dim(values)[3])) values <- array(values, c(M, N, 1))
  E <- dim(values)[3]
  
  # setup of storage matrix
  samplecov_array <- array(numeric(N**2), c(N, N, E))
  
  # loop through experiments
  for (e in 1:E) {
    
    # loop through matrix
    for (i in 1:(N - 1)) {
      samplecov_array[i,i,e] <- var(values[ ,i,e])
      for (j in (i + 1):N) {
        samplecov_array[i,j,e] <- cov(values[ ,i,e], values[ ,j,e])
      }
    }
    samplecov_array[N,N,e] <- var(values[ ,N,e])
  }
  
  # get average for multiple realizations
  out <- Symmetrize(apply(samplecov_array, c(1, 2), mean))
  
  return(out)
}



#-------------------------------------------------------------------------------   
#                         Covariance Functions
#
# Each of the following covariance functions has the same structure:
#
# No dependencies.
#
# INPUT
#   alpha:  angle between two points on the sphere
#  param1:  hyperparamater 1
#  param2:  hyperparameter 2
#
# OUTPUT
#    Corresponding covariance as a number
#-------------------------------------------------------------------------------  

CovExp <- function(alpha,
                   param1 = def_powexp_p1,
                   param2 = def_powexp_p2) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 0 < param1
  # pararm2 unused
  #-----------------------------------------------------------------------------
  return(exp(-(param1 * alpha)))
}

CovPowexp <- function(alpha,
                      param1 = def_powexp_p1,
                      param2 = def_powexp_p2) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 0 < param1
  # 0 < param2  <= 1
  #-----------------------------------------------------------------------------
  return(exp(-(param1 * alpha)**param2))
}

CovSinepow <- function(alpha,
                       param1 = def_sinepow_p1,
                       param2 = def_sinepow_p2) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # param1 unused
  # 0 < param2 < 2
  #-----------------------------------------------------------------------------
  return(1 - (sin(alpha / 2))**param2)
}

CovSpherical <- function(alpha,
                         param1 = def_spherical_p1,
                         param2 = def_spherical_p2) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 0 < param1
  # param2 unused
  #-----------------------------------------------------------------------------
  return((1 + (param1 * alpha)/2) * pmax(0, (1 - param1 * alpha))**2)
}

CovAskey <- function(alpha,
                     param1 = def_askey_p1,
                     param2 = def_askey_p2) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 0 < param1
  # 2 <= param2
  #-----------------------------------------------------------------------------
  return(pmax(0, (1 - param1 * alpha))**param2 * 1**alpha)
}

CovWendc2 <- function(alpha,
                      param1 = def_wendc2_p1,
                      param2 = def_wendc2_p2) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 1/pi < param1
  # 4 <= param2
  #-----------------------------------------------------------------------------
  return((1 + param2 * param1 * alpha) *
         (pmax(0, (1 - param1 * alpha))**param2))
}

CovWendc4 <- function(alpha,
                      param1 = def_wendc4_p1,
                      param2 = def_wendc4_p2) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 1 < param1
  # 6 <= param2
  #-----------------------------------------------------------------------------
  return((1 + param2 * param1 * alpha + (param2**2 - 1) / 3 *
          (param1 * alpha)**2) * (pmax(0, (1 - param1 * alpha)))**param2)
}

CovCmatern <- function(alpha,
                       param1 = def_cmatern_p1,
                       param2 = def_cmatern_p2) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 0 < param1
  # 0 < param2
  # If parameters are too low, covariance is too high to get exact results.
  # Recommended parameters
  # param1 > param2, param2 <= 2
  #-----------------------------------------------------------------------------
  return(sign(alpha)*(1 / (2**(param2 - 1) * gamma(param2)) *
           (param1 * 2 * sin(alpha / 2))**param2 * 
           besselK(param1 * 2 * sin(alpha / 2)**sign(alpha), param2)) + 
           1 - sign(alpha))
}


CovGeometric <- function(alpha,
                         param1 = def_spect_geometric,
                         param2 = NULL) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 0 < param1 < 1
  #-----------------------------------------------------------------------------
  return((1 - param1) / sqrt(1 - 2 * param1 * cos(alpha) +
                                 param1**2))
}

CovPoisson <- function(alpha,
                       param1 = def_spect_poisson,
                       param2 = NULL) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 0 < param1
  #-----------------------------------------------------------------------------
  return(exp(param1 * (cos(alpha) - 1)) * 
         besselJ(param1 * sin(alpha), nu = 0))
}

CovBessel <- function(alpha,
                      param1 = def_spect_bessel,
                      param2 = NULL) {
  #-----------------------------------------------------------------------------
  # hyperparameter ranges:
  # 0 < param1
  #-----------------------------------------------------------------------------
  return(exp(param1 * (cos(alpha) - 1)))
}

CovLinear <- function(alpha,
                      param1 = def_spect_linear,
                      param2 = NULL) {
  #-----------------------------------------------------------------------------
  # no hyperparameters
  #-----------------------------------------------------------------------------
  return(1 - (2 * alpha) / pi)
}

