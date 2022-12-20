#-------------------------------------------------------------------------------   
#     Functions for computing line process covariances (for ATBM on sphere)
#------------------------------------------------------------------------------- 
# File:    linecovariance.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# DEPENDENCIES
#  - default parameters from "defaultparameters.R"
#  - Symmetrize() from "helpfunctions.R"
#
# CONTENT
# Main functions
#   LineCov()          Computes covariance values for distances on line
#   LineCovMat()       Computes covariance matrix for line process
# 
# Covariance functions
#   LineCovExp()       Exponential
#   LineCovPowexp()    Powered Exponential
#   LineCovAskey()     Askey
#   LineCovSinepow()   Sine power
#   LineCovCmatern()   Chordal Matèrn
#-------------------------------------------------------------------------------

# Sources
source("R/simGRFS/models/defaultparameters.R")



LineCov <- function(r, grfmodel) {
  #-----------------------------------------------------------------------------
  # Computes covariance of distance r on line process.
  #
  # DEPENDENCIES
  # - Linecovariance functions
  #
  # INPUT
  #         r:  distance on line
  #  grfmodel:  covariance model (object of class GRFModel)
  #
  # OUTPUT
  #  covariance values
  #-----------------------------------------------------------------------------
  
  # handle wrong input models
  if (grfmodel$modeltype != "covar" ||
      grfmodel$covmodel %in% c("wendc2", "wendc4", "spherical",
                               "geometric", "bessel", "poisson", "linear")) {
    stop("Incompatible GRFmodel provided")
  }
  
  # get model parameters
  covmodel <- grfmodel$covmodel
  param1 <- grfmodel$param1
  param2 <- grfmodel$param2
  
  # compute covariance
  if (covmodel == "exp")            out <- LineCovExp(r, param1, param2)
  else if (covmodel == "powexp")    out <- LineCovPowexp(r, param1, param2)
  else if (covmodel == "sinepow")   out <- LineCovSinepow(r, param1, param2)
  else if (covmodel == "askey")     out <- LineCovAskey(r, param1, param2)
  else if (covmodel == "cmatern")   out <- LineCovCmatern(r, param1, param2)
  
  return(out) 
}



LineCovMat <- function(grfmodel, P, variance = 1) {
  #-----------------------------------------------------------------------------
  # Computes covariance matrix for line process using P points
  #
  # DEPENDENCIES
  # - Linecovariance functions
  # - Symmetrize() from "helpfunctions.R"
  #
  # INPUT
  #  grfmodel     covariance model (object of class GRFModel)
  #  P            number of points to simulate on the line
  #
  # OUTPUT
  #  covariance matrix of line process
  #-----------------------------------------------------------------------------
  
  # handle wrong input models
  if (grfmodel$modeltype != "covar" ||
      grfmodel$covmodel %in% c("wendc2", "wendc4", "spherical",
                               "geometric", "bessel", "poisson", "linear")) {
    stop("Incompatible GRFmodel provided")
  }
  
  # get model parameters
  covar <- grfmodel$covmodel
  param1 <- grfmodel$param1
  param2 <- grfmodel$param2
  
  # get default distance
  distance <- 2 / P
  
  # compute covariances
  if (covar == "exp") {
    covmat <- toeplitz(LineCovExp((0:(P - 1)) * distance, param1 = param1))
  } else if (covar == "powexp") {
    covmat <- toeplitz(LineCovPowexp((0:(P - 1)) * distance,
                                     param1 = param1,
                                     param2 = param2))
  } else if (covar == "sinepow") {
    covmat <- toeplitz(LineCovSinepow((0:(P - 1)) * distance,
                                      param1 = param1,
                                      param2 = param2))
  } else if (covar == "askey") {
    covmat <- toeplitz(LineCovAskey((0:(P - 1)) * distance,
                                    param1 = param1,
                                    param2 = param2))
  } else if (covar == "cmatern") {
    covmat <- toeplitz(LineCovCmatern((0:(P - 1)) * distance,
                                      param1 = param1,
                                      param2 = param2))
  }
  
  # if non-unit variance, multiply with variance
  if (variance != 1) covmat <- variance * covmat
  
  return(covmat)
}



#------------------------------------------------------------------------------- 
#               Covariance Functions for line processes
#
# Each of the following covariance functions has the same structure:
#
# No dependencies.
#
# INPUT
#       r:  distance(s) between two points on the line 
#  param1:  hyperparamater 1
#  param2:  hyperparameter 2
#
# OUTPUT
#    Corresponding covariance(s) as a number
#------------------------------------------------------------------------------- 

LineCovExp <- function(r,
                       param1 = def_exp_p1,
                       param2 = def_exp_p2) {
  #-----------------------------------------------------------------------------
  # Parameter range: param1 > 0
  #
  # Corresponds to spherical covariance:
  # C(alpha) = exp(-param1*alpha)
  #-----------------------------------------------------------------------------
  return(exp(-2 * param1 * asin(r / 2)) * (1 - 2 * param1 * r / sqrt(4 - r**2)))
}

LineCovPowexp <- function(r,
                          param1 = def_powexp_p1,
                          param2 = def_powexp_p2) {
  #-----------------------------------------------------------------------------
  # Parameter range: param1 > 0, 0 < param2 <= 1
  #
  # Corresponds to spherical covariance:
  # C(alpha) = exp(-(param1*alpha)**param2)
  #-----------------------------------------------------------------------------
  return(exp(-(2 * param1 * asin(r / 2))**param2) *
         (1 - r * param2 * (2 * param1)**param2 * (asin(r / 2))**(param2 - 1) /
          sqrt(4 - r**2)))
}

LineCovAskey <- function(r,
                         param1 = def_askey_p1,
                         param2 = def_askey_p2) {
  #-----------------------------------------------------------------------------
  # Parameter range: param1 > 0, param2 >= 2
  #
  # Corresponds to spherical covariance:
  # C(alpha) = (1-param1*alpha)**param2_+
  #-----------------------------------------------------------------------------
  return(pmax(0, (1 - param1 * 2 * asin(r / 2)))**param2 -
         2 * r * param1 * param2 / sqrt(4 - r**2) * 
         pmax(0, (1 - param1 * 2 * asin(r / 2)))**(param2 - 1))
}

LineCovSinepow <- function(r,
                           param1 = def_sinepow_p1,
                           param2 = def_sinepow_p2) {
  #-----------------------------------------------------------------------------
  # Parameter range: 0 < param2 < 1
  #
  # Corresponds to spherical covariance:
  # C(alpha) = 1-sin(alpha/2)**param2
  #-----------------------------------------------------------------------------
  return(1 - (r / 2)**param2 * (1 + param2))
}

LineCovCmatern <- function(r,
                           param1=def_cmatern_p1,
                           param2=def_cmatern_p2) {
  #-----------------------------------------------------------------------------
  # Parameter range: param1, param2 > 0
  #
  # Corresponds to spherical covariance:
  # C(alpha) =(param1*2*sin(alpha/2))**param2 * K_param2(param1*2*sin(alpha/2))
  #-----------------------------------------------------------------------------
  return(sign(r) * (param1 * r)**param2 *
         (1 / (2**(param2 - 1) * gamma(param2))) *
         ((1 + param2) * besselK(param1 * r**sign(r), nu = param2) -
         param1 * r / 2 * (besselK(param1 * r**sign(r), nu = param2 - 1) +
         besselK(param1 * r**sign(r), nu = param2 + 1))) + 1 - sign(r))
}
