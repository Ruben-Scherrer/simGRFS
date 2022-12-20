#-------------------------------------------------------------------------------  
#          Method 6: Approximate TBM - "atbm"
#------------------------------------------------------------------------------- 
# File:    method6atbm.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9

# COMMENTS
#  - This method requires a grfmodel of type "covar".
#  - Different simulation approached are available for this method. "approx" 
#    uses a suitable amount of points on the line processes, with the risk that
#    very close points on the sphere will have the same simulated value.
#    "exact" computes the minimal distance between two points and then 
#    adapts the number of points P to be simulated on the line to this distance.
#    "nugget" allows to define a small nugget effect, which will eliminate the
#    problem of duplicate values.
#  - At the moment, implemented covariance models for this method are: "exp", 
#    "powexp", "sinepow", "askey", "cmatern"
#
# DEPENDENCIES
#  - functions from "linecovariance.R"
#  - functions from "helpfunctions.R"
#  - RandomSimGRF() from "main.R"
#
# CONTENT
#   .SimGRF_atbm()         Main simulation function
#   Precompute_atbm()   Computes cholesky triangle of covmat for line pr.
#-------------------------------------------------------------------------------

# Sources
source("R/simGRFS/helpers/helpfunctions.R")
source("R/simGRFS/models/linecovariance.R")



.SimGRF_atbm <- function(dummyvar,
                         grfmodel,
                         pts,
                         Q = 30,
                         P = NULL,
                         approach = "approx",
                         nugget = 0,
                         correct_variance = TRUE,
                         precomputations = NULL,
                         silent = FALSE,
                         contrast_only = FALSE) {
  #-----------------------------------------------------------------------------
  # Simulates values of a GRF on the sphere given the points and a grfmodel.
  #
  # Precomputation option: The cholesky factor of the line process can be 
  #                        passed as an argument.
  #
  # DEPENDENCIES
  #  - RandomSimGRF() from "main.R"
  #  - functions from "helpfunctions.R"
  #  - LineCov() from "linecovariance.R"
  #  - Precompute_atbm()
  #
  # INPUT
  #        dummyvar:  dummy variable for parallel computing
  #             pts:  (3xN) matrix of points to simulate in format c(x,y,z)
  #        grfmodel:  instance of class GRFmodel, specifying covariance model
  #               Q:  number of frequencies to use
  #        approach:  either "approx" for ignoring points close to each other
  #                   and allowing to use a specified or computed P,
  #                   or "exact" for evaluating the minimal distance and 
  #                   adapting P to it,
  #                   or "nugget" to allow for a given nugget or a nugget to be
  #                   assessed by the algorithm
  #               P:  (optional) number of points to simulate on lines
  #          nugget:  (optional) nugget effect to use on GRF
  #correct_variance:  whether to correct for variance for "nugget" option
  #          silent:  If true ignores warning messages about "approx" approach.
  #   contrast_only:  if true, doesn't simulate a GRF. instead uses uniform
  #                   values on a suitable scale to detect in values for grids
  #
  # OUTPUT
  #   vector (N) with simulated values in the same order as input
  #-----------------------------------------------------------------------------
  N <- dim(pts)[2]
  
  #------------------- PREPARATIONS:  without precomputations ------------------
  if (is.null(precomputations)) {
    if (approach == "exact") {
      md <- MinDistSphere(pts)
      approx_interval_length <- sqrt((md**2) / 2)
    } else approx_interval_length <- exp(-log(N, base = 10))
    if (is.null(P)) P <- ceiling(2 / approx_interval_length)
    
    # get time computation for large simulations (or too close points)
    if (P >= 1000 && is.null(precomputations)) {
      
      # get reference time (computer specific) and compute eta
      reftimes <- numeric(6)
      for (i in 1:6) {
        reftimes[i] <- system.time(RandomSimGRF(100, method = "md",
                                                disp_info = FALSE))[3]
      }
      eta <- P**log(P, base = 400) * mean(reftimes[3:6]) +
             N * mean(reftimes) / 100
      
      # display information
      if (!silent) {
        cat(" >> information:  large simulation, ", 
            "estimated waiting time: ", eta, " seconds \n")
      }
    }
    
    # info about P
    if (!silent) cat(" P =", P, "\n") 
    
    # get cholesky factor of line process
    cholesky_factor <- Precompute_atbm(grfmodel, P)
  }
  
  #-------------------- PREPARATIONS:  with precomputations --------------------
  else {
    cholesky_factor <- precomputations
    P <- dim(cholesky_factor)[1]
  }
  
  #-------------------- FURTHER PREPARATIONS -----------------------------------
  
  # compute nugget if not provided
  if (approach == "nugget" && is.null(nugget)) {
    nugget <- 2 - 2 * LineCov(1 / P, grfmodel)
  }
  
  # get directions (evenly distributed here)
  directions <- UniformSpherePoints(Q)

  #----------------------- MAIN COMPUTATIONS -----------------------------------
  
  # contrast only option
  if (contrast_only) line_values <- matrix(runif(P * Q, min = 0, max = P), P, Q)
  
  # simulation using MD-method
  else line_values <- cholesky_factor %*% matrix(rnorm(P * Q), P, Q)

  # compute simulation values of GRF (using projections)
  indices <- round((t(pts) %*% directions + 1) * P / 2 - 0.5) + 1
  indices[indices == P + 1] <- P  # correct case of pts=directions
  mat_indices <- indices + rep(P * (0:(Q - 1)), each = N)
  zsim <- colSums(matrix(line_values[t(mat_indices)], nrow = Q)) / sqrt(Q)
  
  # correct for nugget 
  if (nugget != 0) {
    if (correct_variance) zsim <- zsim * sqrt(1 - nugget) +
                                 rnorm(N, mean = 0, sd = sqrt(nugget))
    else zsim <- zsim + rnorm(N, mean = 0, sd = sqrt(nugget))
  }

  #---------------------- DISPLAY INFORMATION ----------------------------------
  if (approach == "approx" && nugget == 0 && !silent) {
    n_duplicates <- N - length(unique(zsim))
    if (n_duplicates > 0) {
      cat(" >> information: Approx approach without nugget,",
        n_duplicates / N * 100, "% of points are duplicates. \n")
    }
  }
  
  return(zsim)
}



Precompute_atbm <- function(grfmodel, P) {
  #-----------------------------------------------------------
  # Precompute the cholesky decomposition of the covariance matrix for
  # P points on the line for some line process for ATBM method.
  #
  # DEPENDENCIES
  # - LineCovMat() from "linecovariance.R"
  #
  # INPUT
  #  grfmodel   specified covariance model, object of class "GRFModel"
  #  P          number of points to simulate on the line
  #
  # OUTPUT
  #  cholesky factor of the covariance matrix of line process
  #-----------------------------------------------------------
  
  # get covariance matrix and cholesky factorization
  covmat <- LineCovMat(grfmodel, P)
  cholesky_factor <- t(chol(covmat))
  
  return(cholesky_factor)
}
