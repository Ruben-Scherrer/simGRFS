#-------------------------------------------------------------------------------   
#     Setup of GRFModel class for using it in GRF simulations
#-------------------------------------------------------------------------------  
# File:    grfmodels.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - The parameters param1, param2 correspond to rho_1, rho_2 in the report.
#
# DEPENDENCIES
#  - default parameters from "defaultparameters.R"
#  - covariance functions from "spectralcovariance.R"
#
# CONTENT
#   GRFModel()      creates an instance of class GRFModel 
#-------------------------------------------------------------------------------   

# Sources
source("R/simGRFS/models/defaultparameters.R")
source("R/simGRFS/models/spectralcovariance.R")



GRFModel <- function(modeltype = NULL, 
                     covmodel = NULL,
                     param1 = NULL,
                     param2 = NULL,
                     spectmodel = NULL,
                     D = 30,
                     coefs = NULL,
                     spectparam = NULL) {
  #-----------------------------------------------------------------------------
  # Creates an object of S3 class GRFModel using either covariance functions
  # (modeltype "covar") or the angular power spectrum (modeltype "spectral"), 
  # also called Schoenberg coefficients.
  #
  # If nothing is specified, per default a covariance powered exponential model
  # with default parameters (2,1) gets created.
  #
  # Comments modeltype "covar":
  #     - Implemented covmodels: "exp", "powexp", "sinepow", "spherical",
  #                              "askey", "wendc2", "wendc4", "cmatern"
  # 
  # Comments modeltype "spectral":
  #     - Implemented spectmodels: "exp", "geometric", "poisson", "bessel", 
  #                                  "linear"
  #     - Schoenberg coefficients can either be provided directly via the
  #       argument coefs (option A), or they can be generated using
  #       the arguments spectmodel (and spectparam) to choose the desired
  #       model (option B).
  #
  # DEPENDENCIES
  # - default parameters from "defaultparameters.R"
  # - functions from "spectralcovariance.R"
  #
  # INPUT
  #       modeltype:   either "covar" or "spectral"
  #
  # INPUT additionally for modeltype "covar":
  #        covmodel:   model of covariance function
  #          param1:   (optional) hyperparameter 1 of covariance model 
  #          param2:   (optional) hyperparamter 2 of covariance model
  #
  # INPUT additionally for modeltype "spectral":
  #           coefs:   (option A) Schoenberg coefficients of power spectrum
  #      spectmodel:   (option B) model of Schoenberg coefficients
  #               D:   (option B) maximal degree of Schoenberg coefficients
  #      spectparam:   (option B) (optional) hyperparameter of Schoenberg model
  #
  # OUTPUT
  # instance of class "GRFModel", which is a standardized list with entries:
  #      $modeltype:   either "covar" or "spectral"
  #
  # (if modeltype is "covar":)
  #       $covmodel:   model of covariance function
  #         $param1:   hyperparameter 1 of covariance model 
  #         $param2:   hyperparamter 2 of covariance model
  #
  # (if modeltype is "spectral":)
  #     $spectmodel:   model of Schoenberg coefficients
  #              $D:   maximal degree of Schoenberg coefficients computed
  #     $spectparam:   hyperparameter of spectral model
  #          $coefs:   Schoenberg coefficients of power spectrum
  #-----------------------------------------------------------------------------
  
  #--------------------------- PREPARATIONS ------------------------------------
  # if no modeltype specified, set modeltype
  if (is.null(modeltype) && is.null(spectmodel) &&
      is.null(coefs)) modeltype <- "covar"
  else if (!is.null(spectmodel) && is.null(covmodel) ||
           !is.null(coefs) && is.null(covmodel)) modeltype <- "spectral"
  else if (!is.null(spectmodel) && !is.null(covmodel) ||
           !is.null(coefs) && !is.null(covmodel)) {
    stop("Error: Incompatible model specifications provided.") 
  }
  
  #-------------------------- ERROR HANDLING -----------------------------------
  
  # PROBLEM 1: incompatible parameters and modeltype
  if (modeltype == "covar" && !is.null(coefs) ||
      modeltype == "covar" && !is.null(spectparam) ||
      modeltype == "covar" && !is.null(spectmodel)) {
    warning("Warning: Incompatible parameters to modeltype.
             Assuming spectral modeltype")
    modeltype <- "spectral"
  } else if (modeltype == "spectral" && !is.null(covmodel) ||
             modeltype == "spectral" && !is.null(param1) ||
             modeltype == "spectral" && !is.null(param2)) {
    warning("Warning: Incompatible specifications to modeltype.")
  }
  
  #--------------------- MAIN COMPUTATION covar-models -------------------------
  if (modeltype == "covar") {
    
    # if no model provided: set default model
    if (is.null(covmodel)) covmodel <- "exp"
    
    # set default parameters for parameters that are not provided
    if (is.null(param1)) {
      if (covmodel == "exp") {
        param1 = def_exp_p1
      } else if (covmodel == "powexp") {
        param1 = def_powexp_p1
      } else if (covmodel == "sinepow") {
        param1 = def_sinepow_p1
      } else if (covmodel == "spherical") {
        param1 = def_spherical_p1
      } else if (covmodel == "askey") {
        param1 = def_askey_p1
      } else if (covmodel == "wendc2") {
        param1 = def_wendc2_p1
      } else if (covmodel == "wendc4") {
        param1 = def_wendc4_p1
      } else if (covmodel == "cmatern") {
        param1 = def_cmatern_p1
      } else if (covmodel == "geometric") {
        param1 = def_spect_geometric
      } else if (covmodel == "poisson") {
        param1 = def_spect_poisson
      } else if (covmodel == "bessel") {
        param1 = def_spect_bessel
      } else if (covmodel == "linear") {
        param1 = def_spect_linear
      } else stop("Error: Incompatible model provided")
    }
    
    if (is.null(param2)) {
      if (covmodel == "exp") {
        param2 = def_exp_p2
      } else if (covmodel == "powexp") {
        param2 = def_powexp_p2
      } else if (covmodel == "sinepow") {
        param2 = def_sinepow_p2
      } else if (covmodel == "spherical") {
        param2 = def_spherical_p2
      } else if (covmodel == "askey") {
        param2 = def_askey_p2
      } else if (covmodel == "wendc2") {
        param2 = def_wendc2_p2
      } else if (covmodel == "wendc4") {
        param2 = def_wendc4_p2
      } else if (covmodel == "cmatern") {
        param2 = def_cmatern_p2
      } else if (covmodel == "geometric") {
        param2 = NULL
      } else if (covmodel == "poisson") {
        param2 = NULL
      } else if (covmodel == "bessel") {
        param2 = NULL
      } else if (covmodel == "linear") {
        param2 = NULL
      } else stop("Error: Incompatible model provided")
    }
    
    # create model
    model_out <- list(modeltype = "covar",
                      covmodel = covmodel,
                      param1 = param1, param2 = param2)

    # define class
    class(model_out) <- "GRFModel"
    
    return(model_out)
  }
  
  #--------------------- MAIN COMPUTATION spectral-models ----------------------
  else if (modeltype == "spectral") {
    
    # if no model provided: set default model
    if (is.null(spectmodel)) spectmodel <- "exp"
    
    # if coefficients provided directly, create model
    if (!is.null(coefs)) {
      D <- length(coefs) - 1
      model_out <- list(modeltype = "spectral", D = D, coefs = coefs)
      class(model_out) <- "GRFModel"
      return(model_out)
    }
    
    # if no parameters provided, set default values
    else if (is.null(spectparam)) {
      if (spectmodel == "exp") {
        spectparam <- def_spect_exp
      } else if (spectmodel == "geometric") {
        spectparam <- def_spect_geometric
      } else if (spectmodel == "poisson") {
        spectparam <- def_spect_poisson
      } else if (spectmodel == "bessel") {
        spectparam <- def_spect_bessel
      } else if (spectmodel == "linear") {
        spectparam <- def_spect_linear
      } else stop("Error: incompatible model provided")
    }
    
    # get Schoenberg coefficients
    if (spectmodel == "exp") {
      coefs <- SpectCovExp(D = D, spectparam = spectparam)
    } else if (spectmodel == "geometric") {
      coefs <- SpectCovGeometric(D = D, spectparam = spectparam)
    } else if (spectmodel == "poisson") {
      coefs <- SpectCovPoisson(D = D, spectparam = spectparam)
    } else if (spectmodel == "bessel") {
      coefs <- SpectCovBessel(D = D, spectparam = spectparam)
    } else if (spectmodel == "linear") {
      coefs <- SpectCovLinear(D = D, spectparam = spectparam)
    }
    
    # create model
    model_out <- list(modeltype = "spectral",
                      spectmodel = spectmodel, D = D,
                      spectparam = spectparam, coefs = coefs)
    
    # assign class
    class(model_out) <- "GRFModel"
    
    return(model_out)
  }
  
  else {
    stop("Error: incompatible modeltype provided.")
  }
}
