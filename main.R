#------------------------------------------------------------------------------- 
#
#           MAIN FILE
#           Simulating GRFs on the sphere
#
#------------------------------------------------------------------------------- 
# File:    main.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - See the README for basic usage.
#  - There are 2 standard methods ("md", "atbm") and 4 spectral methods ("kle",
#    "ls", "etbm", "etbm2") available.
#  - For MD method, all models are available.
#  - For spectral methods the following models are avaliable:
#           "exp", "geometric", "poisson", "linear", "bessel"
#  - For ATBM the following models are available: 
#           "exp", "powexp", "askey", "sinepow", "cmatern"
#
# DEPENDENCIES
#  - functions from "grfmodels.R", "helpfunctions.R", "method1md.R",
#                   "method2kle.R", "method3ls.R", "method4etbm.R",
#                   "method5etbm2.R" and "method6atbm.R"
#  - for displaying: packages "rgl", "plotrix", functions from "displaysphere.R"
#  - for fast performance: package "gsl"
#  - indirectly: functions from "covariance.R", "spectralcovariance.R",
#                   "linecovariance.R", parameters from "defaultparameters.R",
#                   and either functions from
#                   "legendre.R" and "sphericalharmonics.R" or
#                   "legendre_no_gsl.R" and "sphericalharmonics_no_gsl.R"
#
# CONTENT
#   SimGRF()            Main function: simulate isotropic GRF on the sphere
#   RandomSimGRF()      Wrapper: Simulate GRF with random points on the sphere
#   MultiSimGRF()       Compute many realizations of the same points and GRFs
#   EvalGRFComponents() Evaluates KL-components for given covariance model
#------------------------------------------------------------------------------- 

# Sources
source("R/simGRFS/helpers/helpfunctions.R")
source("R/simGRFS/displaysphere.R")
source("R/simGRFS/models/grfmodels.R")
source("R/simGRFS/methods/method1md.R")
source("R/simGRFS/methods/method2kle.R")
source("R/simGRFS/methods/method3ls.R")
source("R/simGRFS/methods/method4etbm.R")
source("R/simGRFS/methods/method5etbm2.R")
source("R/simGRFS/methods/method6atbm.R")



SimGRF <- function(pts,
                   grfmodel = NULL,
                   method = "md",
                   precomputations = NULL,
                   display = FALSE,
                   in_polar_coord = FALSE,
                   out_format = "default",
                   disp_info = TRUE,
                   seed = NULL,
                   L = NULL,                    # ls, etbm
                   Q = NULL,                    # etbm, etbm2, atbm
                   P = NULL,                    # atbm
                   atbm_approach = "approx",    # atbm
                   atbm_nugget = 0,             # atbm
                   atbm_nugget_cv = TRUE,       # atbm
                   etbm_progress = FALSE,       # etbm, etbm2
                   slim = FALSE) {     
  #-----------------------------------------------------------------------------
  # Main Function for simulating GRFs on the sphere
  #
  # Either returns values of simulated GRF (if display set to FALSE, default),
  # or displays the simulated points on a sphere with radius 1 and color
  # coded values.
  #
  # DEPENDENCIES
  #  - functions from "helpfunctions.R"
  #  - functions from and file "displaysphere.R"
  #  - GRFModel() from "grfmodels.R"
  #  - functions from "method1md.R", "method2kle.R", "method3ls.R"
  #                   "method4etbm.R", "method5etbm2.R" and "method6atbm.R"
  #  - package "rgl" for displaying options
  #
  # INPUT
  #              pts:  (3xN) matrix of points to simulate in col format c(x,y,z)
  #                    or (3xN) matrix in column format c(lat,long,r)
  #                    or (2xN) matrix in column format c(lat,long)
  #                    (if input in polar coord, set in_polar_coord=TRUE)
  #         grfmodel:  instance of class GRFModel, specifying covariance model
  #                    - methods "md" and "atbm" require modeltype "covar"
  #                    - all other methods require modeltype "spectral"
  #           method:  string, specifies simulation method. options: "md", 
  #                    "kle", "ls", "etbm", "etbm2", "atbm"
  #  precomputations:  (optional) precomputations to pass
  #          display:  boolean, if TRUE, displays points, no return value
  #   in_polar_coord:  boolean, if input is given in polar coord., set TRUE
  #       out_format:  string, either "default" for returning a matrix with
  #                    columns in format c(x,y,z,values)
  #                    or "values" for a vector of only the simulated values
  #                    or "components" for a matrix of simulated KL-components
  #        disp_info:  boolean, whether to display information about simulation
  #             seed:  (optional) number used as seed to allow reproducibility
  #                L:  number of ind. realizations for "ls", "etbm", "etbm2" 
  #                Q:  number of frequencies to use for "etbm", "etbm2", "atbm"
  #                P:  number of points to simulate on line processes for "atbm"
  #    atbm_approach:  string, either "approx", "exact" or "nugget"
  #      atbm_nugget:  number, size of nugget for "nugget" approach
  #   atbm_nugget_cv:  boolean, whether to correct variances for "atbm"
  #    etbm_progress:  boolean, set TRUE to get progress reports for "etbm"
  #             slim:  boolean, option to deactivate coordinate transforms and
  #                    error handlers for better performance (coordinates must
  #                    be in correct format and are returned in the same format
  #                    as the input)
  #
  # OUTPUT (if display = FALSE, slim = FALSE)
  #   matrix (4xN) with data points in column format: c(x,y,z,values)
  #   or
  #   vector (N) with simulated values in the same order as input
  #   or
  #   matrix ((D+1)xN) with simulated KL-components (independent of model)
  #-----------------------------------------------------------------------------
  
  # --------------------------- PREPARATIONS -----------------------------------
  if (!is.null(seed)) set.seed(seed)
  N <- dim(pts)[2]
  
  # set components flag
  if (out_format == "components") {
    components_out = TRUE
  } else components_out = FALSE
  
  # if no method parameters provided: set default values
  if (method == "ls" && is.null(L)) L <- 30
  else if (method == "etbm") {
    if (is.null(L)) L <- 30
    if (is.null(Q)) Q <- 30
  } else if (method == "etbm2") {
    if (is.null(L)) L <- 1
    if (is.null(Q)) Q <- 30
  } else if (method == "atbm" && is.null(Q)) Q <- 30
  
  # if no model is provided: set up default model
  if (is.null(grfmodel)) {
    if (method == "md" || method == "atbm") grfmodel <- GRFModel()
    else grfmodel <- GRFModel(modeltype = "spectral")
  }
  
  # compute correct coordinates 
  if (!slim) {
    if (in_polar_coord) {
      if (dim(pts)[1] == 2) pts <- rbind(pts, rep(1, N))
      polar_pts <- pts
      cart_pts <- Spher2Cart(pts)
    } else if (method == "kle" || method == "ls") {
      cart_pts <- pts
      polar_pts <- Cart2Spher(pts)
    } else cart_pts <- pts
  } else if (method == "kle" || method == "ls") {
    polar_pts <- pts
    cart_pts <- pts # still polar coordinates, no transform for slim mode
  } else {
    cart_pts <- pts
  }
  
  # ------------------------- ERROR HANDLING -----------------------------------
  if (!slim) {
    # PROBLEM 1: points not on unit sphere
    cosums <- colSums(cart_pts**2)
    if (sum(abs(cosums - rep(1, N)) > 1e-15) > 0) {
      stop("Error: Input points not on unit sphere.")
    }
    
    # PROBLEM 2: wrong grfmodel for required method
    if ((grfmodel$modeltype == "covar" && method != "md" && method != "atbm") ||
        (grfmodel$modeltype == "spectral" && method == "md") ||
        (grfmodel$modeltype == "spectral" && method == "atbm")) {
      warning("Warning: Incompatible grfmodel provided. Using default model.")
      if (method %in% c("md", "atbm")) grfmodel <- GRFModel()
      else grfmodel <- GRFModel("spectral")
    }
    
    # PROBLEM 3: wrong method for components output
    if (components_out && method %in% c("md", "ls", "atbm")) {
      warning("Warning: Invalid output format. Returning default format.")
      components_out <- FALSE
      out_format <- "default"
    }
    
    # PROBLEM 4: unsuitable model for ATBM
    if (method == "atbm") {
      if (TrueCov(pi, grfmodel) > 0.01) {
        warning("Warning: Unsuitable model for ATBM method provided.")
      }
    }
    
    # PROBLEM 5: unknown method
    if (!(method %in% c("md", "ls", "kle", "etbm", "etbm2", "atbm"))) {
      warning("Warning: Method unknown. Switching to default method.")
      method <- "md"
    }
  }
  # ----------------------- DISPLAYING INFORMATION -----------------------------
  if (disp_info) {
    cat("\n---------- GRF SIMULATION ----------  ",
        "\n N =", N, " points", "\n Method: ", method)
    
    # covariance information
    if (method == "md" || method == "atbm") {
      cat("\n Covariance model: ", grfmodel$covmodel,
          " ( p1 =", grfmodel$param1, ", p2 =", grfmodel$param2, ")\n")
    } else {
      cat("\n Covariance model: ", grfmodel$spectmodel, 
          " ( p =", grfmodel$spectparam, ")\n",
          "D =", grfmodel$D, " Schoenberg coefficients\n")
    }
    
    # method information
    if (method == "ls") cat(" Method parameter: L =", L, "\n")
    else if (method %in% c("etbm", "etbm2")) {
      cat(" Method parameters: L =", L, " Q =", Q, "\n")
    } else if (method == "atbm") {
      cat(" Method parameter: Q =", Q, "\n", "ATBM approach: ", atbm_approach)
      if (atbm_nugget != 0) cat(" \n Nugget =", atbm_nugget)
      cat("\n")
    }
  }
                     
  # ----------------------- MAIN COMPUTATIONS ----------------------------------
  
  # switch to correct method for simulating and pass arguments
  if (method == "md") {
    out <- .SimGRF_md(pts = cart_pts, grfmodel = grfmodel)
  } else if (method == "kle") {
    out <- .SimGRF_kle(pts = polar_pts,
                       grfmodel = grfmodel,
                       components_out = components_out,
                       precomputations = precomputations)
  } else if (method == "ls") {
    out <- .SimGRF_ls(pts = polar_pts,
                      grfmodel = grfmodel,
                      L = L,
                      precomputations = precomputations)
  } else if (method == "etbm") {
    out <- .SimGRF_etbm(pts = cart_pts,
                        grfmodel = grfmodel,
                        L = L,
                        Q = Q,
                        disp_progress = etbm_progress,
                        components_out = components_out)
  } else if (method == "etbm2") {
    out <- .SimGRF_etbm2(pts = cart_pts,
                         grfmodel = grfmodel,
                         L = L,
                         Q = Q,
                         disp_progress = etbm_progress,
                         components_out = components_out,
                         precomputations = precomputations)
  } else if (method == "atbm") {
    out <- .SimGRF_atbm(pts = cart_pts,
                        grfmodel = grfmodel,
                        Q = Q,
                        approach = atbm_approach,
                        nugget = atbm_nugget,
                        P = P,
                        correct_variance = atbm_nugget_cv,
                        precomputations = precomputations,
                        silent = !disp_info)
  }
  
  if (disp_info) cat("------ END OF GRF SIMULATION -------  \n\n")
  
  
  # ------------------ DISPLAYING POINTS ON SPHERE -----------------------------
  if (display) {
    
    # for older rgl-versions set transparency to 0
    if (as.numeric(substr(packageVersion("rgl"), 3, 50)) < 109) {
      sphere_transparency <- 0
    } else sphere_transparency <- 0.2
    
    source("R/simGRFS/displaysphere.R")    # source for updated resolutions
    if (rgl::cur3d() != 0) rgl::close3d()  # close active device
    
    DisplaySpherePoints(cart_pts,
                        values = out,
                        sphere_transparency = sphere_transparency,
                        title = "")
  }
  
  # ---------------------------- RETURN ----------------------------------------
  else {
    if (out_format == "values" || out_format == "components") return(out)
    else return(rbind(cart_pts, t(out)))
  }
}



RandomSimGRF <- function(N = 20,
                         seed = NULL,
                         grfmodel = NULL,
                         method = "md",
                         display = FALSE,
                         out_format = "default",
                         disp_info = TRUE,
                         L = NULL,                   # ls, etbm
                         Q = NULL,                   # etbm, etbm2, atbm
                         P = NULL,                   # atbm
                         atbm_approach = "approx",   # atbm
                         atbm_nugget = 0,            # atbm
                         atbm_nugget_cv = TRUE,      # atbm
                         etbm_progress = TRUE) {     # etbm, etbm2
  #-----------------------------------------------------------------------------
  # -Wrapper- Simulate a GRF on the sphere using random points on the sphere.
  #
  # DEPENDENCIES
  #  - RandomSpherePoints() from "helpfunctions.R"
  #  - SimGRF()
  #
  # INPUT
  #        N:  number of points to be simulated
  #   ....     (same as for SimGRF, but without pts, precomputations and
  #             in_polar_coord) 
  #
  # OUTPUT (if display = FALSE)
  #   matrix (4xN) with data points in column format: c(x,y,z,values)
  #   or
  #   vector (N) with simulated values in the same order as input
  #   or
  #   matrix ((D+1)xN) with simulated KL-components (independent of model)
  #-----------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  
  # Error handler for unknown methods
  if (!(method %in% c("md", "kle", "ls", "etbm", "etbm2", "atbm"))) {
    warning("Warning: Method unknown. Switching to default method.")
    method <- "md"
  }
  
  # distinguish methods for simulating directly in correct coordinates
  if (method == "md") {
    pts <- RandomSpherePoints(N)
    out <- SimGRF(pts = pts,
                  grfmodel = grfmodel,
                  method = "md",
                  display = display,
                  out_format = out_format,
                  disp_info = disp_info)
  } else if (method == "kle") {
    pts <- RandomSpherePoints(N, out_polar = TRUE)
    out <- SimGRF(pts = pts,
                  grfmodel = grfmodel,
                  method = "kle",
                  in_polar_coord = TRUE,
                  display = display,
                  out_format = out_format,
                  disp_info = disp_info)
  } else if (method == "ls") {
    pts <- RandomSpherePoints(N, out_polar = TRUE)
    out <- SimGRF(pts = pts,
                  grfmodel = grfmodel,
                  method = "ls",
                  in_polar_coord = TRUE,
                  display = display,
                  L = L,
                  out_format = out_format,
                  disp_info = disp_info)
  } else if (method == "etbm") {
    pts <- RandomSpherePoints(N, out_polar = FALSE)
    out <- SimGRF(pts = pts,
                  grfmodel = grfmodel,
                  method = "etbm",
                  in_polar_coord = FALSE, 
                  display = display,
                  L = L,
                  Q = Q,
                  out_format = out_format,
                  disp_info = disp_info,
                  etbm_progress = etbm_progress)
  } else if (method == "etbm2") {
    pts <- RandomSpherePoints(N, out_polar = FALSE)
    out <- SimGRF(pts = pts,
                  grfmodel = grfmodel,
                  method = "etbm2",
                  in_polar_coord = FALSE,
                  display = display,
                  L = L,
                  Q = Q,
                  out_format = out_format,
                  disp_info = disp_info,
                  etbm_progress = etbm_progress)
  } else if (method == "atbm") {
    pts <- RandomSpherePoints(N, out_polar = FALSE)
    out <- SimGRF(pts = pts,
                  grfmodel = grfmodel,
                  method = "atbm",
                  in_polar_coord = FALSE,
                  display = display,
                  Q = Q,
                  atbm_approach = atbm_approach,
                  atbm_nugget = atbm_nugget, 
                  atbm_nugget_cv = atbm_nugget_cv,
                  P = P,
                  out_format = out_format,
                  disp_info = disp_info,
                  etbm_progress = etbm_progress)
  } 
  
  return(out)
}



MultiSimGRF <- function(pts,
                        grfmodel = NULL,
                        method = "md",
                        M = 50,
                        precomputations = NULL,
                        parallel = FALSE,
                        in_polar_coord = FALSE,
                        disp_info = TRUE,
                        disp_progress = TRUE,
                        seed = NULL,
                        L = NULL,                   # ls, etbm
                        Q = NULL,                   # etbm, etbm2, atbm
                        P = NULL,                   # atbm
                        atbm_approach = "approx",   # atbm
                        atbm_nugget = 0,            # atbm
                        atbm_nugget_cv = TRUE,      # atbm
                        slim = FALSE) {
  #-----------------------------------------------------------------------------
  # Simulates multiple realizations of a GRF using optimized methods.
  # 
  # DEPENDENCIES
  #  - GRFModel() from "grfmodels.R"
  #  - functions from "helpfunctions.R"
  #  - functions from "method1md.R", "method2kle.R", "method3ls.R"
  #                   "method4etbm.R", "method5etbm2.R" and "method6atbm.R"
  #  - package "parallel" for parallel computation
  #
  # INPUT
  #             M:  number of realizations of GRF to simulate
  #      parallel:  boolean, whether to use parallel computations
  # disp_progress:  get simulation progress (not working for parallel option)
  #     ....   (same as SimGRF, without display and out_format) 
  #
  # OUTPUT
  #   matrix (MxN) with rows the realizations of simulated values
  #-----------------------------------------------------------------------------
  
  # --------------------------- PREPARATIONS -----------------------------------
  if (!is.null(seed)) set.seed(seed)
  N <- dim(pts)[2]
  
  # prepare progress reports
  if (disp_progress) progresslist <- floor(seq(0, 100, length = M))
  
  # prepare for parallelization
  if (parallel) num_cores <- parallel::detectCores()
  
  # if no method parameters provided: set default values
  if (method == "ls" && is.null(L)) L <- 30
  else if (method == "etbm") {
    if (is.null(L)) L <- 30
    if (is.null(Q)) Q <- 30
  } else if (method == "etbm2") {
    if (is.null(L)) L <- 1
    if (is.null(Q)) Q <- 30
  } else if (method == "atbm" && is.null(Q)) Q <- 30
  
  # if no grfmodel is provided: set default model
  if (is.null(grfmodel)) {
    if (method == "md" || method == "atbm") grfmodel <- GRFModel()
    else grfmodel <- GRFModel(modeltype = "spectral")
  }
  
  # compute both coordinate types
  if (!slim) {
    if (in_polar_coord) {
      if (dim(pts)[1] == 2) pts <- rbind(pts, rep(1, N))
      polar_pts <- pts
      cart_pts <- Spher2Cart(pts)
    } else {
      cart_pts <- pts
      polar_pts <- Cart2Spher(pts)
    }
  } else if (method == "kle" || method == "ls") {
    polar_pts <- pts
    cart_pts <- pts
  } else {
    cart_pts <- pts
  }
  
  # set up output matrix
  out <- matrix(numeric(M * N), M, N)
  
  # ------------------------- ERROR HANDLING -----------------------------------
  if (!slim) {
    # PROBLEM 1: points not on unit sphere
    cosums <- colSums(cart_pts**2)
    if (sum(abs(cosums - rep(1, N)) > 1e-15) > 0) {
      stop("Error: Input points not on unit sphere.")
    }
    
    # PROBLEM 2: wrong grfmodel for required method
    if ((grfmodel$modeltype == "covar" && method != "md" && method != "atbm") ||
        (grfmodel$modeltype == "spectral" && method == "md") ||
        (grfmodel$modeltype == "spectral" && method == "atbm")) {
      warning("Warning: Incompatible grfmodel provided. Using default model.")
      if (method %in% c("md", "atbm")) grfmodel <- GRFModel()
      else grfmodel <- GRFModel("spectral")
    }

    # PROBLEM 4: unsuitable model for ATBM
    if (method == "atbm") {
      if (TrueCov(pi, grfmodel) > 0.01) {
        warning("Warning: Unsuitable model for ATBM method provided.")
      }
    }
    
    # PROBLEM 5: unknown method
    if (!(method %in% c("md", "ls", "kle", "etbm", "etbm2", "atbm"))) {
      warning("Warning: Method unknown. Switching to default method.")
      method <- "md"
    }
  }
  
  # ----------------------- DISPLAYING INFORMATION -----------------------------
  if (disp_info) {
    cat("\n---------- GRF SIMULATION ----------  ",
        "\n M =", M, " simulations",
        "\n N =", N, " points", "\n Method: ", method)
    
    # covariance information
    if (method == "md" || method == "atbm") {
      cat("\n Covariance model: ", grfmodel$covmodel,
          " ( p1 =", grfmodel$param1, ", p2 =", grfmodel$param2, ")\n")
    } else {
      cat("\n Covariance model: ", grfmodel$spectmodel, 
          " ( p =", grfmodel$spectparam, ")\n",
          "D =", grfmodel$D, " Schoenberg coefficients\n")
    }
    
    # method information
    if (method == "ls") cat(" Method parameter: L =", L, "\n")
    else if (method %in% c("etbm", "etbm2")) {
      cat(" Method parameters: L =", L, " Q =", Q, "\n")
    } else if (method == "atbm") {
      cat(" Method parameter: Q =", Q, "\n", "ATBM approach: ", atbm_approach)
      if (atbm_nugget != 0) cat(" \n Nugget =", atbm_nugget)
      cat("\n")
    }
  }
  
  # ----------------------- MAIN COMPUTATIONS ----------------------------------
  # ----------------------  Method 1: MD  --------------------------------------
  if (method == "md") {
    
    # precompute cholesky factor
    if (disp_progress) cat("\r Running precomputations ...")
    if (!is.null(precomputations)) cholesky_factor <- precomputations
    else cholesky_factor <- Precompute_md(pts = pts, grfmodel = grfmodel)
    
    # parallel computations
    if (parallel) {
      out_list <- parallel::mclapply(numeric(M), .SimGRF_md,
                           pts = cart_pts,
                           precomputations = cholesky_factor,
                           mc.cores = num_cores)
      if (disp_info) cat(" Parallel simulation")
    }
    
    # normal computations 
    else {
      for (m in 1:M) {
        out[m, ] <- .SimGRF_md(pts = cart_pts,
                               precomputations = cholesky_factor)
        
        # display progress (and eta for long computations)
        if (disp_progress) {
          if  (m == 1) {
            eta <- "..."
            tic <- NA   
            times <- rep(NA, M)
            eta_flag <- FALSE
          }
          if (!is.na(tic)) {
            times[m] <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
            eta <- round(mean(times, na.rm = TRUE) * (M - m))
            if (m == 2 && eta > 30) eta_flag <- TRUE
          }
          tic <- Sys.time()
          if ((m - 1) %% ceiling(M / 100) == 0 &&
              !is.na(eta) && eta_flag == TRUE) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100), "% ",
                "(estimated time remaining: ", eta,
                " seconds)                             ")
          } else if ((m - 1) %% ceiling(M / 100) == 0) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100),
                "%                                     ")
          }
        }
      }
    }
  }
  
  # ----------------------  Method 2: KLE  -------------------------------------
  else if (method == "kle") {
    
    # precompute spherical harmonics
    if (disp_progress) cat("\r Running precomputations ...")
    if (!is.null(precomputations)) spher_harm <- precomputations
    else spher_harm <- Precompute_kle(pts = polar_pts, D = grfmodel$D)
    
    # parallel computations
    if (parallel) {
      out_list <- parallel::mclapply(numeric(M), .SimGRF_kle,
                                     pts = polar_pts,
                                     grfmodel = grfmodel,
                                     precomputations = spher_harm,
                                     mc.cores = num_cores)
      if (disp_info) cat(" Parallel simulation")
    }
    
    # normal computations
    else {
      for (m in 1:M) {
        out[m, ] <- .SimGRF_kle(pts = polar_pts,
                                grfmodel = grfmodel,
                                precomputations = spher_harm)
        
        # display progress (and eta for long computations)
        if (disp_progress) {
          if  (m == 1) {
            eta <- "..."
            tic <- NA   
            times <- rep(NA, M)
            eta_flag <- FALSE
          }
          if (!is.na(tic)) {
            times[m] <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
            eta <- round(mean(times, na.rm = TRUE) * (M - m))
            if (m == 2 && eta > 30) eta_flag <- TRUE
          }
          tic <- Sys.time()
          if ((m - 1) %% ceiling(M / 100) == 0 &&
              !is.na(eta) && eta_flag == TRUE) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100), "% ",
                "(estimated time remaining: ", eta,
                " seconds)                             ")
          } else if ((m - 1) %% ceiling(M / 100) == 0) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100),
                "%                                     ")
          }
        }
      }
    }
  }
  
  # -------------------------  Method 3: LS  -----------------------------------
  else if (method == "ls") {
    
    # precompute spherical harmonics
    if (disp_progress) cat("\r Running precomputations ...")
    if (!is.null(precomputations)) spher_harm <- precomputations
    else spher_harm <- Precompute_ls(pts = polar_pts, D = grfmodel$D)
    
    # parallel computations
    if (parallel) {
      out_list <- parallel::mclapply(numeric(M), .SimGRF_ls,
                                     pts = polar_pts,
                                     grfmodel = grfmodel,
                                     L = L,
                                     precomputations = spher_harm,
                                     mc.cores = num_cores)
      if (disp_info) cat(" Parallel simulation")
    }
    
    # normal computations
    else {
      for (m in 1:M) {
        out[m, ] <- .SimGRF_ls(pts = polar_pts,
                               grfmodel = grfmodel,
                               L = L,
                               precomputations = spher_harm)
        
        # display progress (and eta for long computations)
        if (disp_progress) {
          if  (m == 1) {
            eta <- "..."
            tic <- NA   
            times <- rep(NA, M)
            eta_flag <- FALSE
          }
          if (!is.na(tic)) {
            times[m] <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
            eta <- round(mean(times, na.rm = TRUE) * (M - m))
            if (m == 2 && eta > 30) eta_flag <- TRUE
          }
          tic <- Sys.time()
          if ((m - 1) %% ceiling(M / 100) == 0 &&
              !is.na(eta) && eta_flag == TRUE) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100), "% ",
                "(estimated time remaining: ", eta,
                " seconds)                             ")
          } else if ((m - 1) %% ceiling(M / 100) == 0) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100),
                "%                                     ")
          }
        }
      }
    }
  }
  
  # ----------------------  Method 4: ETBM  ------------------------------------
  else if (method == "etbm") {
    
    # parallel computations
    if (parallel) {
      out_list <- parallel::mclapply(numeric(M), .SimGRF_etbm,
                           pts = cart_pts,
                           grfmodel = grfmodel,
                           L = L,
                           Q = Q,
                           disp_progress = FALSE,
                           mc.cores = num_cores)
      if (disp_info) cat(" Parallel simulation")
    }
    
    # normal computations
    else {
      for (m in 1:M) {
        out[m, ] <- .SimGRF_etbm(pts = cart_pts,
                                 grfmodel = grfmodel,
                                 L = L,
                                 Q = Q,
                                 disp_progress = FALSE)
        
        # display progress (and eta for long computations)
        if (disp_progress) {
          if  (m == 1) {
            eta <- "..."
            tic <- NA   
            times <- rep(NA, M)
            eta_flag <- FALSE
          }
          if (!is.na(tic)) {
            times[m] <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
            eta <- round(mean(times, na.rm = TRUE) * (M - m))
            if (m == 2 && eta > 30) eta_flag <- TRUE
          }
          tic <- Sys.time()
          if ((m - 1) %% ceiling(M / 100) == 0 &&
              !is.na(eta) && eta_flag == TRUE) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100), "% ",
                "(estimated time remaining: ", eta,
                " seconds)                             ")
          } else if ((m - 1) %% ceiling(M / 100) == 0) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100),
                "%                                     ")
          }
        }
      }
    }
  }
  
  # ----------------------  Method 5: ETBM2  -----------------------------------
  else if (method == "etbm2") {
    
    # precomputations
    if (disp_progress) cat("\r Running precomputations ...")
    if (!is.null(precomputations)) a_array <- precomputations
    else a_array <- Precompute_etbm(Q = Q, D = grfmodel$D)
    
    # parallel computations
    if (parallel) {
      out_list <- parallel::mclapply(numeric(M), .SimGRF_etbm2,
                                     pts = cart_pts,
                                     grfmodel = grfmodel,
                                     L = L,
                                     Q = Q,
                                     disp_progress = FALSE,
                                     precomputations = a_array,
                                     mc.cores = num_cores)
      if (disp_info) cat(" Parallel simulation")
    }
    
    # normal computations
    else {
      for (m in 1:M) {
        out[m, ] <- .SimGRF_etbm2(pts = cart_pts,
                                  grfmodel = grfmodel,
                                  L = L,
                                  Q = Q,
                                  disp_progress = FALSE,
                                  precomputations = a_array)
        
        # display progress (and eta for long computations)
        if (disp_progress) {
          if  (m == 1) {
            eta <- "..."
            tic <- NA   
            times <- rep(NA, M)
            eta_flag <- FALSE
          }
          if (!is.na(tic)) {
            times[m] <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
            eta <- round(mean(times, na.rm = TRUE) * (M - m))
            if (m == 2 && eta > 30) eta_flag <- TRUE
          }
          tic <- Sys.time()
          if ((m - 1) %% ceiling(M / 100) == 0 &&
              !is.na(eta) && eta_flag == TRUE) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100), "% ",
                "(estimated time remaining: ", eta,
                " seconds)                             ")
          } else if ((m - 1) %% ceiling(M / 100) == 0) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100),
                "%                                     ")
          }
        }
      }
    }
  }
  
  # ----------------------  Method 6: ATBM  ------------------------------------
  else if (method == "atbm") {
    
    # precompute line process values (simulation)
    if (disp_progress) cat("\r Running precomputations ...")
    if (!is.null(precomputations)) cholesky_factor <- precomputations
    else {
      
      # compute P (number of points to simulate on lines)
      if (atbm_approach == "exact") {
        md <- MinDistSphere(pts)
        approx_interval_length <- sqrt((md**2) / 2)
      } else approx_interval_length <- exp(-log(N, base = 10))
      if (is.null(P)) P <- ceiling(2 / approx_interval_length)
      
      # simulate line processes
      cholesky_factor <- Precompute_atbm(grfmodel = grfmodel, P = P)
      if (disp_info) cat("\r P =", P, "                     \n") 
    }
    
    # parallel computations
    if (parallel) {
      out_list <- parallel::mclapply(numeric(M), .SimGRF_atbm,
                                     pts = cart_pts,
                                     grfmodel = grfmodel,
                                     Q = Q,
                                     approach = atbm_approach,
                                     P = P,
                                     nugget = atbm_nugget,
                                     correct_variance = atbm_nugget_cv,    
                                     precomputations = cholesky_factor,
                                     silent = TRUE,
                                     mc.cores = num_cores)
      if (disp_info) cat(" Parallel simulation")
    }
    
    # normal computations
    else {
      for (m in 1:M) {
        out[m, ] <- .SimGRF_atbm(pts = cart_pts,
                                 grfmodel = grfmodel,
                                 Q = Q,
                                 approach = atbm_approach,
                                 nugget = atbm_nugget,
                                 correct_variance = atbm_nugget_cv,
                                 P = P,
                                 precomputations = cholesky_factor,
                                 silent = TRUE)
        
        # display progress (and eta for long computations)
        if (disp_progress) {
          if  (m == 1) {
            eta <- "..."
            tic <- NA   
            times <- rep(NA, M)
            eta_flag <- FALSE
          }
          if (!is.na(tic)) {
            times[m] <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
            eta <- round(mean(times, na.rm = TRUE) * (M - m))
            if (m == 2 && eta > 30) eta_flag <- TRUE
          }
          tic <- Sys.time()
          if ((m - 1) %% ceiling(M / 100) == 0 &&
              !is.na(eta) && eta_flag == TRUE) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100), "% ",
                "(estimated time remaining: ", eta,
                " seconds)                             ")
          } else if ((m - 1) %% ceiling(M / 100) == 0) {
            cat("\r Simulation Progress: ", round((m - 1) / M * 100),
                "%                                     ")
          }
        }
      }
    }
  }
  
  # ---------------------------- RETURN ----------------------------------------
  if (disp_progress) cat("\r Simulation Progress: 100%              ",
                         "                                        \n")
  if (disp_info) cat("------ END OF GRF SIMULATION -------  \n\n")
  
  # if parallel computed: convert list to matrix
  if (parallel) out <- matrix(unlist(out_list), nrow = M)
  
  return(out)
}



EvalGRFComponents <- function(components,
                              grfmodel = NULL,
                              coefs = NULL,
                              pts = NULL,
                              in_polar_coord = FALSE,
                              out_format = "values",
                              display = FALSE) {
  #-----------------------------------------------------------------------------
  # Evaluates KL-components for given model (Schoenberg coefficients).
  #
  # Model can be provided either via grfmodel argument or by directly
  # providing a vector of Schoenberg coefficients.
  #
  # Points can be optionally provided to return default format output,
  # or to use the option to display the values on the sphere.
  #
  # DEPENDENCIES
  #  - Spher2Cart() from "helpfunctions.R"
  #  - GRFModel() from "grfmodels.R"
  #  - package "rgl" for displaying options
  #  - functions from "displaysphere.R" for displaying options
  # 
  # INPUT
  #     components:  (D+1)xN matrix of KL-components
  #       grfmodel:  instance of class GRFmodel, option to provide cov model.
  #          coefs:  (D+1) vector of Schoenberg coefs, option to prov. cov.
  #            pts:  optional argument to provide points in all formats
  #                  (necessary for using display option or "default" format)
  # in_polar_coord:  boolean, set TRUE if input is given in poler coordinates
  #     out_format:  string, either "values" or "default", as in main function
  #                  option "default" requires input points
  #        display:  boolean, whether to display the results, requires points
  #
  # OUTPUT (if display=FALSE)
  #   matrix (4xN) with data points in column format: c(x,y,z,values)
  #   or
  #   vector (N) with simulated values in the same order as input
  #-----------------------------------------------------------------------------
  N <- dim(pts)[2]
  
  # if points are provided, compute correct coordinates
  if (!is.null(pts)) {
    if (in_polar_coord) {
      if (dim(pts)[1] == 2) pts <- rbind(pts, rep(1, N))
      cart_pts <- Spher2Cart(pts)
    } else cart_pts <- pts
  }
  
  # if no model specified: get default model coefs
  if (is.null(grfmodel) && is.null(coefs)) {
    coefs <- GRFModel("spectral")$coefs
  } else if (!is.null(grfmodel) && is.null(coefs)) {
    coefs <- grfmodel$coefs
  } else if (!is.null(grfmodel) && !is.null(coefs) &&
           grfmodel$coefs != coefs) {
    warning("Error: coefs provided do not coincide with provided model")
  }
  
  # main computation
  out <- as.vector(sqrt(coefs) %*% components)
  
  # display if requested (and points provided)
  if (display && !is.null(pts)) {
    source("R/simGRFS/displaysphere.R")         # source for updated resolution
    if (rgl::cur3d() != 0) rgl::close3d()
    DisplaySpherePoints(cart_pts, values = out)
  } else if (display && is.null(pts)) {
    warning("No points provided. Values cannot be displayed.")
  }
  
  # else return output
  else {
    if (out_format == "default" && !is.null(pts)) {
      return(rbind(cart_pts, t(out)))
    } else return(out) 
  }
}
  