#------------------------------------------------------------------------------- 
#          Helping functions for simulating GRFs on the sphere
#------------------------------------------------------------------------------- 
# File:    helpfunctions.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# DEPENDENCIES
#  none
# 
# CONTENT
# point generators:
#   RandomSpherePoints()   sample random points on the sphere
#   UniformSpherePoints()  generate points on the sphere uniformly
#   GridSpherePoints()     generate points in a regular grid on the sphere
# 
# coordinate transforms:
#   Spher2Cart()           convert spherical polar coord. to cartesian coord
#   Cart2Spher()           convert cartesian coord. to spherical polar coord
#   Lat2Colat()            converts spherical coordinates to globe coordinates
#
# distance transforms:
#   Geod2Eucl()            converts geodesic distance to euclidean distance
#   Eucl2Geod()            converts euclidean distance to geodesic distance 
#
# further helper functions:
#   MinDistSphere()        computes min. distance between points on the sphere
#   .MinDistSphere_small() Wrapper for dist(), only works for few points
#   Symmetrize()           takes an upper triangular matrix and adds lower tri
#   .CleanMatrix()         cleans a matrix for potentially numerical errors
#-------------------------------------------------------------------------------



RandomSpherePoints <- function(N,
                               radius = 1,
                               out_polar = FALSE,
                               no_poles = FALSE,
                               seed = NULL) {
  #-----------------------------------------------------------------------------
  # Sample random points on the sphere 
  #
  # DEPENDENCIES
  # - Cart2Spher()
  # 
  # INPUT
  #          N:  number of points
  #     radius:  radius of sphere    
  #  out_polar:  if output should be in polar coordinates, set TRUE
  #   no_poles:  if TRUE, samples no points very close to poles
  #       seed:  shortcut to set seed
  #
  # OUTPUT
  #   matrix (3xN) with data points in columns
  #          in format c(x,y,z)
  #       or in format c(lat,long,r) using out_polar=TRUE
  #-----------------------------------------------------------------------------
  
  # set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # sample gaussians
  rpts <- matrix(rnorm(N * 3), 3, N)
  
  # project on sphere
  if (N > 10000) {
    out <- radius * (sweep(rpts, 2, 1 / sqrt(colSums(rpts**2)), FUN = "*"))
  } else {
    fa <- 1 / sqrt(colSums(rpts**2))
    out <- radius * (rpts * rep(fa, rep(3, N)))
  }
  
  # if required, recompute points too close to poles (recursively)
  if (no_poles) {
    npols <- length(out[abs(out[3, ]) > 0.999999])
    if (npols > 0) {
      out[abs(out[3, ]) > 0.999999] <- RandomSpherePoints(npols,
                                                          no_poles = TRUE)
    }
  }
  
  # if required, convert to polar coordinates, then return
  if (out_polar) return(Cart2Spher(out))
  else return(out)
}



UniformSpherePoints <- function(N, radius = 1, out_polar = FALSE) {
  #-----------------------------------------------------------------------------
  # Generate uniformly distributed points on the sphere (not random!)
  # (According to fibonacci method)
  # 
  # (Method by user Fnord, stackoverflow.com/questions/9600801/evenly-
  # distributing-n-points-on-a-sphere)
  #
  # DEPENDENCIES
  # - Cart2Spher()
  #
  # INPUT
  #         N:  number of points
  #    radius:  radius of sphere
  # out_polar:  if output should be in polar coordinates, set TRUE
  #
  # OUTPUT
  #   matrix (3xN) with data points in columns
  #          in format c(x,y,z)
  #       or in format c(lat,long,r) using out_polar=TRUE
  #-----------------------------------------------------------------------------
  
  # setup output matrix
  out <- matrix(numeric(3 * N), 3, N)
  
  # case N = 1 separately
  if (N == 1) out[ ,1] <- c(0, 1, 0)
  else {
    
    # constant golden ratio
    golden <- pi * (3 - sqrt(5))
    
    # loop over all points
    for (n in 0:(N - 1)) {
      y  <- 1 - 2 * (n / (N - 1))
      fa <- sqrt(1 - y**2)
      theta <- golden * n
      x <- cos(theta) * fa
      z <- sin(theta) * fa
      out[ ,n + 1] <- radius * c(x, y, z)
    }
  }
  
  # if required, convert to polar coordinates, then return
  if (out_polar) return(Cart2Spher(out))
  else return(out)
}



GridSpherePoints <- function(N = 100,
                             lat_span = 0.5,
                             long_span = 0.5,
                             location = c(0, -1, 0)) {
  #-----------------------------------------------------------------------------
  # Creates a regular NxN-grid on the sphere at location (c(x,y,z)), using 
  # the distances specified by lat_span and long_span (both arclengths).
  #
  # DEPENDENCIES
  #  - Cart2Spher()
  #  - Spher2Cart()
  #
  # INPUT
  #         N: number of points per direction
  #  lat_span: maximal latitude distance of grid
  # long_span: maximal longtitude distance of grid
  #  location: center points of the grid in col format c(x,y,z)
  #
  # OUTPUT
  #   matrix (3xN^2) with grid points in column format c(x,y,z)
  #-----------------------------------------------------------------------------
  
  # transform location to polar coordinates and extract components
  location_pol <- Cart2Spher(location)
  locx <- location_pol[1]
  locy <- location_pol[2]
  
  # set up components
  z <- rep(1, N**2)
  x <- seq(locx + lat_span / 2, locx - lat_span / 2, len = N)
  y <- seq(locy + long_span / 2, locy - long_span / 2, len = N)
  
  return(Spher2Cart(t(cbind(as.matrix(expand.grid(x = x, y = y)), z))))
}



Spher2Cart <- function(pts_polar){
  #-----------------------------------------------------------------------------
  # Convert spherical polar coordinates to cartesian coordinates.
  # 
  # If only lat and long is provided, the radius is assumed to be 1.
  #
  # No dependencies.
  #
  # INPUT
  # pts_polar:    matrix (3xN) with data points in column format c(lat,long,r)
  #            or matrix (2xN) with data points in column format c(lat,long)
  #
  # OUTPUT
  #   matrix (3xN) with data points in column format c(x,y,z)
  #-----------------------------------------------------------------------------
  
  # handle (incorrect) input as vector
  if (is.vector(pts_polar)) {
    pts_polar <- matrix(pts_polar, length(pts_polar), 1)
  }
  
  N <- dim(pts_polar)[2]
  
  # set up output
  out <- matrix(nrow = 3, ncol = N)
  
  # handle case where only two variables provided 
  if (dim(pts_polar)[1] == 2) pts_polar <- rbind(pts_polar, rep(1, N))
  
  # loop over all points and convert 
  for (n in 1:N) {
    out[1,n] <- pts_polar[3,n] * cos(pts_polar[1,n]) * cos(pts_polar[2,n])
    out[2,n] <- pts_polar[3,n] * cos(pts_polar[1,n]) * sin(pts_polar[2,n])
    out[3,n] <- pts_polar[3,n] * sin(pts_polar[1,n])
  }
  
  return(out)
}



Cart2Spher <- function(pts_cart) {
  #-----------------------------------------------------------------------------
  # Convert cartesian coordinates to spherical polar coordinates.
  #
  # No dependencies.
  #
  # INPUT
  #  pts_cart:  matrix (3xN) with data points in column format c(x,y,z)
  #
  # OUTPUT
  #   matrix (3xN) with data points in column format c(lat,long,r)
  #-----------------------------------------------------------------------------

  # handle (incorrect) input as vector
  if (is.vector(pts_cart)) pts_cart <- matrix(pts_cart, length(pts_cart), 1)
  
  # set up output
  N <- dim(pts_cart)[2]
  out <- matrix(nrow = 3, ncol = N)

  # loop over all points and convert 
  for (n in 1:N) {
    
    # special case c(0,0,0)
    if (all(pts_cart[ ,n] == c(0, 0, 0))) {
      out[ ,n] <- c(0, 0, 0)
      next()
    }
    
    # regular case
    out[3,n] <- sqrt(sum(pts_cart[ ,n]**2))
    out[1,n] <- asin(pts_cart[3,n] / out[3,n])
    out[2,n] <- atan2(pts_cart[2,n], pts_cart[1,n])
  }
  
  return(out)
}



Lat2Colat <- function(pts_polar) {
  #-----------------------------------------------------------------------------
  # Takes polar coordinates using latitudes and converts to polar coordinates
  # using colatitudes.
  #
  # No dependencies.
  #
  # INPUT
  #  pts_polar:  (3xN) matrix of points in format c(lat,long,r)
  #
  # OUTPUT
  #  (3xN) matrix of points in format c(colat,long,r)
  #-----------------------------------------------------------------------------
  
  return(matrix(rep(c(pi / 2, 0, 0), dim(pts_polar)[2]), dim(pts_polar)) +
         diag(c(-1, 1, 1)) %*% pts_polar)
}



Geod2Eucl <- function(d) {
  #-----------------------------------------------------------------------------
  # Takes in a geodesic distance(s) of points on the sphere and returns euclidean
  # distance.
  #
  # No dependencies.
  #
  # INPUT
  #  d:  geodesic distance(s)
  #
  # OUTPUT
  #  number or vector of euclidean distances
  #-----------------------------------------------------------------------------
  return(2 * asin(d / 2))
}



Eucl2Geod <- function(d) {
  #-----------------------------------------------------------------------------
  # Takes in a euclidean distance of points on the sphere and returns geodesic
  # distance.
  #
  # No dependencies.
  #
  # INPUT
  #  d:  euclidean distances
  #
  # OUTPUT
  #  number or vector of geodesic distances
  #-----------------------------------------------------------------------------
  return(2 * sin(d / 2))
}



.MinDistSphere_small <- function(pts) {
  #-----------------------------------------------------------------------------
  # Computes minimal distance between a small amount of points on the sphere.
  #
  # The number of points this function can handle depends on the computer,
  # usually 10'000 points should be feasible.
  #
  # DEPENDENCIES
  #  - Eucl2Geod()
  #
  # INPUT
  #  pts: points on the sphere in col format c(x,y,z)
  #
  # OUTPUT
  #  minimal distance between points
  #-----------------------------------------------------------------------------
  
  # handle erroneous inputs
  if (!is.matrix(pts)) return(NA)
  if (dim(pts)[2] == 0) return(NA)

  return(Eucl2Geod(min(dist(t(pts)))))
}



MinDistSphere <- function(pts) {
  #-----------------------------------------------------------------------------
  # Computes the minimal distance between points on the sphere.
  #
  # For small amount of points, the implemented function dist() is used. 
  # However, for larger amount of points, this method leads to way too large
  # matrices. Therefore, a divide and conquer method has been implemented.
  #
  # DEPENDENCIES
  #  - .MinDistSphere_small()
  #
  # INPUT
  #  pts: points on the sphere in col format c(x,y,z)
  #
  # OUTPUT
  #  minimal distance between points
  #-----------------------------------------------------------------------------
  
  #---------------------------- PREPARATIONS -----------------------------------
  # get number of points and apply base function for small number of points
  N <- dim(pts)[2]
  if (N <= 10000) return(.MinDistSphere_small(pts))
  
  # set starting value value
  minim <- 1
  
  # get parameters
  steps <- ceiling(sqrt(N / 500))
  hdelta <- 2 / steps
  overlap <- 10**(-(log(N, base = 10) - 3))

  # get parameters for boundary cases
  h <- 1 - hdelta
  d <- sqrt(4 * (1 - h**2))
  steps2 <- ceiling(steps * d)
  hdelta2 <- d / steps2
  
  #---------------------- MAIN COMPUTATION: bottom boundary --------------------
  pts_sel <- pts[ ,pts[3, ] <= -1 + hdelta]
  pts_group1 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 + hdelta2]
  dist <- .MinDistSphere_small(pts_group1)
  if (dist < minim && !is.na(dist)) minim <- dist
  pts_group2 <- pts_sel[ ,pts_sel[2, ] > d / 2 - hdelta2]
  dist <- .MinDistSphere_small(pts_group2)
  if (dist < minim && !is.na(dist)) minim <- dist
  for (j in 1:(steps2 - 1)) {
    pts_sel2 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 + (j + 1 + overlap) * hdelta2 &
                          pts_sel[2, ] > -d / 2 + (j - overlap) * hdelta2 &
                          pts_sel[1, ] > 0]
    dist <- .MinDistSphere_small(pts_sel2)
    if (dist < minim && !is.na(dist)) minim <- dist
    pts_sel2 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 + (j + 1 + overlap) * hdelta2 &
                          pts_sel[2,] > -d / 2 + (j - overlap) * hdelta2 &
                          pts_sel[1,] < 0]
    dist <- .MinDistSphere_small(pts_sel2)
    if (dist < minim && !is.na(dist)) minim <- dist
  }
  
  #---------------------- MAIN COMPUTATION: top boundary -----------------------
  pts_sel <- pts[ ,pts[3, ] > 1 - hdelta]
  pts_group3 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 + hdelta2]
  dist <- .MinDistSphere_small(pts_group3)
  if (dist < minim && !is.na(dist)) minim <- dist
  pts_group4 <- pts_sel[ ,pts_sel[2, ] > d / 2 - hdelta2]
  dist <- .MinDistSphere_small(pts_group4)
  if (dist < minim && !is.na(dist)) minim <- dist
  for (j in 1:(steps2 - 1)) {
    pts_sel2 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 + (j + 1 + overlap) * hdelta2 &
                          pts_sel[2, ] > -d / 2 + (j - overlap) * hdelta2 &
                          pts_sel[1, ] > 0]
    dist <- .MinDistSphere_small(pts_sel2)
    if (dist < minim && !is.na(dist)) minim <- dist
    pts_sel2 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 + (j + overlap) * hdelta2 &
                          pts_sel[2, ] > -d / 2 + (j - overlap) * hdelta2 &
                          pts_sel[1, ] < 0]
    dist <- .MinDistSphere_small(pts_sel2)
    if (dist < minim && !is.na(dist)) minim <- dist
  }
  
  
  #---------------------- MAIN COMPUTATION: middle sections --------------------
  for (i in 1:(steps - 1)) {
    pts_sel <- pts[ ,pts[3, ] <= -1 + (i + 1 + overlap) * hdelta &
                     pts[3, ] > -1 + (i - overlap) * hdelta]
    h <- max(1 - hdelta * i, -1 - hdelta * i)
    d <- sqrt(4 * (1 - h**2)) #(length of max diameter)
    steps2 <- ceiling(steps * d)
    hdelta2 <- d / steps2
    pts_sel2 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 + hdelta2]
    dist <- .MinDistSphere_small(pts_sel2)
    if (dist < minim && !is.na(dist)) minim <- dist
    pts_sel2 <- pts_sel[ ,pts_sel[2, ] > d / 2 - hdelta2]
    dist <- .MinDistSphere_small(pts_sel2)
    if (dist < minim && !is.na(dist)) minim <- dist
    
    # loop through y-sections
    for (j in 1:(steps2 - 1)) {
      pts_sel2 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 +
                                            (j + 1 + overlap) * hdelta2 &
                            pts_sel[2, ] > -d / 2 + (j - overlap) * hdelta2 &
                            pts_sel[1, ] > 0]
      dist <- .MinDistSphere_small(pts_sel2)
      if (dist < minim && !is.na(dist)) minim <- dist
      pts_sel2 <- pts_sel[ ,pts_sel[2, ] <= -d / 2 +
                                            (j + 1 + overlap) * hdelta2 &
                            pts_sel[2, ] > -d / 2 + (j - overlap) * hdelta2 &
                            pts_sel[1, ] < 0]
      dist <- .MinDistSphere_small(pts_sel2)
      if (dist < minim && !is.na(dist)) minim <- dist
    }
  }
  
  return(minim)
}



Symmetrize <- function(mat) {
  #-----------------------------------------------------------------------------
  # Takes an upper triangular matrix and adds the missing lower triangle part.
  # (No substractions take part, thus there is no cancellation errors.)
  #
  # Method by user3318600
  # (stackoverflow.com/questions/18165320/creating-a-symmetric-matrix-in-r)
  #
  # No dependencies.
  # 
  # INPUT
  #  mat:  upper triangular matrix
  #
  # OUTPUT
  #  symmetric matrix of the same dimensions as input matrix
  #-----------------------------------------------------------------------------
  
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  return(mat)
}



.CleanMatrix <- function(mat, target_values = 0, precision = 1e-15) {
  #-----------------------------------------------------------------------------
  # (DUBIOUS METHOD, APPLY WITH CARE)
  # Takes in a matrix with values that are potentially slightly erroneous
  # due to numerical errors and returns a cleaned matrix for the values
  # that are provided via target_values.
  #
  # No dependencies.
  # 
  # INPUT
  #           mat:  any matrix
  # target_values:  values for which the matrix should be cleaned
  #     precision:  deviation that should be considered numerical error
  #
  # OUTPUT
  #  the same matrix cleaned for values that are virtually target_values
  #-----------------------------------------------------------------------------

  # loop over all target values and replace
  for (tv in target_values) {
    mat[abs(mat - tv) < precision] <- tv
  }
  
  return(mat)
}
