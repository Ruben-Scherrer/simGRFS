#-------------------------------------------------------------------------------
#           Functions for evaluating (associated) Legendre functions 
#------------------------------------------------------------------------------- 
# File:    legendre.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - The following functions evaluate (associated) Legendre functions using
#    the package "gsl".
#  - For ALP functions the keyword "csphase" can be used to indicate whether
#    or not the Condon-Shortley phase factor of (-1)**k should be included.
#  - Numerical instability in computation of unnormalized ALP starts at around
#    degree 150.
#
# DEPENDENCIES
#  - functions from package "gsl"
#
# CONTENT
#   Legendre()        Evaluates Legendre polynomial (LP)
#   Legendres()       Evaluates all LPs up to degree deg
#   AssocLegendre()   Evaluates associated Legendre polynomial (ALP) 
#   AssocLegendres()  Evaluates all ALPs up to some degree
#   .NextLegendre()   Evaluates next LP based on previous two
#-------------------------------------------------------------------------------



Legendre <- function(n, x) {
  #-----------------------------------------------------------------------------
  # Evaluates Legendre polynomial at degree deg.
  #
  # (Wrapper for gsl::legendre_Pl())
  #
  # DEPENDENCIES
  #  - gsl::legendre_Pl()
  #
  # INPUT
  #   n:   degree of Legendre polynomial to be computed 
  #   x:   data point to be evaluated
  #
  # OUPUT
  #  value of Legendre polynomial of degree n at point x
  #-----------------------------------------------------------------------------
  return(gsl::legendre_Pl(n, x))
}



Legendres <- function(D, x) {
  #-----------------------------------------------------------------------------
  # Evaluates all Legendre polynomials up to degree D at data point x.
  # 
  # (Wrapper for GSL function legendre_Pl_array())
  #
  # DEPENDENCIES
  #  - gsl::legendre_Pl_array()
  #
  # INPUT
  #      D:   maximal degree of Legendre polynomial to be computed 
  #      x:   data point to be evaluated
  #
  # OUPUT
  #  (D+1) vector of values of Legendre polynomials at point x
  #-----------------------------------------------------------------------------
  return(as.vector(gsl::legendre_Pl_array(D, x)))
}



AssocLegendre <- function(n, k, x, csphase = TRUE) {
  #-----------------------------------------------------------------------------
  # Computes the Associated Legendre function of degree n and order k
  # at data point x.
  #
  # DEPENDENCIES
  #  - gsl::legendre_Plm()
  #
  # INPUT
  #       n:   degree of ALP to be computed 
  #       k:   order of ALP to be computed (-n <= k <= n)
  #       x:   data point to be evaluated
  # csphase:   set FALSE for convention not using Condon-Shortley phase
  #
  # OUPUT
  #  value of ALP(n,k) at point x
  #-----------------------------------------------------------------------------
  return(gsl::legendre_Plm(n, abs(k), x) * (-1)**(k * (!csphase)) *
         ((-1)**k * factorial(n + k) /
         factorial(n - k))**(1 - sign(sign(k) + 1)))
}



AssocLegendres <- function(D, x, csphase = TRUE) {
  #-----------------------------------------------------------------------------
  # Evaluates all Associated Legendre functions up to degree D at data point x.
  #
  # DEPENDENCIES
  #  - gsl::legendre_array()
  #
  # INPUT
  #       D:  maximal degree of ALP to be computed 
  #       x:  data point to be evaluated
  # csphase:  set FALSE for convention not using Condon-Shortley phase
  #
  # OUPUT   
  #  ((D+1)x(2D+1)) matrix with generic entry [i,j] the evaluated ALP
  #         of degree i-1 and order j-D-1
  # (that is: column D+1 is order 0, column 1 is order -D and
  #           column 2D+1 is order D)
  #-----------------------------------------------------------------------------
  
  # translate csphase to argument of GSL-function
  if (!csphase) csphase_arg <- 1
  else csphase_arg <- -1
  
  # prepare output
  out <- matrix(numeric((D + 1) * (2 * D + 1)), D + 1, 2 * D + 1)
  
  # get positive order values with GSL function legendre_array()
  out[ ,(D + 1):(2 * D + 1)] <- gsl::legendre_array(x, D, csphase = csphase_arg)
  
  # prepare factor matrix
  factormat <- matrix(numeric((D + 1) * D), D + 1, D)
  factormat[2:(D + 1),D] <- -1 / (1:D) / ((1:D) + 1)
  
  # other orders
  for (k in 2:D) {
    factormat[ ,D + 1 - k] <- c(rep(0, k),
                                -factormat[(k + 1):(D + 1),D + 2 - k] / 
                                ((k:D) + k) / ((k:D) - k + 1))
  }
  out[ ,1:D] <- out[ ,(2 * D + 1):(D + 2)] * factormat
  
  return(out)
}



.NextLegendre <- function(n, x, prev1, prev2) {
  #-----------------------------------------------------------------------------
  # Helper - Evaluates the next order Legendre polynomial based on the previous
  #          two orders
  # 
  # INPUT
  #     n:   degree of Legendre polynomial to be evaluated
  #     x:   argument to be evaluated
  # prev1:   value of evaluated Legendre polynomial of previous order
  # prev2:   value of evaluated Legendre polynomial of second previous order
  #
  # OUPUT
  #  value of next order ALP as numeric
  #-----------------------------------------------------------------------------
  return(((2 * n - 1) * x * prev1 - (n - 1) * prev2) / n)
}

