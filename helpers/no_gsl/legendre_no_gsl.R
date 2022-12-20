#-------------------------------------------------------------------------------   
#   Complete functions for evaluating (associated) Legendre polynomials
#   without package "gsl"
#-------------------------------------------------------------------------------  
# File:    legendre_no_gsl.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - The following functions evaluate (associated) Legendre polynomials not
#    using any further packages.
#  - Using this file instead of legendre.R with "gsl" results in 
#                 x0.3 slower performance for Legendre()
#                 x0.1 slower performance for Legendres()
#                 x0.3 slower performance for AssocLegendre()
#                 x0.05 slower performance for AssocLegendres()
#  - The algorithms for computing assoc. Legendre polynomials are stable for:
#                 range x in (-0.9999999999,0.9999999999) up to degree 50
#                 range x in (-0.9999,0.9999)             up to degree 85
#                 range x in (-0.99,0.99)                 up to degree 100 
#                 range x in (-0.5,0.5)                   up to degree 150
#  - For the computation of the Legendre polynomials no instability issues
#    have been detected.
#  - There are different conventions about the (-1) factor for odd orders. Per
#    default, the Condon-Shortley convention is used, else use the optional
#    argument csphase = FALSE.
#
# DEPENDENCIES
#  none
#
# CONTENT
# Main functions:
#   Legendre()        Evaluates Legendre polynomial 
#   Legendres()       Evaluates all Legendre polynomials up to some degree
#   AssocLegendre()   Evaluates associated Legendre polynomial 
#   AssocLegendres()  Evaluates all assoc. Leg. polys. up to some degree
#
# Helper functions:  
#   .NextLegendre()
#   .NextALP()
#   .PreviousALP()
#-------------------------------------------------------------------------------  



Legendre <- function(n, x) {
  #-----------------------------------------------------------------------------
  # Evaluates recursively the Legendre polynomial of degree d at points x
  # x range: -1 <= x <= 1
  #
  # DEPENDENCIES
  #  - .NextLegendre() 
  #
  # INPUT
  #  n:  degree of Legendre polynomial to be computed 
  #  x:  data point to be evaluated
  #
  # OUPUT
  #  value of Legendre polynomial of degree d at point x
  #-----------------------------------------------------------------------------
  
  # handle special cases n = 0 and n = 1
  if (n == 0) return(1**x)
  else if (n == 1) return(x)
  
  # initial cases
  prev2 <- 1
  prev1 <- x
  
  # loop over degrees and assign
  for (i in 2:n) {
    new   <- .NextLegendre(i, x, prev1, prev2)
    prev2 <- prev1
    prev1 <- new
  }
  
  return(new)
}



Legendres <- function(D, x) {
  #-----------------------------------------------------------------------------
  # Evaluates all Legendre polynomials up to  degree D at point x
  # x range: -1 <= x <= 1
  #
  # DEPENDENCIES
  #  - .NextLegendre() 
  #
  # INPUT
  #  D:   maximal degree of Legendre polynomial to be computed 
  #  x:   data point to be evaluated
  #
  # OUPUT
  #  (D+1) vector of values of Legendre polynomials at point x
  #-----------------------------------------------------------------------------
  
  # handle special cases D = 0 and D = 1
  if (D == 0) return(1)
  else if (D == 1) return(c(1, x))
  
  # initial cases
  out <- numeric(D + 1)
  out[1] <- 1
  out[2] <- x
  
  # loop over degrees and assign
  for (i in 2:D) {
    out[i + 1] <- .NextLegendre(i, x, out[i], out[i - 1])
  }
  
  return(out)
}



AssocLegendre <- function(n, k, x, csphase = TRUE) {
  #-----------------------------------------------------------------------------
  # Evaluates recursively the Associated Legendre Polynomial (ALP)
  # of degree n and order k at point x
  # parameter range: n >= k
  #
  # The algorithm is stable for: 
  #             range x in (-0.9999999999,0.9999999999) up to degree 50
  #             range x in (-0.9999,0.9999)             up to degree 85
  #             range x in (-0.99,0.99)                 up to degree 100 
  #             range x in (-0.5,0.5)                   up to degree 150
  #
  # Algorithm inspired by Algorithm from Selezneva et al.
  # http://www.scielo.org.co/pdf/racefn/v37n145/v37n145a09.pdf
  # But is expanded in both directions.
  #
  # DEPENDENCIES
  #  - .PreviousALP() 
  #  - .NextALP()
  #
  # INPUT
  #       n:   degree of ALP to be computed 
  #       k:   order of ALP to be computed (-n <= k <= n)
  #       x:   data point to be evaluated
  # csphase:   set TRUE for convention using (-1)**k
  #
  # OUPUT
  #  value of ALP(n,k) at point x
  #-----------------------------------------------------------------------------
  
  # handle base cases
  if (n == 0 & k == 0) return(1)
  if (n == 1 & k == 0) return(x)
  
  # handle special cases when abs(x)=1
  if (x == 1) {
    if (k == 0) return(1) else return(0)
  } else if (x == -1) {
    if (k == 0) return((-1)**n) else return(0)
  }
  
  #--------------------- MAIN COMPUTATION: non-positive order ------------------
  if (k < 1) {
    
    # get initial values
    prev2 <- (-1)**n / (2**n * factorial(n)) * (1 - x**2)**(n / 2)
    
    # stop in boundary case
    if (k == (-n)) {
      out <- prev2
      if (csphase) out <- out * (-1)**k
      return(out)
    }
    
    # get initial values 2
    prev1 <-  -2 * n * x * (1 - x**2)**(-0.5) * prev2
    
    # stop in boundary case
    if (k == (-n + 1)) {
      out <- prev1
      if (csphase) out <- out * (-1)**k
      return(out)
    }
    
    # use reccurence relation up to required order
    for (i in (-n + 2):(k)) {
      nextord <- .NextALP(n, i, prev1, prev2, x)
      prev2 <- prev1
      prev1 <- nextord
    }
    out <- nextord
  }
  
  #--------------------- MAIN COMPUTATION: positive order ----------------------
  else {
    
    # get initial values
    next2 <- prod((n + 1):(2 * n)) / (2**n) * (1 - x**2)**(n / 2)
    
    # stop in boundary case
    if (k == n) {
      out <- next2
      if (csphase) out <- out * (-1)**k
      return(out)
    }
    
    # get initial values 2
    next1 <-  x * (1 - x**2)**(-0.5) * next2
    
    # stop in boundary case
    if (k == (n - 1)) {
      out <- next1
      if (csphase) out <- out * (-1)**k
      return(out)
    }
    
    # use reccurence relation up to required order
    for (i in (n - 2):k) {
      prevord <- .PreviousALP(n, i, next1, next2, x)
      next2 <- next1
      next1 <- prevord
    }
    out <- prevord
  }

  # correcting for (-1)**k convention
  if (csphase) out <- out * (-1)**k
  
  return(out)
}



AssocLegendres <- function(D, x, csphase = TRUE) {
  #-----------------------------------------------------------------------------
  # Evaluates recursively the Associated Legendre Polynomial (ALP)
  # for all degrees and orders up to degree n.
  #
  # The algorithm is stable for: 
  #             range x in (-0.9999999999,0.9999999999) up to degree 50
  #             range x in (-0.9999,0.9999)             up to degree 85
  #             range x in (-0.99,0.99)                 up to degree 100 
  #             range x in (-0.5,0.5)                   up to degree 150
  #
  # Algorithm inspired by Algorithm from Selezneva et al.
  # http://www.scielo.org.co/pdf/racefn/v37n145/v37n145a09.pdf
  # But is expanded in both directions.
  #
  # DEPENDENCIES
  #  - .PreviousALP() 
  #  - .NextALP()
  #  - AssocLegendre()
  #
  # INPUT
  #  D:   maximal degree of ALP to be computed 
  #  x:   data point to be evaluated
  #
  # OUPUT   
  #  ((D+1)x(2D+1)) matrix with generic entry [i,j] the evaluated ALP
  #         of degree i-1 and order j-D-1
  # (that is: column D+1 is order 0, column 1 is order -D and
  #           column 2D+1 is order D)
  #-----------------------------------------------------------------------------
  
  # setup output matrix
  out <- matrix(numeric((D + 1) * (2 * D + 1)), D + 1, 2 * D + 1)
  out[1,D + 1] <- 1 
  
  # handle special cases n < 3
  if (D == 0) {
    out[1,1] <- 1
    return(out)
  } else if (D == 1) {
    out[1,2] <- 1
    out[2,2] <- x
    out[2,1] <- AssocLegendre(1, -1, x, csphase = csphase)
    out[2,3] <- AssocLegendre(1, 1, x , csphase = csphase)
    return(out)
  } else if (D == 2) {
    out[1,3] <- 1
    out[2,3] <- x
    out[2,2] <- AssocLegendre(1, -1, x, csphase = csphase)
    out[2,4] <- AssocLegendre(1, 1, x, csphase = csphase)
    for (k in (-2):2) {
      out[3,3 + k] <- AssocLegendre(2, k, x, csphase = csphase)
    }
    return(out)
  }
  
  # handle special cases when abs(x)=1
  if (x == 1) {
    out[2:(D + 1),D + 1] <- rep(1, D)
    return(out)
  } else if (x == -1) {
    if (D %% 2 == 0) {
      out[2:(D + 1),D + 1] <- rep(c(-1, 1), D / 2)
    } else {
      out[2:D,D + 1] <- rep(c(-1, 1), (D - 1) / 2)
      out[D + 1,D + 1] <- -1
    }
    return(out)
  }
  
  # get square root
  sq <- sqrt(1 - x**2)
  
  # set up first degree initial values
  out[2,D]     <- -0.5 * sq
  out[2,D + 1] <- x
  out[2,D + 2] <- sq
  
  # recursively compute initial values for recursion
  for (d in 2:D) {
    factor1 <- sq * (2 * d - 1)
    out[d + 1,D + 1 - d] <- -out[d,D + 2 - d] * sq / (2 * d)
    out[d + 1,D + 2 - d] <- -out[d,D + 3 - d] * sq / (2 * (d - 1))
    out[d + 1,D + 1 + d] <- out[d,D + d] * factor1
    out[d + 1,D + d]     <- out[d,D + d - 1] * factor1
  }
  
  # set up degree = 2 (special case, only order 0 is missing)
  out[3,D + 1] <- .NextALP(2, 0, out[3,D], out[3,D - 1], x)
  
  # recursively compute all other orders for each degree based on initial values
  for (d in 3:D) {
    for (k in (-d + 2):0) {
      out[d + 1,D + 1 + k] <- .NextALP(d, k, out[d + 1,D + k],
                                       out[d + 1,D + k - 1], x)
    }
    for (k in (d - 2):1) {
      out[d + 1,D + 1 + k] <- .PreviousALP(d, k, out[d + 1,D + 2 + k],
                                           out[d + 1,D + 3 + k], x)
    }
  }
  
  # correcting for (-1)**k convention
  if (csphase) {
    if (D %% 2 == 0) {
      for (column in seq(2, 2 * D, by = 2)) {
        out[ ,column] <- -out[ ,column]
      }
    }
    else {
      for (column in seq(1, 2 * D + 1, by = 2)) {
        out[ ,column] <- -out[ ,column]
      }
    }
  }
  
  return(out)
}



#-------------------------------------------------------------------------------   
#                           Helper Functions
#-------------------------------------------------------------------------------



.NextALP <- function(d, k, prev1, prev2, x) {
  #-----------------------------------------------------------------------------
  # Helper - Evaluates the next order ALP based on the previous
  #          two orders
  #
  # No dependencies.
  # 
  # INPUT
  #     d:   degree of ALP to be evaluated
  #     k:   order of ALP to be evaluated
  # prev1:   value of evaluated ALP of previous order
  # prev2:   value of evaluated ALP or second previous order
  #     x:   argument to be evaluated
  #
  # OUPUT
  #  value of next order ALP as numeric
  #-----------------------------------------------------------------------------
  return(2 * (k - 1) * x / sqrt(1 - x**2) * prev1 -
         (d - k + 2) * (d + k - 1) * prev2)
}



.PreviousALP <- function(d, k, next1, next2, x) {
  #-----------------------------------------------------------------------------
  # Helper - Evaluates the previous order ALP based on the next
  #          two orders
  #
  # No dependencies.
  # 
  # INPUT
  #     d:   degree of ALP to be evaluated
  #     k:   order of ALP to be evaluated
  # next1:   value of evaluated ALP of next order
  # next2:   value of evaluated ALP or second next order
  #     x:   argument to be evaluated
  #
  # OUPUT
  #  value of previous order ALP as numeric
  #-----------------------------------------------------------------------------
  return((-next2 + (2 * (k + 1) * x) / sqrt(1 - x**2) * next1) /
         ((d - k) * (d + k + 1)))
}



.NextLegendre <- function(d, x, prev1, prev2) {
  #-----------------------------------------------------------------------------
  # Helper - Evaluates the next order Legendre polynomial based on the previous
  #          two orders
  # 
  # No dependencies.
  #
  # INPUT
  #     d:   degree of Legendre polynomial to be evaluated
  #     x:   argument to be evaluated
  # prev1:   value of evaluated Legendre polynomial of previous order
  # prev2:   value of evaluated Legendre polynomial or second previous order
  #
  # OUPUT
  #  value of next order ALP as numeric
  #-----------------------------------------------------------------------------
  return(((2 * d - 1) * x * prev1 - (d - 1) * prev2) / d)
}
