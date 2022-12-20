#------------------------------------------------------------------------------   
#     Default parameters for covariance functions
#------------------------------------------------------------------------------  
# File:    defaultparameters.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - Contains only default parameters to be loaded.
#
# DEPENDENCIES
#  none
#------------------------------------------------------------------------------  

def_exp_p1 <- 5               # range: 0 < p1  (>0.5 for atbm)
def_exp_p2 <- NULL            # unused
def_powexp_p1 <- 5            # range: 0 < p1
def_powexp_p2 <- 1            # range: 0 < p2 <= 1
def_sinepow_p1 <- NULL        # unused
def_sinepow_p2 <- 0.3         # range: 0 < p2 < 2 (atbm <=1)
def_spherical_p1 <- 3         # range: 0 < p1
def_spherical_p2 <- NULL      # unused
def_askey_p1 <- 2           # range: 0 < p1
def_askey_p2 <- 2             # range: 2 <= p2
def_wendc2_p1 <- 3            # range: 1/pi < p1
def_wendc2_p2 <- 4            # range: 4 <= p2
def_wendc4_p1 <- 2            # range: 1 <= p1
def_wendc4_p2 <- 6            # range: 6 <= p2
def_cmatern_p1 <- 5           # range: 0 < p1 (>> p2)
def_cmatern_p2 <- 0.5           # range: 0 < p2 (<= 2)

def_spect_exp <- 5            # range: 0 < p
def_spect_geometric <- 0.93    # range: 0 < p < 1
def_spect_poisson <- 12        # range: 0 < p
def_spect_bessel <- 100        # range: 0 < p
def_spect_linear <- NULL      # no parameter used



