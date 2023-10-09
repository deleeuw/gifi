source ("coneRegression.R")
source ("splineBasis.R")
source ("gifiStructures.R")

homals <- function (data, 
                    ndim = 2, 
                    degrees = rep (0, ncol (data)), 
                    knots = knotsD (data), 
                    ordinal = rep (FALSE, ncol (data)), 
                    copies = rep (ndim, ncol (data)), 
                    active = rep (TRUE, ncol (data)),
                    ties = rep ("S", ncol (data)),
                    missing = rep,
                    itmax = 1000, 
                    eps = 1e-10, 
                    verbose = FALSE) {
  
}

