#-------------------------------------------------------------------------------   
#           Functions for visualizing spheres and points on spheres 
#-------------------------------------------------------------------------------  
# File:    displaysphere.R
# Author:  Ruben Scherrer (ruben.scherrer@uzh.ch)
# Date:    August 2022
# Status:  Delivery
# Version: 0.9
#
# COMMENTS
#  - After using visualization methods from this file always close
#    the rgl-widget with 'rgl::close3d()'.
#
# DEPENDENCIES
#  - packages "rgl", "plotrix"
#  - Spher2Cart() from "helpfunctions.R"
#
# CONTENT
#   DisplaySphereGrid()    display a spherical grid
#   DisplaySpherePoints()  display points on the sphere
#   .ColorCode()           function for colour coding values for points
#-------------------------------------------------------------------------------   

# Sources
source("R/simGRFS/helpers/helpfunctions.R")

# get session specific default resolutions when loading
default_res <- dev.size("px")
default_ratio <- default_res[1] / default_res[2]



DisplaySphereGrid <- function(radius = 1,
                              resol = 100,
                              grid_nlat = 17, grid_nlong = 18,   
                              title = "sphere",
                              color = "grey",
                              box = TRUE,
                              legend = FALSE,
                              display_res_h = default_res[2],
                              display_res_w = default_res[1]) { 
  #-----------------------------------------------------------------------------
  # Displays a spherical grid
  #
  # REMARK: resol and grid_nlong should be even, grid_nlat should be odd
  #
  # DEPENDENCIES
  #  - package "rgl"
  #  - Spher2Cart() from helpfunctions.R
  #
  # INPUT
  #        radius:    radius of sphere
  #         resol:    resolution of sphere (in lines between poles)
  #     grid_nlat:    number of latitude gridlines (default 8+8+1=17)
  #    grid_nlong:    number of longitude gridlines (default 18)
  #         title:    title of visualization
  #         color:    color of gridlines
  #           box:    whether reference box is displayed
  #        legend:    whether displaying a legend is allowed
  # display_res_h:    resolution of display window height
  # display_res_w:    resolution of display window width
  #
  # OUTPUT
  #  rglwidget display
  #-----------------------------------------------------------------------------
  
  # prepare sphere resolution and gridlines
  grid_nlat  <- grid_nlat + 2                     # (add poles)
  grid_nlong <- grid_nlong + 1                    # (double circle)
  lat        <- seq(-pi / 2, pi / 2, len = resol)
  long       <- seq(-pi,     pi,     len = resol)
  grid_lat   <- seq(-pi / 2, pi / 2, len = grid_nlat)
  grid_long  <- seq(-pi,     pi,     len = grid_nlong)
  
  # set up rgl visualizing tool
  rgl::open3d(windowRect = c(0, 0, display_res_w, display_res_h))
  if (legend == TRUE) rgl::layout3d(matrix(1:2, 1, 2), c(0.7, 0.3), 1)
  if (box == TRUE) rgl::axes3d()
  rgl::title3d(main = title, pos = c(-1, 0, 1.4), floating = NA)
  
  # add circles to grid
  for (i in 1:grid_nlat) {
    circles      <- rbind(rep(grid_lat[i], resol), long, rep(radius, resol))
    circles_cart <- Spher2Cart(circles)
    rgl::lines3d(circles_cart[1, ],
                 circles_cart[2, ],
                 circles_cart[3, ],
                 color = color)
  }
  
  # add meridians to grid
  for (i in 1:grid_nlong) {
    meridians <- rbind(lat, rep(grid_long[i], resol), rep(radius, resol))
    meridians_cart <- Spher2Cart(meridians)
    rgl::lines3d(meridians_cart[1, ],
                 meridians_cart[2, ],
                 meridians_cart[3, ],
                 color = color)
  }
  
  # visualize
  rgl::rglwidget()
}



DisplaySpherePoints <- function(pts,
                                values = NULL,
                                colpts = "blue",
                                add = FALSE,
                                in_polar = FALSE,
                                ratio = default_ratio,
                                title = "",
                                col_span = 100,
                                col_ratio = 0.5,
                                col_reference = NULL,
                                sphere_transparency = 0.2,
                                high_contrast = FALSE,
                                box = TRUE) {
  #-----------------------------------------------------------------------------
  # Display given points on the sphere. If values are provided, color codes
  # these values accordingly.
  #
  # Comments
  #  - If add=FALSE, creates a grid by default.
  #  - Points can only be added to an existing plot if no legend is active.
  #  - If values are provided, add will always be set to FALSE to allow legend.
  #  - Older rgl-versions have trouble with sphere transparency. In these cases
  #    sphere_transparency is set automatically to 0. 
  #
  # DEPENDENCIES
  # - packages "rgl", "plotrix"
  # - Spher2Cart() from helpfunctions.R
  # - DisplaySphereGrid()
  # - .ColorCode()
  #
  # INPUT
  #                 pts:  matrix (3xN) or matrix (4xN) with data points
  #                       in columns in format c(x,y,z) or c(lat,long,r)
  #                       or c(x,y,z,value) or c(lat,long,r,values)
  #              values:  (optional) if provided, points get color coded
  #              colpts:  color of points to be plotted, if no values provided
  #                 add:  whether points should be added to existing plot
  #            in_polar:  set TRUE if input in spherical polar coordinates
  #               ratio:  ratio of resolution for displaying legend
  #               title:  title to display above sphere
  #            col_span:  span of colors to use
  #           col_ratio:  ratio of colors (between 0 and 1)
  #       col_reference:  reference values for colors (for simulations)
  # sphere_transparency:  transparency of sphere
  #       high_contrast:  instead of using linear colors, use random colors
  #
  # OUTPUT
  #  rglwidget display
  #-----------------------------------------------------------------------------

  # if values are provided via pts: assign variables
  if (dim(pts)[1] == 4) {
    values <- pts[4, ]
    pts    <- pts[1:3, ]
  }
  
  # if values are provided: setting add=FALSE and error handling 
  if (!is.null(values)) {
    add <- FALSE
    if (length(values) != dim(pts)[2]) {
      stop("Error: Amount of values and amount of points do not coincide.")
    }
  }
  
  # set radius of sphere
  radius <- sqrt(sum(pts[ ,1]**2))
  
  # if input in polar coordinates: convert coordinates
  if (in_polar) pts <- Spher2Cart(pts)
  
  # if no active rgl-decive or add=FALSE: create new grid 
  if (rgl::cur3d() == 0 || !add) {
    if (rgl::cur3d() != 0 && !add) rgl::close3d()
    if (!is.null(values)) legend = TRUE
    else legend = FALSE
    DisplaySphereGrid(radius = radius, title = title, legend = legend,
                      box = box)
  }
  
  # prepare colors
  if (!is.null(values)) {
    if (!is.null(col_reference)) {
      values <- c(values, col_reference)
      values[values > max(col_reference)] <- NA
      values[values < min(col_reference)] <- NA
    }
    colcode_list <- .ColorCode(values,
                               breaks = ceiling(col_ratio * col_span),
                               ncols = col_span,
                               random_colors = high_contrast)
    breaks    <- as.integer(colcode_list[1])
    ncols     <- as.integer(colcode_list[2])
    ceil_val  <- as.numeric(colcode_list[3])
    floor_val <- as.numeric(colcode_list[4])
    colpts    <- colcode_list[5:length(colcode_list)]
  }
  
  # add sphere to improve contrast
  rgl::shade3d(rgl::ellipse3d(diag(3), centre = c(0, 0, 0), t = 1),
               color = "lightgray",
               alpha = 1 - sphere_transparency,
               lit = FALSE,
               add = TRUE)

  # plot points
  rgl::plot3d(pts[1, ], pts[2, ], pts[3, ], col = colpts, add = TRUE, size = 5)
  
  # if values are provided: display legend
  if (!is.null(values) & length(rgl::subsceneList()) > 1 ) {
    rgl::rglwidget()
    resize <- 1 / ratio
    rgl::next3d()
    rgl::bgplot3d({
      plot.new()
      plotrix::color.legend(0.5 - resize / 10, 0.15 * resize,
                            0.5 + resize / 10, 1 - 0.15 * resize,
                            rect.col = rainbow(ncols)[1:breaks],
                            legend = round(seq(floor_val, ceil_val, len = 5),
                                           8),
                            gradient = "y",
                            cex = 1.5 + ratio / 10)
    })
  }
  
  # visualize
  rgl::rglwidget()
}



.ColorCode <- function(values,
                       breaks = 50,
                       ncols = 100,
                       random_colors = FALSE) {
  #-----------------------------------------------------------------------------
  # Color codes values for displaying points on the sphere
  #
  # No dependencies.
  #
  # INPUT
  #        values:  values to be displayed
  #        breaks:  color smoothness for values (< ncols)
  #         ncols:  total amount of colors to be used
  # random_colors:  set TRUE to use random (non-linear) colors
  #
  # OUTPUT 
  #  list of strings in the following order
  #    1. number of breaks
  #    2. number of total colours used
  #    3. ceiled maximal value
  #    4. floored minimal value
  #    5. color code for values given
  #
  # REMARK: if ratio of ncols/breaks is about 2, good visibility!
  #-----------------------------------------------------------------------------
  
  # compute spread and order of values, and optimal increment for color code
  spread <- max(values, na.rm = TRUE) - min(values, na.rm = TRUE)
  order <- round(log(spread, base = 10))
  incre <- 10**(order - 1) *
           c(0.5, 1, 2.5)[ceiling(3 * (0.5 - (order - log(spread, base = 10))))]
  
  # compute sequence and reference values
  floor_val <- floor(min(values, na.rm = TRUE) / incre) * incre
  ceil_val <- ceiling(max(values, na.rm = TRUE) / incre) * incre
  val_seq <- seq(floor_val, ceil_val, by = incre)

  # compute color index and color code
  col_index <- as.integer(cut(c(values, floor_val, ceil_val), breaks = breaks))
  
  # set colors (either randomly or linearly)
  if (random_colors) out_colors <- sample(rainbow(n = ncols))[col_index]
  else out_colors <- rainbow(n = ncols)[col_index]
  
  # set NA to grey
  out_colors[is.na(out_colors)] <-  "#FDFEFE"

  
  return(c(breaks, ncols, ceil_val, floor_val, out_colors))
}

