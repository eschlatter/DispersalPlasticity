#' Create a map with regular shapes or a Gaussian random field
#' 
#' This function returns a map with an approximate proportion of the area in idealized shapes (circles, squares, or fractal).
#' 
#' @param row integer. Number of rows in the map.
#' @param col integer. Number of columns in the map.
#' @param p numeric. Proportion of the map to be in cover type "1" (or a vector of proportions of multiple cover types if type=="fractal").
#' @param nshape integer. Number of shapes to create (maximum of 5).
#' @param type string. The type of map to create: 'fractal', 'circle', 'square', 'gradient', 'stripes', 'checkerboard'
#' @param maxval numeric. The maximum value of the map.
#' @param minval numeric. The minimum value of the map.
#' @param otherparm list. Other parameters passed to be extracted.
#' @param stripe.size numeric. Width of stripes if type == 'stripes'.
#' @param range numeric. The level of aggregation if type == 'fractal'.
#' @param binmap logical. If TRUE, a binary map will be returned.
#' @param gradtype string. If 'row', gradient will cross rows, if 'col' it will go across columns.
#' @param rasterflag logical. If TRUE, will create a georeferenced map.
#' @param minx numeric. The minimum easting or longitude value for the map.
#' @param miny numeric. The minimum northing or latitude value for the map.
#' @param cellsize numeric. The cellsize (in map units defined in projstr).
#' @param projstr string. The projection string that defines the map.
#' @author James D. Forester
#' @export
ideal.map <- function(row, col, p = 0.1, nshape = 1, type = "fractal", maxval = 1, minval = 0, otherparm = NULL, stripe.size = 0.5, range = row/10, 
                      binmap = TRUE, gradtype = "col", rasterflag = FALSE, minx = 562569, 
                      miny = 5230469, cellsize = 30, projstr = "+proj=utm +zone=15 +datum=NAD83", plotflag=FALSE) {
  ## returns a map with an approximate proportion of the area in idealized shapes (circles or squares) The maximum number of shapes is 5.
  ## maplist<-list(row=50,col=50,p=0.05,nshape=5) testmap<-ideal.map(100,100,otherparm=maplist);image(1:nrow(testmap),1:ncol(testmap),testmap)
  ## parameters in otherparm will OVERRIDE those listed explicitly in the function call.
  
  ## stripe.size must be < 1 (i.e. what proportion of the rows should be in a given stripe or checker width)
  
  ## type may be 'circle', 'square', 'gradient', 'stripes', 'checkerboard'
  
  ## extract any parameters passed as a list via otherparm
  if (is.list(otherparm)) {
    for (i in 1:length(otherparm)) {
      ## first extract the object value
      tempobj <- otherparm[[i]]
      ## now create a new variable with the original name of the list item
      eval(parse(text = paste(names(otherparm)[[i]], "= tempobj")))
    }
  }
  
  
  shapearea <- round(row * col * p/nshape)
  
  squareside <- round(sqrt(shapearea))
  radius <- round(sqrt(shapearea/pi))
  center.row <- round(0.5 * row)
  center.col <- round(0.5 * col)
  
  map <- matrix(minval, row, col)
  
  if (type %in% c("square", "circle")) {
    
    ## this hard-wires the location of the (maximum five) shapes. Could be smarter, but this works.
    # shapelocs <- list(shape1 = cbind(center.row, center.col), 
    #                   shape2 = cbind(round(center.row/2), round(center.col/2)), 
    #                   shape3 = cbind(round(center.row * 1.5), round(center.col * 1.5)), 
    #                   shape4 = cbind(round(center.row * 1.5), round(center.col/2)), 
    #                   shape5 = cbind(round(center.row/2), round(center.col * 1.5)))
    shapelocs <- list(shape1 = c(round(runif(n=1,min=1,max=row)),round(runif(n=1,min=1,max=col))),
                      shape2 = c(round(runif(n=1,min=1,max=row)),round(runif(n=1,min=1,max=col))),
                      shape3 = c(round(runif(n=1,min=1,max=row)),round(runif(n=1,min=1,max=col))),
                      shape4 = c(round(runif(n=1,min=1,max=row)),round(runif(n=1,min=1,max=col))),
                      shape5 = c(round(runif(n=1,min=1,max=row)),round(runif(n=1,min=1,max=col))))
    
    for (shape in 1:min(nshape, 5)) {
      shaperow <- shapelocs[[shape]][1]
      shapecol <- shapelocs[[shape]][2]
      
      if (type == "square") {
        ## draw the squares
        startrow <- round(shaperow - 0.5 * squareside)
        startcol <- round(shapecol - 0.5 * squareside)
        map[startrow:(startrow + squareside), startcol:(startcol + squareside)] <- maxval
      } else {
        ## draw the circles
        startrow <- max(round(shaperow - radius), 1)
        startcol <- max(round(shapecol - radius), 1)
        for (i in startrow:min((startrow + (2 * radius)), row)) {
          for (j in startcol:min((startcol + (2 * radius)), col)) {
            if (sqrt((shaperow - i)^2 + (shapecol - j)^2) <= radius) {
              map[i, j] <- maxval
            }
          }
        }
      }
    }
  } else if (type == "gradient") {
    if (gradtype == "row") {
      slope <- (minval - maxval)/row
      for (i in 1:row) {
        map[i, ] <- slope * (i - 1) + maxval
      }
    } else {
      slope <- (minval - maxval)/col
      for (i in 1:col) {
        map[, i] <- slope * (i - 1) + maxval
      }
    }
  } else if (type %in% c("stripes", "checkerboard")) {
    width <- round(row * stripe.size)
    start.stripe <- 1
    while (start.stripe < row) {
      map[start.stripe:min((start.stripe + width - 1), row), ] <- maxval
      start.stripe <- start.stripe + 2 * width
    }
    if (type == "checkerboard") {
      map <- abs(map + (-1 * t(map)))
    }
  } else if (type == "fractal") {
    map <- sim.grf2(row - 1, col - 1, range = range, binmap = binmap, p1 = p)
    
    
    if(!binmap){
      ##move map to desired range of values
      minmap=min(map)
      map=map-minmap
      
      maxmap=max(map)
      map = map/maxmap
      map=map*(maxval-minval) + minval
    }
    #map[map < 0] <- 0
    
  }
  
  if (rasterflag) {
    library(raster)
    mapdat <- list()
    mapdat$x <- seq(minx+0.5*cellsize, by = cellsize, len = ncol(map))
    mapdat$y <- seq(miny+0.5*cellsize, by = cellsize, len = nrow(map))
    mapdat$z <- map
    map <- raster(mapdat$z, xmn = range(mapdat$x)[1]-0.5*cellsize, xmx = range(mapdat$x)[2]+0.5*cellsize, ymn = range(mapdat$y)[1]-0.5*cellsize, ymx = range(mapdat$y)[2]+0.5*cellsize, crs = CRS(projstr))
    
  }
  if(plotflag) plot(map) 
  
  return(map)
  
}



#' An exponential covariance function
#' 
fncova <- function(x, r) {
  ## x distance; r range
  exp(-x/r)
}

#' Simulate a Gaussian random field
#' 
#' Simulation of a Gaussian random field using periodic embedding and an exponential covariance function
#' 
#' @param m integer. Number of rows in the map.
#' @param n integer. Number of columns in the map.
#' @param range numeric. The range of the exponential covariance function.
#' @param plotit logical. If TRUE, a map will be plotted.
#' @param binmap logical. If TRUE, a binary map will be returned.
#' @param p1 numeric. Proportion of covertype "1" if binmap == TRUE. If length(p1) > 1 and sum(p1) == 1, length(p1) cover types will be created (if sum(p1) < 1, length(p1) + 1 cover types will be created, and cover-type "0" will consist of 1-sum(p1) of the map).
#' @author Hae Kyung Im, James Forester
sim.grf2 <- function(m, n, range, plotit = F, binmap = F, p1 = 0.5) {
  ## range in cell units
  ncov=length(p1)
  if(sum(p1) >1){
   print("Cover-type proportions sum to > 1.")
   break 
  }else if(sum(p1) <1) {
    ##This will assign the residual area to class "0"
    p1=c(1-sum(p1), p1)
    ncov=ncov+1
  }
  
  N <- 4 * n * m
  h1 <- 1
  h2 <- 1
  ## latdist distances from (0,0) to (i,j) points i=0...m j=0...n
  latdist <- matrix(NA, m + 1, n + 1)
  rowmat <- row(latdist)
  colmat <- col(latdist)
  latdist <- sqrt((rowmat - 1)^2 * h1^2 + (colmat - 1)^2 * h2^2)
  ## Mr is matrix with first column of covariance matrix ordered by column, ie first col of Mr is (c00,c01,c02,...)'
  Mr <- fncova(latdist, range)
  ## Ms is matrix containing the first column of embedded cov matrix S.  S is circulant
  Ms <- matrix(NA, 2 * m, 2 * n)
  Ms[1:(m + 1), 1:(n + 1)] <- Mr
  Ms[(m + 2):(2 * m), 1:(n + 1)] <- Ms[m:2, 1:(n + 1)]
  Ms[1:(2 * m), (n + 2):(2 * n)] <- Ms[1:(2 * m), n:2]
  ## s <- as.vector(Ms) ##=c(col1,col2,...)
  st <- fft(Ms)/N
  ep <- complex(real = rnorm(N), imaginary = rnorm(N))
  ep <- matrix(ep, nrow = 2 * m)
  et <- sqrt(st) * ep
  e <- fft(et)
  Z <- e[1:(m + 1), 1:(n + 1)]
  X <- Re(Z)
  Y <- Im(Z)
  if (binmap) {
    X=X*0+(ncov-1)
    pres=1
    cut.val.old <- quantile(Re(Z), pres)
    for(i in 1:(ncov-1)){
      pres= pres - p1[i]
      cut.val <- quantile(Re(Z), pres)
      #print(paste(cut.val, cut.val.old))
      X[which(Re(Z) > cut.val & Re(Z) <= cut.val.old)] <- (i-1)
      cut.val.old=cut.val
    }
  }
  
  
  if (plotit) 
    image(1:(n + 1), 1:(m + 1), X,asp=1)
  return(X)
  
}





#' Create neutral landscape maps
#' 
#' Use standard methods to generate fractal maps. Binary and continuous surfaces may be produced.
#' 
#' @param k integer. The extent of the map (2^k+1)^2 pixels
#' @param h numeric. Level of aggregation in the map.
#' @param p numeric (0,1). The proportion of map in habitat 1
#' @param binary logical. If TRUE, a 0/1 categorical landscape is produced.
#' @param plotflag logical. If TRUE, the map will be plotted.
#' @param rasterflag logical. If TRUE, a spatially-referenced raster will be returned.
#' @param minx numeric. The minimum x coordinate of the raster.
#' @param miny numeric. The minimum y coordinate of the raster.
#' @param cellsize integer. The number of spatial units per raster pixel. 
#' @param projstr string. The proj4string describing the spatial reference of the raster.
#' @author Shannon Pittman, James Forester
#' @export
#' @example examples/neutral.landscape_example.R
fracland <- function(k, h, p, binary = TRUE, plotflag = FALSE, rasterflag = FALSE, minx = 562569, 
                     miny = 5230469, cellsize = 30, projstr = "+proj=utm +zone=15 +datum=NAD83",  range= NULL, ...) {
  ## Function for creating neutral landscapes Shannon Pittman University of Minnesota May, 2013 k = the extent of the map (2^k+1)^2 pixels h =
  ## how clumped the map should be (ranging from ?? to ??) -- weird behavior at higher values p = proportion of map in habitat 1 binary =
  ## plotflag == if TRUE will plot a filled contour version of the matrix
  
  ## function call: testmap=land(6,1,.5,FALSE,TRUE)
  library(sp)
  
  # k <- 6 # Scalar-determines size of landscape
  A <- 2^k + 1  # Scalar-determines length of landscape matrix
  
  # if(algorithm=="Gaussian"){
  #   if(is.null(range)) range = A*h
  #   B=ideal.map(A,A,p=p, range=range, binmap=binary, ...)
  #   
  # }else{
    
    B <- matrix(0, A, A)  # Creates landscape matrix
    #HabValue <- 0.5  # Determines value assigned to 'Habitat' points
    # h <- 0.9 # how clumped the landscape is
    #PixelsPerMeter <- 1  #set scale of landscape
    # p <- 0.5 #percentage of landscape that is 'Habitat'
    B[1, 1] <- 0
    B[1, A] <- 0
    B[A, 1] <- 0
    B[A, A] <- 0
    
    
    iter <- 1
    for (iter in 1:k) {
      scalef <- (0.5 + (1 - h)/2)^(iter)
      
      d <- 2^(k - iter)
      
      # ALL SQUARE STEPS#
      for (i in seq(d + 1, A - d, 2 * d)) {
        for (j in seq(d + 1, A - d, 2 * d)) {
          B[i, j] <- mean(c(B[i - d, j - d], B[i - d, j + d], B[i + d, j - d], B[i + d, j + d])) + scalef * rnorm(n = 1)
        }
      }
      
      # OUTSIDE DIAMOND STEP#
      for (j in seq(d + 1, A - d, 2 * d)) {
        B[1, j] <- mean(c(B[1, j - d], B[1, j + d], B[1 + d, j])) + scalef * rnorm(n = 1)
        B[A, j] <- mean(c(B[A, j - d], B[A, j + d], B[A - d, j])) + scalef * rnorm(n = 1)
      }
      
      for (i in seq(d + 1, A - d, 2 * d)) {
        B[i, 1] <- mean(c(B[i - d, 1], B[i + d, 1], B[i, 1 + d])) + scalef * rnorm(n = 1)
        B[i, A] <- mean(c(B[i - d, A], B[i + d, A], B[i, A - d])) + scalef * rnorm(n = 1)
      }
      
      # INSIDE DIAMOND STEP#
      if (2 * d + 1 <= A - 2 * d) {
        for (i in seq(d + 1, A - d, 2 * d)) {
          for (j in seq(2 * d + 1, A - 2 * d, 2 * d)) {
            B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
          }
        }
        
        for (i in seq(2 * d + 1, A - 2 * d, 2 * d)) {
          for (j in seq(d + 1, A - d, 2 * d)) {
            B[i, j] <- mean(c(B[i - d, j], B[i + d, j], B[i, j - d], B[i, j + d])) + scalef * rnorm(n = 1)
          }
        }
      }
      
      iter <- iter + 1
    }
    
    if (binary == T) {
      R <- sort(B)
      PosR <- (1 - p) * length(R)  #larger values become habitat, designated as 0
      pval <- R[PosR]
      T1 <- which(B > pval)
      T2 <- which(B <= pval)
      B[T1] <- 0  #habitat is 0
      B[T2] <- 1
      if (plotflag) 
        filled.contour(B, levels = c(0, 0.5, 1), col = c("black", "white"))
    } else {
      if (plotflag) 
        filled.contour(B)
    }
#  }
  
  if (rasterflag) {
    library(raster)
    mapdat <- list()
    mapdat$x <- seq(minx+0.5*cellsize, by = cellsize, len = ncol(B))
    mapdat$y <- seq(miny+0.5*cellsize, by = cellsize, len = nrow(B))
    mapdat$z <- B
    B <- raster(mapdat$z, xmn = range(mapdat$x)[1]-0.5*cellsize, xmx = range(mapdat$x)[2]+0.5*cellsize, ymn = range(mapdat$y)[1]-0.5*cellsize, ymx = range(mapdat$y)[2]+0.5*cellsize, crs = CRS(projstr))
    
  }
  return(B)
  
}
