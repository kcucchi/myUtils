

# Define functions for splines ---------------------------------

#' Definition of the cubic B-spline (cardinal B-spline of order 4)
#'
#' @param x a numeric vector where to evaluate the cubic spline
#' @return evaluation of the cubic spline at values in \code{x}
#' @examples
#' N_4_original(2)
N_4_original <- function(x){

  res <- rep(0,length(x))

  for(i in 1:length(x)){
    if(x[i] >=0 & x[i] < 1){
      res[i] <- 1/6*x[i]^3
    }
    if(x[i] >= 1 & x[i] <= 2){
      res[i] <- (-1/2*x[i]^3 + 2 * x[i]^2 - 2*x[i] + 2/3)
    }
    if(x[i] >= 2 & x[i] <= 3){
      res[i] <- (1/2*x[i]^3 - 4 * x[i]^2 + 10*x[i] - 22/3)
    }
    if(x[i] >= 3 & x[i] <= 4){
      res[i] <- (-1/6*x[i]^3 + 2 * x[i]^2 - 8*x[i] + 32/3)
    }
  }

  return(res)

}


#' Definition of the unit cubic B-spline (centered at 0, width 1, amplitude 1)
#'
#' @param x a numeric vector where to evaluate the unit cubic spline
#' @return evaluation of the unit cubic spline at values in \code{x}
#' @examples
#' N_4(0)
#' @export
N_4_unit <- function(x){
  return(N_4_original((x*4)+2) / (2/3))
}

#' Definition of the unit cubic B-spline (centered at x_0, width sigma, amplitude A)
#'
#' @param x a numeric vector where to evaluate the unit cubic spline
#' @return evaluation of the unit cubic spline at values in \code{x}
#' @examples
#' x <- seq(-5,5,by = 0.1)
#' plot(x,N_4_original(x),ylim = c(0,1))
#' lines(x,N_4_unit(x),col='red')
#' lines(x,N_4(x,x_0 = 1,A = 0.4,sigma = 2),col='green')
#' @export
N_4 <- function(x,x_0,A,sigma){
  return(A * N_4_unit((x-x_0)/sigma))
}


# Define functions for optimization ---------------------------------


#' Definition of the distance function used in optimisation
#'
#' @param param the vector of parameters to be optimized (amplitude, timing and width)
#' @param t numeric vector of time values
#' @param y numeric vector of observed values at time \code{t}
#' @return the residual sum of squares between N_4(t,x_0,A,sigma) and y
#' @export
dist <- function(param,t,y){
  x_0 <- param[1]
  A <- param[2]
  sigma <- param[3]
  val_spline <- N_4(x = t,x_0 = x_0,A = A,sigma = sigma)
  val_diff <- (y - val_spline)
  return(sum(val_diff^2))
}

#'Function optimizing the timing, spread and amplitude of peak
#'
#' @param t numeric vector of time values
#' @param y numeric vector of observed values at time \code{t}
#' @return the residual sum of squares between N_4(t,x_0,A,sigma) and y
#' @export
spline_optim <- function(t,y){
  param_optim <- optim(par = c(8,1,1),fn = function(param){dist(param,t = t,y = y)})$par
  names(param_optim) <- c('x_0','A','sigma')
  return(param_optim)
}




