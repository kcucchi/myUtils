
# Define functions for wavelet analysis ---------------------------------
# sources:
# Christopher Torrence and Gilbert P. Compo
# https://stackoverflow.com/questions/22131080/wavelet-reconstruction-of-time-series

#' Definition of the Morlet wavelet
#'
#' @param t a numeric vector corresponding to times
#' where to evaluate the Morlet wavelet
#' @param t_0 (optional) a time indicating the timing translation
#' @param j (optional) an index corresponding to the calculated scale s_j
#' (the corresponding scale is s_0  2^(j*dj))
#' @return evaluation of the Morlet wavelet at values in \code{t}
#' @examples
#' N_4_original(2)
#' @export
# Lets first write a function for Wavelet decomposition as in formula (1):
mo<-function(t,t_0=0,omega=6,j=0){
  dial<-2*2^(j*.125)
  sqrt((1/dial))*pi^(-1/4)*exp(1i*omega*((t-t_0)/dial))*exp(-((t-t_0)/dial)^2/2)
}

#'
#'Perform Morlet Wavelet transform
#'@param y timeseries to transform
#'@param J determines the largest scale s_J = s_0 2^(J dj)
#'@param vect_time (optional) vector ot time points corresponding to y
#'@param dt (optional) the time step
#'@return the wavelet transform
#'@export
wt <- function(y,J=110,vect_time=NULL,dt=1){
  mat_res <- matrix(data = NA,nrow = length(y),ncol=(J+1))
  y.madj <- y - mean(y)
  for(j in 0:J){
    for(k in 1:length(y.madj)){
      mat_res[k,j+1]<-mo(t=1:(length(y.madj)),j=j,t_0=k)%*%y.madj
    }
  }
  # name the dimensions
  # rows are time points
  rownames(mat_res) <- format(vect_time)
  # columns are frequencies
  colnames(mat_res) <- dt*2*2^(0:J*.125)
  return(mat_res)
}

#'
#'Extract components between 2 scales
#'@param mat_wt matrix containing the wavelet transform
#'@param range_scale vector of length 2 containing scales to extract
#'@param fill variable to fill matrix outside of scale range
#'@return wt matrix with scales in range_scale and zero otherwise
#'@export
extract_scales <- function(mat_wt,range_scale,fill=0){
  res_wt <- matrix(data = fill,nrow = nrow(mat_wt),ncol = ncol(mat_wt))
  idx_cols <- as.numeric(colnames(mat_wt)) >= min(range_scale) &
    as.numeric(colnames(mat_wt)) <= max(range_scale)
  res_wt[,idx_cols] <- mat_wt[,idx_cols]
  return(res_wt)
}


#'
#' Reconstruct from wavelet transform
#'@param mat_wt matrix containing the wavelet transform
#'@param y.m (optional) mean of original time series
#'@param J determines the largest scale s_J = s_0 2^(J dj)
#'@return the inverse of the wavelet decomposition
#'@export
wt_inverse <- function(mat_wt,y.m=0,J=110){

  dial<-2*2^(0:J*.125)
  rec<-rep(NA,(length(y.madj)))
  for(l in 1:(length(y.madj))){
    rec[l]<-0.2144548*sum(mat_wt[l,]/sqrt(dial))
  }

  return(rec + y.m)

}

#'
#' extract phase from mat_wt
#'@param mat_wt matrix containing the wavelet transform
#'@return the time series of phase angles
#'@export
phase_angle <- function(mat_wt){

  return(atan2(y = Im(mat_wt),x = Re(mat_wt)))

}

#'plot real part of the wavelet decomposition
#'@param x vector containing values of x
#'@param f_x vector containing values of functions
#'@param bool_inc boolean indicating whether
#'returning zeros with positive slopes
#'@param bool_dec boolean indicating whether
#'returning zeros with negative slopes
#'@return the zeros of the function
#'@export
find_zeros_lin <- function(x,f_x,bool_inc=T,bool_dec=T){

  idx_change <- numeric(0)

  if(bool_inc){idx_change <- c(idx_change,which(diff(sign(f_x)) == 2))}
  if(bool_dec){idx_change <- c(idx_change,which(diff(sign(f_x)) == -2))}

  x_0 <- rep(NA,length(idx_change))

  for(i in 1:length(idx_change)){
    x_i <- x[c(idx_change[i] , idx_change[i]+1)]
    f_x_i <- f_x[c(idx_change[i] , idx_change[i]+1)]
    x_0[i] <- x_i[1] + (x_i[2] - x_i[1]) *
      (abs(f_x_i[1])/(abs(f_x_i[2])+abs(f_x_i[1])))
  }

  return(sort(x_0))

}

#'plot real part of the wavelet decomposition
#'@param mat_wt matrix containing the wavelet transform
#'@return ggplot
#'@export
plot_re <- function(mat_wt){

  wt.r <- Re(mat_wt)

  # check plot
  df_wt <- as.data.frame(wt.r)
  df_wt$date <- rownames(wt.r)
  df_wt_long <- tidyr::gather(data = df_wt,
                              key=period,value=value,-date)

  # str(df_wt_long)

  df_wt_long$period <- as.numeric(df_wt_long$period)
  df_wt_long$date <- as.Date(df_wt_long$date)

  g_res <-
    ggplot() +
    geom_tile(data = df_wt_long,
              mapping = aes(x=date,y=period,fill=Re(value))) +
    scale_y_continuous(trans='log10',
                       breaks = c(0.5,1,2)) +
    scale_fill_distiller(palette = "Spectral") +
    theme_bw()

  return(g_res)

}



# # test functions
# y = cos(2*pi/12*1:100) + cos(pi/4+2*pi/(2*12)*1:100)
# y = y + 0.1*rnorm(100) # add some noise
# t <- seq.Date(from = as.Date('2000-01-01'),by = "1 month",length.out = 100)
# plot(x = t,y = y,type='l')
#
# # wavelet analysis
# mat_wt <- wt(y = y,vect_time = t,dt=1/12)
# str(mat_wt)
#
# rec_wt <- wt_inverse(mat_wt = mat_wt,dt=1/12,y.m = mean(y))
# plot(x = t,y = rec_wt,type='l')
# lines(x = t,y = y,col='red')
#
# #
# # Plot the real part
# #
#
# print(plot_re(mat_wt) +
#         scale_y_continuous(trans='log10',
#                            breaks = c(0.5,1,2),
#                            limits = c(2/12,3)))
#
#
# #
# # keep only periodicity around 1 year
# #
#
# mat_wt_1y <-
#   extract_scales(mat_wt = mat_wt,
#                  range_scale = c(0.8,1.2))
# rec_wt_1y <-
#   wt_inverse(mat_wt = mat_wt_1y,y.m = mean(y))
# plot(x = t,y = rec_wt,type='l')
# lines(x = t,y = y,col='red')
# lines(x = t,y = rec_wt_1y,col='green')
#
# #
# # plot phase angles corresponding to c(0.8,1.2)
# #
#
# phase_1y <-
#   rowMeans(x = phase_angle(extract_scales(mat_wt = mat_wt,
#                                           range_scale = c(0.9,1.2),
#                                           fill = NA)),na.rm = T)
#
#
#
# plot(x = 1:length(t),y = phase_1y,type='l')
# t_0 <- find_zeros_lin(x = 1:length(t),f_x = phase_1y,bool_dec = F)
# points(x = t_0,
#        y = rep(0,length(t_0)),
#        col = 'red')
#
# # extract timings with phases 0
# plot(x = 1:length(t), y = rec_wt_1y,type='l')
# abline(v=t_0,col="red")

