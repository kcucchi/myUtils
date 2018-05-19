

#' Fit yearly and bi-yearly sinusoids
#'
#' @param t_months a numeric vector indicating time (in months)
#' @param prop a numeric vector indicating the proportion of cases
#' in year reported during month
#' @return a list containing seasonal coefficients and fitted values
#' @examples
#' t_months = 1:12
#' prop = c(0.1,0.1,0.2,0.2,0.3,0.3,0.3,0.2,0.3,0.2,0.2,0.1)
#' @export
fit_seasonal <- function(t_months,prop){

  # setup dataframe for harmonic regression
  dat_model <- data.frame(prop = prop)
  dat_model$term_b <- cos(2*pi/12*t_months)
  dat_model$term_c <- sin(2*pi/12*t_months)
  dat_model$term_d <- cos(4*pi/12*t_months)
  dat_model$term_e <- sin(4*pi/12*t_months)

  # fit model
  model_fit <- lm(formula = prop ~ term_b + term_c + term_d + term_e,
                  data = dat_model)

  # name resulting coefficients
  a <- model_fit$coefficients[['(Intercept)']]
  b <- model_fit$coefficients[['term_b']]
  c <- model_fit$coefficients[['term_c']]
  d <- model_fit$coefficients[['term_d']]
  e <- model_fit$coefficients[['term_e']]

  # output annual and seasonal components
  df_components <-
    data.frame(t_months = t_months,
               ts_yearly = a + b * dat_model$term_b +
                 c * dat_model$term_c,
               ts_biyearly = a + d * dat_model$term_d +
                 e * dat_model$term_e,
               ts_fit = model_fit$fitted.values,
               ts_dat = prop)

  # create results vector
  res <- rep(NA,5)
  names(res) <- c('A','SA','P','P_fit','R')

  # amplitude of annual component
  res['A'] <- sqrt(b^2 + c^2)

  # amplitude of semiannual component
  res['SA'] <- sqrt(d^2 + e^2)

  # timing of peak
  if(b>0 & c>0){
    res['P'] <- atan(c/b)
  }else if (b <0 & c<0 ){
    res['P'] <- atan(c/b ) + pi
  }else if(c > 0 & b < 0){
    res['P'] <- (-1* atan(c/b)  + pi/2)
  }else{res['P'] <- (-1*atan(c/b ) + 3*pi/2) }

  # alternative :
  # choose timing of peak as month with highest fitted counts

  res['P_fit'] <-
    df_components$t_months[which(df_components$ts_fit == max(df_components$ts_fit))]


  # relative contribution of semi-annual component
  res['R'] <- res['SA'] / (res['A'] + res['SA'])


  return(list(res=res,
              df_components=df_components))

}
