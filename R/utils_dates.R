#' Calculate number of days since start of epidemics
#'
#' @param date_0 a Date
#' @return the number of days
#' between date in date_0
#' and start of the corresponding epidemic season
#' (1st Feb if date is in Feb-Jul,
#' 1st Aug if date is in Aug-Jan)
#' @export
date2nbDaysSinceEpid <- function(date_0){

  spring <- (lubridate::month(date_0) %in% 2:7)
  year <- lubridate::year(date_0)

  if(spring){
    return(as.numeric(date_0 -
                        as.Date(paste0('02/01/',as.numeric(year)),
                                format='%m/%d/%Y')))
  }else{
    res <- as.numeric(date_0 -
                        as.Date(paste0('08/01/',as.numeric(year)),
                                format='%m/%d/%Y'))
    if(res < 0){
      # if the peak is in january,
      # the start of the season is in the previous year
      res <- as.numeric(date_0 -
                          as.Date(paste0('08/01/',as.numeric(year)-1),
                                  format='%m/%d/%Y'))
    }
    return(res)
  }
}
