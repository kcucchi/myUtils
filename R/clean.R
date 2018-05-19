

#' remove '_cnt' or '_pref' from GBcodes
#'
#' @param vect_names a vector with names to clean
#' @return vact_names with 'GBcode' instead of 'GBcode_cnt' or 'GBcode_pref'
#' @export
trim_GBcode <- function(vect_names){

  vect_names[which(grepl(pattern = 'GBcode',x = vect_names))] <-
    'GBcode'

  return(vect_names)

}
