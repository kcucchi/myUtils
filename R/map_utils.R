

#' Transforms a sp object into a dataframe to plot into ggplot
#'
#' @param x_sp the sp object to transform into dataframe
#' @return the corresponding dataframe to be plotted using geom_polygon
#' @export
sp2df <- function(x_sp){

  x_sp@data$id <- rownames(x_sp@data)
  x_df <- fortify(x_sp, region='id')
  x_df <- dplyr::left_join(x = x_df,y = x_sp@data,by = 'id')

  return(x_df)

}

