# Documenting datasets
# R packages by Hadley Wickham (http://r-pkgs.had.co.nz/data.html)
# use devtools::load_all() to load during package development

#' Plasmodium falciparum parasite rate in 2-10 year olds in 2000-2017 in Sub-Saharan Africa
#'
#' A dataset containing longitude, latitude, point estimates, and credible intervals bounds
#'
#' @format A data frame with 144000 rows and 5 variables:
#' \describe{
#'   \item{x}{Longitude in decimal degrees}
#'   \item{y}{Latitude in decimal degrees}
#'   \item{Median}{The median posterior predictive rate of Plasmodium falciparum in 2-10 year olds, 2000-2017}
#'   \item{UCI}{The lower extreme of the credible interval around the median posterior predictive rate of Plasmodium falciparum in 2-10 year olds, 2000-2017}
#'   \item{UCI}{The upper extreme of the credible interval around the median posterior predictive rate of Plasmodium falciparum in 2-10 year olds, 2000-2017}
#' }
#' @source \url{https://map.ox.ac.uk/}
"SubSaharanAfrica_Pf_incidence"


#' Shape files for Sub-Saharan Africa
#'
#' A SpatialPolygonsDataFrame for Sub-Saharan Africa
#'
#' @format A SpatialPolygonsDataFrame with dimensions 59 and 16
#' \describe{
#'   \item{iso}{To-be-completed}
#'   \item{admn_level}{To-be-completed}
#' }
#' @source Obtained using malariaAtlas::getShp (see data-raw/get_shape_files.R)
"SubSaharanAfrica_shp"
