# Documenting data sets
# R packages by Hadley Wickham (http://r-pkgs.had.co.nz/data.html)
# use devtools::load_all() to load during package development

#' \emph{Plasmodium falciparum} predicted incidence
#'
#' \emph{P. falciparum} predicted all-age incidence (clinical cases per 1,000 population per annum)
#' in 2017 for sub-Saharan Africa at 2.5 arcminute (approximately 5km) resolution [1].
#'
#' The median and credible interval were computed using samples from a posterior
#' predictive simulation that approximated the joint posterior predictive
#' distribution thereby accounting for spatial covariance [1,2].
#'
#' @format A data frame with 1794240 observations and four variables:
#' \describe{
#'   \item{x}{Longitude in decimal degrees}
#'   \item{y}{Latitude in decimal degrees}
#'   \item{z}{Median predicted incidence at location x y}
#'   \item{u}{Width of the 95\% predicted incidence credible interval at location x y}
#' }
#'
#' @source These data are available at the Malaria Atlas Project (MAP) website \url{https://map.ox.ac.uk/}.
#'    Specifically, they were obtained by selecting
#'    'ANNUAL MEAN OF PF INCIDENCE' at \url{https://map.ox.ac.uk/malaria-burden-data-download/}.
#' @references
#' \describe{
#'   \item{[1]}{Weiss DJ, Lucas TCD, Nguyen M, et al. Mapping the global
#'   prevalence, incidence, and mortality of Plasmodium falciparum, 2000--17: a
#'   spatial and temporal modelling study. Lancet 2019.}
#'   \item{[2]}{Gething PW, Patil AP, and Hay SI. Quantifying aggregated
#'   uncertainty in Plasmodium falciparum malaria prevalence and populations at
#'   risk via efficient space-time geostatistical joint simulation. PLoS
#'   Computational Biology 2010.}
#'}
#' @examples
#' str(SubSaharanAfrica_Pf_incidence)
#' head(SubSaharanAfrica_Pf_incidence)
"SubSaharanAfrica_Pf_incidence"


#' \emph{Plasmodium falciparum} predicted incidence
#'
#' \emph{P. falciparum} predicted all-age incidence (clinical cases per 1,000 population per annum)
#' in 2017 for central Africa at 2.5 arcminute (approximately 5km) resolution [1].
#'
#' The median and credible interval were computed using samples from a posterior
#' predictive simulation that approximated the joint posterior predictive
#' distribution thereby accounting for spatial covariance [1,2].
#'
#' @format An data frame with 270083 observations and four variables:
#' \describe{
#'   \item{x}{Longitude in decimal degrees}
#'   \item{y}{Latitude in decimal degrees}
#'   \item{z}{Median predicted incidence at location x y}
#'   \item{u}{Width of the 95\% predicted incidence credible interval at location x y}
#' }
#'
#' @source These data are available at the Malaria Atlas Project (MAP) website \url{https://map.ox.ac.uk/}.
#'    Specifically, they were obtained by selecting
#'    'ANNUAL MEAN OF PF INCIDENCE' at \url{https://map.ox.ac.uk/malaria-burden-data-download/}.
#' @references
#' \describe{
#'   \item{[1]}{Weiss DJ, Lucas TCD, Nguyen M, et al. Mapping the global
#'   prevalence, incidence, and mortality of Plasmodium falciparum, 2000--17: a
#'   spatial and temporal modelling study. Lancet 2019.}
#'   \item{[2]}{Gething PW, Patil AP, and Hay SI. Quantifying aggregated
#'   uncertainty in Plasmodium falciparum malaria prevalence and populations at
#'   risk via efficient space-time geostatistical joint simulation. PLoS
#'   Computational Biology 2010.}
#'}
#' @examples
#' str(CentralAfrica_Pf_incidence)
#' head(CentralAfrica_Pf_incidence)
"CentralAfrica_Pf_incidence"


#' Shape files for sub-Saharan Africa
#'
#' An object of class SpatialPolygonsDataFrame from the R package \strong{sp} v1.3-1 (see reference) containing shape file data for sub-Saharan Africa.
#'
#' @format
#' @source Obtained using malariaAtlas::getShp; see \url{https://github.com/artaylor85/pixelate/blob/master/data-raw/get_shape_files.R}.
#' @seealso \url{https://www.rdocumentation.org/packages/sp/versions/1.3-1/topics/SpatialPolygonsDataFrame-class}
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support
#'   for Spatial Vector Data. The R Journal 10 (1), 439-446,
#'   https://doi.org/10.32614/RJ-2018-009
"SubSaharanAfrica_shp"


#' Shape files for central Africa
#'
#' An object of class SpatialPolygonsDataFrame from the R package \strong{sp} v1.3-1 (see reference) containing shape file data for central Africa.
#'
#' @format
#' @source Obtained using malariaAtlas::getShp; see \url{https://github.com/artaylor85/pixelate/blob/master/data-raw/get_shape_files.R}.
#' @seealso \url{https://www.rdocumentation.org/packages/sp/versions/1.3-1/topics/SpatialPolygonsDataFrame-class}
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support
#'   for Spatial Vector Data. The R Journal 10 (1), 439-446,
#'   https://doi.org/10.32614/RJ-2018-009
"CentralAfrica_shp"

