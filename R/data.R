# Documenting datasets
# R packages by Hadley Wickham (http://r-pkgs.had.co.nz/data.html)
# use devtools::load_all() to load during package development

#' \emph{Plasmodium falciparum} predicted incidence rate
#'
#' \emph{P. falciparum} predicted all-age incidence rate (clinical cases per 1,000 population per annum)
#' in 2017 for sub-Saharan Africa at 2.5 arcminute (approximately 5km) resolution.
#'
#' @format An data frame with a row per prediction (obs) and four variables:
#' \describe{
#'   \item{x}{Longitude in decimal degrees}
#'   \item{y}{Latitude in decimal degrees}
#'   \item{z}{The median of a ‘per-pixel’ posterior predictive simulation of \emph{P. falciparum} incidence rate in 2017}
#'   \item{u}{The width of the 95\% credible interval of a ‘per-pixel’ posterior predictive simulation of \emph{P. falciparum}
#'   incidence rate in 2017.}
#' }
#' @source These data are available at the Malaria Atlas Project (MAP) website \url{https://map.ox.ac.uk/}.
#'    Specifically, they were obtained by selecting
#'    'ANNUAL MEAN OF PF INCIDENCE' at \url{https://map.ox.ac.uk/malaria-burden-data-download/}.
#' @references Weiss DJ, Lucas TCD, Nguyen M, et al.
#'    Mapping the global prevalence, incidence, and mortality of Plasmodium falciparum, 2000--17:
#'    a spatial and temporal modelling study. Lancet 2019; published online June 19.
#' @examples
#' str(SubSaharanAfrica_Pf_incidence)
#' head(SubSaharanAfrica_Pf_incidence)
"SubSaharanAfrica_Pf_incidence"


#' Shape files for sub-Saharan Africa
#'
#' An object of class SpatialPolygonsDataFrame from the R package \strong{sp} v1.3-1 containing shape file data for sub-Saharan Africa.
#'
#' @format
#' @source Obtained using malariaAtlas::getShp; see \url{https://github.com/artaylor85/pixelate/blob/master/data-raw/get_shape_files.R}.
#' @seealso \url{https://www.rdocumentation.org/packages/sp/versions/1.3-1/topics/SpatialPolygonsDataFrame-class}
"SubSaharanAfrica_shp"
