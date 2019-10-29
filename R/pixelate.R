#==============================================================================
#' Pixelate as per average uncertainty
#'
#' Pixelate spatially continuous predictions as per their average uncertainties.
#'
#' This is a wrapper function which, given a dot matrix and several arguments,
#' pixelates as follows. Let a single dot denote a set containing a prediction,
#' its coordinates, and its uncertainty measure represented by a single value,
#' e.g. 95\% credible interval width. Let a pixel refer to a square or rectangle
#' comprising one or more predictions and thus dots. By default, pixels are
#' square. Uncertainties are averaged over a limited number of large pixels
#' (pixels of the bigk-th size). We specify a lower bound on the number of large
#' pixels. The function pixelate internally calculates the smallest number of
#' large pixels greater than or equal to the specified lower bound. By default
#' the lower bound is the same in both the x and y dimension. The actual number
#' of large pixels in either dimension is be determined by whichever dimension
#' is smaller (this behaviour is hard-coded within pixelate). Average
#' uncertainties are classified as high, medium (with bigk-2 subdivisions), or
#' low, according to the quantile interval they fall into, where the number of
#' quantile intervals is equal to a specified number of different pixel sizes (k
#' = 1,...,bigk) and the quantiles are based on the empirical distribution of
#' average uncertainties (calculated internally). For a given k, each pixel has
#' the same number of dots per pixel (dpp) in the x and y direction. However, we
#' do not specify dpps directly; they are calculated internally to best match
#' the specified parameters. The arguments scale and scale_factor determine the
#' rates at which dpps scale. If the average uncertainty is high (falls within
#' the top quantile interval), predictions within the large pixel are averaged.
#' If the average uncertainty is intermediate (falls with an intermediate
#' quantile interval), predictions within the large pixel are averaged over
#' smaller nested pixels. If the average uncertainty is low (falls within the
#' bottom quantile interval), predictions are not averaged. Importantly, missing
#' predictions and predictions that are zero with certainty are excluded from
#' the entire pixelation process (i.e. computation and classification of average
#' uncertainty, and computation of average prediction across large or nested
#' pixel sizes). This behaviour is hard-coded within a function that pixelate
#' calls internally.
#'
#' @param dot_matrix Data frame. Contains a row per dot with at least
#'   four variables: longitude, x; latitude, y; prediction, z; and uncertainty measure u.
#' @param num_pix_xy_bigk Integer or integer vector length two. Specifies a lower bound on the number of large
#'   pixels in the x and y direction.
#' @param bigk Integer. Specifies the number of average uncertainty quantile intervals and
#'   thus different pixel sizes to pixelate into.
#' @param scale String equal to either "linear" or "exponential". Specifies
#'   whether to scale pixel sizes (in units of dots) from class k = 2,...,bigk
#'   linearly or exponentially.
#' @param scale_factor Integer. Specifies a factor by which to either multiply
#'   or exponentiate pixel sizes (in units of dots), depending on the scale.
#' @return List.
#' \describe{
#'   \item{pix_matrix}{The original dot matrix with additional variables:
#'   average uncertainty, u_bigk; the average uncertainty quantile interval
#'   allocation, bins; and averaged predictions, pix_z.}
#'   \item{pix_matrix_expanded}{A spatially expanded dot matrix with additional
#'   variables: the average uncertainty, u_bigk; average uncertainty quantile
#'   interval allocation, bins; and averaged predictions, pix_z.}
#'   \item{uncertainty_breaks}{The values of average uncertainty at the bigk+1
#'   quantiles of the empirical distribution of average uncertainties.}
#'   \item{dpp}{The dots per pixel (dpp) for k = 1,...,bigk pixel sizes in the x
#'   and y direction.}
#'   \item{dot_matrix_dim}{The dimensions (in dots) of the
#'   original dot matrix.}
#'   \item{dot_mem}{A matrix of dot memberships, where each membership specifies
#'   which quantile interval that the large pixel containing the specified dot
#'   falls into.}
#'   \item{arguments}{The arguments passed to pixelate when it was called.}
#' }
#' @examples
#' #=================================================
#' # Use pixelate and inspect its output
#' #=================================================
#' # Pixelate using default parameters
#' px_def <- pixelate(SubSaharanAfrica_Pf_incidence)
#'
#' # Inspect list returned by pixelate
#' str(px_def)
#'
#' # Inspect a sample of uncertain pixelated predictions
#' uncertain_ind = which(px_def$pix_matrix$u > 0)
#' head(px_def$pix_matrix[uncertain_ind, ])
#'
#' # Pixelate using alternative parameters
#' px_alt <- pixelate(SubSaharanAfrica_Pf_incidence,
#'                    num_pix_xy_bigk = 25, bigk = 5)
#'
#'
#' #=================================================
#' # Plotting pixelate's output
#' #=================================================
#' # Load and attach ggplot2
#' if (!require("ggplot2")){
#'    stop("Package ggplot2 needed for the following code. Please install it.")
#' }
#'
#' # Define a plotting function
#' plot_sp_pred <- function(sp_pred){
#'
#'  ggplot(sp_pred) +
#'
#'    # Add raster surface
#'    geom_raster(mapping = aes(x = x, y = y, fill = pix_z)) +
#'
#'    # Add gradient
#'    scale_fill_gradientn(name = "Median incidence rate",
#'                         colors = c("seashell", "tomato", "darkred"),
#'                         na.value = 'lightblue') +
#'
#'    # Add axis labels
#'    ylab('Latitude (degrees)') +
#'    xlab('Longitude (degrees)') +
#'
#'    # Ensure the plotting space is not expanded
#'    coord_fixed(expand = FALSE) +
#'
#'    # Modify the legend and add a plot border:
#'    theme(legend.justification = c(0, 0),
#'          legend.position = c(0.02, 0.01),
#'          legend.background = element_rect(fill = NA),
#'          legend.title = element_text(size = 8),
#'          legend.text = element_text(size = 8),
#'          panel.border = element_rect(fill = NA))
#'
#' }
#'
#' # Plot default and alternative
#' plot_sp_pred(px_def$pix_matrix)
#' plot_sp_pred(px_alt$pix_matrix)
#'
#' @export
pixelate <- function(dot_matrix,
                     num_pix_xy_bigk = 15,
                     bigk = 6L,
                     scale = "linear",
                     scale_factor = 1L) {

  # Record arguments for reference
  arguments <- as.list(environment())

  warning("\n
        Please be aware, pixelate works by averaging uncertainty across predictions
        within large pixels. Averaging uncertainty is only valid when predictions are
        generated jointly. When they are not, e.g. by simulating from a 'per-pixel'
        posterior predictive distribution, uncertainty due to spatial covariance
        between predictions is not accounted for. Regards pixelation, this omission
        will be consequential if and only if it changes the quantile interval to
        which the average uncertainty of a given large pixel belongs. Quantile
        interval allocation is relatively robust, since it depends on the average
        uncertainty rank and can accommodate some variation around ranks within
        quantile intervals. Allocation will change if covariance is both
        non-negligible in relation to other sources of uncertainty and spatially
        non-uniform.\n")

  # Set dot_matrix to data.frame if not already
  if (class(dot_matrix) != "data.frame") {
    dot_matrix <- as.data.frame(dot_matrix)
  }

  # Initial check
  if (!setequal(c("u", "x", "y", "z"), ls(dot_matrix))) {
    stop("The dot matrix must contain variables 'u', 'x', 'y', and 'z' only.")
  }

  # Set num_pix_xy_bigk in both x and y direction if not already
  if (is.na(num_pix_xy_bigk[2])) {num_pix_xy_bigk[2] <- num_pix_xy_bigk[1]}

  # Calculate dot dimensions of dot matrix
  dot_matrix_dim <- apply(dot_matrix[, c("x", "y")], 2, function(j) {length(unique(j))})

  # Additional checks
  errors = c(dmat = ifelse (all(dot_matrix_dim >= c(2,4)) | all(dot_matrix_dim >= c(4,2)), FALSE, TRUE),
             npix = ifelse (all(num_pix_xy_bigk < 2) | nz_remainder(num_pix_xy_bigk), TRUE, FALSE),
             bigk = ifelse (bigk < 2 | nz_remainder(bigk), TRUE, FALSE),
             scale = ifelse (!scale %in% c("linear", "exponential"), TRUE, FALSE),
             scale_factor = ifelse(scale_factor < 1 | nz_remainder(scale_factor), TRUE, FALSE))

  error_mgs = c(dmat = "\nThe dot matrix must have at least two and four predictions in the x and y direction or vice versa.",
                npix = "\nThe number of bigK pixels must be integer and exceed one in a least one direction x or y.",
                bigk = "\nThe number of pixel sizes, bigK, must be an integer and exceed two.",
                scale = "\nUnrecognised scale: please provide either 'linear' or 'exponential'.",
                scale_factor = "\nThe scale_factor needs to be an integer and exceed zero.")

  if (any(errors)) stop(error_mgs[errors])

  # Calculate size in dots for pixel k = 2
  dpp_2 <- compute_dpp_2(num_pix_xy_bigk, bigk, scale, scale_factor, dot_matrix_dim)

  # Dot matrix compatibility check
  dpp_min <- compute_dpp(2, bigk, scale, scale_factor)
  dot_req <- dpp_min[bigk,] * num_pix_xy_bigk

  if (any(dpp_2 < 2) | any(dot_matrix_dim < dot_req)) {
    stop(sprintf("\n
  Together, arguments num_pix_xy_bigk, bigk, scale and scale_factor are
  incompatible with the dot matrix dimensions. At least %s
  spatial predictions are required in the x and y direction for the arguments
  as currently specified. The dot matrix has %s in the x and y direction.
  Consider reducing num_pix_xy_bigk, bigk, scale_factor and/or using a linear
  scale.\n", paste0(signif(dot_req), collapse = ' and '),
                 paste0(dot_matrix_dim, collapse = ' and ')))
  }

  # Calculate size in dots for pixel k = 1,...,bigk
  dpp <- compute_dpp(min(dpp_2), bigk, scale, scale_factor)

  # Expand the dot matrix to enable vectorisation
  expanded_dot_matrix <- expand_dot_matrix(dpp, dot_matrix)

  # Calculate dot dimensions of dot matrix
  expanded_dot_matrix_dim <- apply(expanded_dot_matrix[, c("x", "y")], 2, function(j) {length(unique(j))})

  # Calculate dot pixel membership
  dot_mem <- allocate_dot_mem(dpp, expanded_dot_matrix_dim)

  # Pixelate
  pix_output <- pixelate_by_u(expanded_dot_matrix, dot_mem, dpp)

  # Unpackage results
  pix_matrix <- pix_output$pix_matrix
  uncertainty_breaks <- pix_output$uncertainty_breaks

  # Crop Pix_matrix back to its original size
  pix_matrix_trimmed <- pix_matrix[(pix_matrix$x <= max(dot_matrix$x)) & (pix_matrix$y <= max(dot_matrix$y)), ]

  # Group together output
  to_return <- list(pix_matrix = pix_matrix_trimmed,
                    pix_matrix_expanded = pix_matrix,
                    uncertainty_breaks = uncertainty_breaks,
                    dpp = dpp,
                    dot_matrix_dim = dot_matrix_dim,
                    dot_mem = dot_mem,
                    arguments = arguments)

  return(to_return)
}
