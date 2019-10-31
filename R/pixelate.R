#==============================================================================
#' Pixelate as per average uncertainty
#'
#' Pixelate spatially continuous predictions according to the uncertainty that
#' surrounds them.
#'
#' This is a wrapper function which, given a data frame of observations and
#' several arguments, pixelates as follows.
#'
#' Let a single observation denote a set containing a prediction, its
#' coordinates, and its uncertainty represented by a single value, e.g. 95\%
#' credible interval width. Let a pixel refer to a square or rectangle
#' comprising one or more observations and thus predictions. By default, pixels
#' are square.
#'
#' Uncertainties are averaged over a limited number of large pixels (pixels of
#' the bigk-th size). We specify a lower bound on the number of large pixels.
#' The function pixelate internally calculates the smallest number of large
#' pixels greater than or equal to the specified lower bound. The lower bound
#' can either be an integer or integer vector length two. If a single integer is
#' specified, all pixels will be square and the actual number of pixels will be
#' calculated relative to the lower bound in the smallest dimension. This is the
#' default. If an integer vector of length two is specified, pixels will likely
#' be rectangular and calculated relative to the lower bounds in both
#' directions x then y.
#'
#' Average uncertainties are classified as high, intermediate (with bigk-2
#' subdivisions), or low, according to the quantile interval they fall into,
#' where the number of quantile intervals is equal to a specified number of
#' different pixel sizes (k = 1,...,bigk) and the quantiles are based on the
#' empirical distribution of average uncertainties.
#'
#' The k-th pixel size is defined by a count of observations per pixel (opp) in
#' the x and y direction. We do not specify opps directly; they are calculated
#' internally to best match the specified parameters. Arguments scale and
#' scale_factor determine the rate at which opps scale. There are two scales,
#' imult and iexpn. Both scale over k = 3,...,bigk for bigk > 2, because
#' \eqn{opp_1 = 1} always, and \eqn{opp_2} is calculated internally to best
#' match the specified parameters. imult specifies scaling by iterative
#' multiplication (i.e. a geometric series): \deqn{opp_k = opp_2 * (2 *
#' scale_factor)^(bigk-2).} iexpn specifies scaling by iterative exponentiation:
#' \deqn{opp_k = opp_2 ^ ((2 * scale_factor)^(bigk-2)).} The factor 2 is
#' necessary to ensure pixels nest within one another.
#'
#' If the average uncertainty is high (falls within the top quantile interval),
#' predictions within the large pixel are averaged. If the average uncertainty
#' is intermediate (falls with an intermediate quantile interval), predictions
#' are averaged across smaller pixels nested within the large pixel. If the
#' average uncertainty is low (falls within the bottom quantile interval),
#' predictions are not averaged (\eqn{opp_1 = 1}).
#'
#' Importantly, observations containing missing predictions and predictions that
#' are zero with certainty are excluded from the entire pixelation process (i.e.
#' computation and classification of average uncertainty, and computation of
#' average prediction across large or nested pixel sizes).
#'
#' @param obs_df Data frame. Contains a row per observation with four variables:
#'   longitude, x; latitude, y; prediction, z; and uncertainty measure u.
#' @param num_bigk_pix Integer or integer vector length two. Specifies a lower
#'   bound on the number of large pixels (pixels of the bigk-th size) in the x
#'   and y direction. If a single integer is specified, all pixels will be
#'   square and the actual number of pixels will be calculated relative to the
#'   lower bound in the smallest dimension. If an integer vector of length two
#'   is specified, pixels will likely be rectangular and calculated relative to
#'   the lower bounds in both directions x then y.
#' @param bigk Integer. Specifies the number of average uncertainty quantile
#'   intervals and thus different pixel sizes.
#' @param scale Character equal to either "imult" or "iexpn". Specifies whether
#'   to scale pixel sizes (in units of observations) from class k = 3,...,bigk
#'   by iterative multiplication or iterative exponentiation (see Details).
#' @param scale_factor Integer. Specifies a factor (in units of observations)
#'   that features in either iterative multiplication or iterative
#'   exponentiation (see Details)
#' @return pixelate returns a list.
#' \describe{
#'   \item{pix_df}{The original observation data frame with additional
#'   variables: average uncertainty, u_bigk; the average uncertainty quantile
#'   interval allocation, bins; and averaged predictions, pix_z.}
#'   \item{pix_df_expanded}{A spatially expanded observation data frame with
#'   additional variables: the average uncertainty, u_bigk; average uncertainty
#'   quantile interval allocation, bins; and averaged predictions, pix_z. All
#'   variables besides x and y are NA in spatially expanded observations.}
#'   \item{uncertainty_breaks}{The values of average uncertainty at the bigk+1
#'   quantiles of the empirical distribution of average uncertainties.}
#'   \item{opp}{The observations per pixel (opp) for k = 1,...,bigk pixel sizes
#'   in the x and y direction.}
#'   \item{obs_df_dim}{The dimensions (in units of observations) of the original
#'   observation data frame.}
#'   \item{obs_mem}{A data frame of observation memberships, where each
#'   membership specifies the quantile interval that the large pixel containing
#'   the specified observation falls into.}
#'   \item{arguments}{The arguments passed to pixelate when it was called.}
#'}
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
#' uncertain_ind = which(px_def$pix_df$u > 0)
#' head(px_def$pix_df[uncertain_ind, ])
#'
#' # Pixelate using alternative parameters
#' px_alt <- pixelate(SubSaharanAfrica_Pf_incidence,
#'                    num_bigk_pix = 25, bigk = 5)
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
#' plot_sp_pred(px_def$pix_df)
#' plot_sp_pred(px_alt$pix_df)
#'
#' @export
pixelate <- function(obs_df,
                     num_bigk_pix = 15,
                     bigk = 6,
                     scale = "imult",
                     scale_factor = 1) {

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

  # Set obs_df to data.frame if not already
  if (class(obs_df) != "data.frame") {
    obs_df <- as.data.frame(obs_df)
  }

  # Initial check
  if (!setequal(c("u", "x", "y", "z"), ls(obs_df))) {
    stop("The observation data frame must contain variables 'u', 'x', 'y', and 'z' only.")
  }

  # Set num_bigk_pix in both x and y direction if not already
  # Create a variable to ensure pixels are square
  if (is.na(num_bigk_pix[2])) {
    num_bigk_pix[2] <- num_bigk_pix[1]
    square_pix = TRUE
  } else {
    square_pix = FALSE
  }

  # Calculate observation dimensions of observation data frame
  obs_df_dim <- apply(obs_df[, c("x", "y")], 2, function(j) {length(unique(j))})

  # Additional checks
  errors = c(dmat = ifelse (all(obs_df_dim >= c(2,4)) | all(obs_df_dim >= c(4,2)), FALSE, TRUE),
             npix = ifelse (all(num_bigk_pix < 2) | nz_remainder(num_bigk_pix), TRUE, FALSE),
             bigk = ifelse (bigk < 2 | nz_remainder(bigk), TRUE, FALSE),
             scale = ifelse (!scale %in% c("imult", "iexpn"), TRUE, FALSE),
             scale_factor = ifelse(scale_factor < 1 | nz_remainder(scale_factor), TRUE, FALSE))

  error_mgs = c(dmat = "\nThe observation data frame must have at least two and four predictions in the x and y direction or vice versa.",
                npix = "\nThe number of bigK pixels must be integer and exceed one in a least one direction x or y.",
                bigk = "\nThe number of pixel sizes, bigK, must be an integer and exceed two.",
                scale = "\nUnrecognised scale: please provide either 'imult' or 'iexpn'.",
                scale_factor = "\nThe scale_factor needs to be an integer and exceed zero.")

  if (any(errors)) stop(error_mgs[errors])

  # Calculate size in observations for pixel k = 2
  opp_2 <- compute_opp_2(num_bigk_pix, bigk, scale, scale_factor, obs_df_dim)

  # Observation data frame compatibility check
  opp_min <- compute_opp(2, bigk, scale, scale_factor)
  obs_req <- opp_min[bigk,] * num_bigk_pix

  if (any(opp_2 < 2) | any(obs_df_dim < obs_req)) {
    stop(sprintf("\n
    Together, arguments num_bigk_pix, bigk, scale and scale_factor are
    incompatible with the observation data frame dimensions.
    At least %s spatial predictions are required in
    the x and y direction for the arguments as currently specified. The
    observation data frame has %s in the x and y
    direction. Consider reducing num_bigk_pix, bigk, scale_factor and/or using an
    imult scale.\n", paste0(signif(obs_req), collapse = ' and '),
                 paste0(signif(obs_df_dim), collapse = ' and ')))
  }

  # Determine whether to use one or both elements of opp
  if (square_pix) {
    opp_2_to_use = min(opp_2)
  } else {
    opp_2_to_use = opp_2
  }

  # Calculate size in observations for pixel k = 1,...,bigk
  opp <- compute_opp(opp_2_to_use, bigk, scale, scale_factor)

  # Expand the observation data frame to enable vectorisation
  expanded_obs_df <- expand_obs_df(opp, obs_df)

  # Calculate observation dimensions of observation data frame
  expanded_obs_df_dim <- apply(expanded_obs_df[, c("x", "y")], 2, function(j) {length(unique(j))})

  # Calculate observation pixel membership
  obs_mem <- allocate_obs_mem(opp, expanded_obs_df_dim)

  # Pixelate
  pix_output <- pixelate_by_u(expanded_obs_df, obs_mem, opp)

  # Unpackage results
  pix_df <- pix_output$pix_df
  uncertainty_breaks <- pix_output$uncertainty_breaks

  # Crop Pix_df back to its original size
  pix_df_trimmed <- pix_df[(pix_df$x <= max(obs_df$x)) & (pix_df$y <= max(obs_df$y)), ]

  # Group together output
  to_return <- list(pix_df = pix_df_trimmed,
                    pix_df_expanded = pix_df,
                    uncertainty_breaks = uncertainty_breaks,
                    opp = opp,
                    obs_df_dim = obs_df_dim,
                    obs_mem = obs_mem,
                    arguments = arguments)

  return(to_return)
}
