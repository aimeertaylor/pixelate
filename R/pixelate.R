#==============================================================================
#' Pixelate as per average uncertainty
#'
#' Pixelate point estimates as per their average uncertainties.
#'
#' This is a wrapper function which, given a dot matrix and several arguements,
#' pixelates as follows. Note that, if there are few point estimates, many
#' combinations of arguments (if not all) will fail. Let a dot denote a set of
#' coordinates associated with a point estimate and its uncertainty, which is
#' represented by a single value, e.g. the width of the estimate's credible
#' interval. Let a pixel refer to a square or rectangle comprising one or more
#' dots. Uncertainties are averaged over a limited number of large pixels (bigk-th
#' pixel size, many dots per pixel). Average uncertainties are classified as
#' high, medium (with one or more subdivisions if argument bigk exceeds three), or
#' low. If the average uncertainty is high, point estimates within the large
#' pixel are averaged. If the average uncertainty is intermediate, point
#' estimates within the large pixel are averaged over smaller pixels (fewer dots
#' per pixel) that are nested within the large pixel. There are bigk-2 nested pixel
#' sizes in total. If the average uncertainty is low, point estimates are not
#' averaged.
#'
#' @param dot_matrix Data frame. Supplies dot coordinates, point estimates and
#'   credible interval bounds; see example data (once added).
#' @param num_pix_xy_bigk Integer vector length two. Specifies the number of large
#'   pixels (bigk-th pixel size) in the x and y direction
#' @param bigk Integer. Specifies the number of different uncertainity classes and
#'   thus pixel sizes.
#' @param scale String equal to either "linear" or "exponential". Specifies
#'   whether to scale pixel sizes (in units of dots) from class k = 2,...,bigk
#'   linearly or exponentially
#' @param scale_factor Integer. Specifies a factor by which to either multiply
#'   or exponentiate pixel sizes (in units of dots), depending on the scale.
#' @param verbose Logical. Set to FALSE to suppress messages
#' @return List. The list contains the original and expanded dot matrices
#'   complete with pixelated point estimates (variable Pix_z); bigk+1 quantiles of
#'   the distribution of average uncertainties used to classify average
#'   uncertainties; Dots per pixel for k = 1,...,bigk pixel sizes in the x and y
#'   direction; the dot matrix dimensions (in dots); and a matrix of dot
#'   memberships (the index of the pixel each dot belongs to for k = 1,...,bigk
#'   pixel sizes).
#' @export
#=============================================================================
pixelate <- function(dot_matrix,
                     num_pix_xy_bigk = 25,
                     bigk = 5L,
                     scale = "linear",
                     scale_factor = 1L,
                     verbose = TRUE) {

  warning("
 Please be aware, pixelate works by averaging uncertainty across predictions
 and ranking average uncertainty. Averaging uncertainty is only strictly valid
 when predictions are generated jointly. When they are not, e.g. by simulating
 from a 'per-pixel' posterior predictive distribution, uncertainty due to
 covariance between predictions is not accounted for. Regards pixelation, this
 omission will be consequential if and only if it changes the ranks of average
 uncertainty. The ranks of average uncertainty will change if covariance is both
 non-negligable in relation to other sources of uncertainty and spatially
 non-uniform.\n")

  # Set num_pix_xy_bigk in both x and y direction if not already
  if (is.na(num_pix_xy_bigk[2])) {num_pix_xy_bigk[2] <- num_pix_xy_bigk[1]}

  # Calculate dot dimensions of dot matrix -------------------------------------
  dot_matrix_dim <- apply(dot_matrix[, c("x", "y")], 2, function(j) {length(unique(j))})

  # Initial checks ------------------------------------------------------------
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

  # Calculate size in dots for pixel k = 2 -------------------------------------
  dpp_2 <- compute_dpp_2(num_pix_xy_bigk, bigk, scale, scale_factor, dot_matrix_dim)

  # Dot matrix compatibility check ---------------------------------------------
  dpp_min <- compute_dpp(2, bigk, scale, scale_factor)
  dot_req <- dpp_min[bigk,] * num_pix_xy_bigk

  if (any(dpp_2 < 2) | any(dot_matrix_dim < dot_req)) {
    stop(sprintf("
      Together, arguments num_pix_xy_bigk, bigk, scale and scale_factor are
      incompatible with the dot matrix dimensions. At least %s spatial predictions
      are required in the x and y direction for the arguments as currently specified.
      The dot matrix has %s in the x and y direction. Consider reducing
      num_pix_xy_bigk, bigk, scale_factor and/or using a linear scale.",
                 paste0(dot_req, collapse = ' and '), paste0(dot_matrix_dim, collapse = ' and ')))
  }

  # Calculate size in dots for pixel k = 1,...,bigk
  dpp <- compute_dpp(min(dpp_2), bigk, scale, scale_factor)

  # Expand the dot matrix to enable vectorisation
  expanded_dot_matrix <- expand_dot_matrix(dpp, dot_matrix, verbose)

  # Calculate dot dimensions of dot matrix
  expanded_dot_matrix_dim <- apply(expanded_dot_matrix[, c("x", "y")], 2, function(j) {length(unique(j))})

  # Calculate dot pixel membership
  dot_mem <- allocate_dot_mem(dpp, expanded_dot_matrix_dim)

  # Pixelate
  pix_output <- pixelate_by_u(expanded_dot_matrix, dot_mem, dpp, verbose)

  # Unpackage results
  pix_matrix <- pix_output$pix_matrix
  uncertainty_breaks <- pix_output$uncertainty_breaks

  # Crop Pix_matrix back to its original size
  pix_matrix_trimmed <- pix_matrix[(pix_matrix$x <= max(dot_matrix$x)) & (pix_matrix$y <= max(dot_matrix$y)), ]

  # Group together output
  to_return <- list(pix_matrix = pix_matrix_trimmed,
                    uncertainty_breaks = uncertainty_breaks,
                    pix_matrix_expanded = pix_matrix,
                    dpp = dpp,
                    dot_matrix_dim = dot_matrix_dim,
                    dot_mem = dot_mem)

  return(to_return)
}



