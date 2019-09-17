#============================================================
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
#' dots. Uncertainties are averaged over a limited number of large pixels (Kth
#' pixel size, many dots per pixel). Average uncertainties are classified as
#' high, medium (with one or more subdivisions if argument K exceeds three), or
#' low. If the average uncertainty is high, point estimates within the large
#' pixel are averaged. If the average uncertainty is intermediate, point
#' estimates within the large pixel are averaged over smaller pixels (fewer dots
#' per pixel) that are nested within the large pixel. There are K-2 nested pixel
#' sizes in total. If the average uncertainty is low, point estimates are not
#' averaged.
#'
#' @param dot_matrix Data frame. Supplies dot coordinates, point estimates and
#'   credible interval bounds; see example data (once added).
#' @param num_pix_xy_K Integer vector length two. Specifies the number of large
#'   pixels (Kth pixel size) in the x and y direction
#' @param K Integer. Specifies the number of different uncertainity classes and
#'   thus pixel sizes.
#' @param scale String equal to either "linear" or "exponential". Specifies
#'   whether to scale pixel sizes (in units of dots) from class k = 2,...,K
#'   linearly or exponentially
#' @param scale_factor Integer. Specifies a factor by which to either multiply
#'   or exponentiate pixel sizes (in units of dots), depending on the scale.
#' @return List. The list contains the original and expanded dot matrices
#'   complete with pixelated point estimates (variable Pix_z); K+1 quantiles of
#'   the distribution of average uncertainties used to classify average
#'   uncertainties; Dots per pixel for k = 1,...,K pixel sizes in the x and y
#'   direction; the dot matrix dimensions (in dots); and a matrix of dot
#'   memberships (the index of the pixel each dot belongs to for k = 1,...,K
#'   pixel sizes).
#' @export
#============================================================
pixelate = function(dot_matrix, num_pix_xy_K = c(10,10), K = 3L, scale = 'linear',
                    scale_factor = 1L, no.free_cores = 2){

  # Calculate dot dimensions of dot matrix
  dot_matrix_dim = apply(dot_matrix[,c('x','y')],2,function(j){length(unique(j))})

  # Calculate size in dots for pixel k = 2
  dpp_2 = compute_dpp_2(num_pix_xy_K, K, scale, scale_factor, dot_matrix_dim)

  # Calculate size in dots for pixel k = 1,...,K
  dpp = compute_dpp(min(dpp_2), K, scale, scale_factor)

  # Expand the dot matrix to enable vectorisation
  expanded_dot_matrix = expand_dot_matrix(dpp, dot_matrix)

  # Calculate dot dimensions of dot matrix
  expanded_dot_matrix_dim = apply(expanded_dot_matrix[,c('x','y')],2,function(j){length(unique(j))})

  # Calculate dot pixel membership
  dot_mem = allocate_dot_mem(dpp, expanded_dot_matrix_dim)

  # Pixelate
  Pix_output = pixelate_by_CI_width(expanded_dot_matrix, dot_mem, dpp, no.free_cores)

  # Unpackage results
  Pix_matrix = Pix_output$Pix_matrix
  Breaks = Pix_output$Breaks

  # Crop Pix_matrix back to its original size
  Pix_matrix_trimmed = Pix_matrix[(Pix_matrix$x <= max(dot_matrix$x)) & (Pix_matrix$y <= max(dot_matrix$y)),]

  # Group together output
  To_return = list(Pix_matrix = Pix_matrix_trimmed, Breaks = Breaks,
                   Pix_matrix_expanded = Pix_matrix,
                   dpp = dpp, dot_matrix_dim = dot_matrix_dim, dot_mem = dot_mem)

  return(To_return)
}
