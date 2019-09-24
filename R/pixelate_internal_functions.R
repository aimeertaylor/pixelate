#=================================================================================
# Compute dots per pixel for pixel size k = 2 (smallest, most-resolved, multidot
# pixel size)
#
# This function essentially back engineers the dots per pixel for pixel size k =
# 2 given a target number of pixels of the K-th size (largest, least resolved
# pixels) in the x and y direction, which is specified by num_pix_xy_bigk. Note
# that we use bigk rather than K in the code to avoid problems due to case insensitivity.
#
# Rationale: for a given dot matrix, it is easy to envision a target number of
# least-resolved large pixels. We thus specifiy this target, then back
# calculate the dpp for pixel size k = 2, which is the smallest, most resolved,
# multi-dot pixel size.
#=================================================================================
compute_dpp_2 <- function(num_pix_xy_bigk, # bigk pixels in x and y direction
                          bigk, scale, scale_factor, # see compute_ddp
                          dot_matrix_dim) { # dot matrix dimensions

  # Compute the dpp for k = bigk in the x and y direction
  # Use floor to allow for some partial pixels at edges
  dpp_bigk <- ceiling(dot_matrix_dim / num_pix_xy_bigk)

  # Compute the dpp for k = 2 in the x and y direction
  if (scale == "linear") {
    dpp_2 <- sapply(floor(dpp_bigk / ( (2 * scale_factor) ^ (bigk - 2))), as.integer)
  } else if (scale == "exponential") {
    dpp_2 <- sapply(floor(dpp_bigk ^ (1 / ( (2 * scale_factor) ^ (bigk - 2)))), as.integer)
  } else {
    stop("Scale must either be 'linear' or 'exponential'")
  }

  # Compatibility check
  if (any(dpp_2 < 2)) {
    stop("Inputs (num_pix_xy_bigk, bigk, scale, scale_factor) are incompatible with the dot matrix.
          Consider a linear scale if not already, and/or reducing num_pix_xy_bigk, bigk, scale_factor.")
  }

  return(sapply(dpp_2, as.integer))
}



#=================================================================================
# Compute dots per pixel
#
# This function computes dots per pixel (dpp) for k = 3,...,bigk differently sized
# pixels and returns the dpps for all k = 1,...,bigk.
#
# The dpp for pixel size k = 1 is one in both the x and y direction by default.
# The dpp for pixel size k = 2 is specified by the user using dpp_2, which is
# either an integer or an integer vector length two. If the former, pixels will
# be square. This function doesn't depend on the input dot matrix. To ensure
# that all dpps are compatible with a specified dot matrix, use compute_dpp_2 to
# compute the dpp for pixel size k = 2. Given the dpp for pixel size k = 2, the
# dpps for k = 3,...,bigk are computed.
#=================================================================================
compute_dpp <- function(dpp_2, # Dots per pixel for k = 2
                        bigk, # Number of different pixel sizes, e.g. 3L
                        scale, # scale i.e. 'exponential' or 'linear'
                        scale_factor) { # scale_factor e.g. 1L

  if (is.null(dpp_2[2])) {dpp_2[2] <- dpp_2[1]} # Set dpp_2 in y direction

  #----------------------------------------------------------
  # Checks
  #----------------------------------------------------------
  if (!scale %in% c("linear", "exponential")) {
    stop("Unrecognised scale: please provide either 'linear' or 'exponential'")
  }
  if (bigk < 2 | !is.integer(bigk)) {
    stop("bigk needs to be an integer greater than or equal to two, e.g. 2L")
  }
  if (any(dpp_2 < 2) | any(!sapply(dpp_2, is.integer))) {
    stop("dots per pixel for k = 2 in the x and y direction need to be integers
         greater than or equal to two, e.g. c(2L,2L)")
  }
  if (scale_factor < 1 | !is.integer(scale_factor)) {
    stop("scale_factor needs to be an integer greater than or equal to one, e.g. lL")
  }


  #----------------------------------------------------------
  # Compute dpp for k = 1,...,bigk
  #----------------------------------------------------------
  dpp <- array(NA, dim = c(bigk, 2), dimnames = list(1:bigk, c("x", "y"))) # Create matrix
  dpp[1, ] <- c(1, 1) # Always set dpp for k = 1 to one
  dpp[2, ] <- dpp_2 # Set dpp for k = 2 (smallest, most resolved, multi-dot pixel size)

  # Set subsequent multi-dot pixel levels, k = 3,..,bigk
  if (scale == "linear") {
    if (bigk > 2) {
      for (i in 3:bigk) {
        dpp[i, ] <- dpp[i - 1, ] * 2 * scale_factor # inc. 2 s.t. always divisible by 2
      }
    }
  } else if (scale == "linear") {
    if (bigk > 2) {
      for (i in 3:bigk) {
        dpp[i, ] <- dpp[i - 1, ] ^ ( 2 * scale_factor)  # inc. 2 s.t. always divisible by 2
      }
    }
  } else {
    stop("scale should either be 'linear' or 'exponential'")
  }

  #----------------------------------------------------------
  # Check all rows are nested
  #----------------------------------------------------------
  for (i in bigk:2) {
    if (!all(dpp[i, ] %% dpp[i - 1, ] == 0)) {stop("Pixels are not nested")}
  }

  #----------------------------------------------------------
  # End of function
  #----------------------------------------------------------
  return(data.frame(dpp))
}




#=================================================================================
# Compute number of pixels of in the x and y direction
#
# This function computes the number of pixels in the x and y direction for pixel
# sizes k = 1,...,bigk depending on the dot dimensions of the dot matrix
#
# Partial pixels are include (i.e.) non-integer pixel counts are rounded up.
#
# dpp specifies dots per pixel for k = 1,...,bigk pixel sizes in the x and y
# direction dot_matrix_dim specifies the dimensions of the dot matrix in the x
# and y direction.
#=================================================================================
compute_num_pix_xy <- function(dpp, dot_matrix_dim) {

  # Note that dot_matrix_dim/dpp returns the wrong ordering
  # due to R's treatment of rows and columns
  num_pix_xy <- ceiling(data.frame(x = dot_matrix_dim["x"] / dpp$x,
                                   y = dot_matrix_dim["y"] / dpp$y))
  return(num_pix_xy)
}


#=================================================================================
# This function expands the dimension of the dot matrix with NA dots
#
# The dot matrix is expanded such that, for pixel sizes k = 1,...,bigk, an integer
# number of pixels fit along both the x and y directions. Despite enlarging the
# data set, expansion enables vectorised pixelation (since the expanded data set
# has no partial pixels on the borders, whereas the original dot matrix may
# have), which is considerably faster than looping over pixels.
#=================================================================================
expand_dot_matrix <- function(dpp, dot_matrix) {

  bigk <- nrow(dpp)

  # Calculate the dimensions of the original dot matrix
  dot_matrix_dim <- apply(dot_matrix[, c("x", "y")], 2, function(j) {length(unique(j))})

  # Calculate the number of dots required to round up to complete pixels
  to_expand_xy <- compute_num_pix_xy(dpp, dot_matrix_dim)[bigk, ] * dpp[bigk, ] - dot_matrix_dim

  # Check expansion does not a single pixel
  if (!all(to_expand_xy < dpp[bigk, ])) {stop("Expansion exceeds one x or y pixel")}

  # Extract the current x and y coordinates and the maximums
  unq_x <- sort(unique(dot_matrix$x))
  unq_y <- sort(unique(dot_matrix$y))
  max_x <- max(unq_x)
  max_y <- max(unq_y)

  # State expansion if necessary
  if (any(to_expand_xy > 0)) {writeLines("\nExpanding dot matrix...")}

  # Calculate the extra x and y coordinates
  # Note that it doesn't really matter what they are as we will delete
  # them before returning the final pixelated matrix. However, it is necessary
  # to ensure the additional coordinates are distinct from the existing ones.
  # For this we add 1:to_expand_xy to max x or y.
  if (all(to_expand_xy > 0)) {
    square <- expand.grid(x = max_x + 1:to_expand_xy$x, y = max_y + 1:to_expand_xy$y)
    lower_left_rectangle <- expand.grid(x = unq_x, y = max_y + 1:to_expand_xy$y)
    upper_right_rectangle <- expand.grid(x = max_x + 1:to_expand_xy$x, y = unq_y)
    extra <- rbind(square, lower_left_rectangle, upper_right_rectangle)
  } else if (to_expand_xy$x > 0 & to_expand_xy$y == 0) {
    extra <- expand.grid(x = max_x + 1:to_expand_xy$x, y = unq_y)
  } else if (to_expand_xy$y > 0 & to_expand_xy$x == 0) {
    extra <- expand.grid(x = unq_x, y = max_y + 1:to_expand_xy$y)
  } else {
    return(dot_matrix)
  }

  # Add extra colunms: this step assumes dot matrix has these columns
  extra$z <- NA
  extra$u <- NA

  # Return expanded matrix
  expanded_dot_matrix <- rbind(dot_matrix, extra)

  # Reorder for accordance with membership allocation
  expanded_dot_matrix <- dplyr::arrange(expanded_dot_matrix, y, x)

  # End of function
  return(expanded_dot_matrix)
}



#=================================================================================
# Allocate dot pixel membership
#
# For each pixel size k = 1,...,bigk, allocate a pixel membership to each dot
#
# This function assumes coordinates in the dot matrix are sorted by y and then
# by x, e.g. dot_matrix = cbind(x = c(1,2,3,1,2,3), y = c(1,1,1,2,2,2)). The
# function trims the allocated memberships to match the input dot_matrix_dim.
# However, trim is redundant if dot_matrix_dim is based on an expanded dot
# matrix.
#
# dpp specifies dots per pixel for k = 1,...,bigk pixel sizes in the x and y
# direction dot_matrix_dim specifies the dimensions of the dot matrix in the x
# and y direction. It returns a N by bigk matrix of pixel memberships where N is
# the number of dots in the dot matrix (equal to prod(dot_matrix_dim)).
#=================================================================================
allocate_dot_mem <- function(dpp, dot_matrix_dim) {

  bigk <- nrow(dpp)
  num_pix_xy <- compute_num_pix_xy(dpp, dot_matrix_dim)

  # Membership vector
  dot_mem <- apply(array(1:bigk, dim = c(bigk, 1)), 1, function(k) {
    # For each k:
    dot_mem_k <- as.vector(sapply(1:num_pix_xy$y[k], function(i) {
      # Create x row trimmed to data set width
      row_i <- rep(1:num_pix_xy$x[k] + num_pix_xy$x[k] * (i - 1), each = dpp$x[k])[1:dot_matrix_dim["x"]]
      rep(row_i, dpp$y[k]) # duplicate row_i in y-direction
    }))
    dot_mem_k <- dot_mem_k[1:prod(dot_matrix_dim)] # trim to data set height
  })

  # Check pixel memberships are correct
  if (!all(apply(dot_mem, 2, max, na.rm = T) == apply(num_pix_xy, 1, prod))) {stop("Something is wrong")}

  return(dot_mem)
}


#============================================================
# Function to compute CI widths averaged over over k = bigk pixels
# and pixelate.
#
# To mitigate differring extent of CLT over different dpp
# and to ensure pixel variation, we determine breaks using the
# distribution of CI widths averaged at the bigk-th pixel size and
# no smaller.
#============================================================
pixelate_by_u <- function(dot_matrix, dot_mem, dpp) {

  # Arrange s.t. compatible with membership allocation
  dot_matrix <- dplyr::arrange(dot_matrix, y, x)

  # Extract numbers of pixels for k = 1,...,bigk
  num_pix <- apply(dot_mem, 2, function(p) {length(unique(p))})
  bigk <- length(num_pix)

  # Temporarily remove dots where the median value is zero with certainty
  # to prevent biases to the pixelation process
  # Aside, typically there are no z > 0 with u == 0; check by:
  # any(dot_matrix$z > 0 & dot_matrix$u == 0, na.rm = T)
  dot_certain_zero <- which(dot_matrix$z == 0 & dot_matrix$u == 0)
  dot_matrix$z[dot_certain_zero] <- NA
  dot_matrix$u[dot_certain_zero] <- NA

  # Compute CI widths for pixels where k = bigk
  writeLines("\nComputing CI widths averaged over dpp for k = bigk")
  dot_matrix$u_bigk <- NA
  # Sort by dot_mem: allows vectorisation via temp_matrix since all pixels have the same dpp
  sorted_dot_mem_bigk <- sort.int(dot_mem[, bigk], index.return = T)
  # Create a temporary matrix with pixel entries per column
  temp_matrix <- matrix(dot_matrix$u[sorted_dot_mem_bigk$ix], nrow = prod(dpp[bigk, ]))
  # temp_matrix enables calculation of us_bigk using colMeans(), which calls C code directly
  us_bigk <- colMeans(temp_matrix, na.rm = T)
  # Expand us_bigk such that there is one per dot and allocate using indices
  dot_matrix$u_bigk[sorted_dot_mem_bigk$ix] <- rep(us_bigk, each = prod(dpp[bigk, ]))

  # Compute breaks based on quantiles of CI width for k = bigk
  uncertainty_breaks <- quantile(dot_matrix$u_bigk,
                                 probs = seq(0, 1, length.out = bigk + 1),
                                 na.rm = TRUE)

  # Allocate pixel types (.bincode is fast for assigning numeric bins)
  dot_matrix$bins <- .bincode(dot_matrix$u_bigk, breaks = uncertainty_breaks, include.lowest = T)

  # First populate with NA
  dot_matrix$pix_z <- NA

  # For k = 1, use dot estimate. Note that which(.) excludes NA bins
  dot_matrix$pix_z[which(dot_matrix$bins == 1)] <- dot_matrix$z[which(dot_matrix$bins == 1)]

  writeLines("\nComputing point estimates for k = 2,...,bigk")

  #Obsolete code for parallelising over k = 1,...,bigk
  #Parallisation did not speed code and introduced a bug
  #Dependency: library(doParallel)
  #registerDoParallel(cores = detectCores()-no.free_cores)
  #foreach(k = 1:bigk) %dopar%

  for (k in 1:bigk) {
    # Sort by dot_mem: allows vectorisation via temp_matrix since all pixels have the same dpp
    dot_inds_k <- which(dot_matrix$bins == k) # which discards NA
    sorted_dot_mem_k <- sort.int(dot_mem[dot_inds_k, k], index.return = T)
    #----------------------------------------
    # Aside, uncomment to check ordering is by pixel:
    # temp_matrix = matrix(dot_mem[dot_inds_k,k][sorted_dot_mem_k$ix], nrow = prod(dpp[k,]))
    # all(colMeans(temp_matrix, na.rm = T) == unique(dot_mem[,k][which(dot_matrix$bins == k)]))
    #----------------------------------------
    # Create a temporary matrix with pixel entries per column
    temp_matrix <- matrix(dot_matrix$z[dot_inds_k][sorted_dot_mem_k$ix], nrow = prod(dpp[k, ]))
    # temp_matrix enables calculation of pixelated point estimates using colMeans(), which calls C code directly
    pix_z_temp <- colMeans(temp_matrix, na.rm = T)
    # Expand us_bigk such that there is one per dot and allocate using indices
    dot_matrix$pix_z[dot_inds_k][sorted_dot_mem_k$ix] <- rep(pix_z_temp, each = prod(dpp[k, ]))
  }

  # Reset NAs: some NA dots have non NA values due to mean(us, na.rm = T)
  dot_matrix$pix_z[is.na(dot_matrix$z)] <- NA

  # Re-addition of certain zeros
  dot_matrix$z[dot_certain_zero] <- 0
  dot_matrix$u[dot_certain_zero] <- 0
  dot_matrix$pix_z[dot_certain_zero] <- 0

  # Return unercetainty breaks as well as pixelated matrix
  to_return <- list(pix_matrix = dot_matrix,
                    uncertainty_breaks = uncertainty_breaks)

  return(to_return)
}
