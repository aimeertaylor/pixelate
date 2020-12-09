#=================================================================================
# Compute observations per pixel (opp) for pixel size k = 2 (smallest, most
# resolved, multi-observation pixel size)
#
# This function essentially back calculates the observations per pixel for pixel
# size k = 2 given a target number of complete pixels of the bigk-th size
# (largest, least resolved pixel size) in the x and y direction, which is
# specified by num_bigk_pix. Note that we use bigk rather than K to avoid
# problems with case insensitivity.
#
# Rationale: for a given data frame of observations, it is easy to envision a
# target number of least resolved large pixels. We thus specify this target,
# then back calculate the opp for pixel size k = 2, which is the smallest, most
# resolved, multi-observation pixel size.
#=================================================================================
compute_opp_2 <- function(num_bigk_pix, # bigk pixels in x and y direction
                          bigk, scale, scale_factor, # pixelate arguments
                          obs_df_dim) { # observation data frame dimensions

  # Compute the opp for k = bigk in the x and y direction
  opp_bigk <- obs_df_dim / num_bigk_pix # floor not needed since features below

  # Compute the opp for k = 2 in the x and y direction
  if (scale == "imult") {
    opp_2 <- sapply(floor(opp_bigk / ( (2 * scale_factor) ^ (bigk - 2))), as.integer)
  } else if (scale == "iexpn") {
    opp_2 <- sapply(floor(opp_bigk ^ (1 / ( (2 * scale_factor) ^ (bigk - 2)))), as.integer)
  } else {
    stop("Scale must either be 'imult' or 'iexpn'")
  }

  return(sapply(opp_2, as.integer))
}



#=================================================================================
# Compute observations per pixel (opp)
#
# This function computes opp for k = 3,...,bigk different pixel sizes and
# returns opp for k = 1,...,bigk.
#
# The opp for pixel size k = 1 is one in both the x and y direction by design.
# The opp for pixel size k = 2 is specified by opp_2. This function does not
# depend on the input observation data frame. To ensure that all opps are
# compatible with a specified observation data frame, use compute_opp_2 to
# compute the opp for pixel size k = 2. Given the opp for pixel size k = 2, the
# opps for k = 3,...,bigk are computed.
#=================================================================================
compute_opp <- function(opp_2, # Observations per pixel for k = 2
                        bigk, # Number of different pixel sizes, e.g. 3L
                        scale, # scale i.e. 'imult' or 'iexpn'
                        scale_factor) { # scale_factor e.g. 1L

  # Compute opp for k = 1,...,bigk
  opp <- array(NA, dim = c(bigk, 2), dimnames = list(1:bigk, c("x", "y"))) # Create array
  opp[1, ] <- c(1, 1) # Always set opp for k = 1 to one
  opp[2, ] <- opp_2 # Set opp for k = 2 (smallest, most resolved, multi-observation pixel size)

  # Set subsequent multi-observation pixel sizes, k = 3,..,bigk
  if (scale == "imult") {
    if (bigk > 2) {
      for (i in 3:bigk) {
        opp[i, ] <- opp[i - 1, ] * 2 * scale_factor # inc. 2 s.t. always divisible by 2
      }
    }
  } else if (scale == "iexpn") {
    if (bigk > 2) {
      for (i in 3:bigk) {
        opp[i, ] <- opp[i - 1, ] ^ ( 2 * scale_factor)  # inc. 2 s.t. always divisible by 2
      }
    }
  } else {
    stop("scale should either be 'imult' or 'iexpn'")
  }

  # Check all rows are nested
  for (i in bigk:2) {
    if (!all(opp[i, ] %% opp[i - 1, ] == 0)) {stop("Pixels are not nested")}
  }

  return(data.frame(opp))
}




#=================================================================================
# Compute the number of pixels in the x and y direction
#
# This function computes the number of pixels in the x and y direction for pixel
# sizes k = 1,...,bigk depending on the observation dimensions of the observation
# data frame
#
# Partial pixels are included (i.e. we take the ceiling of non-integer pixel counts).
#=================================================================================
compute_num_pix_xy <- function(opp, obs_df_dim) {

  # Note that obs_df_dim/opp returns the wrong ordering
  # due to R's treatment of rows and columns hence explicit division
  num_pix_xy <- ceiling(data.frame(x = obs_df_dim["x"] / opp$x,
                                   y = obs_df_dim["y"] / opp$y))
  return(num_pix_xy)
}


#=================================================================================
# This function expands the dimension of the observation data frame with NA
# observations
#
# The observation data frame is expanded such that, for pixel sizes k =
# 1,...,bigk, an integer number of pixels fit along both the x and y directions.
# Despite enlarging the data set, expansion increases the speed of pixelation.
# This is because it enables vectorised pixelation (since the expanded data set
# has no partial pixels on the borders, whereas the original observation data
# frame may have), which is considerably faster than looping over pixels.
#=================================================================================
expand_obs_df <- function(opp, obs_df) {

  # Extract the number of pixel sizes
  bigk <- nrow(opp)

  # Calculate the dimensions of the original observation data frame
  obs_df_dim <- apply(obs_df[, c("x", "y")], 2, function(j) {length(unique(j))})

  # Calculate the number of observations required to round up to complete large pixels
  to_expand_xy <- compute_num_pix_xy(opp, obs_df_dim)[bigk, ] * opp[bigk, ] - obs_df_dim

  # Check expansion does not exceed a single large pixel on either edge
  if (!all(to_expand_xy < opp[bigk, ])) {stop("Expansion exceeds one x or y pixel")}

  # Extract the current x and y coordinates and their maximums
  unq_x <- sort(unique(obs_df$x))
  unq_y <- sort(unique(obs_df$y))
  max_x <- max(unq_x)
  max_y <- max(unq_y)

  # State expansion if necessary
  if (any(to_expand_xy > 0)) {message("\nExpanding the input observation data frame...")}

  # Calculate the extra x and y coordinates
  # Note that it does not really matter what they are as we will delete
  # them before returning the final pixelated data frame. However, it is necessary
  # to ensure the additional coordinates are distinct from the existing ones.
  # To ensure they are distinct, we add 1:to_expand_xy to max x or y.
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
    return(obs_df)
  }

  # Add extra columns
  # (this step assumes the oberservation data frame has only variables x, y, z and u)
  extra$z <- NA
  extra$u <- NA

  # Return expanded data frame
  expanded_obs_df <- rbind(obs_df, extra)

  # Reorder for accordance with membership allocation
  expanded_obs_df <- dplyr::arrange(expanded_obs_df, y, x)

  # End of function
  return(expanded_obs_df)
}



#=================================================================================
# Allocate observation pixel membership
#
# For each pixel size k = 1,...,bigk allocate a pixel membership to each observation
#
# This function assumes coordinates in the observation data frame are sorted by y and then
# by x, e.g. obs_df[,("x","y")] = cbind(x = c(1,2,3,1,2,3), y = c(1,1,1,2,2,2)). The
# function trims the allocated memberships to match the input obs_df_dim.
# Trim is redundant if obs_df_dim is based on an expanded observation data frame.
#
# opp specifies observations per pixel for k = 1,...,bigk pixel sizes in the x and y
# direction. obs_df_dim specifies the dimensions of the observation data frame in the x
# and y direction. allocate_obs_mem returns a N by bigk data frame of pixel memberships
# where N is the number of observations in the observation data frame (equal to prod(obs_df_dim)).
#=================================================================================
allocate_obs_mem <- function(opp, obs_df_dim) {

  bigk <- nrow(opp)
  num_pix_xy <- compute_num_pix_xy(opp, obs_df_dim)

  # Membership vector
  obs_mem <- apply(array(1:bigk, dim = c(bigk, 1)), 1, function(k) {
    obs_mem_k <- as.vector(sapply(1:num_pix_xy$y[k], function(i) {
      # Create x row trimmed to data set width
      row_i <- rep(1:num_pix_xy$x[k] + num_pix_xy$x[k] * (i - 1), each = opp$x[k])[1:obs_df_dim["x"]]
      rep(row_i, opp$y[k]) # duplicate row_i in y-direction
    }))
    obs_mem_k <- obs_mem_k[1:prod(obs_df_dim)] # trim to data set height
  })

  # Check pixel memberships are correct
  if (!all(apply(obs_mem, 2, max, na.rm = T) == apply(num_pix_xy, 1, prod))) {stop("Something is wrong")}

  return(obs_mem)
}

#============================================================
# Function to compute uncertainty breaks
#============================================================
compute_breaks <- function(sample, interval_count, type = "equally_occupied")  {
  if (type == "equally_occupied") {

    breaks <- quantile(sample, na.rm = TRUE,
                       probs = seq(0, 1, length.out = interval_count + 1))

  } else if (type == "equally_spaced") {

    breaks <- seq(min(sample, na.rm = TRUE),
                  max(sample, na.rm = TRUE), length.out = interval_count + 1)

    # Add probability names (for consistency with equally_occupied)
    names(breaks) <- paste0(round(stats::ecdf(sample)(breaks)*100,5),0)

  } else {

    stop("type must be either 'equally_occupied' or 'equally_spaced'")

  }

  return(breaks)
}


#============================================================
# Function to pixelate as per average uncertainty averaged over
# large pixels (pixels of the bigk-th size)
#============================================================
pixelate_by_u <- function(obs_df, obs_mem, opp, interval_type) {

  # Arrange s.t. compatible with membership allocation
  obs_df <- dplyr::arrange(obs_df, y, x)

  # Extract numbers of pixels for k = 1,...,bigk
  num_pix <- apply(obs_mem, 2, function(p) {length(unique(p))})
  bigk <- length(num_pix)

  # Temporarily remove observations that are zero with certainty
  # Aside, typically there are no z > 0 with u == 0; check by:
  if (any(obs_df$z > 0 & obs_df$u == 0, na.rm = T)) {
    message("Note that there are some certain yet non-zero predictions")
  }
  obs_certain_zero <- which(obs_df$z == 0 & obs_df$u == 0)
  obs_df$z[obs_certain_zero] <- NA
  obs_df$u[obs_certain_zero] <- NA

  # Compute average uncertainty
  message("Averaging uncertainty over large pixels...")
  # Allocate a vector for average uncertainty
  obs_df$u_bigk <- NA
  # Sort by obs_mem: allows vectorisation via temp_df since all pixels have the same opp
  sorted_obs_mem_bigk <- sort.int(obs_mem[, bigk], index.return = T)
  # Create a temporary data frame with pixel entries per column
  temp_df <- matrix(obs_df$u[sorted_obs_mem_bigk$ix], nrow = prod(opp[bigk, ]))
  # temp_df enables calculation of us_bigk using colMeans, which is fast
  us_bigk <- colMeans(temp_df, na.rm = T)
  # Expand us_bigk such that there is one per observation and allocate using indices
  obs_df$u_bigk[sorted_obs_mem_bigk$ix] <- rep(us_bigk, each = prod(opp[bigk, ]))

  # Compute breaks based on quantiles of average uncertainty
  uncertainty_breaks <- compute_breaks(sample = obs_df$u_bigk,
                                       interval_count = bigk,
                                       type = interval_type)

  message("Allocating observations to different pixel sizes...")
  # Allocate pixel types (.bincode is fast for assigning numeric bins)
  obs_df$bins <- .bincode(obs_df$u_bigk, breaks = uncertainty_breaks, include.lowest = T)

  # First populate with NA
  obs_df$pix_z <- NA
  obs_df$pix_u <- NA # For uncertainty averaged over nested pixels

  # For k = 1, use observation. Note that which(.) excludes NA bins
  obs_df$pix_z[which(obs_df$bins == 1)] <- obs_df$z[which(obs_df$bins == 1)]
  obs_df$pix_u[which(obs_df$bins == 1)] <- obs_df$u[which(obs_df$bins == 1)]

  message("Averaging predictions over differently sized pixels...\n")

  for (k in 1:bigk) {
    # Sort by obs_mem: allows vectorisation via temp_df since all pixels have the same opp
    obs_inds_k <- which(obs_df$bins == k) # which() discards NA
    sorted_obs_mem_k <- sort.int(obs_mem[obs_inds_k, k], index.return = T)

    # Create a temporary data frame with pixel entries per column
    temp_df <- matrix(obs_df$z[obs_inds_k][sorted_obs_mem_k$ix], nrow = prod(opp[k, ]))
    # temp_df enables calculation of averaged predictions using colMeans, which is fast
    pix_z_temp <- colMeans(temp_df, na.rm = T)
    # Expand us_bigk such that there is one per observation and allocate using indices
    obs_df$pix_z[obs_inds_k][sorted_obs_mem_k$ix] <- rep(pix_z_temp, each = prod(opp[k, ]))

    # Repeat three steps above but for uncertainty
    temp_df <- matrix(obs_df$u[obs_inds_k][sorted_obs_mem_k$ix], nrow = prod(opp[k, ]))
    pix_u_temp <- colMeans(temp_df, na.rm = T)
    obs_df$pix_u[obs_inds_k][sorted_obs_mem_k$ix] <- rep(pix_u_temp, each = prod(opp[k, ]))
  }

  # Reset NAs. Otherwise, some NA observations will have non NA values due to mean(us, na.rm = T)
  obs_df$pix_z[is.na(obs_df$z)] <- NA
  obs_df$pix_u[is.na(obs_df$u)] <- NA

  # Re-addition of certain zeros
  obs_df$z[obs_certain_zero] <- 0
  obs_df$u[obs_certain_zero] <- 0
  obs_df$pix_z[obs_certain_zero] <- 0
  obs_df$pix_u[obs_certain_zero] <- 0

  # Return uncertainty breaks as well as pixelated data frame
  to_return <- list(pix_df = obs_df,
                    uncertainty_breaks = uncertainty_breaks)

  return(to_return)
}


#============================================================
# Function to check arguments are integers
#
# We use nz_remainder instead of is.integer because most people
# are likely to specify an integer as a numeric rather than integer
# e.g. bigk = 3 rather than bigk = 3L
#============================================================
nz_remainder = function(x){
  y = !all(x %% sapply(x, floor) == 0)
  return(y)
}
