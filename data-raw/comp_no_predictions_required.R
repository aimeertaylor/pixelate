# To access compute_opp in development model use devtools::load_all()

# Minimum number of spatial predictions required in either direction
num_pix_xy_bigk_min = c(1L,2L); bigk_min = 2L
opp_min <- compute_opp(opp_2 = 2L, bigk_min, scale = "linear", scale_factor = 1L)
opp_min[bigk_min,] * num_pix_xy_bigk_min

# Number of spatial predictions required given default parameters
num_pix_xy_bigk = formals(pixelate)$num_pix_xy_bigk
if (is.na(num_pix_xy_bigk[2])) {num_pix_xy_bigk[2] <- num_pix_xy_bigk[1]}
bigk = formals(pixelate)$bigk
scale = formals(pixelate)$scale
scale_factor = formals(pixelate)$scale_factor
opp_dft <- compute_opp(opp_2 = 2L, bigk, scale, scale_factor)
opp_dft[bigk,] * num_pix_xy_bigk
