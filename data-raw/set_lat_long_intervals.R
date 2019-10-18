# Set the lat-log intervals of interest and save
lat_long_intervals <- data.frame(lon_min = c(-18,11.7),
                                 lon_max = c(52, 35.5),
                                 lat_min = c(-26.5,-14),
                                 lat_max = c(18, 5.7))
save(lat_long_intervals, file = "data-raw/lat_long_intervals.RData")


