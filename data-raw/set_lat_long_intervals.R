# Set the lat-log intervals of interest and save
lat_long_intervals <- data.frame(lon_min = -17, lon_max = 52, lat_min = -26, lat_max = 21)
save(lat_long_intervals, file = "data-raw/lat_long_intervals.RData")
