# This code is adapted from the following link:
# https://www.r-bloggers.com/how-to-quickly-enrich-a-map-with-natural-and-anthropic-details/

library(sf)
rm(list = ls())

# Load the limits of the spatial predictions
load(file = "data-raw/lat_long_intervals.RData")

# Set a Coordinate Reference System using an integer EPSG code.
# 4326 is the EPSG code for the World Geodetic System 1984 (WGS 84)
EPSG_code = 4326

# Make an enclosed spatial polygon and
plygn <- matrix(c(lat_long_intervals$lon_min,
                  lat_long_intervals$lon_min,
                  lat_long_intervals$lon_max,
                  lat_long_intervals$lon_max,
                  lat_long_intervals$lon_min,
                  lat_long_intervals$lat_max,
                  lat_long_intervals$lat_min,
                  lat_long_intervals$lat_min,
                  lat_long_intervals$lat_max,
                  lat_long_intervals$lat_max), ncol = 2)

bbox.sf <- sf::st_geometry(st_polygon(x = list(plygn)))
sf::st_crs(bbox.sf) = EPSG_code

# Save example bbox
usethis::use_data(bbox.sf, overwrite = TRUE)
