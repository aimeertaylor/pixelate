# Create an example simple feature (sf) bounding box
# for plots with country borders cropped to extent of
# spatial predictaions.

# This code is adapted from the following link:
# https://www.r-bloggers.com/how-to-quickly-enrich-a-map-with-natural-and-anthropic-details/
# By Francesco Bailo, August 9, 2018

library(sf) # for simple feature spatial objects
rm(list = ls())

# Load the limits of the spatial predictions
load(file = "data-raw/lat_long_intervals.RData")

# Set a Coordinate Reference System (CRS) using an integer EPSG code.
# 4326 is the EPSG code for the World Geodetic System 1984 (WGS 84) in degrees
# (EPSG:32632 is the EPSQ CRS of WGS 84 in meters)
EPSG_code = 4326

# Matrix of coordinates defining an enclosed bounding box
bbox <- matrix(c(lat_long_intervals$lon_min,lat_long_intervals$lon_min,
                 lat_long_intervals$lon_max,lat_long_intervals$lon_max,
                 lat_long_intervals$lon_min,lat_long_intervals$lat_max,
                 lat_long_intervals$lat_min,lat_long_intervals$lat_min,
                 lat_long_intervals$lat_max,lat_long_intervals$lat_max), ncol = 2)

# Create a sf polygon with sf::st_polgon and get its geometry with sf::st_geometry
bbox.sf <- sf::st_geometry(sf::st_polygon(x = list(bbox)))

# Set the crs using sf::st_crs
sf::st_crs(bbox.sf) = EPSG_code

# Save the example simple feature bounding box
usethis::use_data(bbox.sf, overwrite = TRUE)
