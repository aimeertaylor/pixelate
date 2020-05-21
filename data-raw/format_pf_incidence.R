library(raster)
library(dplyr)
rm(list=ls())

load(file = "data-raw/lat_long_intervals.RData")

# Import predictions and lower and upper credible interval bounds
pf_incid_med <- raster::raster('data-raw/Pf_Incidence/Raster Data/Pf_incidence_rate_median/incidence_rate_median_Global_admin0_2017.tif')
pf_incid_LCI <- raster::raster('data-raw/Pf_Incidence/Raster Data/Pf_incidence_rate_LCI/pf_incidence_rate_LCI_Global_admin0_2017.tif')
pf_incid_UCI <- raster::raster('data-raw/Pf_Incidence/Raster Data/Pf_incidence_rate_UCI/pf_incidence_rate_UCI_Global_admin0_2017.tif')

# Convert to data.frame
pf_incid_med <- raster::as.data.frame(pf_incid_med, xy = TRUE)
pf_incid_LCI <- raster::as.data.frame(pf_incid_LCI, xy = TRUE)
pf_incid_UCI <- raster::as.data.frame(pf_incid_UCI, xy = TRUE)

# Take the convention that the outcome variable is denoted z
# (manually check column names first)
colnames(pf_incid_med)[3] <- 'z'
colnames(pf_incid_LCI)[3] <- 'z'
colnames(pf_incid_UCI)[3] <- 'z'


# SubSaharanAfrica_Pf_incidence -------------------------------
lat_long = lat_long_intervals['CentralAfrica',]

# Filter by lat long intervals
med_filtered = dplyr::filter(pf_incid_med,
                                x > lat_long$lon_min,
                                x < lat_long$lon_max,
                                y > lat_long$lat_min,
                                y < lat_long$lat_max)
LCI_filtered = dplyr::filter(pf_incid_LCI,
                             x> lat_long$lon_min,
                             x< lat_long$lon_max,
                             y> lat_long$lat_min,
                             y< lat_long$lat_max)
UCI_filtered = dplyr::filter(pf_incid_UCI,
                             x> lat_long$lon_min,
                             x< lat_long$lon_max,
                             y> lat_long$lat_min,
                             y< lat_long$lat_max)

# Check they all match
writeLines('check all are zeros..')
any((med_filtered$x != LCI_filtered$x) & (med_filtered$x != UCI_filtered$x))
any((med_filtered$y != LCI_filtered$y) & (med_filtered$y != UCI_filtered$y))

# Create and save data frame
CentralAfrica_Pf_incidence = data.frame(x = med_filtered$x,
                                           y = med_filtered$y,
                                           z = med_filtered$z,
                                           u = UCI_filtered$z - LCI_filtered$z)


# CentralAfrica_Pf_incidence -------------------------------
lat_long = lat_long_intervals['SubSaharanAfrica',]

# Filter by lat long intervals
med_filtered = dplyr::filter(pf_incid_med,
                             x > lat_long$lon_min,
                             x < lat_long$lon_max,
                             y > lat_long$lat_min,
                             y < lat_long$lat_max)
LCI_filtered = dplyr::filter(pf_incid_LCI,
                             x> lat_long$lon_min,
                             x< lat_long$lon_max,
                             y> lat_long$lat_min,
                             y< lat_long$lat_max)
UCI_filtered = dplyr::filter(pf_incid_UCI,
                             x> lat_long$lon_min,
                             x< lat_long$lon_max,
                             y> lat_long$lat_min,
                             y< lat_long$lat_max)

# Check they all match
writeLines('check all are zeros..')
any((med_filtered$x != LCI_filtered$x) & (med_filtered$x != UCI_filtered$x))
any((med_filtered$y != LCI_filtered$y) & (med_filtered$y != UCI_filtered$y))

# Create and save data frame
SubSaharanAfrica_Pf_incidence = data.frame(x = med_filtered$x,
                                           y = med_filtered$y,
                                           z = med_filtered$z,
                                           u = UCI_filtered$z - LCI_filtered$z)


# Save example data
usethis::use_data(CentralAfrica_Pf_incidence, SubSaharanAfrica_Pf_incidence, overwrite = TRUE)


