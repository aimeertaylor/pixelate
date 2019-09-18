library(raster) # Suggests
library(dplyr) # Imports
# library(rgdal) # James, is this really needed?
# library(malariaAtlas) # James, is this really needed?
rm(list=ls())

# Import point estimates and lower and upper credible interval bounds
pf_incid_median = raster::raster('incidence_rate_median_Global_admin0_2017.tif')
pf_incid_LCI = raster::raster('pf_incidence_rate_LCI_Global_admin0_2017.tif')
pf_incid_UCI = raster::raster('pf_incidence_rate_UCI_Global_admin0_2017.tif')

# Convert to data.frame
pf_incid_median <- as.data.frame(pf_incid_median, xy = TRUE)
pf_incid_LCI <- as.data.frame(pf_incid_LCI, xy = TRUE)
pf_incid_UCI <- as.data.frame(pf_incid_UCI, xy = TRUE)

# Take the convention that the outcome variable is denoted z
colnames(pf_incid_median)[3] = 'z'
colnames(pf_incid_LCI)[3] = 'z'
colnames(pf_incid_UCI)[3] = 'z'


# Central Africa ---------------------------------------

# The lat-log intervals of interest
lat_long_intervals = c(15,35,-7.5,5)

# Filter by lat long intervals
pf_incid_median = dplyr::filter(pf_incid_median,
                              x > lat_long_intervals[1],
                              x < lat_long_intervals[2],
                              y > lat_long_intervals[3],
                              y < lat_long_intervals[4])
pf_incid_LCI = dplyr::filter(pf_incid_LCI,
                                 x> lat_long_intervals[1],
                                 x< lat_long_intervals[2],
                                 y> lat_long_intervals[3],
                                 y< lat_long_intervals[4])
pf_incid_UCI = dplyr::filter(pf_incid_UCI,
                      x> lat_long_intervals[1],
                      x< lat_long_intervals[2],
                      y> lat_long_intervals[3],
                      y< lat_long_intervals[4])

# Check they all match
writeLines('check all are zeros..')
sum((pf_incid_median$x != pf_incid_LCI$x) & (pf_incid_median$x != pf_incid_UCI$x))
sum((pf_incid_median$y != pf_incid_LCI$y) & (pf_incid_median$y != pf_incid_UCI$y))

# Create and save data frame
CentralAfrica_Pf_incidence = data.frame(x=pf_incid_median$x,
                        y = pf_incid_median$y,
                        Median = pf_incid_median$z,
                        UCI = pf_incid_UCI$z,
                        LCI = pf_incid_LCI$z)

save(CentralAfrica_Pf_incidence, file = '~/Documents/Pixel/RData/CentralAfrica_Pf_incidence.RData')

# Sub-Saharan Africa ----------------------------------

# The lat-log intervals of interest
lat_long_intervals = c(-17,52,-26,21)

# Filter by lat long intervals
pf_incid_median = dplyr::filter(pf_incid_median,
                         x > lat_long_intervals[1],
                         x < lat_long_intervals[2],
                         y > lat_long_intervals[3],
                         y < lat_long_intervals[4])
pf_incid_LCI = dplyr::filter(pf_incid_LCI,
                      x> lat_long_intervals[1],
                      x< lat_long_intervals[2],
                      y> lat_long_intervals[3],
                      y< lat_long_intervals[4])
pf_incid_UCI = dplyr::filter(pf_incid_UCI,
                      x> lat_long_intervals[1],
                      x< lat_long_intervals[2],
                      y> lat_long_intervals[3],
                      y< lat_long_intervals[4])

# Check they all match
writeLines('check all are zeros..')
sum((pf_incid_median$x != pf_incid_LCI$x) & (pf_incid_median$x != pf_incid_UCI$x))
sum((pf_incid_median$y != pf_incid_LCI$y) & (pf_incid_median$y != pf_incid_UCI$y))

# Create and save data frame
SubSaharanAfrica_Pf_incidence = data.frame(x=pf_incid_median$x,
                        y = pf_incid_median$y,
                        Median = pf_incid_median$z,
                        UCI = pf_incid_UCI$z,
                        LCI = pf_incid_LCI$z)

save(SubSaharanAfrica_Pf_incidence, file = '../RData/SubSaharanAfrica_Pf_incidence.RData')
