library(malariaAtlas)

load(file = "data-raw/lat_long_intervals.RData")

lat_lon = lat_long_intervals['SubSaharanAfrica',]
SubSaharanAfrica_shp <- malariaAtlas::getShp(lat = c(lat_lon$lat_min, lat_lon$lat_max),
                                             long = c(lat_lon$lon_min, lat_lon$lon_max),
                                             admin_level = "admin0")

lat_lon = lat_long_intervals['CentralAfrica',]
CentralAfrica_shp <- malariaAtlas::getShp(lat = c(lat_lon$lat_min, lat_lon$lat_max),
                                          long = c(lat_lon$lon_min, lat_lon$lon_max),
                                          admin_level = "admin0")

usethis::use_data(SubSaharanAfrica_shp, CentralAfrica_shp, overwrite = TRUE)

