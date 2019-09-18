library(malariaAtlas)

CentralAfrica_shp <- malariaAtlas::getShp(lat = c(-7.5,5), long = c(15,35), admin_level = "admin0")
SubSaharanAfrica_shp <- malariaAtlas::getShp(lat = c(-26,21), long = c(-17,52), admin_level = "admin0")
usethis::use_data(CentralAfrica_shp, SubSaharanAfrica_shp, overwrite = TRUE)

