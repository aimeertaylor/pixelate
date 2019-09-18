library(malariaAtlas)

SubSaharanAfrica_shp <- malariaAtlas::getShp(lat = c(-26,21), long = c(-17,52), admin_level = "admin0")
usethis::use_data(SubSaharanAfrica_shp, overwrite = TRUE)

