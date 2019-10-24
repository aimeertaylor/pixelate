#################################################
# All code in this script is copied and pasted
# from the vignette of Vizumap [1], which is the
# R-package associated with [2].
#
# [1] Petra Kuhnert and Lydia Lucchesi (2017). Vizumap:
#     Visualizing Uncertainty in spatial data. R package version 1.1.0.
# [2] Lucchesi, Lydia R., and Christopher K. Wikle.
#     "Visualizing uncertainty in areal data with bivariate choropleth maps,
#      map pixelation and glyph rotation." Stat 6.1 (2017): 292-302.
#################################################

# install.packages("devtools")
if(!requireNamespace("Vizumap")){
  devtools::install_github(repo = "pkuhnert/Vizumap", build_vignettes = TRUE, force = TRUE)
}

library('Vizumap')
vignette("Vizumap")

# The code
data(us_data)
us_data$GEO.id2 <- as.numeric(us_data$GEO.id2)
ca_data <- subset(us_data, us_data$GEO.id2 > 6000 & us_data$GEO.id2 < 7000)
ca_data <- read.uv(data = ca_data, estimate = "pov_rate", error = "pov_moe")
row.names(ca_data) <- seq(1, nrow(ca_data), 1)

data(us_geo)
ca_geo <- subset(us_geo, us_geo@data$STATE == "06")
pix <- pixelate(ca_geo, id = "region")

df <- data.frame(region = sapply(slot(ca_geo, "polygons"), function(x) slot(x, "ID")), name = unique(ca_geo@data$GEO_ID))
ca_data$region <- df[match(ca_data$GEO_ID, df$name), 1]
ca_data$region <- as.character(ca_data$region)

#check that values in shared column match
ca_data$region %in% pix$region

#uniform distribution
u_m <- build_pmap(data = ca_data, distribution = "uniform", pixelGeo = pix, id = "region")

#normal distribution
ca_data$se <- ca_data$pov_moe / 1.645
ca_data <- read.uv(data = ca_data, estimate = "pov_rate", error = "se")

n_m <- build_pmap(data = ca_data, distribution = "normal", pixelGeo = pix, id = "region")

#experiment with discrete distribution
#exponential - example for q argument
ca_data.q <- with(ca_data, data.frame(p0.05 = qexp(0.05, 1/pov_rate), p0.25 = qexp(0.25, 1/pov_rate), p0.5 = qexp(0.5, 1/pov_rate), p0.75 = qexp(0.75, 1/pov_rate), p0.95 = qexp(0.95, 1/pov_rate)))

head(ca_data.q)

d_m <- build_pmap(data = ca_data, distribution = "discrete",
                  pixelGeo = pix, id = "region", q = ca_data.q)


#animate the normal distribution map
a <- animate(n_m, aniLength = 30)

# Plot pixelated and annimated
view(u_m)
view(n_m)
view(d_m)
view(a)
view(d_m)
