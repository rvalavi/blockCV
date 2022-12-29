australia <- rnaturalearth::ne_countries(country = "Australia",
                                         scale = 50,
                                         returnclass = "sf") %>%
  sf::st_crop(c(xmin = 110, xmax = 155, ymin = -45, ymax = -5))


download_path <- "~/Desktop/blockCV_update/"

# get the 19 bio-climatic variables
raster_covars <- geodata::worldclim_global(var = "bio",
                                           res = 5,
                                           path = download_path) %>%
  terra::crop(australia) %>%
  terra::mask(australia)
names(raster_covars) <- substr(names(raster_covars), 10, 100)
plot(raster_covars)

# get the 19 bio-climatic variables
raster_covars10 <- geodata::worldclim_global(var = "bio",
                                           res = 10,
                                           path = download_path) %>%
  terra::crop(australia) %>%
  terra::mask(australia)
names(raster_covars10) <- substr(names(raster_covars10), 11, 100)
plot(raster_covars10)

# GDA 2020 Lambert CRS
output_crs <- 7845
rasters <- terra::project(raster_covars, paste0("epsg:", output_crs))
plot(rasters[[1]], main = names(rasters)[1])
rasters[[1]]
ncell(rasters[[1]])

path <- "~/Desktop/blockCV_update/rasters/"

terra::writeRaster(
  x = rasters,
  filename = paste0(path, names(rasters), ".tif"),
  overwrite = TRUE
)






occ <- spocc::occ(query = "Pseudocheirus peregrinus",
                  from = "gbif",
                  date = c("2000-01-01", "2020-12-31"),
                  gbifopts = list(country = "AU"),
                  has_coords = TRUE,
                  limit = 1e4)

occ_df2 <- spocc::occ2df(occ)
head(occ_df2)
nrow(occ_df2)

occ_sf2 <- sf::st_as_sf(occ_df2, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = output_crs) %>%
  st_coordinates() %>%
  as.data.frame()

a <- occ_sf2[occ_sf2$X < 0, ]
plot(rasters[[1]])
points(a)

occ_df <- read.csv("~/Desktop/blockCV_update/records-2022-12-29.csv") %>%
  dplyr::select(occurrenceStatus, year, decimalLatitude,	decimalLongitude) %>%
  dplyr::filter(occurrenceStatus == "PRESENT",
                year > 2000) %>%
  tidyr::drop_na()
head(occ_df)
nrow(occ_df)

occ_sf <- sf::st_as_sf(occ_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs = output_crs) %>%
  st_coordinates() %>%
  as.data.frame()

b <- occ_sf[occ_sf$Y > -35e5, ]

xy <- rbind(a, b)

rcells <- terra::cellFromXY(rasters[[1]], xy)
a_xy <- xy[!duplicated(rcells), ]
nrow(a_xy)

w <- occ_sf[occ_sf$Y < -35e5, ]
rcells <- terra::cellFromXY(rasters[[1]], w)
w <- w[!duplicated(rcells), ]
w <- w[sample(nrow(w), 500), ]
nrow(w)

final_xy <- rbind(a_xy, w)
nrow(final_xy)

con <- which(final_xy$X < 200000 & final_xy$Y > -35e5)
final_xy[con, ]

xy_cleaned <- final_xy[-con, ]
nrow(xy_cleaned)

xy_filtered <- xy_cleaned[sample(nrow(xy_cleaned), 450), ]

plot(rasters[[1]])
points(xy_cleaned, pch = 16, col='red', cex=0.5)

xy_bg <- terra::spatSample(rasters[[1]], 10000, na.rm=TRUE, xy=TRUE) %>%
  dplyr::select(x, y) %>%
  setNames(c("X", "Y"))
head(xy_bg)

plot(rasters[[1]])
points(xy_bg, pch = 16, col='blue', cex=0.5)
points(xy_cleaned, pch = 16, col='red', cex=0.5)

library(tidyverse)

all_xy <- as.data.frame(xy_cleaned) %>%
  mutate(occ = 1) %>%
  bind_rows(mutate(xy_bg, occ = 0))
head(all_xy)
table(all_xy$occ)

vars <- c("bio_4", "bio_5", "bio_12", "bio_15")
extr <- terra::extract(rasters[[vars]], all_xy[, 1:2], ID=FALSE)
head(extr)
nrow(extr)

extr$occ <- all_xy$occ
training <- extr %>%
  drop_na()
head(training)


# loading the package
library(randomForest)

# convert the response to factor for producing class relative likelihood
training$occ <- as.factor(training$occ)

prNum <- as.numeric(table(training$occ)["1"]) # number of presences
# the sample size in each class; the same as presence number
smpsize <- c("0" = prNum, "1" = prNum)

rf_downsample <- randomForest(formula = occ ~.,
                              data = training,
                              ntree = 1000,
                              sampsize = smpsize,
                              replace = TRUE)

preds <- predict(rasters[[vars]], rf_downsample)
plot(preds)
preds <- predict(rasters[[vars]], rf_downsample, type="prob")[[2]]
plot(preds)

rnd_points <- dismo::randomPoints(raster::raster(preds), 5000, prob = TRUE)
# plot(log(preds))
plot(preds)
points(rnd_points, pch = 16, col='blue', cex=0.5)

absences <- rnd_points %>%
  as.data.frame() %>%
  mutate(probs = terra::extract(preds, ., ID=FALSE)[,1],
         occ = rbinom(nrow(.), 1, prob = probs)) %>%
  dplyr::filter(occ == 0)
nrow(absences)

plot(preds)
points(absences[sample(nrow(absences), 1000), 1:2], pch = 16, col='blue', cex=0.5)
head(absences)

nrow(absences)
nrow(xy_cleaned)

full_occ <- absences %>%
  sample_n(1000) %>%
  dplyr::select(x, y, occ) %>%
  setNames(c("X", "Y", "occ")) %>%
  bind_rows(mutate(xy_cleaned, occ = 1)) %>%
  sample_n(1000)

table(full_occ$occ)

plot(preds)
full_occ %>%
  filter(occ == 0) %>%
  dplyr::select(X, Y) %>%
  points(pch = 16, col='blue', cex=0.5)
full_occ %>%
  filter(occ == 1) %>%
  dplyr::select(X, Y) %>%
  points(pch = 16, col='red', cex=0.5)


# write.csv(full_occ, "inst/extdata/species.csv", row.names = FALSE)

species_data <- read.csv("inst/extdata/species.csv")
head(species_data)

species_sf <- sf::st_as_sf(species_data, coords = c("X", "Y"), crs = 7845)
plot(species_sf)
species_sf

my_blocks <- cv_spatial(x = species_sf,
                        column = "occ3",
                        r = rasters[[2]],
                        k = 5,
                        hexagon = FALSE,
                        # selection = "systematic",
                        iteration = 50,
                        # cell_size = 400000,
                        rows_cols = c(10, 5))

cv_ggplot(my_blocks, x = species_sf, num_plots = 3:4)
cv_ggplot(my_blocks, x = species_sf, r = rasters[[1]])
