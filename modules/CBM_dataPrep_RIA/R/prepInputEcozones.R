prepInputsEcozones <- function(url, dPath, masterRaster) {
    ecozones <- prepInputs(
    # this website https://sis.agr.gc.ca/cansis/index.html is hosted by the Canadian Government
    url = url,
    alsoExtract = "similar",
    destinationPath = dPath,
    # rasterToMatch = masterRaster,
    # studyArea = sim$studyArea,
    # useSAcrs = TRUE,
    overwrite = TRUE,
    fun = "sf::st_read",
    filename2 = TRUE
  )
  ecozones <- st_transform(ecozones, st_crs(masterRaster))
  ecozones <- st_crop(ecozones, st_as_sf(as(extent(masterRaster), "SpatialPolygons")))
  ecozones <- ecozones[!ecozones$ZONE_NAME %in% "Pacific Maritime",]
  ecozoneRaster <- fasterize::fasterize(ecozones, masterRaster, field = "ECOZONE")
  ecozoneRaster[is.na(masterRaster[])] <- NA
  ecozoneRaster
}
