prepInputsVRI <- function(VRIurl, dPath, rasterToMatch){
  VRIin <- prepInputs(url = VRIurl,
                      fun = "sf::st_read",
                      destinationPath = dPath)
  RIA_VRI <- st_transform(VRIin, crs = st_crs(rasterToMatch))
  gcIDRaster <- fasterize::fasterize(RIA_VRI, rasterToMatch, field = "curve2")
  ageRaster <- fasterize::fasterize(RIA_VRI, rasterToMatch, field = "PROJ_AGE_1")
  gcIDRaster[] <- as.integer(gcIDRaster[])
  ageRaster[] <- as.integer(ageRaster[])
  VRIraster <- raster::stack(gcIDRaster, ageRaster)
  return(VRIraster)
}
