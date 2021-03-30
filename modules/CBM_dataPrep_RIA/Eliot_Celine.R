
Require::Require(c("sf", "data.table", "reproducible"))
# Get a template raster for StudyArea -- res 250, extent = ..
# b <- prepInputs(url = "https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing")

# 2.5G file -- too big for R's unzip -- this will get it, but not unzip on R
# prepInputsVRI <- function(VRIurl, dPath, crsWanted){
#   VRIin <- prepInputs(url = VRIurl,
#                       fun = "sf::st_read",
#                       destinationPath = dPath)
#   RIA_VRI <- st_transform(VRIin, crs = crsWanted)
#   gcIDRaster <- fasterize::fasterize(RIA_VRI, b, field = "curve2")
#   ageRaster <- fasterize::fasterize(RIA_VRI, b, field = "PROJ_AGE_1")
#   gcIDRaster[] <- as.integer(gcIDRaster[])
#   ageRaster[] <- as.integer(ageRaster[])
#   VRIraster <- raster::stack(gcIDRaster, ageRaster)
#   return(VRIraster)
# }

RIA_VRIstack <- Cache(prepInputsVRI,VRIurl = "https://drive.google.com/file/d/1LXSX8M46EnsTCM3wGhkiMgqWcqTubC12",
                          dPath = "C:/Celine/github/spadesCBM_RIA/modules/CBM_dataPrep_RIA/data",
                          crsWanted = st_crs(b),
                 )

 RIA_VRI <- st_transform(origRIA_VRI, crs = st_crs(b))

VRI3cols <- prepInputs(url = "https://drive.google.com/file/d/1LXSX8M46EnsTCM3wGhkiMgqWcqTubC12/view?usp=sharing",
                fun = "sf::st_read",
                destinationPath = "C:/Celine/github/spadesCBM_RIA/modules/CBM_dataPrep_RIA/data")



# Now again, from local file named **.shp
origRIA_VRI <- prepInputs(targetFile = "ria_vri-final.shp",
                          destinationPath = "data", fun = "sf::st_read")

# Convert it to same crs as b
RIA_VRI <- st_transform(origRIA_VRI, crs = st_crs(b))

gcIDRaster <- fasterize::fasterize(RIA_VRI, b, field = "curve2")
ageRaster <- fasterize::fasterize(RIA_VRI, b, field = "PROJ_AGE_1")
ageWData <- !is.na(ageRaster[])
gcIDWData <- !is.na(gcIDRaster[])

masterRaster <- raster::raster(gcIDRaster)
masterRaster[ageWData] <- gcIDRaster[ageWData]
gcIDRaster[!ageWData] <- NA


dPath <- "data"


prepInputsEcozones <- function(url, dPath, masterRaster) {
  ecozones <- prepInputs(
    # this website https://sis.agr.gc.ca/cansis/index.html is hosted by the Canadian Government
    url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
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

ecozone <- Cache(prepInputsEcozones, url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
                 dPath = dPath,
                 masterRaster = masterRaster)


cbmAdmin <- prepInputs(url= "https://drive.google.com/file/d/1NwVPqMho8fOcf9mVVzm5SlfDbbyQmQlh/view?usp=sharing",
                       fun = "data.table::fread")
provs <- getData("GADM", country = "CAN", level = 1)

bc <- provs[provs$NAME_1 == "British Columbia",]

cbmAdminThisSA <- cbmAdmin[adminName == "British Columbia", ]

rows <- match(ecozone[], cbmAdminThisSA$EcoBoundaryID)
spatialUnitID <- cbmAdminThisSA[rows,"SpatialUnitID"]

# Assertion
table(spatialUnitID)


dt <- data.table(gcID = gcIDRaster[], age = ageRaster[], ecozone = ecozone[], spatialUnitID = spatialUnitID)

# assertion -- if there are both NAs or both have data, then the colums with be the same, so sum is either 0 or 2
bbb <- apply(dt, 1, function(x) sum(is.na(x)))
if (!all(names(table(bbb)) %in% c("0", "4")))
  stop("should be only 0 or 4s")

VRI_3Cols <- origRIA_VRI[, c("curve1", "curve2", "PROJ_AGE_1")]
st_write(VRI_3Cols, "VRI_3Cols.shp")
files <- dir(pattern = "VRI_3Cols")
zip(zipfile = "VRI_3Cols.zip", files = files)


object.size(VRI_3Cols)
object.size(origRIA_VRI)
