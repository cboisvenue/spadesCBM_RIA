

# old rasters from Greg 2015 I think.
rastersOLD <- prepInputs(
  url = "https://drive.google.com/file/d/1DN31xcXh97u6v8NaVcy0O3vzKpLpld69/view?usp=sharing",
  fun = "raster::stack",
  #rasterToMatch = masterRaster, useGDAL = FALSE) # this was Eliot's
  destinationPath = 'inputs')

stackIan <- prepInputs(url = "https://drive.google.com/file/d/1keh-3KIvgy8G5nBBxolH_LsoC53w3qqN/view?usp=sharing",
                       fun = "raster::stack",
                       destinationPath = 'inputs')
# current masterRaster
fireRetMasterRaster <- spadesCBMout$masterRaster

# trying things
rastersOLDt1 <- prepInputs(
  url = "https://drive.google.com/file/d/1DN31xcXh97u6v8NaVcy0O3vzKpLpld69/view?usp=sharing",
  fun = "raster::stack",
  rasterToMatch = spadesCBMout$masterRaster,
  useGDAL = FALSE,
  destinationPath = 'inputs')

age <- postProcess(x = stackIan$RIAlandscape.5,
                   rasterToMatch = spadesCBMout$masterRaster,
                   filename2 = "inputs/pixelAges.tif")

aValues <- age[]
table(aValues, useNA = 'ifany')
mValues <- spadesCBMout$masterRaster[]
table(mValues)
# there are this many more NAs in the age raster provided by Ian
# 5435395 - 5605119
# [1] -169724

# and here they are
agesDT <- data.table(spadesCBMout$allPixDT[,.(ages, pixelIndex)],aValues)
setnames(agesDT,c("ages","aValues"), c("ages2020", "ages2015"))
agesNA <- agesDT[is.na(ages2015) & !is.na(ages2020),]

# need to check if these pixels are disturbed
# are they burned?
presentDayBurns <- prepInputs(url = 'https://drive.google.com/file/d/1MjQ5y9Txr1ezv0BatG4b_6FpM6wti1b5/view?usp=sharing',
                              destinationPath = 'inputs',
                              overwrite = TRUE,
                              fun = 'readRDS')
# there are 3426/169724 that get burned between 1985-2015
agesNA[pixelIndex %in% presentDayBurns$pixelID,]

# are they harvested?
presentDayHarvest <- prepInputs(url = 'https://drive.google.com/file/d/1Ca-kPun7_VF2xV8s40IJ9hri3Ok-6qx5/view?usp=sharing',
                              destinationPath = 'inputs',
                              overwrite = TRUE,
                              fun = 'readRDS')
# there are 14969/169724 that get burned between 1985-2015
agesNA[pixelIndex %in% presentDayHarvest$pixelID,]

(14969 + 3426)/169724
#only 11%


#read-in the VRI2015
# this is an ESRI file
#https://pub.data.gov.bc.ca/datasets/02dba161-fdb7-48ae-a4bb-bd6ef017c36d/2015/VEG_COMP_LYR_L1_POLY_2015.gdb.zip

# back-tracking VRI2015 age raster to 1985 based on the disturbance info
presentDayBurns[, events := 1L]
presentDayHarvest[, events := 2L]
allDist <- rbind(presentDayBurns, presentDayHarvest)

# there needs to be a part here to make the right age raster
RIArtm <- prepInputs(url = "https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing",
                     destinationPath = 'inputs')
VRI2015ageRaster <- Cache(prepInputsVRIage, VRIurl = "https://drive.google.com/file/d/1spcJq_4r4TRA2jfjclQcwBcIfnbniEXZ",
                            dPath = "C:/Celine/github/spadesCBM_RIA/modules/CBM_dataPrep_RIA/data",
                            rasterToMatch = RIArtm,
                            cacheRepo = cacheDir)
# make a data table from the raster


# this is the wrong raster but waiting for the VRI2015
ages2015dt <- agesDT[,.(pixelIndex, ages2015)]




prepInputsVRIage <- function(VRIurl, dPath, rasterToMatch, field = "PROJ_AGE_1"){
  VRIin <- prepInputs(url = VRIurl,
                      fun = "sf::st_read",
                      destinationPath = dPath)
  RIA_VRI <- st_transform(VRIin, crs = st_crs(rasterToMatch))
  #gcIDRaster <- fasterize::fasterize(RIA_VRI, rasterToMatch, field = "curve2")
  ageRaster <- fasterize::fasterize(RIA_VRI, rasterToMatch, field = "PROJ_AGE_1")
  #gcIDRaster[] <- as.integer(gcIDRaster[])
  ageRaster[] <- as.integer(ageRaster[])
  VRIageRaster <- raster::raster(ageRaster)
  return(VRIageRaster)
}

# Eliot modified reproducible so prepInputs can deal with .gbd files
remotes::install_github("PredictiveEcology/reproducible@gdb_archiveNA")
# I am using "dummy" masterRaster --> use your real one
masterRaster <- raster(extent(-1963750, -1321250, 7407500, 8239000),
                       crs = "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m", res = 250)
# Make a studyArea --> use your real one
sa <- as(extent(masterRaster), "SpatialPolygons")
crs(sa) <- crs(masterRaster)
loadAge <- function(x, field = "PROJ_AGE_1") {
  # a <- Cache(sf::st_read, x) # I used Cache during my development because this takes 37 minutes to run -- I was sick of running it again and again
  a <- sf::st_read(x)
  a1 <- a[, field]
  return(a1)
}
a <- prepInputs(url =
                  "https://pub.data.gov.bc.ca/datasets/02dba161-fdb7-48ae-a4bb-bd6ef017c36d/2015/VEG_COMP_LYR_L1_POLY_2015.gdb.zip",
                fun = quote(loadAge(x = targetFilePath,
                                    field = "PROJ_AGE_1")),
                targetFile = "VEG_COMP_LYR_L1_POLY_2015.gdb.zip", archive = NA,
                studyArea = sa
                #rasterTomatch = masterRaster
)
# starting from the read-in age raster in 2015 and keeping the same masterRaster
# as the fireReturnInterval runs

## NOTEs: there are 104427 pixels with NAs in 2015 and ages in 2020. 15018 of
## those are disturbed during the 1985-2015 period in the Wulder and White data.
## There are no pixels that have no ages in 2015 and have ages in 2020. Options:

## figure out what ages the NAs have in 2020.
## Note that some ages in the 2020 that have no ages in the 2015 raster AND that
## are disturbed (so in the sim$disturbanceRasters data table) are old...older
## then the age from the dist...but most are below 50.

ageDTboth <- data.table(ageDT2015, ageDT2020)
ageNoMathc <- ageDTboth[is.na(ageDT2015) & !is.na(ageDT2020),]
library(ggplot2)
distPixWageIn2020 <- qplot(ageNoMathc[pixelIndex %in% allDist$pixelID]$ageDT2020, geom="histogram")

## check the age distribution of the ageNoMathc for the remaining 89409 pixels
## (8046.81 ha over 280118.2)

noAge2015 <- ageNoMathc[!(pixelIndex %in% allDist$pixelID)]
noAge2015hist <- qplot(noAge2015$ageDT2020, geom="histogram")

## have to bring those to 1985...how many are over 35 (2020-1985)?
noAge2015[ageDT2020>=35,]
#81250. These I can do -5.
# look at the others? Can I buffer around them??
probPix <- noAge2015[ageDT2020<35,]
# double checking if any of these are in the disturbed pixels
probPix[pixelIndex %in% allDist$pixelID,]
# none
# trying to see if I can find pixels around those pixels that are disturbed.
# need to create a raster of these pixels and give them a value
probPixRaster <- raster::raster(masterRaster)
probPixRaster[] <- NA
probPixRaster[probPix$pixelIndex] <- 1

## trying to find the closest pixel with a disturbance in the 1985-2015
## disturbance data.table.
# are there pixels that are disturbed twice in this data?
countDist <- allDist[, .N, by = pixelID]
# yes
# are those in the proPix?
probPix[pixelIndex %in% countDist[N>1,]$pixelID,]
# no
# create a raster with all the disturbances
allDistRaster <- raster::raster(masterRaster)
allDistRaster[] <- NA
allDistRaster[countDist$pixelID] <- 1
# trying focal. This takes a raster and gives me back a raster
f3 <- function(x){
  theNAs <- is.na(x)
  if (all(theNAs))
    NA
  else
    x[sample(which(!theNAs),1)]
}
f1 <- focal(allDistRaster,
            w=matrix(1,nrow=5,ncol=5),
            fun = f3)

set(ageDTboth, NULL, "f25", f1[])
checkF25 <- ageDTboth[pixelIndex %in% probPix$pixelIndex,]
checkF25[!is.na(f25),]
#ageDT2015 pixelIndex ageDT2020 pixelIndex f25
# 1:        NA    1073507        15          1   1
# 2:        NA    1078649        15          1   1
# 3:        NA    1081227        15          1   1
# 4:        NA    1083788        15          2   1
# 5:        NA    1083795        15          1   1
# ---
#   6159:        NA    8453863        31          1   1
# 6160:        NA    8456433        31          1   1
# 6161:        NA    8456435        31          1   1
# 6162:        NA    8466726        33          2   1
# 6163:        NA    8469298        33          2   1

## trying a little bigger
f49 <- focal(allDistRaster,
             w=matrix(1,nrow=7,ncol=7),
             fun = f3)

set(ageDTboth, NULL, "f49", f49[])

checkf49 <- ageDTboth[pixelIndex %in% probPix$pixelIndex,]
checkf49[!is.na(f49),]
checkf49[is.na(f49),]

## next to do repeat this but with one disturbance at a time since I have to
## know which disturbance to attribute (1 =fire and 2 = harvest)

## what is up with these last pixels??
lastPix <- checkf49[is.na(f49),]$pixelIndex

sim$allPixDT[pixelIndex %in% lastPix,]
range(sim$allPixDT[pixelIndex %in% lastPix,]$pixelIndex)
lastBad <- NA
allPix <- data.table(sim$allPixDT,lastBad)
allPix$lastBad[which(allPix$pixelIndex %in% lastPix)] <- 1
# Now what? take the values of the ones in the row close by??






