# sratch for dataPrep_RIA
library(raster)
library(SpaDES)
library(data.table)

dataPath <- file.path(getwd(),"modules/CBM_dataPrep_RIA/data")

### Raster messin' ########################################################################
# some raster tasks I am trying to figure out
## THIS RASTER LIST IS WRONG
raslist1 <- "https://drive.google.com/file/d/1O6Laf-y7s-N_WSUtbeTHALRLvffZG-Bp"

rasstack <- raster::stack(raslist1)

masterRaster <- Cache(prepInputs,url = "https://drive.google.com/file/d/1ceodWoiKHyK1_fJGDlUMRID3HsneWMJf")

# this is a test to see if the pixelAge2015.tif can be read and match the masterRaster
ageRaster2015 <- Cache(prepInputs,url = "https://drive.google.com/file/d/1Kbwdee_JkpCSWwaZUzOqILXUXZUq7luq")
age2015rtm <- postProcess(x = ageRaster2015, rasterToMatch = masterRaster, filename2 = "data/RIA2019/pixelAge2015.tif")

## RIAlandscape.tif from RIA2019
rasters <- Cache(prepInputs,url = "https://drive.google.com/file/d/1O6Laf-y7s-N_WSUtbeTHALRLvffZG-Bp",
                 fun = "raster::stack", rasterToMatch = masterRaster, useGDAL = FALSE)

names(rasters) <- c("TSA","THLB","AU","blockId","age")

quickPlot::dev.useRSGD(FALSE)
dev()
clearPlot()

Plot(rasters$age)
clearPlot()
Plot(rasters$TSA)
clearPlot()
Plot(rasters$AU)

getwd()


writeRaster(rasters$age,file.path(dataPath,"pixelAge2015.tif"))
writeRaster(rasters$TSA,file.path(dataPath,"pixelTSA.tif"))
writeRaster(rasters$AU,file.path(dataPath,"pixelAU.tif"))

# what are the raster ids for the fire (scfm) runs?
## DON'T USE THIS FILE
## Ian produced a raster stack instead since it did not match my masterRaster
# scfmRuns <- readRDS("https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M")
# scfmRuns <- prepInputs(url = "https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M",
#                            fun = "readRDS",
#                            destinationPath = dataPath,
#                            #purge = 7,
#                            filename2 = "scfmAnnualBurn_2015-2540.rds")
#
#
# ## rebuilding the burn raster from the scfmRuns data.table
# # function in 04 - scfmPostHoc.R
# burns2015 <-rstCurrentBurnFromDT(rasterToMatch = masterRaster,scfmAnnualBurns = scfmRuns, 2015)
#
# # checking so I understand
# vals <- values(burns2015)


## Read in scfm fires raster stack

# using brick post Eliot reproducible adjustment
# This should now work as desired
options(reproducible.useGDAL = FALSE)
scfmFires <- Cache(prepInputs,
                   url = "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8",
                   fun = "raster::brick",
                   rasterToMatch = masterRaster,
                   datatype = "INT1U",
                   useGDAL = FALSE)

options(reproducible.useGDAL = FALSE)
scfmFires <- Cache(prepInputs,url = "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8",
                   fun = "raster::brick", rasterToMatch = masterRaster, useGDAL = FALSE)
#NOT working
scfmURL <- "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8/"
scfmStack <- raster::stack(scfmURL)

#downloaded the file manually
scfmStack <- raster::stack(file.path(dataPath,"annualFires525yrs.tif"))
names(scfmStack)
scfmStack$annualFires525yrs.1

fires2015 <- postProcess(x = scfmStack$annualFires525yrs.1, rasterToMatch = masterRaster, filename2 = file.path(dataPath,"scfm2015.tif"))


## From Eliot:
# Pseudo code:
# There will be 2 scales -- the "large" one with the fire data, and
# the small one that is masterRaster (I refer below "small" and "large")
# 1. Create
# a RasterLayer with 1:ncell(ras) as the values in the raster, i..e, 1, 2, 3,
# 4... 8547820.
# 2.Run prepInputs with rasterToMatch on that layer. It will, by
# default, use "nearest neighbour", which is what you want. You will get a new,
# smaller, raster indexMasterRaster with values like c(3, 4, 5, 8, 9, 10...)
# which will be the index on the bigger raster, its position on the smaller
# raster.
# 3.For each year you need disturbance: take the data.table of fire
# locations (likely columns: "year", "pixelId"... make sure it is
# as.integer(pixelId) to keep memory low), indexMasterRaster , and an empty
# "large" raster and empty small raster;
# emptyLarge[data.tableFiresYearXPixelId]
# <- 1 emptySmall[] <- emptyLarge[indexMasterRaster[]]

# from Ian:I've been pushing these scripts to WBI_RIA, FYI, but here it is.
library(Require)
Require(c("data.table", "raster", "reproducible"))

####get data####
# this is the RTM that was used to create the data.table
RIA_RTM <- prepInputs(url = 'https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing',
                      destinationPath = 'inputs') #you probably already have this raster - RIA_RTM.tif on google
# this is the RTM_THBL, the raster the spadesCBM sims and ws3 sims will be
# completed on
tempTHLB <- prepInputs(url = "https://drive.google.com/file/d/1LlE5NXDPKS6Ljv2SvwQnw2A_Kze2uun-",
                       destinationPath = 'inputs')
# this is the data.table
scfmAnnualBurns <- prepInputs(url = 'https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M',
                              destinationPath = 'inputs',
                              overwrite = TRUE,
                              fun = 'readRDS')

#####steps that only happen once
#replace pixel values with indices
IndexRTM <- setValues(RIA_RTM, 1:ncell(RIA_RTM))

#postProcess to match tempTHLB
IndexTHLB <- setValues(tempTHLB, 1:ncell(tempTHLB))

#postProcess the RTM
IndexRTM <- postProcess(IndexRTM, rasterToMatch = tempTHLB)
#build matching data.table

indexDT <- data.table(scfmIndex = getValues(IndexRTM),
                      thlbIndex = getValues(IndexTHLB))

#the NAs in scfmIndex are pixels that are not in THLB (but inside the landscape) - we can remove them
indexDT <- indexDT[!is.na(scfmIndex)]
   # In this file, scfmIndex is the pixel index on the scfm map,
   # THLB is the pixel index on the THLB map that matches that location

###function for annualSubsets
indexAnnualFire <- function(indexDT, annualBurns, time) {
  #get only fires from current year
  annualBurns <- annualBurns[year == time]
  annualBurns <- indexDT[annualBurns, on = c("scfmIndex" = "pixels")]
  annualBurns <- annualBurns[!is.na(thlbIndex)]
  return(annualBurns$thlbIndex)

}
##it returns a vector of pixelIDs, all of which burned

thisYearBurned <- indexAnnualFire(indexDT = indexDT, annualBurns = scfmAnnualBurns, time = 2016)

someRaster[thisYearBurned] <- 0 #if this was age and you wanted burn ages to be 0
someDT[pixelID %in% thisYearBurned, age := 0] # if you were using a data.table instead

# to run annually, youll just have to make sure those objects indexDT and
# scfmAnnualBurns are in the simList somewhere, and change 2016 to time(sim)

indexAnnualFire(indexDT = indexDT, annualBurns = scfmAnnualBurns, time = 2016)


#######################
## Wulder and WHite dist - present day runs
options(reproducible.useGDAL = FALSE)
## these rasters were removed? No, looks like they were zipped...histFires2
histFires <- Cache(prepInputs,
                   url = "https://drive.google.com/file/d/1cM7Ube-CBncqIol1KMWAByh0Oa-k1HXV",
                   fun = "raster::stack",
                   rasterToMatch = masterRaster,
                   datatype = "INT1U",
                   useGDAL = FALSE)

histFires2 <- prepInputs(
                  url = "https://drive.google.com/file/d/1kxCL-i311yd3cS7QDQ2GwHHtyQFiiXoo/view?usp=sharing",
                  destinationPath = 'inputs')

### TRYING ALL THREE DTs Ian provided############################

# see the IndexAnnualDisturbances.R script
# library(Require)
# Require(c("data.table", "raster", "reproducible"))

####get data####

RIA_RTM <- prepInputs(url = 'https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing',
                      destinationPath = 'inputs') #you probably already have this raster - RIA_RTM.tif on google


### Trying to create my own masterRaster #################################
# files obtained from greg's drive on March 8th 2021
# using only the ria_vri_final_shp

## Mine
# This is the RIA_RTM you gave me Ian
RIA_RTM <- prepInputs(url = 'https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing',
                      destinationPath = 'inputs') #you probably already have this raster - RIA_RTM.tif on google

masterShape <- prepInputs(
  url = "https://drive.google.com/file/d/1j8P780gwsSiEekANMm7qIphoVMVm74BJ",
  destinationPath = dataPath,
  rasterToMatch = RIA_RTM,
  overwrite = TRUE,
  fun = "raster::shapefile",
  filename2 = TRUE)

masterRaster <- fasterize::fasterize(masterShape, raster = RIA_RTM, field = "thlb2")

## Ian's with comments
###1. We need a template raster so when we rasterize the VRI shapefile, it will have a resolution, CRS, etc. So, we will use the RIA_RTM -
#this is a bare-bones prepInputs call that will get the RIA_RTM
RIA_RTM <- prepInputs(url = 'https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing',
                      destinationPath = 'inputs')
####2. Then we get the VRI shapefile - the VRI is an enormous shapefile, reading in shapefiles like that will take along time unless we use sf::st_read
#we also want to crop and reproject it as soon as possible, to speed things up. So we pass our RIA_RIA as the rasterToMatch in prepInputs
cacheDir <- reproducible::checkPath("cache", create = TRUE)
masterShape <- prepInputs(
  url = "https://drive.google.com/file/d/1j8P780gwsSiEekANMm7qIphoVMVm74BJ",
  destinationPath = dataPath,
  rasterToMatch = RIA_RTM,
  useCache = TRUE,
  overwrite = TRUE,
  fun = "sf::st_read", #sf::st_read will load shapefiles into R as an 'sf' objects  - this is important for speed!
  filename2 = TRUE)
#SF and SPDF objects are both GIS 'vectors' datasets but rely on a different set of packages/functions for GIS operations in R.
# fortunately prepInputs is aware of this and should pick correctly.
####3. now we will make your THLB raster, the true rasterToMatch - by rasterizing masterShape.
# the raster will have the spatial characteristics (CRS, res, origin, etc) of RIA_RTM, but values determined by the field 'thlb2'.
#I believe it will be NA where thlb2 is NA, and where RIA_RTM is NA, but I'm not 100% about the latter
masterRaster <- fasterize::fasterize(masterShape, raster = RIA_RTM, field = "thlb2")
#fasterize requires an sf object for the shapefile (masterShape), which is another reason we pass sf::st_read instead of raster::shapefile


masterShape <- prepInputs(
  url = "https://drive.google.com/file/d/1j8P780gwsSiEekANMm7qIphoVMVm74BJ",
  destinationPath = dataPath,
  targetFile = "ria_vri-final.shp",
  alsoExtract = c("ria_vri-final.shx", "ria_vri-final.prj", "ria_vri-final.dbf"),
  rasterToMatch = RIA_RTM,
  useCache = TRUE,
  purge = 7,
  overwrite = TRUE,
  fun = "sf::st_read", #sf::st_read will load shapefiles into R as an 'sf' objects  - this is important for speed!
  filename2 = TRUE)
### CONCLUSION: can't, Ian will build the masterRaster

## age raster
ageRaster2020 <- prepInputs(
  url = "https://drive.google.com/file/d/16pXx9pufiWY45H8rgCiFl6aFbEcC9GeD",
  fun = "raster::raster",
  destinationPath = dataPath
)
## au rasrter
gcIdRaster <- prepInputs(
  url = "https://drive.google.com/file/d/1LKblpqSZTVlaNske-C_LzNFEtubn7Ms7",
  fun = "raster::raster",
  purge = 7,
  destinationPath = dataPath
)

## age raster? ### CONCLUSION: this .dbf is 16GB, my laptop can't handle it.
# check out what is in the .dbf
library(foreign)
gregDbf <- as.data.table(read.dbf(file.path(getwd(),"data/ria_vri-final.dbf")))

gregDbf <- as.data.table(read.dbf(
  "C:/Celine/github/spadesCBM_RIA/modules/CBM_dataPrep_RIA/data/GregParadis/data/spadescbm_bundle/ria_tmp-shp/ria_tmp.dbf"))

### END OF raster messin' ########################################################################

## reading in userDist #############################

userDist <- prepInputs(url = "https://drive.google.com/file/d/1Gr_oIfxR11G1ahynZ5LhjVekOIr2uH8X",
                       fun = "data.table::fread",
                       destinationPath = dataPath,
                       #purge = 7,
                       filename2 = "mySpuDmids.csv")
## End user dist####################################


### growth and yield###################################################

## read-in the three files
au_table <- prepInputs(url = "https://drive.google.com/file/d/1YmQ6sNucpEmF8gYkRMocPoeKt2P26ZiX/view?usp=sharing",
                       fun = "data.table::fread",
                       destinationPath = dataPath,
                       overwrite = TRUE,
                       purge = 7,
                       filename2 = "au_table.csv")

curve_points_table <- prepInputs(url = "https://drive.google.com/file/d/1BYHhuuhSGIILV1gmoo9sNjAfMaxs7qAj",
                       fun = "data.table::fread",
                       destinationPath = dataPath,
                       overwrite = TRUE,
                       #purge = 7,
                       filename2 = "curve_points_table.csv")

curve_table <- prepInputs(url = "https://drive.google.com/file/d/1SJJCdrmVmAnCLeEahfdAY5lgTws0oilJ",
                       fun = "data.table::fread",
                       destinationPath = dataPath,
                       overwrite = TRUE,
                       #purge = 7,
                       filename2 = "curve_table.csv")
library(ggplot2)
volCurves <- ggplot(data = curve_points_table, aes(x = x, y = y,
                                                   group = curve_id,
                                                   colour = curve_id)) +
                    geom_line() ## TODO: move to plotInit event

quickPlot::dev.useRSGD(FALSE)
dev()
clearPlot()

curve_points_table[, V1 := NULL]
names(curve_points_table) <- c("GrowthCurveComponentID", "Age", "MerchVolume")

#checks on curve_point_table###################################

## read-in an SK run to check on things
SKsim <- readRDS(file = "C:/Celine/github/temp/CBMsimList1993.rds")


# curves that have a min age >10
  ten <- which(curve_points_table[,.(min = min(Age)), by = .(GrowthCurveComponentID)]$min > 10)
  moreThan10startCurves <- curve_points_table[,.(min = min(Age)), by = .(GrowthCurveComponentID)][ten,]
  ## there are 40 curves with min age >10...one STARTS at 96...
  names(moreThan10startCurves)
  #[1] "curve_id" "min"
  setnames(moreThan10startCurves,"min","Age")
  names(moreThan10startCurves)
  cols <- c("GrowthCurveComponentID", "Age")
  checkVolOnLateStart10 <- merge(moreThan10startCurves,curve_points_table, by = cols)
  # all MerchVol == 0 at these min Ages

  checkVolOnLateStart10b <- curve_points_table[GrowthCurveComponentID %in% moreThan10startCurves$GrowthCurveComponentID,]
  # they all seem to have 0 vol at the start age. This one min age is 14 but all 0s until 23.
  curve_points_table[curve_id == 801003,][y > 0,]
  ## still 9 years before there is volume there...

  # what age does the vol start?
  lateStartNon0 <- checkVolOnLateStart10b[MerchVolume > 0,]
  range(lateStartNon0[, .(minAgeWVol = min(Age)), by = .(GrowthCurveComponentID)]$minAgeWVol)
  #[1] 19 96


# how many have age 0? and how many of those have vol at age 0? any?
  curveWithAge0 <- curve_points_table[ x == 0,]
  # 80 - they all also have y == 0
  ## This means I can get rid of the 0s ##
  curve_points_table[ x == 0 & y > 0,]
  # none ## GOOD

  curveNo0 <- curve_points_table[Age > 0 & MerchVolume > 0, ]

  #Plot these
  volCurvesNo0 <- ggplot(data = curveNo0, aes(x = Age, y = MerchVolume,
                                                     group = GrowthCurveComponentID,
                                                     colour = GrowthCurveComponentID)) +
    geom_line()


  # in the curves with age 0, at what age does vol start?
  curveWithAge0a <- curve_points_table[curve_id %in% curveWithAge0$curve_id,]
  # nope can't get this to work...
  #setDT(curveWithAge0a)[,firstVol := first(y), curve_id]
  curveWithAge0a[y > 0,][,.(min = min(y)), by = .(curve_id)]
  ## @Greg ### All start small except for curve 823006 7.5

  ## some curve don't have vole until they are 50
  firstVolAllcurves <- curve_points_table[y > 0,]
  ## min age at firstVol
  range(firstVolAllcurves[,.(minAge = min(x)), by = .(curve_id)]$minAge)
  firstVolAllcurves[,.(minAge = min(x)), by = .(curve_id)][which(firstVolAllcurves[,.(minAge = min(x)), by = .(curve_id)]$minAge > 20),]
  ## 103 of the curves have their first volume at age > 20


  ### TYPSY VS VDYP
  ## TYPSY
  ## are there differences between the VDYP and the TYPSY curves?
  cols <- c("au_id", "managed_curve_id")
  ## the VDPY curve ids are repeated in the typsy curve id, it was to avoid NAs
  typsy <- au_table[which(managed_curve_id != unmanaged_curve_id),]
  typsy[,.(au_id,managed_curve_id)]

  typsyVol <- curve_points_table[GrowthCurveComponentID %in% typsy$managed_curve_id,]

  # Plot typsy curves only
  volTypsyPlots <- ggplot(data = typsyVol, aes(x = Age, y = MerchVolume,
                                              group = GrowthCurveComponentID,
                                              colour = GrowthCurveComponentID)) +
    geom_line()

  clearPlot()
  volTypsyPlots

  typsyVol[Age == 0 & MerchVolume != 0,]
  typsyCheck1 <- typsyVol[MerchVolume > 0, .(minAgeWvol = min(Age)), by = .(GrowthCurveComponentID)]
  range(typsyCheck1$minAgeWvol)
  #[1] 18 70
  ## OBSERVATION: smallest age the typsy curves have vol is 18 and some is up to 70
  ## what vols at those ages?
  tempTypsy <- typsyVol[typsyCheck1, on = .(GrowthCurveComponentID)]
  tempTypsy[Age == minAgeWvol,]
  ## that looks good

  ## VDYP
  vdypVol <- curve_points_table[GrowthCurveComponentID %in% au_table$unmanaged_curve_id,]
  volVdyplots <- ggplot(data = vdypVol, aes(x = Age, y = MerchVolume,
                                               group = GrowthCurveComponentID,
                                               colour = GrowthCurveComponentID)) +
    geom_line()

  clearPlot()
  volVdyplots

  vdypVol[Age == 0 & MerchVolume != 0,] # none: good!

  vdypCheck1 <- vdypVol[MerchVolume > 0, .(minAgeWvol = min(Age)), by = .(GrowthCurveComponentID)]
  range(vdypCheck1$minAgeWvol)
  # 1 96
  tempVdyp <- typsyVol[typsyCheck1, on = .(GrowthCurveComponentID)]
  tempVdyp[Age == minAgeWvol,]

  ## VDYP curves


## read-in the portion of the VRI table
vri_sample <- prepInputs(url = "https://drive.google.com/file/d/1GPnavp64rN36-5L2xdVq3r7ooxKF1XGL",
                        fun = "data.table::fread",
                        destinationPath = dataPath,
                        #purge = 7,
                        filename2 = "samplVRI.csv")

vriNames <-  prepInputs(url = "https://drive.google.com/file/d/1Pkh8Xo5o_aAufGCbukxeaWuNg_5U6IUQ",
                        fun = "data.table::fread",
                        destinationPath = dataPath,
                        #purge = 7,
                        filename2 = "samplVRI.csv")
# there are biomass calculations for each of the stands...chance for data assimilation??
grep("BIO",vriNames$V2)
#[1] 185 187

auCurve1Curve2 <- prepInputs(url = "https://drive.google.com/file/d/103z0_z3ACjEf70d-x8m8Zv1uxbZk5RXd",
                             fun = "data.table::fread",
                             destinationPath = dataPath,
                             #purge = 7,
                             filename2 = "samplVRI.csv")

