# sratch for dataPrep_RIA
library(raster)
library(SpaDES)
library(data.table)

dataPath <- file.path(getwd(),"modules/CBM_dataPrep_RIA/data")

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
tempTHLB <- prepInputs(url = 'https://drive.google.com/file/d/1ceodWoiKHyK1_fJGDlUMRID3HsneWMJf/view?usp=sharing',
                       destinationPath = 'inputs')
# this is the data.table
scfmAnnualBurns <- prepInputs(url = 'https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M/view?usp=sharing',
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

### TRYING ALL THREE DTs Ian provided


library(Require)
Require(c("data.table", "raster", "reproducible"))
####get data####
RIA_RTM <- prepInputs(url = 'https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing',
                      destinationPath = 'inputs') #you probably already have this raster - RIA_RTM.tif on google


### growth and yield

## read-in the three files
au_table <- prepInputs(url = "https://drive.google.com/file/d/1gnFkqDyIvF9LN8f-K_9i95_MaWjXhgfc",
                       fun = "data.table::fread",
                       destinationPath = dataPath,
                       #purge = 7,
                       filename2 = "au_table.csv")

curve_points_table <- prepInputs(url = "https://drive.google.com/file/d/1MIH4roV7lQdGFFmdLCN8bNZ0eNoakmlL",
                       fun = "data.table::fread",
                       destinationPath = dataPath,
                       #purge = 7,
                       filename2 = "curve_points_table.csv")

curve_table <- prepInputs(url = "https://drive.google.com/file/d/18wOyyf4QCw6XoSH2yAGDG4jhtXST074l",
                       fun = "data.table::fread",
                       destinationPath = dataPath,
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

## reading in userDist

userDist <- prepInputs(url = "https://drive.google.com/file/d/1Gr_oIfxR11G1ahynZ5LhjVekOIr2uH8X",
                          fun = "data.table::fread",
                          destinationPath = dataPath,
                          #purge = 7,
                          filename2 = "mySpuDmids.csv")

#checks on curve_point_table###################################

# curves that have a min age >10
  ten <- which(curve_points_table[,.(min = min(x)), by = .(curve_id)]$min > 10)
  moreThan10startCurves <- curve_points_table[,.(min = min(x)), by = .(curve_id)][ten,]
  names(moreThan10startCurves)
  #[1] "curve_id" "min"
  setnames(moreThan10startCurves,"min","x")
  names(moreThan10startCurves)
  cols <- c("curve_id","x")
  checkVolOnLateStart10 <- merge(moreThan10startCurves,curve_points_table, by = cols)
  # they all seem to have 0s at the start age. This one min age is 14 but all 0s until 25.
  curve_points_table[curve_id == 801003,][y > 0,]
  ## still 11 years before there is volume there...

# how many have age 0? and how many of those have vol at age 0? any?
  curveWithAge0 <- curve_points_table[ x == 0,]
  # 68
  curve_points_table[ x == 0 & y > 0,]
  # none
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

