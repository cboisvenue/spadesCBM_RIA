# sratch for dataPrep_RIA
library(raster)
library(SpaDES)
library(data.table)

# some raster tasks I am trying to figure out
raslist1 <- "https://drive.google.com/file/d/1O6Laf-y7s-N_WSUtbeTHALRLvffZG-Bp"

rasstack <- raster::stack(raslist1)

masterRaster <- Cache(prepInputs,url = "https://drive.google.com/file/d/1ceodWoiKHyK1_fJGDlUMRID3HsneWMJf")

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
dataPath <- file.path(getwd(),"modules/CBM_dataPrep_RIA/data")



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

#NOT working
scfmFires <- Cache(prepInputs,url = "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8/",
                   fun = "raster::stack", rasterToMatch = masterRaster, useGDAL = FALSE)
#NOT working
scfmURL <- "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8/"
scfmStack <- raster::stack(scfmURL)

#downloaded the file mannually

scfmStack <- raster::stack(file.path(dataPath,"annualFires525yrs.tif"))
names(scfmStack)
scfmStack$annualFires525yrs.1
