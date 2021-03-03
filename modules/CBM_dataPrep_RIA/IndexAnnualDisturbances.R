library(Require)
Require(c("data.table", "raster", "reproducible"))

####get data####
RIA_RTM <- prepInputs(url = 'https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing',
                      destinationPath = 'inputs') #you probably already have this raster - RIA_RTM.tif on google


#Make fire and harvest into RDS objects - Note this product is Txomin Hermosilla's Analysis from
#Landsat Normalized Burn Ratio, different from the White Wulder product (but very similar)
#This does not need to be run by Céline - code is here for reference

# options("reproducible.cachePath" = 'cache')
# harvestRas <- prepInputs(url = 'https://drive.google.com/file/d/1m7mjcx5Sz--RB7x4N3cPYpGkfmxX8KPB/view?usp=sharing',
#                          destinationPath = 'inputs',
#                          fun = 'raster::brick',
#                          useCache = 'overwrite',
#                          filename2 = "CA_Harvest_1985-2015.tif",
#                          # rasterToMatch = RIA_RTM, #they were made with this rtm, but this confirms it
#                          overwrite = TRUE)
# fireRas <- prepInputs(url = 'https://drive.google.com/file/d/1kxCL-i311yd3cS7QDQ2GwHHtyQFiiXoo/view?usp=sharing',
#                       destinationPath = 'inputs',
#                       fun = 'raster::brick',
#                       filename2 = "CA_Fire_1985-2015.tif",
#                       # rasterToMatch = RIA_RTM,
#                       overwrite = TRUE)
#
# makeAnnual <- function(distRaster) {
#   years <- 1985:2015
#   names(distRaster) <- paste0("year", years) #it screws up if numeric only
#   outDT <- lapply(years, FUN = function(year, ras = distRaster){
#     ras <- ras[[paste0("year", year)]]
#     dt <- data.table(pixelID = 1:ncell(ras), vals = getValues(ras), year = year)
#     dt <- dt[!is.na(vals)]
#     set(dt, NULL, 'vals', NULL) #drop vals - only need to know year and pixels
#     return(dt)
#   })
#   outDT <- rbindlist(outDT)
#   return(outDT)
# }
#
# #note Txomin is
# CAharvestAnnual <- makeAnnual(distRaster = harvestRas)
# CAfireAnnual <- makeAnnual(distRaster = fireR)
# saveRDS(CAharvestAnnual, "inputs/CA_harvest_Annual.rds")
# saveRDS(CAfireAnnual, "inputs/CA_fire_Annual.rds")

#####steps that only happen once
#replace pixel values with indices

#tempTHLB because this is not the final THLB
tempTHLB <- prepInputs(url = 'https://drive.google.com/file/d/1ceodWoiKHyK1_fJGDlUMRID3HsneWMJf/view?usp=sharing',
                       destinationPath = 'inputs')

scfmAnnualBurns <- prepInputs(url = 'https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M/view?usp=sharing',
                              destinationPath = 'inputs',
                              overwrite = TRUE,
                              fun = 'readRDS')
#sparse data.table representign Hermosilla et al 30-year change analysis reprojected to RTM
CAharvestAnnual <- prepInputs(url = 'https://drive.google.com/file/d/1Ca-kPun7_VF2xV8s40IJ9hri3Ok-6qx5/view?usp=sharing',
                              destinationPath = 'inputs')
CAfireAnnual <- prepInputs(url = 'https://drive.google.com/file/d/1MjQ5y9Txr1ezv0BatG4b_6FpM6wti1b5/view?usp=sharing',
                           destinationPath = 'inputs')

IndexRTM <- setValues(RIA_RTM, 1:ncell(RIA_RTM))

#postProcess to match tempTHLB
IndexTHLB <- setValues(tempTHLB, 1:ncell(tempTHLB))

#postProcess the RTM
IndexRTM <- postProcess(IndexRTM, rasterToMatch = tempTHLB)
#build matching data.table

indexDT <- data.table(rtmIndex = getValues(IndexRTM),
                      thlbIndex = getValues(IndexTHLB))

#the NAs in scfmIndex are pixels that are not in THLB (but inside the landscape) - we can remove them
indexDT <- indexDT[!is.na(scfmIndex)]


###function for annualSubsets
#if we make this a function
#@param indexDT the data.table with corresponding cell indices, created above
#@param sparseDT a sparse data.table representation of a raster object,
  #with a years column if it represents annual time-series
#@param year at which to subset the data.table
#@param pixelColName the name of the column containing pixelIDs in sparseDT

indexAnnualFire <- function(indexDT, sparseDT, time = NULL, pixelColName = "pixelID") {
  #get only fires from current year
  if (!is.null(time)) {
   sparseDT <- sparseDT[year == time]
  }
  sparseDT <- indexDT[sparseDT, on = c("rtmIndex" = pixelColName)]
  sparseDT <- sparseDT[!is.na(thlbIndex)]
  return(sparseDT$thlbIndex)
}

scfm2018 <- indexAnnualFire(indexDT = indexDT, sparseDT = scfmAnnualBurns, pixelColName = 'pixels', time = 2018)
CAfire2013 <- indexAnnualFire(indexDT = indexDT, sparseDT = CAfireAnnual, pixelColName = "pixelID", time = 2013)
CAharvest2013 <- indexAnnualFire(indexDT = indexDT, sparseDT = CAharvestAnnual, pixelColName = "pixelID", time = 2013)
