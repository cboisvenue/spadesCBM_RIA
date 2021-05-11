## Figure for RIA runs

# the CBMutils::spatialPlot function plot rasters directly, but does not return
# rasters. Here I modify it to return rasters for more flexibility in presenting
# results


spatialRaster <- function(pixelkeep, cbmPools, poolsToPlot, years, masterRaster) {
  cbmPools[is.na(cbmPools)] <- 0
  colnames(cbmPools)[c(1,3,4)] <- c("simYear", "pixelGroup", "age")
  if ("totalCarbon" %in% poolsToPlot) {
    totalCarbon <- apply(cbmPools[, SoftwoodMerch:HardwoodBranchSnag], 1, "sum")
    cbmPools <- cbind(cbmPools, totalCarbon)
  }

  if (any(!poolsToPlot %in% colnames(cbmPools))) {
    stop("The carbon pool you specified for plotting is not contained in the pool definitions")
  }

  cbmPools <- as.data.table( cbmPools)
  #Build rasters for every year and pool
  carbonStacks <- vector(mode = "list", length = length(poolsToPlot))
  names(carbonStacks) <- poolsToPlot

  for (pool in poolsToPlot) {
    carbonStacks[[pool]] <- lapply(years, FUN = function(x, poolsDT = cbmPools,
                                                         var = pool,
                                                         pixelKeep =  pixelkeep) {

      poolsDT <- poolsDT[order(pixelGroup)] %>% #order by stand index
        .[simYear == x, .(pixelGroup, "var" =  get(pool))] #subset by year
      #subset  pixelKeep
      colsToKeep <- c("pixelIndex", paste0("pixelGroup", x))

      pixelKeep <- pixelKeep[, colsToKeep, with = FALSE] %>%
        setnames(., c("pixelIndex", "pixelGroup"))
      # with=FALSE tells data.table colsToKeep isn't a column name until here it
      # is ok...then in 1993 - an extra line gets added from the merge below
      # Keep <- poolsDT[pixelKeep, on = c("pixelGroup")]
      pixelKeep <- merge(pixelKeep, poolsDT, by = "pixelGroup", all.x = TRUE) %>% #join with pixelKeep
        .[order(pixelIndex)] #order by rowOrder for raster prep
      #pixelKeep <- pixelKeep[poolsDT, on = c("pixelGroup")] %>% #join with pixelKeep
      pixels <- getValues(masterRaster)

      plotMaster <- raster(masterRaster)
      plotMaster[] <- 0
      plotMaster[pixelKeep$pixelIndex] <- pixelKeep$var
      # masterRaster[masterRaster == 0] <- NA #Species has zeroes instead of NA. Revisit if masterRaster changes
      # masterRaster[!is.na(masterRaster)] <- pixelKeep$var

      #name will begin with x if no character assigned
      return(plotMaster)
    })
  }

  names(carbonStacks) <- paste0(poolsToPlot)
  unlist(carbonStacks)
  return(carbonStacks)
}

## FRIruns
FRIresultRasters <- spatialRaster(
  pixelkeep = RIAfriRuns$pixelKeep,
  cbmPools = RIAfriRuns$cbmPools,
  poolsToPlot = "totalCarbon",
  years = c(2020, 2540),
  masterRaster = RIAfriRuns$masterRaster)

writeRaster(FRIresultRasters$totalCarbon[[1]], filename = file.path(outputDir,"FRI","TotalCarbon2020.tif"))
writeRaster(FRIresultRasters$totalCarbon[[2]], filename = file.path(outputDir,"FRI","TotalCarbon2540.tif"))

# ## modified NPP
# NPPplot <- function(spatialDT, NPP, masterRaster) {
#   # Calculate the avgNPP (MgC/ha) by pixel group.
#   npp <- as.data.table(copy(NPP))
#   npp[,avgNPP := mean(NPP), by = c("pixelGroup")]
#   cols <- c("simYear", "NPP")
#   avgNPP <- unique(npp[, (cols) := NULL])
#   # link that to the pixels
#   t <- spatialDT[, .(pixelIndex, pixelGroup)]
#   setkey(t,pixelGroup)
#   setkey(avgNPP,pixelGroup)
#   temp <- merge(t, avgNPP, on = "pixelGroup")
#   setkey(temp, pixelIndex)
#   #pixelCount[which(is.na(pixelCount$N)),"N"] <- 0
#   # temp1 <- temp[which(!is.na(temp$simYear)),.(pixelIndex,NPP)]
#   # temp1[order(pixelIndex)]
#   #masterRaster[!masterRaster == 0] <- temp$NPP
#   plotMaster <- raster(masterRaster)
#   plotMaster[] <- 0
#   # instead of tC/ha for each pixel,
#   plotMaster[temp$pixelIndex] <- temp$avgNPP
#   #pixel size in ha
#   pixSize <- prod(res(masterRaster))/10000
#   temp[,pixNPP := avgNPP*pixSize]
#   overallAvgNpp <- sum(temp$pixNPP)/(nrow(temp)*pixSize)
#   return(overallAvgNpp)
# }
#
# quickPlot::Plot(plotMaster, new = TRUE,
#                 title = paste0("Pixel-level average NPP MgC/ha/yr.",
#                                "\n Landscape average: ", round(overallAvgNpp,3), "  MgC/ha/yr."))


NPPplot(
  spatialDT = RIAfriRuns$spatialDT,
  NPP = RIAfriRuns$NPP,
  masterRaster = RIAfriRuns$masterRaster
)

clearPlot()

carbonOutPlot(
  emissionsProducts = RIAfriRuns$emissionsProducts,
  masterRaster = RIAfriRuns$masterRaster
)

FRIageDist2540hist <- qplot(RIAfriRuns$spatialDT$ages, geom = "histogram")
FRIageDist2020hist <- qplot(RIAfriRuns$allPixDT[!is.na(ages),]$ages, geom = "histogram")
clearPlot()
Plot(FRIageDist2540hist)
Plot(FRIageDist2020hist, addTo = TRUE)


### fire runs results and fiddling ############################################

#Outputs for Sam
# getting the individual matrices used in the c-transactions in the sims
fireSimTransferMatrices <- lapply(opMatrices, unique)
fireReturnsDisturbancesMids <- fireSimTransferMatrices$disturbance
fireReturnsDomTurnMids <- fireSimTransferMatrices$domturnover
fireReturnsBioTurnMids <- fireSimTransferMatrices$bioturnover
fireReturnsDomDecayMids <- fireSimTransferMatrices$domDecay
fireReturnsSlowDecayMids <- fireSimTransferMatrices$`slow decay`
fireReturnsSlowMixMids <- fireSimTransferMatrices$`slow mixing`


# looking for old  stands to see what happens passed the oldest age on the growth curves
fireReturnSims <- spadesCBMout
fireReturnSims$cbmPools[simYear == 2540 & ages>400,]
# pixel group 1 is 531 years old
fireReturnSims$cbmPools[1:5,]
# pixelGroup ==1 is in the 1st position in the $allProcesses$Growth1 list
fireReturnSims$allProcesses$Growth1[[1]]

# trying to find the right gcid - fireReturnsSims$growth_increments have id
# column going from 1 to 310. These are factor levels created from  fireReturnSims$curveID
# [1] "growth_curve_component_id" "ecozones". Need to figure out which is matching pixelGroup 1
curveID <- fireReturnSims$curveID
gcidsLevels <- levels(fireReturnSims$level3DT$gcids)
# ecozone and growth_curve_component_id for pixelGroup 1
fireReturnSims$pixelGroupC[pixelGroup == 1,]
# find those
which(gcidsLevels == "4001001_4")
pg1gc <- as.data.table(fireReturnSims$growth_increments[
  fireReturnSims$growth_increments[,1] == which(gcidsLevels == "4001001_4"),])
pg1gc[age == max(age),]
### fire runs results and fiddling END ############################################


### presentDay runs checking

presentDayResultRasters <- spatialRaster(
  pixelkeep = RIApresentDayRuns$pixelKeep,
  cbmPools = RIApresentDayRuns$cbmPools,
  poolsToPlot = "totalCarbon",
  years = c(1985, 2015),
  masterRaster = RIApresentDayRuns$masterRaster)

writeRaster(presentDayResultRasters$totalCarbon[[1]], filename = file.path(outputDir,"presentDay","TotalCarbon1985.tif"))
writeRaster(presentDayResultRasters$totalCarbon[[2]], filename = file.path(outputDir,"presentDay","TotalCarbon2015.tif"))


NPPplot(
  spatialDT = RIApresentDayRuns$spatialDT,
  NPP = RIApresentDayRuns$NPP,
  masterRaster = RIApresentDayRuns$masterRaster
)

clearPlot()

carbonOutPlot(
  emissionsProducts = RIApresentDayRuns$emissionsProducts,
  masterRaster = RIApresentDayRuns$masterRaster
)


presentDayTotCrasters <- CBMutils::plotCarbonRasters(
  pixelkeep = RIApresentDayRuns$pixelKeep,
  cbmPools = RIApresentDayRuns$cbmPools,
  poolsToPlot = "totalCarbon",
  years = c(1985, 2015),
  masterRaster = RIApresentDayRuns$masterRaster
)
presentDayRunsAgeDist2015hist <- qplot(RIApresentDayRuns$spatialDT$ages, geom = "histogram")
presentDayRunsAgeDist1985hist <- qplot(RIApresentDayRuns$allPixDT[!is.na(ages),]$ages, geom = "histogram")
Plot(ageDist1985hist)
Plot(ageDist2015hist, addTo = TRUE)
