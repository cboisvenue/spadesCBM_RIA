## Figure for RIA runs -------------------------------------------------------

library(Require)
Require("magrittr") # this is needed to use "%>%" below
Require("SpaDES.core")

#install_github("PredictiveEcology/CBMutils@development")
#load_all("~/GitHub/PredictiveEcology/CBMutils")
#Require("PredictiveEcology/CBMutils (>= 0.0.6)")
devtools::load_all("C:/Celine/github/CBMutils")

options("reproducible.useRequire" = TRUE)
library(data.table)
library(raster)
library(ggplot2)

# read-in all the results from the paper simulations
RIAfriRuns <- loadSimList("C:/Celine/github/spadesCBM_RIA/outputs/FRI/RIAfriRuns")
RIApresentDayRuns <- loadSimList("C:/Celine/github/spadesCBM_RIA/outputs/presentDay/RIApresentDayRuns")
RIAharvest1Runs <- loadSimList("C:/Celine/github/spadesCBM_RIA/outputs/harvest1/RIAharvest1Runs.rds")
RIAharvest2Runs <- loadSimList("C:/Celine/github/spadesCBM_RIA/outputs/harvest2/RIAharvest2Runs.rds")

#----------------------------------------------------------------------------

# - raster plotting function for paper -----------------------------------------
# the CBMutils::spatialPlot function plot rasters directly, but does not return
# rasters. Here I modify it to return rasters for more flexibility in presenting
# results


spatialRaster <- function(pixelkeep, cbmPools, poolsToPlot, years, masterRaster) {
  cbmPools[is.na(cbmPools)] <- 0
  colnames(cbmPools)[c(1,3,4)] <- c("simYear", "pixelGroup", "age")
  # totalCarbon
  if ("totalCarbon" %in% poolsToPlot) {
    totalCarbon <- apply(cbmPools[, SoftwoodMerch:HardwoodBranchSnag], 1, "sum")
    cbmPools <- cbind(cbmPools, totalCarbon)
  }
## Add AG and BG options here
  if ("aboveGround" %in% poolsToPlot) {
    colsAG <- c("SoftwoodMerch", "SoftwoodFoliage", "SoftwoodOther",
                "HardwoodMerch", "HardwoodFoliage", "HardwoodOther",
                "SoftwoodStemSnag", "SoftwoodBranchSnag",
                "HardwoodStemSnag", "HardwoodBranchSnag",
                "AboveGroundVeryFastSoil", "AboveGroundFastSoil",
                "AboveGroundSlowSoil")
    aboveGround <- apply(cbmPools[, ..colsAG], 1, "sum")
    cbmPools <- cbind(cbmPools, aboveGround)
  }
  ## belowGround
  if ("belowGround" %in% poolsToPlot) {
    colsBG <- c("SoftwoodCoarseRoots", "SoftwoodFineRoots",
                "HardwoodCoarseRoots", "HardwoodFineRoots",
                "BelowGroundVeryFastSoil",
                "BelowGroundFastSoil", "MediumSoil",
                "BelowGroundSlowSoil")
    belowGround <- apply(cbmPools[, ..colsBG], 1, "sum")
    cbmPools <- cbind(cbmPools, belowGround)
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
# END raster plotting function for paper -----------------------------------------

## make rasters, plot rasters and save plots ----------------------

FRIresultRasters <- spatialRaster(
  pixelkeep = RIAfriRuns$pixelKeep,
  cbmPools = RIAfriRuns$cbmPools,
  poolsToPlot = c("totalCarbon","aboveGround"),
  years = c(2020, 2540),
  masterRaster = RIAfriRuns$masterRaster)

presentDayResultRasters <- spatialRaster(
  pixelkeep = RIApresentDayRuns$pixelKeep,
  cbmPools = RIApresentDayRuns$cbmPools,
  poolsToPlot = c("totalCarbon","aboveGround"),
  years = c(1985, 2015),
  masterRaster = RIApresentDayRuns$masterRaster)

harv1baseResultRasters <- spatialRaster(
  pixelkeep = RIAharvest1Runs$pixelKeep,
  cbmPools = RIAharvest1Runs$cbmPools,
  poolsToPlot = c("totalCarbon","aboveGround"),
  years = c(2020, 2099),
  masterRaster = RIAharvest1Runs$masterRaster)

harv2lessResultRasters <- spatialRaster(
  pixelkeep = RIAharvest2Runs$pixelKeep,
  cbmPools = RIAharvest2Runs$cbmPools,
  poolsToPlot = c("totalCarbon","aboveGround"),
  years = c(2020, 2099),
  masterRaster = RIAharvest2Runs$masterRaster)

## Total C raster figures
## start + present day end
clearPlot()
Plot(FRIresultRasters$totalCarbon[[1]], title = "a) Carrying capacity scenario (year 2020)
     Total Carbon t/ha C")
Plot(presentDayResultRasters$totalCarbon[[2]], title = "b) Present day scenario (year 2015)
     Total Carbon t/ha C")
Plot(harv1baseResultRasters$totalCarbon[[1]], title = "c) Base harvest scenario (year 2020)
     Total Carbon t/ha C")
Plot(harv2lessResultRasters$totalCarbon[[1]], title = "d) Less harvest scenario (year 2020)
     Total Carbon t/ha C")
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/rasters/totCstartFor3presenDayEnd", type = "png")
# all end of sims
clearPlot()
Plot(FRIresultRasters$totalCarbon[[2]],  title = "a) Carrying capacity scenario (year 2540)
     Total Carbon t/ha C")
Plot(presentDayResultRasters$totalCarbon[[2]], title = "b) Present scenario (year 2015)
     Total Carbon t/ha C")
Plot(harv1baseResultRasters$totalCarbon[[2]], title = "c) Base harvest scenario (year 2099)
     Total Carbon t/ha C")
Plot(harv2lessResultRasters$totalCarbon[[2]], title = "d) Less harvest scenario (year 2099)
     Total Carbon t/ha C")
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/rasters/totCend", type = "png")
## aboveGround raster figures
## start + present day end
clearPlot()
Plot(FRIresultRasters$aboveGround[[1]], title = "a) Carrying capacity scenario (year 2020)
     Above Ground Carbon t/ha C")
Plot(presentDayResultRasters$aboveGround[[2]], title = "b) Present scenario (year 2015)
     Above Ground Carbon t/ha C")
Plot(harv1baseResultRasters$aboveGround[[1]], title = "c) Base harvest scenario (year 2020)
     Above Ground Carbon t/ha C")
Plot(harv2lessResultRasters$aboveGround[[1]], title = "d) Less harvest scenario (year 2020)
     Above Ground Carbon t/ha C")
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/rasters/AGstartFor3presenDayEnd", type = "png")
# all end of sims
clearPlot()
Plot(FRIresultRasters$aboveGround[[2]],  title = "a) Carrying capacity scenario (year 2540)
     Above Ground Carbon t/ha C")
Plot(presentDayResultRasters$aboveGround[[2]], title = "b) Present scenario (year 2015)
     Above Ground Carbon t/ha C")
Plot(harv1baseResultRasters$aboveGround[[2]], title = "c) Base harvest scenario (year 2099)
     Above Ground Carbon t/ha C")
Plot(harv2lessResultRasters$aboveGround[[2]], title = "d) Less harvest scenario (year 2099)
     Above Ground Carbon t/ha C")
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/rasters/AGend", type = "png")

## make rasters, plot rasters and save plots ----------------------

## plot age class distributions and save plots ----------------------

### Age class distributions
#
# FRIageDist2540hist <- hist(RIAfriRuns$spatialDT$ages, plot=FALSE)
# FRIageDist2020hist <- hist(RIAfriRuns$allPixDT[!is.na(ages),]$ages, plot=FALSE)
# presentDayAgeDist2015hist <- hist(RIApresentDayRuns$spatialDT$ages, plot=FALSE)
# presentDayAgeDist1985hist <- hist(RIApresentDayRuns$allPixDT[!is.na(ages),]$ages, plot=FALSE)
# harv1ageDist2099hist <- hist(RIAharvest1Runs$spatialDT$ages, plot=FALSE)
# harv1ageDist2020hist <- hist(RIAharvest1Runs$allPixDT[!is.na(ages),]$ages, plot=FALSE)
# harv2ageDist2099hist <- hist(RIAharvest2Runs$spatialDT$ages, plot=FALSE)
# harv2ageDist2020hist <- hist(RIAharvest2Runs$allPixDT[!is.na(ages),]$ages, plot=FALSE)
#
# # end 4 scenarios
# clearPlot()
# Plot(FRIageDist2540hist,presentDayAgeDist2015hist,harv1ageDist2099hist,harv2ageDist2099hist)
# Plot(FRIageDist2540hist, title = "a) Carrying capacity scenario year 2540",
#      xlab = "Ages")
# Plot(presentDayAgeDist2015hist, title = "b) Present scenario year 2015",
#      xlab = "Ages", addTo = TRUE)
# Plot(harv1ageDist2099hist, title = "c) Base harvest scenario year 2099",
#      xlab = "Ages", addTo = TRUE)
# Plot(harv2ageDist2099hist, title = "d) Less harvest scenario year 2099",
#      xlab = "Ages", addTo = TRUE)
#
# savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/ageClassDists/endSimsAgeClassDists", type = "png")
# all scenario start

#
#                            main = "a) Carrying capacity scenario year 2540",
#                            xlab = "Ages", plot = FALSE)
# FRIageDist2020hist <- hist(RIAfriRuns$allPixDT[!is.na(ages),]$ages,
#                            main = "a) Carrying capacity scenario year 2020",
#                            xlab = "Ages")
# presentDayAgeDist2015hist <- hist(RIApresentDayRuns$spatialDT$ages,
#                                   main = "b) Present scenario year 2015",
#                                   xlab = "Ages")
# presentDayAgeDist1985hist <- hist(RIApresentDayRuns$allPixDT[!is.na(ages),]$ages,
#                                   main = "b) Present scenario year 1985",
#                                   xlab = "Ages")
# harv1ageDist2099hist <- hist(RIAharvest1Runs$spatialDT$ages,
#                              main =  "c) Base harvest scenario year 2099",
#                              xlab = "Ages")
# harv1ageDist2020hist <- hist(RIAharvest1Runs$allPixDT[!is.na(ages),]$ages,
#                               main =  "c) Base harvest scenario year 2020",
#                               xlab = "Ages")
# harv2ageDist2099hist <- hist(RIAharvest2Runs$spatialDT$ages,
#                              main =  "d) Less harvest scenario year 2099",
#                              xlab = "Ages")
# harv2ageDist2020hist <- hist(RIAharvest2Runs$allPixDT[!is.na(ages),]$ages,
#                              main =  "d) Less harvest scenario year 2020",
#                              xlab = "Ages")

FRIageDist2540hist <- qplot(RIAfriRuns$spatialDT$ages, geom = "histogram",
                            main = "a) Carrying capacity scenario year 2540",
                            xlab = "Ages")
FRIageDist2020hist <- qplot(RIAfriRuns$allPixDT[!is.na(ages),]$ages, geom = "histogram",
                            main = "a) Carrying capacity scenario year 2020",
                            xlab = "Ages")
clearPlot()
Plot(FRIageDist2020hist, FRIageDist2540hist)
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/ageClassDists/AgeClassDistsFRI", type = "png")

presentDayAgeDist2015hist <- qplot(RIApresentDayRuns$spatialDT$ages, geom = "histogram",
                                   main = "b) Present scenario year 2015",
                                   xlab = "Ages")
presentDayAgeDist1985hist <- qplot(RIApresentDayRuns$allPixDT[!is.na(ages),]$ages,
                                   geom = "histogram",
                                   main = "b) Present scenario year 1985",
                                   xlab = "Ages")
clearPlot()
Plot(presentDayAgeDist1985hist, presentDayAgeDist2015hist)
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/ageClassDists/AgeClassDistsPresentDay", type = "png")

harv1ageDist2099hist <- qplot(RIAharvest1Runs$spatialDT$ages, geom = "histogram",
                              main =  "c) Base harvest scenario year 2099",
                              xlab = "Ages")
harv1ageDist2020hist <- qplot(RIAharvest1Runs$allPixDT[!is.na(ages),]$ages,
                              geom = "histogram",
                              main =  "c) Base harvest scenario year 2020",
                              xlab = "Ages")
clearPlot()
Plot(harv1ageDist2020hist, harv1ageDist2099hist)
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/ageClassDists/AgeClassDistsHarv1", type = "png")

harv2ageDist2099hist <- qplot(RIAharvest2Runs$spatialDT$ages, geom = "histogram",
                              main =  "d) Less harvest scenario year 2099",
                              xlab = "Ages")
harv2ageDist2020hist <- qplot(RIAharvest2Runs$allPixDT[!is.na(ages),]$ages,
                              geom = "histogram",
                              main =  "d) Less harvest scenario year 2020",
                              xlab = "Ages")
clearPlot()
Plot(harv2ageDist2020hist, harv2ageDist2099hist)
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/ageClassDists/AgeClassDistsHarv2", type = "png")

clearPlot()
Plot(harv1ageDist2099hist, harv2ageDist2099hist)
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/ageClassDists/AgeClassDistsBothHarv", type = "png")

clearPlot()
Plot(FRIageDist2540hist)
Plot(presentDayAgeDist2015hist)
Plot(harv1ageDist2099hist, addTo = TRUE)
Plot(harv2ageDist2099hist, addTo = TRUE)
savePlot(filename = "C:/Celine/github/spadesCBM_RIA/results/ageClassDists/AgeClassDistsEndAll", type = "png")
## Note that the harv1 and harv2 are bi-modal
## END plot age class distributions and save plots ----------------------

## function to sum carbon for totalCarbon or aboveGround or belowGround -------------
calcC <- function(cbmPools, poolToSum, masterRaster){
  #targetPool <- poolToSum
  # year <- time(RIApresentDayRuns)
  # cbmPools <- RIApresentDayRuns$cbmPools
  # masterRaster <- RIApresentDayRuns$masterRaster
  # calculate total carbon by pixelGroup
  if ("totalCarbon" %in% poolToSum) {
    targetPool <- apply(cbmPools[, SoftwoodMerch:HardwoodBranchSnag], 1, "sum")
    cbmPools <- cbind(cbmPools, targetPool)
  }
  ## Add AG and BG options here
  if ("aboveGround" %in% poolToSum) {
    colsAG <- c("SoftwoodMerch", "SoftwoodFoliage", "SoftwoodOther",
                "HardwoodMerch", "HardwoodFoliage", "HardwoodOther",
                "SoftwoodStemSnag", "SoftwoodBranchSnag",
                "HardwoodStemSnag", "HardwoodBranchSnag",
                "AboveGroundVeryFastSoil", "AboveGroundFastSoil",
                "AboveGroundSlowSoil")
    targetPool <- apply(cbmPools[, ..colsAG], 1, "sum")
    cbmPools <- cbind(cbmPools, targetPool)
  }
  ## belowGround
  if ("belowGround" %in% poolToSum) {
    colsBG <- c("SoftwoodCoarseRoots", "SoftwoodFineRoots",
                "HardwoodCoarseRoots", "HardwoodFineRoots",
                "BelowGroundVeryFastSoil",
                "BelowGroundFastSoil", "MediumSoil",
                "BelowGroundSlowSoil")
    targetPool <- apply(cbmPools[, ..colsBG], 1, "sum")
    cbmPools <- cbind(cbmPools, targetPool)
  }

  sumColsOnly <- cbmPools[,.(simYear,pixelCount, pixelGroup, targetPool)]
  ## check that all is good
  sumColsOnly[,sum(pixelCount), by=simYear]
  # simYear      V1
  # 1:    1985 3112425
  # 2:    1990 3112425
  # 3:    1995 3112425
  # 4:    2000 3112425
  # 5:    2005 3112425
  # 6:    2010 3112425
  # 7:    2013 3112425
  # 8:    2014 3112425
  # 9:    2015 3112425
  resInHa <- res(masterRaster)[1]*res(masterRaster)[2]/10000
  sumColsOnly[, absCarbon := (pixelCount*resInHa*targetPool)]
  landscapeCarbon <- sumColsOnly[,sum(absCarbon)/1000000, by = simYear]
  return(landscapeCarbon)

}
## END function to sum carbon----------------------------------------------

## Same function just for summing total carbon -----------------
calcTotalC <- function(cbmPools, masterRaster){
  # year <- time(RIApresentDayRuns)
  # cbmPools <- RIApresentDayRuns$cbmPools
  # masterRaster <- RIApresentDayRuns$masterRaster
  # calculate total carbon by pixelGroup
  totalCarbon <- apply(cbmPools[, SoftwoodMerch:HardwoodBranchSnag], 1, "sum")
  cbmPools <- cbind(cbmPools, totalCarbon)
  totColsOnly <- cbmPools[,.(simYear,pixelCount, pixelGroup, totalCarbon)]
  ## check that all is good
  totColsOnly[,sum(pixelCount), by=simYear]
  # simYear      V1
  # 1:    1985 3112425
  # 2:    1990 3112425
  # 3:    1995 3112425
  # 4:    2000 3112425
  # 5:    2005 3112425
  # 6:    2010 3112425
  # 7:    2013 3112425
  # 8:    2014 3112425
  # 9:    2015 3112425
  resInHa <- res(masterRaster)[1]*res(masterRaster)[2]/10000
  totColsOnly[, absCarbon := (pixelCount*resInHa*totalCarbon)]
  landscapeCarbon <- totColsOnly[,sum(absCarbon)/1000000, by = simYear]
  return(landscapeCarbon)

}
## END Same function just for summing total carbon -----------------

# these are in Megatonnes of C
presentDayTotalC <- calcTotalC(cbmPools = RIApresentDayRuns$cbmPools,
                                  masterRaster = RIApresentDayRuns$masterRaster)
presentDayAboveGroundC <- calcC(cbmPools = RIApresentDayRuns$cbmPools,
                       poolToSum = "aboveGround",
                       masterRaster = RIApresentDayRuns$masterRaster)
presentDayCresults <- as.data.table(cbind(presentDayTotalC, presentDayAboveGroundC$V1))
setnames(presentDayCresults, names(presentDayCresults),
         c("simYear", "TotalC", "AGC"))
write.csv(presentDayCresults, file = "C:/Celine/github/spadesCBM_RIA/results/presentDayAbsC.csv")


FRITotalC <- calcTotalC(cbmPools = RIAfriRuns$cbmPools,
                               masterRaster = RIAfriRuns$masterRaster)
FRIAGC <- calcC(cbmPools = RIAfriRuns$cbmPools,
                             poolToSum = "aboveGround",
                             masterRaster = RIAfriRuns$masterRaster)
FRICresults <- as.data.table(cbind(FRITotalC,FRIAGC$V1))
setnames(FRICresults, names(FRICresults),
         c("simYear", "TotalC", "AGC"))
write.csv(FRICresults, file = "C:/Celine/github/spadesCBM_RIA/results/FRIabsC.csv")


harv1TotalC <- calcTotalC(cbmPools = RIAharvest1Runs$cbmPools,
                               masterRaster = RIAharvest1Runs$masterRaster)
harv1AGC <- calcC(cbmPools = RIAharvest1Runs$cbmPools,
                  poolToSum = "aboveGround",
                          masterRaster = RIAharvest1Runs$masterRaster)
harv1Cresults <- as.data.table(cbind(harv1TotalC,harv1AGC$V1))
setnames(harv1Cresults, names(harv1Cresults),
         c("simYear", "TotalC", "AGC"))
write.csv(harv1Cresults, file = "C:/Celine/github/spadesCBM_RIA/results/harv1absC.csv")

harv2TotalC <- calcTotalC(cbmPools = RIAharvest2Runs$cbmPools,
                          masterRaster = RIAharvest2Runs$masterRaster)
harv2AGC <- calcC(cbmPools = RIAharvest2Runs$cbmPools,
                  poolToSum = "aboveGround",
                  masterRaster = RIAharvest2Runs$masterRaster)
harv2Cresults <- as.data.table(cbind(harv2TotalC,harv2AGC$V1))
setnames(harv2Cresults, names(harv2Cresults),
         c("simYear", "TotalC", "AGC"))
write.csv(harv2Cresults, file = "C:/Celine/github/spadesCBM_RIA/results/harv2absC.csv")










analyseNPP <- function(NPPDT, cbmPools, masterRaster){
  #NPPDT <- RIApresentDayRuns$NPP
  # adding the pixelCount for each pixelGroup for each year
  #cbmPools <- RIApresentDayRuns$cbmPools
  #masterRaster <- RIApresentDayRuns$masterRaster

  cbmPools <- cbmPools[,.(simYear, pixelCount, pixelGroup)]

  dt1 <- merge.data.table(NPPDT, cbmPools,
                             by = c("pixelGroup", "simYear"), all.x = TRUE)

  avgNPPbyHabyYr <- dt1[, .(avgNPPha = mean(NPP)), by = simYear]

  # calculate the absolute carbon update of the landscape for each year
  resInHa <- res(masterRaster)[1]*res(masterRaster)[2]/10000

  dt1[, absNPP := (pixelCount*resInHa*NPP), by = simYear]
  # this is in Mega tonnes of carbon
  absNPPbyYr <- dt1[, .(absNPP = sum(absNPP)/1000000), by = simYear]

  NPPtable <- as.data.table(cbind(avgNPPbyHabyYr,absNPP=absNPPbyYr$absNPP))

  return(NPPtable)

}

FRINPP <- analyseNPP(RIAfriRuns$NPP, RIAfriRuns$cbmPools, RIAfriRuns$masterRaster)
FRINPP[, scenario := "FRI"]
presentDayNPP <- analyseNPP(RIApresentDayRuns$NPP, RIApresentDayRuns$cbmPools, RIApresentDayRuns$masterRaster)
presentDayNPP[, scenario := "presentDay"]
harv1NPP <- analyseNPP(RIAharvest1Runs$NPP, RIAharvest1Runs$cbmPools, RIAharvest1Runs$masterRaster)
harv1NPP[, scenario := "harvBase"]
harv2NPP <- analyseNPP(RIAharvest2Runs$NPP, RIAharvest2Runs$cbmPools, RIAharvest2Runs$masterRaster)
harv2NPP[, scenario := "harvLess"]

simsNPP <- as.data.table(rbind(FRINPP, presentDayNPP, harv1NPP, harv2NPP))
write.csv(simsNPP, file = "C:/Celine/github/spadesCBM_RIA/results/simsNPP.csv")

# below not performed yet
  avgNPPabs <- mean(absNPPbyYr$absNPP)
  maxNPPbyYr <- max(absNPPbyYr$absNPP)
  minNPPbyYr <- min(absNPPbyYr$absNPP)

  avgNPPha <- mean(avgNPPbyHabyYr$avgNPPha)
  maxNPPha <- max(avgNPPbyHabyYr$avgNPPha)
  minNPPha <- min(avgNPPbyHabyYr$avgNPPha)









### Fiddling below

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


# ploting the spinup resutls
friSpinup <- data.table(RIAfriRuns$spinupResult)
friSpinup[, simYear := 0]
friSpinup[, ages := RIAfriRuns$level3DT$ages]
friSpinup[, pixelGroup := RIAfriRuns$level3DT$pixelGroup]
friSpinup
neworder <- c(27, 1, 29, 28, 2:26)
setcolorder(friSpinup, neworder)

FRIspinupRaster <- spatialRaster(
  pixelkeep = RIAfriRuns$pixelKeep,
  cbmPools = friSpinup,
  poolsToPlot = "totalCarbon",
  years = 0,
  masterRaster = RIAfriRuns$masterRaster)

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


### presentDay runs checking ######################################################################

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
### end presentDay runs checking ######################################################################

### harvest1 runs checking ######################################################################
harv1resultRasters <- spatialRaster(
  pixelkeep = RIAharvest1Runs$pixelKeep,
  cbmPools = RIAharvest1Runs$cbmPools,
  poolsToPlot = "totalCarbon",
  years = c(2020, 2099),
  masterRaster = RIAharvest1Runs$masterRaster)

Plot(harv1resultRasters$totalCarbon[[1]], title = "Total Carbon (t/ha) in 2020")
Plot(harv1resultRasters$totalCarbon[[2]], title = "Total Carbon (t/ha) in 2099")

writeRaster(harv1resultRasters$totalCarbon[[1]], filename = file.path(outputDir,"harvest1","TotalCarbon2020.tif"))
writeRaster(harv1resultRasters$totalCarbon[[2]], filename = file.path(outputDir,"harvest1","TotalCarbon2099.tif"))

harv1ageDist2099hist <- qplot(RIAharvest1Runs$spatialDT$ages, geom = "histogram")
harv1ageDist2020hist <- qplot(RIAharvest1Runs$allPixDT[!is.na(ages),]$ages, geom = "histogram")
clearPlot()
Plot(harv1ageDist2099hist, title = "Age class Distribution 2099")
Plot(harv1ageDist2020hist, title = "Age class Distribution 2020", addTo = TRUE)

# ploting the spinup resutls
harv1Spinup <- data.table(RIAharvest1Runs$spinupResult)
harv1Spinup[, simYear := 0]
harv1Spinup[, ages := RIAharvest1Runss$level3DT$ages]
harv1Spinup[, pixelGroup := RIAharvest1Runs$level3DT$pixelGroup]
harv1Spinup
neworder <- c(27, 1, 29, 28, 2:26)
setcolorder(harv1Spinup, neworder)

harv1spinupRaster <- spatialRaster(
  pixelkeep = RIAharvest1Runs$pixelKeep,
  cbmPools = harv1Spinup,
  poolsToPlot = "totalCarbon",
  years = 0,
  masterRaster = RIAharvest1Runs$masterRaster)
