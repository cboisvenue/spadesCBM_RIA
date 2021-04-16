NPPplot(
  spatialDT = spadesCBMout$spatialDT,
  NPP = spadesCBMout$NPP,
  masterRaster = spadesCBMout$masterRaster
)

carbonOutPlot(
  emissionsProducts = spadesCBMout$emissionsProducts,
  masterRaster = spadesCBMout$masterRaster
)

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
