

# old rasters from Greg 2015 I think.
# rastersOLD <- prepInputs(
#   url = "https://drive.google.com/file/d/1DN31xcXh97u6v8NaVcy0O3vzKpLpld69/view?usp=sharing",
#   fun = "raster::stack",
#   #rasterToMatch = masterRaster, useGDAL = FALSE) # this was Eliot's
#   destinationPath = 'inputs')
#
# stackIan <- prepInputs(url = "https://drive.google.com/file/d/1keh-3KIvgy8G5nBBxolH_LsoC53w3qqN/view?usp=sharing",
#                        fun = "raster::stack",
#                        destinationPath = 'inputs')
# # current masterRaster
# fireRetMasterRaster <- spadesCBMout$masterRaster
#
# # trying things
# rastersOLDt1 <- prepInputs(
#   url = "https://drive.google.com/file/d/1DN31xcXh97u6v8NaVcy0O3vzKpLpld69/view?usp=sharing",
#   fun = "raster::stack",
#   rasterToMatch = spadesCBMout$masterRaster,
#   useGDAL = FALSE,
#   destinationPath = 'inputs')
#
# age <- postProcess(x = stackIan$RIAlandscape.5,
#                    rasterToMatch = spadesCBMout$masterRaster,
#                    filename2 = "inputs/pixelAges.tif")
#
# aValues <- age[]
# table(aValues, useNA = 'ifany')
# mValues <- spadesCBMout$masterRaster[]
# table(mValues)
# # there are this many more NAs in the age raster provided by Ian
# # 5435395 - 5605119
# # [1] -169724

# and here they are
# agesDT <- data.table(spadesCBMout$allPixDT[,.(ages, pixelIndex)],aValues)
# setnames(agesDT,c("ages","aValues"), c("ages2020", "ages2015"))
# agesNA <- agesDT[is.na(ages2015) & !is.na(ages2020),]
#
# # need to check if these pixels are disturbed
# # are they burned?
# presentDayBurns <- prepInputs(url = 'https://drive.google.com/file/d/1MjQ5y9Txr1ezv0BatG4b_6FpM6wti1b5/view?usp=sharing',
#                               destinationPath = 'inputs',
#                               overwrite = TRUE,
#                               fun = 'readRDS')
# # there are 3426/169724 that get burned between 1985-2015
# agesNA[pixelIndex %in% presentDayBurns$pixelID,]
#
# # are they harvested?
# presentDayHarvest <- prepInputs(url = 'https://drive.google.com/file/d/1Ca-kPun7_VF2xV8s40IJ9hri3Ok-6qx5/view?usp=sharing',
#                               destinationPath = 'inputs',
#                               overwrite = TRUE,
#                               fun = 'readRDS')
# # there are 14969/169724 that get burned between 1985-2015
# agesNA[pixelIndex %in% presentDayHarvest$pixelID,]
#
# (14969 + 3426)/169724
# #only 11%


#read-in the VRI2015
# this is an ESRI file
#https://pub.data.gov.bc.ca/datasets/02dba161-fdb7-48ae-a4bb-bd6ef017c36d/2015/VEG_COMP_LYR_L1_POLY_2015.gdb.zip

# # back-tracking VRI2015 age raster to 1985 based on the disturbance info
# presentDayBurns[, events := 1L]
# presentDayHarvest[, events := 2L]
# allDist <- rbind(presentDayBurns, presentDayHarvest)
#
# # there needs to be a part here to make the right age raster
# RIArtm <- prepInputs(url = "https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing",
#                      destinationPath = 'inputs')
# VRI2015ageRaster <- Cache(prepInputsVRIage, VRIurl = "https://drive.google.com/file/d/1spcJq_4r4TRA2jfjclQcwBcIfnbniEXZ",
#                             dPath = "C:/Celine/github/spadesCBM_RIA/modules/CBM_dataPrep_RIA/data",
#                             rasterToMatch = RIArtm,
#                             cacheRepo = cacheDir)
# # make a data table from the raster
#
#
# # this is the wrong raster but waiting for the VRI2015
# ages2015dt <- agesDT[,.(pixelIndex, ages2015)]
#
#
#
#
# prepInputsVRIage <- function(VRIurl, dPath, rasterToMatch, field = "PROJ_AGE_1"){
#   VRIin <- prepInputs(url = VRIurl,
#                       fun = "sf::st_read",
#                       destinationPath = dPath)
#   RIA_VRI <- st_transform(VRIin, crs = st_crs(rasterToMatch))
#   #gcIDRaster <- fasterize::fasterize(RIA_VRI, rasterToMatch, field = "curve2")
#   ageRaster <- fasterize::fasterize(RIA_VRI, rasterToMatch, field = "PROJ_AGE_1")
#   #gcIDRaster[] <- as.integer(gcIDRaster[])
#   ageRaster[] <- as.integer(ageRaster[])
#   VRIageRaster <- raster::raster(ageRaster)
#   return(VRIageRaster)
# }
#
# # Eliot modified reproducible so prepInputs can deal with .gbd files
# remotes::install_github("PredictiveEcology/reproducible@gdb_archiveNA")
# # I am using "dummy" masterRaster --> use your real one
# masterRaster <- raster(extent(-1963750, -1321250, 7407500, 8239000),
#                        crs = "+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77 +x_0=0 +y_0=0 +ellps=GRS80 +units=m", res = 250)
# # Make a studyArea --> use your real one
# sa <- as(extent(masterRaster), "SpatialPolygons")
# crs(sa) <- crs(masterRaster)
# loadAge <- function(x, field = "PROJ_AGE_1") {
#   # a <- Cache(sf::st_read, x) # I used Cache during my development because this takes 37 minutes to run -- I was sick of running it again and again
#   a <- sf::st_read(x)
#   a1 <- a[, field]
#   return(a1)
# }
# a <- prepInputs(url =
#                   "https://pub.data.gov.bc.ca/datasets/02dba161-fdb7-48ae-a4bb-bd6ef017c36d/2015/VEG_COMP_LYR_L1_POLY_2015.gdb.zip",
#                 fun = quote(loadAge(x = targetFilePath,
#                                     field = "PROJ_AGE_1")),
#                 targetFile = "VEG_COMP_LYR_L1_POLY_2015.gdb.zip", archive = NA,
#                 studyArea = sa
#                 #rasterTomatch = masterRaster
# )

# extra tries that failed
# asked for help but got "do it another way" answer

# pickAges <- function(gcID, pickN, allPixDT){
#                       ageDist <- allPixDT[growth_curve_id == as.numeric(gcID) & ages >80,]$ages
#                       return(ageDist[round(runif(
#                         n = as.numeric(pickN),
#                         min = 1,
#                         max = length(ageDist)),digits = 0)])
# }
#
# newVector <- inDT[, lapply(.SD, FUN = pickAges(gcID = growth_curve_id,
#                                                             pickN = N,
#                                                             allPixDT = egPixDT)),
#                                by = 'growth_curve_id']
#
#
# newAgesDistPix <- rndNumBygcID[, lapply(.SD, FUN = pickAges(gcID = growth_curve_id,
#                                                             pickN = N,
#                                                             allPixDT = sim$allPixDT)),
#                                by = 'growth_curve_id']

#newPixDT <- egPixDT[inDT, on = "growth_curve_id", nomatch = 0] # join the 2 objects on "growth_curve_id"
# ageList <- as.list(vector(mode = "list", length = nrow(rndNumBygcID)))
# for(i in 1:length(rndNumBygcID$growth_curve_id)){
#   ageDist <- matchPixDT[growth_curve_id == as.numeric(rndNumBygcID[i,1]),]
#   ageList[[i]] <- sample(ageDist$ages, size = as.numeric(rndNumBygcID[i,2]))
# }
#
# ageNew <- rbindlist(ageList)
#
## another try



############ From the browser post VRI2015 reading-in###########################

# starting from the read-in age raster in 2015 and keeping the same masterRaster
# as the fireReturnInterval runs.

age2015 <- ageRaster2015[]
age2020 <- ageRaster[]
# NAs 2015 5539822
# NAs 2020 5435395


## figure out what ages the NAs have in 2020.
## Note that some ages in the 2020 that have no ages in the 2015 raster AND that
## are disturbed (so in the sim$disturbanceRasters data table) are old...older
## then the age from the dist...but most are below 50.


ageDT <- data.table(pixelIndex = sim$allPixDT$pixelIndex, age2015, age2020)
ageNoMatch <- ageDT[is.na(age2015) & !is.na(age2020),]
## NOTEs: there are 104427 more pixels with NAs in age2015 than in age2020.
ageDT[!is.na(age2015) & is.na(age2020),]
## There are no pixels that have ages in 2015 and are NAs-age in 2020.

# ageNA1985, what is their year of disturbance?

setkeyv(ageNoMatch, "pixelIndex")
setkeyv(allDist,"pixelID")

# do any pixels burn twice?
length(unique(presentDayBurns$pixelID)) == dim(presentDayBurns)[1]
#TRUE
# harvested twice?
length(unique(presentDayHarvest$pixelID)) == dim(presentDayHarvest)[1]
#TRUE

# get only the burnt/harvested pixels in the noMatch
ageNAburns <- presentDayBurns[pixelID %in% ageNoMatch$pixelIndex, ]
setnames(ageNAburns,"pixelID", "pixelIndex")
ageNAharvest <- presentDayHarvest[pixelID %in% ageNoMatch$pixelIndex, ]
setnames(ageNAharvest,"pixelID", "pixelIndex")
# Merge them with ageNoMatch
#dt_a[dt_b, on = .(b = y)]
ageNA1985 <- merge.data.table(ageNoMatch, ageNAburns, by = "pixelIndex", all.x = TRUE)
setnames(ageNA1985,"year", "burnYear")
ageNA1985[, events := NULL]
ageNA1985 <- merge.data.table(ageNA1985, ageNAharvest, by = "pixelIndex", all.x = TRUE)
setnames(ageNA1985,"year", "harvestYear")
ageNA1985[, events := NULL]
ageNA1985$noDist <- 0
ageNA1985[which(is.na(burnYear) & is.na(harvestYear)),]$noDist <- 1

# make a column of "straight substraction"
ageNA1985[, substract := (age2020 - 35)]

## pixels that have no disturbance and are >0 in 1985 get this age ##################
ageNA1985$PixSubtract1985 <- 999
ageNA1985[noDist == 1 & substract >= 0,]$PixSubtract1985 <- ageNA1985[noDist == 1 & substract >= 0,]$substract
ageNA1985[PixSubtract1985 == 999,]
# still 23177 that are not dealt with.

## the ones that have a non-negative age at the time of disturbance will get the
## age in the year prior to disturbance
ageNA1985[, PixFireYrAge := age2020 - (2020-burnYear)]
ageNA1985[, PixCutYrAge := age2020 - (2020-harvestYear)]
## this gives the ages of the burnt pixels in 1985 (values>0)
ageNA1985[, firePix1985age := PixFireYrAge - (burnYear - 1985)]
ageNA1985[, cutPix1985age := PixCutYrAge - (harvestYear - 1985)]

## how many left?
ageNA1 <- ageNA1985[PixSubtract1985 != 999] #81250
ageNA2 <- ageNA1985[ PixSubtract1985 == 999 & (firePix1985age > 0 | cutPix1985age > 0),] #3835

#104427-85085
ageNA1985[pixelIndex %in% ageNA1$pixelIndex | pixelIndex %in% ageNA2$pixelIndex] #85085
probPix <- ageNA1985[!(pixelIndex %in% ageNA1$pixelIndex | pixelIndex %in% ageNA2$pixelIndex)] #19342
# remove these
#colsRemove <- c("substract", "PixSubtract1985", "PixFireYrAge")
#probPix[, c("substract", "PixSubtract1985", "PixFireYrAge") := list(NULL, NULL, NULL)]

# how many of those are disturbed?
table(probPix$noDist)
# 0     1
# 11183  8159

# look at disturbed pixels ages in 1985
probCuts1985 <- qplot(probPix[noDist == 0,]$cutPix1985age, geom = "histogram")
probBurns1985 <- qplot(probPix[noDist == 0,]$firePix1985age, geom = "histogram")

## decision: all the disturbed pixels (noDist == 0) that were cut or burnt will
## have mature ages at the year of disturbance. Those will be selected from the
## age distribution of pixels with the same gcID with ages >80 in age2020.

distPixToMatch <- probPix[noDist == 0,]
length(unique(distPixToMatch$pixelIndex)) # no duplicats
# figure out the gcID for each pixelIndex
# get all the ages>80 at ages2020 in each gcID
# randomly select one for each pixelIndex

#gcIDdistPixNoMatch <- unique(sim$allPixDT[pixelIndex %in% distPixToMatch$pixelIndex,]$growth_curve_id)
#84
rndNumBygcID <- sim$allPixDT[pixelIndex %in% distPixToMatch$pixelIndex, .(pixelIndex, growth_curve_id)][,.N, by = "growth_curve_id"]

matchPixDT <- sim$allPixDT[growth_curve_id %in% rndNumBygcID$growth_curve_id & ages>80,.(growth_curve_id, ages)]

newAgedistPix1 <- matchPixDT[rndNumBygcID, on = "growth_curve_id", nomatch = 0]

#newPixDT <- newPixDT[, .SD[sample(.N, N)], by = "growth_curve_id"] # take a sample of each growth_curve_id, length N
#newAgedistPix <- newAgedistPix1[, lapply(.SD, sample(ages, size = N[1])), by = growth_curve_id]

newAgedistPix <- newAgedistPix1[, .SD[sample(ages,size = N[1])], by = "growth_curve_id"]

# This is just a test -- does each growth curve id have the same # rows as N says
newAgedistPix[, .N == N[1],by = "growth_curve_id"]

# re-attach a pixelIndex
pixIndDistPixToMatch <- sim$allPixDT[pixelIndex %in% distPixToMatch$pixelIndex, .(pixelIndex, growth_curve_id)]
setorder(pixIndDistPixToMatch, growth_curve_id)
setorder(newAgedistPix, growth_curve_id)
set(pixIndDistPixToMatch, NULL, "distPixNegAge1985", newAgedistPix$ages)
pixIndDistPixToMatch[, growth_curve_id := NULL]
ageNA1985 <- merge.data.table(ageNA1985, pixIndDistPixToMatch, on = "pixelIndex", all.x = TRUE)

##
noDistPixToMatch <- probPix[noDist == 1,]

negAgeHist <- hist(noDistPixToMatch$substract, plot = FALSE)
Plot(negAgeHist)

## trying to find the closest pixel with a disturbance in the 1985-2015
## disturbance data.table.
# are there pixels that are disturbed twice in this data?
countDist <- allDist[, .N, by = pixelID]
# yes
# are those in the proPix?
probPix[pixelIndex %in% countDist[N>1,]$pixelID,] #no

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
# maybe too small?
f9 <- focal(allDistRaster,
            w=matrix(1,nrow=3,ncol=3),
            fun = f3)

# agg that column to the all pixels DT
set(ageDT, NULL, "f9", f9[])
# do any of the 3X3 windows cover the last 8159 pixels?
checkF9 <- ageDT[pixelIndex %in% noDistPixToMatch$pixelIndex,]
table(checkF9$f9, useNA = "ifany")
# 1  NaN
# 4410 3749
# figure out what year and what disturbances
# some pixels are disturbed twice but they are not in my probPix
#
# make a raster with the dist year
yrDistRaster <- raster::raster(masterRaster)
setnames(countDist, "pixelID", "pixelIndex")
setnames(allDist, "pixelID", "pixelIndex")
yrDist <- unique(countDist[allDist, on = "pixelIndex", nomatch = 0][,.(pixelIndex, year)])
yrDistRaster[] <- NA
yrDistRaster[yrDist$pixelIndex] <- yrDist$year
yrf9 <- focal(yrDistRaster,
              w=matrix(1,nrow=3,ncol=3),
              fun = f3)
set(ageDT, NULL, "yrf9", yrf9[])

# create a new "distPixels" DT adding the dist to the f9==1
f9pix <- ageDT[f9==1]
f9dist <- merge.data.table(allDist, f9pix, all = TRUE)
# add the growth_curve_id to help the match??
f9dist <- f9dist[sim$allPixDT, on = 'pixelIndex', nomatch = 0][,.(pixelIndex, year, events, age2015, age2020, f9, yrf9, growth_curve_id)]
f9dist[, targetPix := 0L]
f9dist[pixelIndex %in% noDistPixToMatch$pixelIndex]$targetPix <- 1

# are there any pixels with targetPix == 1 and yrf9 !is.na()?
f9dist[targetPix >0 &!is.na(yrf9)]#4217

## NEXT: need to assign the dist years to those pixels and an age at time of dist (as above)
## NEXT 2: try again with a bigger window

## now check if the lag1 adds years to the
# lagNames <- c("yearLag", "eventLag", "growth_curve_idLag")
# f9dist[,.(lagNames) := list(shift(year,1,type="lag"), shift(events,1,type="lag"),shift(year,1,growth_curve_id="lag"))]
f9dist[, yearLag := shift(year,1,type="lag")]
f9dist[!is.na(yearLag) & targetPix>0]
#1886
# would I get more with a "lead"
f9dist[, yearLead := shift(year,1,type="lead")]
f9dist[!is.na(yearLead) & targetPix>0] #1897
# 2lag?
f9dist[, yearLag := shift(year,1,type="lag")]
f9dist[!is.na(yearLag) & targetPix>0]



ageAllDist <- merge.data.table(ageDT, allDist, by = "pixelIndex", all.x = TRUE)
# add a column to indicate the pixels I am looking to fill ages
ageAllDist[, targetPix := 0L]
ageAllDist[pixelIndex %in% noDistPixToMatch$pixelIndex]$targetPix <- 1

f9Pix <- ageAllDist[!is.na(f9)]
thesePix4 <- checkF9[!is.na(f9)]

### HERE

qplot(sim$allPixDT[growth_curve_id == rndNumBygcID[1,1] & ages >80,]$ages, geom = "histogram")


library(ggplot2)
distPixWageIn2020 <- qplot(ageNoMatch[pixelIndex %in% allDist$pixelID,]$age2020, geom="histogram")

## check the age distribution of the ageNoMatch for the remaining 89409 pixels
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






