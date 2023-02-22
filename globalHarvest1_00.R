# creating a based example for spadesCBM_RIA using the harvest1 scenario (see
# RIA paper). This example will be where we start the trnasition to changing the
# C++ core (in CBM_core) to libcbm.
# CBoisvenue
# Feb. 22 2023

## PFC work-around
## this is a work-around for working from PFC...R cannot connect to URL

options("download.file.method" = "wininet")

## This gives me the development branche of SpaDES.core
install.packages("SpaDES.core", repos = "https://predictiveecology.r-universe.dev")

if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("PredictiveEcology/Require@development")
library(Require)

Require("magrittr")
Require("SpaDES.core")

## pulling out all the harvest1 scenario objects from spadesCBMglobal.Rmd needed
## to run a sim
options("reproducible.useRequire" = TRUE)

cacheDir <- reproducible::checkPath("cache", create = TRUE)
moduleDir <- reproducible::checkPath("modules")
inputDir <- reproducible::checkPath("inputs", create = TRUE)
outputDir <- reproducible::checkPath("outputs", create = TRUE)
scratchDir <- file.path(tempdir(), "scratch", "CBM") %>% reproducible::checkPath(create = TRUE)

timesHarvest1 <- list(start = 2020.00, end = 2099.00)

parameters <- list(
  CBM_defaults = list(
    .useCache = TRUE
  ),
  CBM_vol2biomass_RIA = list(
    .useCache = TRUE
  )
)

# harvest1 is the base case
parametersHarvest1 <- parameters
parametersHarvest1$CBM_dataPrep_RIAharvest1 <- list(
  .useCache = TRUE
)

parametersHarvest1$CBM_core <- list(
  #.useCache = TRUE, #"init", #c(".inputObjects", "init")
  # .plotInterval = 5,
  .plotInitialTime = timesHarvest1$start,
  poolsToPlot = c("totalCarbon"),
  spinupDebug = FALSE) ## TODO: temporary

modulesHarvest1 <- list("CBM_defaults", "CBM_dataPrep_RIAharvest1", "CBM_vol2biomass_RIA", "CBM_core")

objects <- list(
  dbPath = file.path(inputDir, "cbm_defaults", "cbm_defaults.db"),
  sqlDir = file.path(inputDir, "cbm_defaults")
)

paths <- list(
  cachePath = cacheDir,
  modulePath = moduleDir,
  inputPath = inputDir,
  rasterPath = scratchDir
)

pathsHarvest1 <- paths
pathsHarvest1$outputPath <- file.path(outputDir,"harvest1")

quickPlot::dev.useRSGD(FALSE)
dev()
clearPlot()
options(spades.moduleCodeChecks = FALSE,
        reproducible.useMemoise = FALSE,
        spades.recoveryMode = FALSE)

whensHarvest1 <- sort(c(timesHarvest1$start, timesHarvest1$start + c(10, 30, 50, 60),
                        timesHarvest1$end - 1:0))
outputsHarvest1 <- as.data.frame(expand.grid(objectName = c("cbmPools", "NPP"), saveTime = whensHarvest1))

opts <- options("spades.useRequire" = FALSE)
options(opts)

##TODO CBMutils does not seem to load - I am connected to the development branch
##of CBMutils
library("devtools")
devtools::load_all("C:/Celine/github/CBMutils")

RIAharvest1Runs <- simInitAndSpades(times = timesHarvest1,
                                    params = parametersHarvest1,
                                    modules = modulesHarvest1,
                                    objects = objects,
                                    paths = pathsHarvest1,
                                    outputs = outputsHarvest1,
                                    loadOrder = unlist(modulesHarvest1),
                                    debug = TRUE)

















##########################NOT USED YET##########################
while (!require("SpaDES.project")) {
  install.packages("SpaDES.project", repos = "https://predictiveecology.r-universe.dev")
  require(SpaDES.project)
}

##TODO current set-up creates input/inputs and input/outputs folder...need to
##figure out which we keep and make sure the "extras" are not created. The
##Spades.project call creates the ones without the "s" at the end but thorughout
##the code the ones with the "s" are used...but this might be only in the
##specific module folders. Modules .Rmd use "outputs" but no code seems to use it.
## Inputs with an "s" seems to be only used here in CBM_core (lines 305-317):
# spinup <- function(sim) {
#   io <- inputObjects(sim, currentModule(sim))
#   objectNamesExpected <- io$objectName
#   available <- objectNamesExpected %in% ls(sim)
#   if (any(!available)) {
#     stop(
#       "The inputObjects for CBM_core are not all available:",
#       "These are missing:", paste(objectNamesExpected[!available], collapse = ", "),
#       ". \n\nHave you run ",
#       paste0("spadesCBM", c("defaults", "inputs", "m3ToBiomass", "userDist"), collapse = ", "),
#       "?"
#     )
#   }
## Need to see if we can change that to "input".

update.packages("SpaDES.core")

