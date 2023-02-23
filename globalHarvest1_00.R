# creating a based example for spadesCBM_RIA using the harvest1 scenario (see
# RIA paper). This example will be where we start the trnasition to changing the
# C++ core (in CBM_core) to libcbm.
# CBoisvenue
# Feb. 22 2023

# project basic setup -------------------------------------------------------------------------

if (file.exists("~/.Renviron")) readRenviron("~/.Renviron") ## GITHUB_PAT
if (file.exists("spadesCBM_RIA.Renviron")) readRenviron("spadesCBM_RIA.Renviron") ## database credentials

.ncores <- min(parallel::detectCores() / 2, 32L) ## default number of CPU cores to use, e.g. for pkg install
.nodename <- Sys.info()[["nodename"]] ## current computer name; used to configure machine-specific settings
.user <- Sys.info()[["user"]] ## current computer username; used to configure user-specific settings

## define project directory - this code expects it is being run from this location
## **do not change the paths defined here**
## if you need to add a machine- or user-specific path, please do so _conditionally_
prjDir <- switch(.user,
                 cboisven = "C:/Celine/github/spadesCBM_RIA",
                 "~/GitHub/spadesCBM_RIA")

## ensure script being run from the project directory
stopifnot(identical(normalizePath(prjDir), normalizePath(getwd())))

options(
  Ncpus = .ncores,
  repos = c(PE = "https://predictiveecology.r-universe.dev",
            CRAN = "https://cloud.r-project.org"),
  Require.RPackageCache = "default" ## will use default package cache directory: `RequirePkgCacheDir()`
)

## work-around for working from PFC...R cannot connect to certain urls
## TODO: improve conditional by only using wininet if *at* PFC, not just on a PFC machine
if ((.Platform$OS.type == "windows") && grepl("[L|W]-VIC", .nodename)) {
  options("download.file.method" = "wininet")
}

# install and load packages -------------------------------------------------------------------

## use project-specific location for packages to avoid conflicts with other projects
pkgDir <- file.path(tools::R_user_dir(basename(prjDir), "data"), "packages",
                    version$platform, getRversion()[, 1:2])
dir.create(pkgDir, recursive = TRUE, showWarnings = FALSE)
.libPaths(pkgDir, include.site = FALSE)
message("Using libPaths:\n", paste(.libPaths(), collapse = "\n"))

## package installation only; do not load module packages until after install
if (!"remotes" %in% rownames(installed.packages(lib.loc = .libPaths()[1]))) {
  install.packages("remotes")
}

if (!"Require" %in% rownames(installed.packages(lib.loc = .libPaths()[1])) ||
    packageVersion("Require", lib.loc = .libPaths()[1]) < "0.2.6") {
  remotes::install_github("PredictiveEcology/Require@development")
}

library(Require)

# setLinuxBinaryRepo() ## setup binary package installation for linux users; currently interferes with remotes installation

## installs development branch of SpaDES.core from https://predictiveecology.r-universe.dev
pkgsToInstall <- c("googledrive", "magrittr", "SpaDES.core")
#Install(pkgsToInstall), upgrade = FALSE, standAlone = TRUE) ## TODO: fails with spatial pkgs + quickPlot
if (!all(pkgsToInstall %in% rownames(installed.packages(lib.loc = .libPaths()[1])))) {
  install.packages(pkgsToInstall)
}

if (.user == "cboisven") {
  ## TODO CBMutils does not seem to load - I am connected to the development branch of CBMutils
  devtools::load_all("C:/Celine/github/CBMutils")
} else {
  remotes::install_github("PredictiveEcology/CBMutils@development", repos = "https://cloud.r-project.org", upgrade = FALSE) ## TODO: Require fails to install
  Require("PredictiveEcology/CBMutils@development", dependencies = TRUE, standAlone = TRUE, upgrade = FALSE)
}

library(magrittr) ## TODO: remove and use base pipe
library(SpaDES.core)

# setting up harvest1 scenarios ---------------------------------------------------------------

## pulling out all the harvest1 scenario objects from spadesCBMglobal.Rmd needed
## to run a sim
options(reproducible.useMemoise = FALSE,
        reproducible.useTerra = TRUE,
        spades.moduleCodeChecks = FALSE,
        spades.recoveryMode = FALSE,
        spades.useRequire = TRUE)

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
pathsHarvest1$outputPath <- file.path(outputDir, "harvest1")

quickPlot::dev.useRSGD(FALSE)
#dev()
clearPlot()

whensHarvest1 <- sort(c(timesHarvest1$start, timesHarvest1$start + c(10, 30, 50, 60),
                        timesHarvest1$end - 1:0))
outputsHarvest1 <- as.data.frame(expand.grid(objectName = c("cbmPools", "NPP"), saveTime = whensHarvest1))

## TODO: need to automatically get cbm_defaults data
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

