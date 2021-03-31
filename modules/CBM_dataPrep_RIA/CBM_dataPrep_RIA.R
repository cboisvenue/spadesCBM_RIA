defineModule(sim, list(
  name = "CBM_dataPrep_RIA",
  description = "A data preparation module to format and prepare user-provided input to the SpaDES forest-carbon modelling familly.",
  keywords = NA,
  authors = c(
    person("Celine", "Boisvenue", email = "Celine.Boisvenue@canada.ca", role = c("aut", "cre"))
  ),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.2", CBM_dataPrep_RIA = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "CBM_dataPrep_RIA.Rmd"),
  reqdPkgs = list(
    "data.table", "fasterize", "magrittr", "raster", "RSQLite", "sf",
    "PredictiveEcology/CBMutils (>= 0.0.6)",
    "PredictiveEcology/LandR@development"
  ),
  parameters = rbind(
    # defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(
      ".plotInitialTime", "numeric", NA, NA, NA,
      "This describes the simulation time at which the first plot event should occur"
    ),
    defineParameter(
      ".plotInterval", "numeric", NA, NA, NA,
      "This describes the simulation time interval between plot events"
    ),
    defineParameter(
      ".saveInitialTime", "numeric", NA, NA, NA,
      "This describes the simulation time at which the first save event should occur"
    ),
    defineParameter(
      ".saveInterval", "numeric", NA, NA, NA,
      "This describes the simulation time interval between save events"
    ),
    defineParameter(
      ".useCache", "logical", FALSE, NA, NA,
      paste(
        "Should this entire module be run with caching activated?",
        "This is generally intended for data-type modules,",
        "where stochasticity and time are not relevant"
      )
    )
  ),
  inputObjects = bindrows(
    expectsInput(
      objectName = "cbmData", objectClass = "dataset",
      desc = "S4 object created from selective reading in of cbm_default.db in CBM_defaults module",
      sourceURL = NA
    ),
    expectsInput(
      objectName = "pooldef", objectClass = "character",
      desc = "Vector of names (characters) for each of the carbon pools, with `Input` being the first one",
      sourceURL = NA
    ),
    expectsInput(
      objectName = "PoolCount", objectClass = "numeric",
      desc = "count of the length of the Vector of names (characters) for each of the carbon pools, with `Input` being the first one",
      sourceURL = NA
    ),
    expectsInput(objectName = "dbPath", objectClass = "character", desc = NA, sourceURL = NA),
    expectsInput(objectName = "sqlDir", objectClass = "character", desc = NA, sourceURL = NA),
    expectsInput(
      objectName = "userDistFile", objectClass = "character", ## TODO: should be a param
      desc = paste("User provided file name that identifies disturbances for simulation",
                   "(key words for searching CBM files, if not there the userDist will be created with defaults"),
      sourceURL = NA
    ),
    ## user input here
    #studyAreaUrl <- "https://drive.google.com/file/d/1LxacDOobTrRUppamkGgVAUFIxNT4iiHU/"
    expectsInput( ## URL RIA CORRECT CHECKED
      objectName = "studyArea", objectClass = "polygon",
      desc = "Limiting the study area to BC",
      sourceURL = "https://drive.google.com/file/d/1LxacDOobTrRUppamkGgVAUFIxNT4iiHU/"
    ),
    expectsInput(
      objectName = "cbmAdmin", objectClass = "dataframe",
      desc = "Provides equivalent between provincial boundaries, CBM-id for provincial boundaries and CBM-spatial unit ids",
      sourceURL = "https://drive.google.com/file/d/1xdQt9JB5KRIw72uaN5m3iOk8e34t9dyz"
    ),

    expectsInput( ## URL RIA CORRECT CHECKED
      objectName = "userDist", objectClass = "data.table",
      desc = "User provided file that identifies disturbances for simulation (distName),
      raster Id if applicable, and wholeStand toggle (1 = whole stand disturbance, 0 = partial disturbance),
      if not there it will use userDistFile",
      sourceURL = "https://drive.google.com/file/d/1Gr_oIfxR11G1ahynZ5LhjVekOIr2uH8X"
    ),
    expectsInput(## URL RIA CORRECT CHECKED
      objectName = "ageRasterURL", objectClass = "character",
      desc = "Pointer to the age raster on a cloud drive",
      sourceURL = "https://drive.google.com/file/d/16pXx9pufiWY45H8rgCiFl6aFbEcC9GeD"
    ),
    expectsInput(
      objectName = "ageRaster", objectClass = "raster",
      desc = "Raster ages for each pixel,ADD URL SPECIFIC TO EACH OF THE RUNS: fireRIs, presentDay, harvest1 and harvest2"
    ),
    expectsInput(
    #   objectName = "gcIndexRasterURL", objectClass = "character", ## TODO: url provided below
    #   desc = "Pointer to URL for spatial link from growth curve to pixels",
    #   sourceURL =  "https://drive.google.com/file/d/1LKblpqSZTVlaNske-C_LzNFEtubn7Ms7"
    # ),
    # expectsInput(## URL RIA CORRECT CHECKED
      objectName = "gcIndexRaster", objectClass = "raster",
      desc = "Raster ages for each pixel",
      sourceURL =  "https://drive.google.com/file/d/1LKblpqSZTVlaNske-C_LzNFEtubn7Ms7"
      #sourceURL = "WILL HAVE TO MAKE THIS FROM GREG'S INFO"
    ),
    expectsInput(## URL RIA CORRECT CHECKED
      objectName = "spuRaster", objectClass = "raster",
      desc = "Raster has spatial units for each pixel"
    ),
    expectsInput(## URL RIA CORRECT CHECKED
      objectName = "ecoRaster", objectClass = "raster",
      desc = "Raster has ecozones for each pixel"
    ),
    expectsInput(
      objectName = "userGcM3File", objectClass = "character",## TODO: should be a param
      desc = paste("User-provided pointer to the file containing: GrowthCurveComponentID,Age,MerchVolume.",
                   "Default name userGcM3"),
      sourceURL = NA
    ),
    expectsInput(## URL RIA CORRECT CHECKED
      objectName = "userGcM3", objectClass = "dataframe",
      desc = "User file containing: GrowthCurveComponentID,Age,MerchVolume. Default name userGcM3",
      sourceURL = "https://drive.google.com/file/d/1BYHhuuhSGIILV1gmoo9sNjAfMaxs7qAj"
    ),
    expectsInput(## URL RIA CORRECT CHECKED
      ## URL RIA FOR scfmFires rasters is this https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8
      ## URL below is for the data table for the 526 years of scfm fire sims
      objectName = "disturbanceRasters", objectClass = "vector",
      desc = "RIA 2020 specific - fires rasters were too big for my laptop. Created a data table for ",
      sourceURL = "https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M"
    ),
    expectsInput( ## URL RIA CORRECT CHECKED
      objectName = "masterRasterURL", objectClass = "character",
      desc = "Pointer to the raster on a cloud drive",
      sourceURL = "https://drive.google.com/file/d/1LlE5NXDPKS6Ljv2SvwQnw2A_Kze2uun-"
    ),
    expectsInput(
      objectName = "masterRaster", objectClass = "raster",
      desc = "Raster that defines the extent of the simulation area for th BC NE RIA. It is used to map results and crop other rasters"#
      #sourceURL = NA
    )
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "pools", objectClass = "matrix", desc = NA),
    createsOutput(
      objectName = "ages", objectClass = "numeric",
      desc = "Ages of the stands from the inventory in 1990"
    ),
    createsOutput(
      objectName = "nStands", objectClass = "numeric",
      desc = "not really the number of stands, but the number of pixel groups"
    ),
    createsOutput(
      objectName = "gcids", objectClass = "numeric",
      desc = "The identification of which growth curves to use on the specific stands provided by..."
    ),
    createsOutput(
      objectName = "historicDMIDs", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating historical disturbance type, linked to the S4 table called cbmData. Only Spinup."
    ),
    createsOutput(
      objectName = "lastPassDMIDS", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating final disturbance type, linked to the S4 table called cbmData. Only Spinup."
    ),
    createsOutput(
      objectName = "delays", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating regeneration delay post disturbance. Only Spinup."
    ),
    createsOutput(
      objectName = "minRotations", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating minimum number of rotations. Only Spinup."
    ),
    createsOutput(
      objectName = "maxRotations", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating maximum number of rotations. Only Spinup."
    ),
    createsOutput(
      objectName = "returnIntervals", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating the fixed fire return interval. Only Spinup."
    ),
    createsOutput(
      objectName = "spatialUnits", objectClass = "numeric",
      desc = "The id given to the intersection of province and ecozones across Canada, linked to the S4 table called cbmData"
    ),
    createsOutput(
      objectName = "ecozones", objectClass = "numeric",
      desc = "Vector, one for each stand, indicating the numeric represenation of the Canadian ecozones, as used in CBM-CFS3"
    ),
    createsOutput(
      objectName = "level3DT", objectClass = "data.table",
      desc = paste("the table linking the spu id, with the disturbance_matrix_id and the events.",
                   "The events are the possible raster values from the disturbance rasters of Wulder and White.")
    ),
    createsOutput(
      objectName = "spatialDT", objectClass = "data.table",
      desc = "the table containing one line per pixel"
    ),
    createsOutput(
      objectName = "mySpuDmids", objectClass = "data.frame",
      desc = "the table containing one line per pixel"
    )
  )
))

doEvent.CBM_dataPrep_RIA <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "CBM_dataPrep_RIA", "save")
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "CBM_dataPrep_RIA", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
      "' in module '", current(sim)[1, "moduleName", with = FALSE], "'",
      sep = ""
    ))
  )
  return(invisible(sim))
}


Init <- function(sim) {
  ## Rasters----------------------------------------------------------------------
  ## user provides raster to match (masterRaster) which is a raster for the
  ## study area, it will define the crs etc, for all other layers. The user also
  ## provides age raster, and a raster linking each growth curve to pixels (gcIndex).
  ## Using the masterRaster, the ecozone raster is made (Canadian ecozones) and the
  ## spatial unit raster. The spatial units are a CBM-CFS3 specific location
  ## that is the intersection of the ecozones and administrative boundaries.
  ## These spatial units (or spu) and the ecozones link the CBM-CFS3 ecological
  ## parameters to the right location (example: decomposition rates).
  ##

  io <- inputObjects(sim, currentModule(sim))
  objectNamesExpected <- io$objectName
  available <- objectNamesExpected %in% ls(sim)

  ## TODO: these aren't required
  omit <- which(objectNamesExpected %in% c("userDistFile", "userGcM3File"))
  available <- available[-omit]
  browser()
  if (any(!available)) {
    stop(
      "The inputObjects for CBM_dataPrep are not all available:",
      "These are missing:", paste(objectNamesExpected[!available], collapse = ", "),
      ". \n\nHave you run ",
      paste0("CBM_", c("defaults"), collapse = ", "),
      "?"
    )
  }

  age <- sim$ageRaster
  gcIndex <- sim$gcIndexRaster
  spuRaster <- sim$spuRaster # made in the .inputObjects
  ecoRaster <- sim$ecoRaster # made in the .inputObjects
  ## End rasters------------------------------------------------------------------
browser()

  ## Create the data table of all pixels and all values for the study area----------------
  level2DT <- data.table(
    spatial_unit_id = spuRaster[], ages = age[], pixelIndex = 1:ncell(age),
    growth_curve_component_id = gcIndex[], growth_curve_id = gcIndex[],
    ecozones = ecoRaster[]
  )
  # keep only the pixels that have all the information: the pixels that will be simulated
  level2DT <- level2DT[!is.na(ages) & !is.na(growth_curve_id)]
  spatialDT <- level2DT
  ## END data.table of all pixels---------------------------------------------------------


  ## Create the pixel groups: groups of pixels with the same attributes ---------------
  spatialDT <- spatialDT[order(pixelIndex), ]
  spatialDT$pixelGroup <- LandR::generatePixelGroups(spatialDT,
    maxPixelGroup = 0,
    columns = c("ages", "spatial_unit_id", "growth_curve_component_id", "ecozones")
  )
  spatialDT <- spatialDT[order(pixelIndex), ]
  spatialDT <- spatialDT[, .(
    ages, spatial_unit_id, pixelIndex,
    growth_curve_component_id, growth_curve_id, ecozones, pixelGroup
  )]
  spatialDT <- spatialDT[order(pixelIndex), ]
  sim$spatialDT <- spatialDT
  # end create pixel groups-------------


  ## Data.table for simulations (one row per pixel group)---------------------
  # this table will be the pixel groups that are used in the spinup procedure in
  # the CBM_core spinup event
  level3DT <- unique(spatialDT[, -("pixelIndex")]) %>% .[order(pixelGroup), ]
  sim$level3DT <- level3DT
  ## End data.table for simulations-------------------------------------------


  ## TODO: problem with ages<=1
  ##################################################### # SK example: can't seem
  #to solve why growth curve id 52 (white birch, good # productivity) will not
  #run with ages= c(0,1,2) it gets stuck in the spinup. Tried ages==1, # and
  #ages==2. Maybe because the first few years of growth are 0 ? (to check) it #
  #does not grow and it does not fill-up the soil pools. # Notes: the GAMs are
  #fit on the cumulative curves of carbon/ha for three # pools. This is to make
  #sure the curves go through 0...but maybe it would # work better for GAMs to
  #be fit on the increments (?). # since all growth curves are for merchantible
  #timber (with diameter limits), it is acceptable to start all increments at
  #the level of year==3.
  #work for this problem for most curves for now: this is from SK runs
  #sim$level3DT[ages==0 & growth_curve_component_id==52,ages:=3]
 ######################################
  ##################### temp fix should

  sim$level3DT[ages <= 1, ages := 3]
  sim$level3DT[order(pixelGroup), ]

  ## Creating all the vectors for the spinup --------------------------------
  sim$ages <- sim$level3DT[, ages]
  sim$nStands <- length(sim$ages)
  sim$pools <- matrix(ncol = sim$PoolCount, nrow = sim$nStands, data = 0)
  colnames(sim$pools) <- sim$pooldef
  sim$pools[, "Input"] <- rep(1.0, nrow(sim$pools))
  sim$gcids <- sim$level3DT[, growth_curve_component_id]
  sim$delays <- rep.int(0, sim$nStands)
  sim$minRotations <- rep.int(10, sim$nStands)
  sim$maxRotations <- rep.int(30, sim$nStands)
  retInt <- merge(sim$level3DT[, ], sim$cbmData@spinupParameters[, c(1, 2)], by = "spatial_unit_id", all.x = TRUE) %>% .[order(pixelGroup)]
  sim$returnIntervals <- retInt[, "return_interval"]
  sim$spatialUnits <- sim$level3DT[, spatial_unit_id]
  sim$ecozones <- sim$level3DT$ecozones

  ################################################################################
  ## matching the disturbances with the Disturbance Matrix IDs in CBM-CFS3 defaults
  ################################################################################
  # Matching disturbances to CBM disturbance matrix id---------------------------------

  ## WBI RIA: skipping the SK building of mySpuDmids.
  ## In this case the user provided complete userDist.


                # # make the disturbance look-up table to the disturbance_matrix_id(s)
                # # making sim$mySpuDmids
                #
                # userDist <- sim$userDist
                #
                # # Most cases will only require fire (wildfire) and a clearcut. There are 426
                # # disturbance matrices identified in the archive of CBM
                # # (sim$cbmData@disturbanceMatrix). Matrices are associated with spatial units
                # # (sim$cbmData@disturbanceMatrixAssociation). User can select any disturbance
                # # they want to represent. Some disturbance matrices are based on data but most
                # # are expert opinion in the CBM-CFS3 archive.
                # # Disturbance Matrices are specific to spatial_units_ids--------------
                # spu <- unique(sim$spatialDT$spatial_unit_id)
                # # what disturbances in those spu(s)?
                # # spuDist() function is in CBMutils package
                # # it lists all the possible disturbances in the CBM-CFS3 archive for that/those
                # # spatial unit with the name of the disturbance in the 3rd colum.
                # listDist <- spuDist(spu, sim$dbPath)
                #
                #
                # ## Example specific for SK (as per Boisvenue et al 2016)
                # # Disturbances are from White and Wulder and provided as yearly rasters
                # # raster values 1 to 5
                # # #C:\Celine\GitHub\CBM_\data\forIan\SK_data\disturbance_Sask\ReadMe.txt
                # # # Fire =  1
                # # # Harvest = 2
                # # # Lcondition = 3
                # # # Road = 4
                # # # Unclass = 5
                # # Whatever number of disturbances identified that will be used in the
                # # simulation, each disturbance has to have one one disturbance matrix id
                # # associated with it.
                # # make mySpuDmids (distNames,rasterId,spatial_unit_id,disturbance_matrix_id)
                # distName <- c(rep(userDist$distName, length(spu)))
                # rasterId <- c(rep(userDist$rasterId, length(spu)))
                # wholeStand <- c(rep(userDist$wholeStand, length(spu)))
                # spatial_unit_id <- c(sort(rep(spu, length(userDist$distName))))
                # mySpuDmids <- data.table(distName, rasterId, spatial_unit_id, wholeStand)
                #
                # dmid <- data.frame(spatial_unit_id = integer(), disturbance_matrix_id = integer())
                #
                # for (i in 1:length(mySpuDmids$distName)) {
                #   ### DANGER HARD CODED FIXES
                #   ## to do: present the user with options that live in listDist for the
                #   ## specific spu or in sim$cbmData@disturbanceMatrix
                #   if (mySpuDmids$distName[i] == "clearcut") {
                #     dmid[i, ] <- cbind(mySpuDmids$spatial_unit_id[i], 409)
                #   } else {
                #     getDist <- listDist[grep(mySpuDmids$distName[i], listDist[, 3], ignore.case = TRUE), 1:2]
                #     getDist <- getDist[getDist$spatial_unit_id == mySpuDmids$spatial_unit_id[i], ]
                #     dmid[i, ] <- getDist[1, ]
                #   }
                # }
                #
                # ## bunch of warnings here...
                # mySpuDmids <- data.table(mySpuDmids, dmid$disturbance_matrix_id)
                # names(mySpuDmids) <- c("distName", "rasterId", "spatial_unit_id", "wholeStand", "disturbance_matrix_id")
                # sim$mySpuDmids <- mySpuDmids
  sim$mySpuDmids <- sim$userDist

  # need to match the historic and last past dist to the spatial unit
  # DECISION: both the last pass and the historic disturbance will be the same
  # for the fire retunr interval runs

  ## TO DO: in Canada historic DMIDs will always be fire, but the last past may
  ## not, it could be harvest. Make this optional and give the user a message
  ## saying these are the defaults.


  mySpuFires <- sim$mySpuDmids[grep("wildfire", sim$mySpuDmids$distName, ignore.case = TRUE), ]

  myFires <- mySpuFires[spatial_unit_id %in% unique(sim$level3DT$spatial_unit_id), ]
  setkey(myFires, spatial_unit_id)
  setkey(sim$level3DT, spatial_unit_id)
  # this is mainly to make them the same length at the number of pixel groups
  histLastDMIDs <- merge(sim$level3DT, myFires)
  sim$historicDMIDs <- histLastDMIDs$disturbance_matrix_id
  ## TO DO: this is where it could be something else then fire
  sim$lastPassDMIDS <- histLastDMIDs$disturbance_matrix_id

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  browser()
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # if we chose to not use the RSQLite library in this module, and extract
  # disturbance matrix id (dmid) from sim$cbmData@disturbanceMatrixAssociation,
  # then $sqlDir and $dbPath are not needed.
  if (!suppliedElsewhere(sim$sqlDir)) {
    sim$sqlDir <- file.path(dPath, "cbm_defaults")
  }
  if (!suppliedElsewhere(sim$dbPath)) {
    sim$dbPath <- file.path(dPath, "cbm_defaults", "cbm_defaults.db")
  }

  if (!suppliedElsewhere(sim$cbmData)) {
    spatialUnitIds <- as.matrix(getTable("spatialUnitIds.sql", sim$dbPath, sim$sqlDir))
    disturbanceMatrix <- as.matrix(getTable("disturbanceMatrix.sql", sim$dbPath, sim$sqlDir))
    sim$cbmData <- new("dataset",
      turnoverRates = as.matrix(getTable("turnoverRates.sql", sim$dbPath, sim$sqlDir)),
      rootParameters = as.matrix(getTable("rootParameters.sql", sim$dbPath, sim$sqlDir)),
      decayParameters = as.matrix(getTable("decayParameters.sql", sim$dbPath, sim$sqlDir)),
      spinupParameters = as.matrix(getTable("spinupParameters.sql", sim$dbPath, sim$sqlDir)),
      climate = as.matrix(getTable("climate.sql", sim$dbPath, sim$sqlDir)),
      spatialUnitIds = spatialUnitIds,
      slowAGtoBGTransferRate = as.matrix(0.006),
      biomassToCarbonRate = as.matrix(0.5),
      stumpParameters = as.matrix(getTable("stumpParameters.sql", sim$dbPath, sim$sqlDir)),
      overmatureDeclineParameters = as.matrix(getTable("overmaturedecline.sql", sim$dbPath, sim$sqlDir)),
      disturbanceMatrix = disturbanceMatrix,
      disturbanceMatrixAssociation = as.matrix(getTable("disturbanceMatrixAssociation.sql", sim$dbPath, sim$sqlDir)),
      disturbanceMatrixValues = as.matrix(getTable("disturbanceMatrixValues.sql", sim$dbPath, sim$sqlDir)),
      landclasses = as.matrix(getTable("landclasses.sql", sim$dbPath, sim$sqlDir)),
      pools = as.matrix(getTable("pools.sql", sim$dbPath, sim$sqlDir)),
      domPools = as.matrix(getTable("domPools.sql", sim$dbPath, sim$sqlDir))
    )
  }
  if (!suppliedElsewhere(sim$pooldef)) {
    sim$pooldef <- CBMutils::.pooldef
    sim$PoolCount <- length(sim$pooldef)
  }

  # user provided data tables------------------------------------------------------

  # 1. growth and yield information
  # userGcM3 and userGcM3File, these files are the m3/ha and age info by growth
  # curve ID, columns should be GrowthCurveComponentID	Age	MerchVolume
  ## TODO add a data manipulation to adjust if the m3 are not given on a yearly basis
  if (!suppliedElsewhere("userGcM3", sim)) {

    if (!suppliedElsewhere("userGcM3File", sim)) {
      sim$userGcM3File <- extractURL("userGcM3")
    }

    sim$userGcM3 <- prepInputs(url = sim$userGcM3File,
                               fun = "data.table::fread",
                               destinationPath = dPath,
                               #purge = 7,
                               filename2 = "curve_points_table.csv")

    # message(
    #   "User has not supplied growth curves (m3 by age or the file name for the growth curves). ",
    #   "The default will be used which is for a region in Saskatchewan."
    # )
    ## RIA 2020 specific
    sim$userGcM3[, V1 := NULL]
    names(sim$userGcM3) <- c("GrowthCurveComponentID", "Age", "MerchVolume")
  }

  # 2. Disturbance information - see disturbance raster below
  # this may be provided by the user, by the defaults or by other modules/family
  # of modules. It is the link between the spatial location of the disturbance
  # (like a raster value) and the disturbance name.
  if (!suppliedElsewhere(sim$userDist)) {
    if (!suppliedElsewhere(sim$userDistFile)) {
      sim$userDistFile <- extractURL("userDist")
    }
    sim$userDist <- prepInputs(url = sim$userDistFile,
                                 fun = "data.table::fread",
                                 destinationPath = dPath,
                                 #purge = 7,
                                 filename2 = "mySpuDmids.csv")
    }

  ## cbmAdmin needed to create spuRaster below (a bit convoluted ## TODO )
  if (!suppliedElsewhere("cbmAdmin", sim)) {
    sim$cbmAdmin <- prepInputs(url = extractURL("cbmAdmin"),
                               fun = "data.table::fread",
                               destinationPath = dPath,
                               #purge = 7,
                               filename2 = "mySpuDmids.csv")
    #fread(file.path(dPath, "cbmAdmin.csv")) ## TODO: use prepInputs with url
  }


  # user provided rasters or spatial information------------------------
  #1. VRI rasters gcID and age
  ## HERE
  browser()
  RIArtm <- prepInputs(url = "https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing")

  RIA_VRIstack <- Cache(prepInputsVRI,VRIurl = "https://drive.google.com/file/d/1LXSX8M46EnsTCM3wGhkiMgqWcqTubC12",
                        dPath = "C:/Celine/github/spadesCBM_RIA/modules/CBM_dataPrep_RIA/data",
                        rasterToMatch = RIArtm
  )

  # # 1. Raster to match (masterRaster). This is the study area.
  #
  # if (!suppliedElsewhere("masterRaster", sim)) {
  #   if (!suppliedElsewhere("masterRasterURL", sim)) {
  #     sim$masterRasterURL <- extractURL("masterRasterURL")
  #     # message(
  #     #   "User has not supplied a masterRaster or a URL for a masterRaster (masterRasterURL object).\n",
  #     #   "masterRaster is going to be read from the default URL given in the inputObjects for ",
  #     #   currentModule(sim)
  #     # )
  #   }
  #
  #   sim$masterRaster <- Cache(
  #     prepInputs,
  #     url = sim$masterRasterURL,
  #     fun = "raster::raster",
  #     destinationPath = dPath
  #   )
  #
  #   #sim$masterRaster[sim$masterRaster == 0] <- NA
  # }
  #
  # # 1.1 studyArea (b/c of a grided problem below when I build the eco and spu rasters)
  # if (!suppliedElsewhere("studyArea", sim)){
  #   sim$studyArea <- prepInputs(url = extractURL("studyArea"),
  #                               destinationPath = dPath) %>%
  #     sf::st_as_sf(.) %>%
  #     .[.$TSA_NUMBER %in% c("40", "08", "41", "24", "16"),] %>%
  #     sf::st_buffer(., 0) %>%
  #     sf::as_Spatial(.) %>%
  #     raster::aggregate(.)
  # }
  #
  # # 2. Age raster from inventory
  #
  # if (!suppliedElsewhere(sim$ageRaster)) {
  #   if (!suppliedElsewhere(sim$ageRasterURL)) {
  #     sim$ageRasterURL <- extractURL("ageRasterURL")
  #   }
  #   sim$ageRaster <- Cache(prepInputs,
  #                          url = sim$ageRasterURL,
  #                          # rasterToMatch = sim$masterRaster,
  #                          fun = "raster::raster",
  #                          destinationPath = dPath
  #   )
  #   ## TODO: put in a message to out pointing out the max age (this has to be
  #   ## sinked to the max age on the growth curve max age for the spinup)
  #   # maxAge <- max(sim$ageRaster)
  #   # message(max age on the raster is XX)
  # }
  #
  # # 3. What growth curve should be applied to what pixels?
  # if (!suppliedElsewhere(sim$gcIndexRaster)) {
  #     sim$gcIndexRaster <- Cache(prepInputs,
  #                              url = extractURL("gcIndexRaster"),
  #                              #rasterToMatch = sim$masterRaster,
  #                              fun = "raster::raster",
  #                              destinationPath = dPath)
  # }
  #
  # # 4. Ecozone raster. This takes the masterRaster (study area) and figures
  # # out what ecozones each pixels are in. This determines some
  # # defaults CBM-parameters across Canada.
  # if (!suppliedElsewhere(sim$ecoRaster)) {
  #   #         ecozones <- prepInputs(
  #   #               # this website https://sis.agr.gc.ca/cansis/index.html is hosted by the Canadian Government
  #   #               url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
  #   #               alsoExtract = "similar",
  #   #               destinationPath = dPath,
  #   #               studyArea = sim$studyArea,
  #   #               useSAcrs = TRUE,
  #   #               overwrite = TRUE,
  #   #               fun = "raster::shapefile",
  #   #               filename2 = TRUE
  #   #               ) %>%
  #   #             cropInputs(., rasterToMatch = sim$masterRaster)
  #   # sim$ecoRaster <- fasterize::fasterize(sf::st_as_sf(ecozones),
  #   #                                       raster = sim$masterRaster,
  #   #                                       field = "ECOZONE"
  #   #                                       )
  #   # }
  #   ### if there are too many momery issues, this raster can be dowloaded here:
  #   ### https://landr-team-group.slack.com/files/UCNUAJ6HK/F01RR6YRR5G/ecozoneraster.tif
  #   sim$ecoRaster <- prepInputs(
  #     url = "https://drive.google.com/file/d/1mnPIi1YfBZ6ej3VnBf_cm5D3nmR0WMBB/view?usp=sharing",
  #     destinationPath = dPath,
  #     targetFile = "ecoRas.tif",
  #     fun = "raster"
  #     #"C:/Celine/github/spadesCBM_RIA/modules/CBM_dataPrep_RIA/data")
  #   )
  # }
  # # 5. Spatial Unit raster. This takes the masterRaster (study area) and figures
  # # out what CBM-specific spatial units each pixels. This determines some
  # # defaults CBM-parameters across Canada.
  # if (!suppliedElsewhere(sim$spuRaster)) {
  #   ## NEW
  #   CanadaAdminras <- prepInputs(url = 'https://drive.google.com/file/d/1mnPIi1YfBZ6ej3VnBf_cm5D3nmR0WMBB/view?usp=sharing',
  #                                 targetFile = "canadaAdminRas.tif",
  #                                 fun = "raster",
  #                                 destinationPath = dPath)
  #
  #   SPUras <- data.table(pixelID = 1:ncell(sim$masterRaster),
  #                        AdminBoundaryID = values(CanadaAdminras),
  #                        ecozone = getValues(sim$ecoRaster))
  #   SPUras <- sim$cbmAdmin[SPUras, on = c("AdminBoundaryID" = "AdminBoundaryID",
  #                                  "EcoBoundaryID" = "ecozone")]
  #
  #   sim$spuRaster <- setValues(sim$masterRaster, SPUras$SpatialUnitID)
  # }
  #   ## NEW to here
  # #   canadaSpu <- prepInputs(targetFile = "spUnit_Locator.shp",
  # #                           url = "https://drive.google.com/file/d/17UH3TDuEA_NQISeavu77HWz_P4tMYb2X",
  # #                           destinationPath = dPath,
  # #                           alsoExtract = "similar")
  # #   spuShp <- postProcess(canadaSpu,
  # #                         rasterToMatch = sim$masterRaster,
  # #                         targetCRS = crs(sim$masterRaster),
  # #                         useCache = FALSE, filename2 = NULL
  # #   )
  # #   sim$spuRaster <- fasterize::fasterize(sf::st_as_sf(spuShp), raster = sim$masterRaster, field = "spu_id")
  # #
  # #   spuShp <- prepInputs(studyArea = sim$studyArea,
  # #                        useSAcrs = TRUE,
  # #                        url = "https://drive.google.com/file/d/17UH3TDuEA_NQISeavu77HWz_P4tMYb2X",
  # #                        destinationPath = dPath)
  # #   spuShp <- spTransform(spuShp, crs(sim$masterRaster))
  # #
  # #   sim$spuRaster <- fasterize::fasterize(sf::st_as_sf(spuShp),
  # #                                         raster = sim$masterRaster, field = "spu_id")
  # # }
  #
  # # # 5. Ecozone raster. This takes the masterRaster (study area) and figures
  # # # out what ecozones each pixels are in. This determines some
  # # # defaults CBM-parameters across Canada.
  # # if (!suppliedElsewhere(sim$ecoRaster)) {
  # #   ecozones <- prepInputs(
  # #     # this website https://sis.agr.gc.ca/cansis/index.html is hosted by the Canadian Government
  # #     url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
  # #     alsoExtract = "similar",
  # #     destinationPath = dPath,
  # #     rasterToMatch = sim$masterRaster,
  # #     overwrite = TRUE,
  # #     fun = "raster::shapefile",
  # #     filename2 = TRUE
  # #   ) %>%
  # #     cropInputs(., rasterToMatch = sim$masterRaster)
  # #   sim$ecoRaster <- fasterize::fasterize(sf::st_as_sf(ecozones),
  # #                                         raster = sim$masterRaster,
  # #                                         field = "ECOZONE"
  # #   )
  # # }
  #
  #
  # # 6. Disturbance rasters. The default example is a list of rasters, one for
  # # each year. But these can be provided by another family of modules in the
  # # annual event.
  ### TO DO: add a message if no information is provided asking the user if
  ### disturbances will be provided on a yearly basis.
  if (!suppliedElsewhere("disturbanceRasters", sim)) {
    ## this case is reading in a sparseDT.
    # RTM that the datatable was created with
    RIA_RTM <- prepInputs(url = 'https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing',
                          destinationPath = dPath) #you probably already have this raster - RIA_RTM.tif on google
    # need the masterRaster (sim$masterRaster)
    # sparseDT
    scfmAnnualBurns <- prepInputs(url = 'https://drive.google.com/file/d/1P41fr5fimmxOTGfNRBgjwXetceW6YS1M/view?usp=sharing',
                                  destinationPath = 'inputs',
                                  overwrite = TRUE,
                                  fun = 'readRDS')

    IndexRTM <- setValues(RIA_RTM, 1:ncell(RIA_RTM))

    #postProcess to match tempTHLB
    IndexTHLB <- setValues(sim$masterRaster, 1:ncell(sim$masterRaster))

    #postProcess the RTM
    IndexRTM <- postProcess(IndexRTM, rasterToMatch = sim$masterRaster)
    #build matching data.table

    indexDT <- data.table(rtmIndex = getValues(IndexRTM),
                          thlbIndex = getValues(IndexTHLB))

    #the NAs in rtmIndex are pixels that are not in THLB (but inside the landscape) - we can remove them
    indexDT <- indexDT[!is.na(rtmIndex)]

    ## the function indexAnnualFire() will have to be used in the annual event
    ## of the CBM_core module



    # options(reproducible.useGDAL = FALSE)
    # sim$disturbanceRasters <- Cache(prepInputs,
    #                    url = "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8",
    #                    fun = "raster::brick",
    #                    rasterToMatch = sim$masterRaster,
    #                    datatype = "INT1U",
    #                    useGDAL = FALSE)

      ## TODO this is a brick, with 526 rasters that low RAM system cannot
      ## handle. Instead a sparseDT was create and is read-in here. We need to
      ## add here the capacity to deal with a bunch of rasters in a folder,
      ## (like in the SK runs), a stack of rasters (below), raster brick (above)
      ## or polygons? OR a sparseDT
    # stack
      # sim$disturbanceRasters <- Cache(prepInputs,
      #                               url = "https://drive.google.com/file/d/1fJIPVMyDu66CopA-YP-xSdP2Zx1Ll_q8",
      #                               fun = "raster::stack",
      #                               rasterToMatch = masterRaster,
      #                               useGDAL = FALSE)
    # rasters in a folder
      # distHere <- extractURL(disturbanceRasters)
      # sim$disturbanceRasters <- list.files(distHere,full.names = TRUE) %>%
      #   grep(., pattern = ".grd$", value = TRUE)
      # # if all fails
    # or just one raster? or polygons?
  }

  #stopifnot(length(sim$disturbanceRasters) > 0)

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}
