defineModule(sim, list(
  name = "CBM_vol2biomass_RIA",
  description = paste("A module to prepare the user-provided growth and yield information for use",
                      "in the family of models spadesCBM - CBM-CFS3-like simulation of forest",
                      "carbon in the platform SpaDES. This module takes in user-provided m3/ha",
                      "and meta data for teh growth curves and returns annual increments for",
                      "the aboveground live c-pools."),
  keywords = "",
  authors = c(
    person("Celine", "Boisvenue", email = "Celine.Boisvenue@canada.ca", role = c("aut", "cre"))
  ),
  childModules = character(0),
  version = list(CBM_vol2biomass_RIA = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.txt", "CBM_vol2biomass_RIA.Rmd")),
  reqdPkgs = list(
    "ggplot2", "ggpubr", "mgcv", "quickPlot", "PredictiveEcology/CBMutils (>= 0.0.6)",
    "robustbase", "ggforce"
  ),
  parameters = rbind(
    # defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(
      ".plotInitialTime", "numeric", NA, NA, NA,
      "Describes the simulation time at which the first plot event should occur."
    ),
    defineParameter(
      ".plotInterval", "numeric", NA, NA, NA,
      "Describes the simulation time interval between plot events."
    ),
    defineParameter(
      ".saveInitialTime", "numeric", NA, NA, NA,
      "Describes the simulation time at which the first save event should occur."
    ),
    defineParameter(
      ".saveInterval", "numeric", NA, NA, NA,
      "This describes the simulation time interval between save events."
    ),
    defineParameter(
      ".useCache", "logical", FALSE, NA, NA,
      paste(
        "Should this entire module be run with caching activated?",
        "This is generally intended for data-type modules, where stochasticity",
        "and time are not relevant"
      )
    )
  ),
  inputObjects = bindrows(
    # expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    # this are variables in inputed data.tables:SpatialUnitID, EcoBoundaryID, juris_id, ecozone, jur, eco, name, GrowthCurveComponentID, plotsRawCumulativeBiomass, checkInc
    expectsInput(objectName = "curveID", objectClass = "character",
                 desc = "Vector of column names that together, uniquely define growth curve id"),
    expectsInput(
      objectName = "table3",
      objectClass = "dataframe",
      desc = "Stem wood biomass model parameters for merchantable-sized trees from Boudewyn et al 2007",
      sourceURL = "https://nfi.nfis.org/resources/biomass_models/appendix2_table3.csv"
    ),
    expectsInput(
      objectName = "table4", objectClass = "dataframe", desc = "Stem wood biomass model parameters for nonmerchantable-sized trees from Boudewyn et al 2007",
      sourceURL = "https://nfi.nfis.org/resources/biomass_models/appendix2_table4.csv"
    ),
    expectsInput(
      objectName = "table5", objectClass = "dataframe", desc = "Stem wood biomass model parameters for sapling-sized trees from Boudewyn et al 2007",
      sourceURL = "https://nfi.nfis.org/resources/biomass_models/appendix2_table5.csv"
    ),
    expectsInput(
      objectName = "table6", objectClass = "dataframe", desc = "Proportion model parameters from Boudewyn et al 2007",
      sourceURL = "https://nfi.nfis.org/resources/biomass_models/appendix2_table6.csv"
    ),
    expectsInput(
      objectName = "table7", objectClass = "dataframe", desc = "Caps on proportion models from Boudewyn et al 2007",
      sourceURL = "https://nfi.nfis.org/resources/biomass_models/appendix2_table7.csv"
    ),
    expectsInput(
      objectName = "TSAspuEco", objectClass = "dataframe",
      desc = "Links TSA, ecozone and spatial unit ids",## TODO: URL doesn't seem to work
      sourceURL = "https://drive.google.com/file/d/1_iRjcK7PbbFjCZWzpOp-jf0rIXmaMHtB"
    ),
    expectsInput(
      objectName = "cbmAdmin", objectClass = "dataframe",
      desc = "Provides equivalent between provincial boundaries, CBM-id for provincial boundaries and CBM-spatial unit ids",
      sourceURL = "https://drive.google.com/file/d/1NwVPqMho8fOcf9mVVzm5SlfDbbyQmQlh/view?usp=sharing"
    ),
    expectsInput(objectName = "gcMeta",
                 objectClass = "dataframe",
                 desc = "Provides equivalent between provincial boundaries,
                 CBM-id for provincial boundaries and CBM-spatial unit ids",
                 sourceURL = NA),
    ## RIA URL LINK CORRECT (au_table.csv)
    expectsInput(objectName = "gcMetaFile",
                 objectClass = "character",
                 desc = "File name and location for the user provided gcMeta dataframe",
                 sourceURL = "https://drive.google.com/file/d/1YmQ6sNucpEmF8gYkRMocPoeKt2P26ZiX"
                 ),
    expectsInput(objectName = "canfi_species",
                 objectClass = "dataframe",
                 desc = "File containing the possible species in the Boudewyn table - note
                 that if Boudewyn et al added species, this should be updated. Also note that such an update is very unlikely",
                 sourceURL = "https://drive.google.com/file/d/1l9b9V7czTZdiCIFX3dsvAsKpQxmN-Epo"),
    expectsInput(
      objectName = "userGcM3File", objectClass = "character",## TODO: should be a param
      desc = paste("Pointer to the user file name for the files containing: GrowthCurveComponentID,Age,MerchVolume.",
                   "Default name userGcM3"),
      sourceURL = NA
    ),
    expectsInput(## URL RIA CORRECT CHECKED
      objectName = "userGcM3", objectClass = "dataframe",
      desc = "User file containing: GrowthCurveComponentID,Age,MerchVolume. Default name userGcM3",
      sourceURL = "https://drive.google.com/file/d/1BYHhuuhSGIILV1gmoo9sNjAfMaxs7qAj"
    ),
    expectsInput(objectName = "ecozones", objectClass = "data.table", desc = "the table linking the spu id, with the
                  disturbance_matrix_id and the events. The events are the possible raster values from the disturbance rasters of Wulder and White"),
    expectsInput(objectName = "gcids", objectClass = "data.table", desc = "the table linking the spu id, with the
                  disturbance_matrix_id and the events. The events are the possible raster values from the disturbance rasters of Wulder and White"),
    expectsInput(objectName = "spatialUnits", objectClass = "data.table", desc = "the table linking the spu id, with the
                  disturbance_matrix_id and the events. The events are the possible raster values from the disturbance rasters of Wulder and White")
  ),
  outputObjects = bindrows(
    # createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = NA, objectClass = NA, desc = NA),
    createsOutput(objectName = "volCurves", objectClass = "plot", desc = "Plot of all the growth curve provided by the user"),
    createsOutput(objectName = "plotsRawCumulativeBiomass", objectClass = "plot", desc = "Plot of cumulative m3/ha curves
                  translated into tonnes of carbon/ha, per AG pool, prior to any smoothing"),
    createsOutput(objectName = "checkInc", objectClass = "plot", desc = "Plot of 1/2 of the increment per AG pool,
                  calculated from the smoothed cumulative tonnes c/ha, derived into increments, per AG pool. "),
    createsOutput(objectName = "growth_increments", objectClass = "matrix", desc = "Matrix of the 1/2 increment that will be used to create the gcHash"),
    createsOutput(objectName = "gcHash", objectClass = "environment", desc = "Environment pointing to each gcID, that is itself an environment,
                  pointing to each year of growth for all AG pools.Hashed matrix of the 1/2 growth increment.
                  This is used in the c++ functions to increment AG pools two times in an annual event (in the spadesCBMcore.R module.")
  )
))

## event types
#   - type `init` is required for initialization

doEvent.CBM_vol2biomass_RIA <- function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "CBM_vol2biomass_RIA", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "CBM_vol2biomass_RIA", "save")
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # plotFun(sim) # uncomment this, replace with object to plot
      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "CBM_vol2biomass_RIA", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "CBM_vol2biomass_RIA", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    event1 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "CBM_vol2biomass_RIA", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
    },
    event2 = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "CBM_vol2biomass_RIA", "templateEvent")

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
      "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'",
      sep = ""
    ))
  )
  return(invisible(sim))
}

## event functions

Init <- function(sim) {
  # user provides userGcM3: incoming cumulative m3/ha
  # plot
  # Test for steps of 1 in the yield curves
  ageJumps <- sim$userGcM3[, list(jumps = unique(diff(as.numeric(Age)))), by = "GrowthCurveComponentID"]
  idsWithJumpGT1 <- ageJumps[jumps > 1]$GrowthCurveComponentID
  if (length(idsWithJumpGT1)) {
    missingAboveMin <- sim$userGcM3[, approx(Age, MerchVolume, xout = setdiff(seq(0, max(Age)), Age)),
                          by = "GrowthCurveComponentID"]
    setnames(missingAboveMin, c("x", "y"), c("Age", "MerchVolume"))
    missingAboveMin <- na.omit(missingAboveMin)
    sim$userGcM3 <- rbindlist(list(sim$userGcM3, missingAboveMin))
    setorderv(sim$userGcM3, c("GrowthCurveComponentID", "Age"))

    # Assertion
    ageJumps <- sim$userGcM3[, list(jumps = unique(diff(as.numeric(Age)))), by = "GrowthCurveComponentID"]
    idsWithJumpGT1 <- ageJumps[jumps > 1]$GrowthCurveComponentID
    if (length(idsWithJumpGT1) > 0)
      stop("There are still yield curves that are not annually resolved")
  }


  sim$volCurves <- ggplot(data = sim$userGcM3, aes(x = Age, y = MerchVolume, group = GrowthCurveComponentID, colour = GrowthCurveComponentID)) +
    geom_line() ## TODO: move to plotInit event
  message("User: please look at the curve you provided via sim$volCurves")

  ## not all curves provided are used in the simulation - and ***FOR NOW*** each
  ## pixels only gets assigned one growth curve (no transition, no change in
  ## productivity).
  ## To run module independently, the gcID used in this translation can be specified here
  # if(!suppliedElsewhere("level3DT",sim)){
  #   userGcM3 <- sim$userGcM3
  # }else{
  userGcM3 <- sim$userGcM3#[GrowthCurveComponentID %in% unique(sim$gcids), ]
  # }
  # START reducing Biomass model parameter tables -----------------------------------------------
  # reducing the parameter tables to the jurisdiction or ecozone we have in the study area
  ## To run module independently, the gcID used in this translation can be specified here
  # if(!suppliedElsewhere("spatialUnits",sim)){
  #   spu  <- ### USER TO PROVIDE SPU FOR EACH gcID###########
  # }else{

  spu <- unique(sim$spatialUnits)
  # }
  # if(!suppliedElsewhere("ecozones",sim)){
  #   eco <- ### USER TO PROVIDE SPU FOR EACH gcID###########
  # }else{

  eco <- unique(sim$ecozones)
  # }
  thisAdmin <- sim$cbmAdmin[sim$cbmAdmin$SpatialUnitID %in% spu & sim$cbmAdmin$EcoBoundaryID %in% eco, ]

  # "s" table for small table3, 4, 5, 6, 7 - tables limited to the targeted
  # ecozones and jurisdictions
  stable3 <- as.data.table(sim$table3[sim$table3$juris_id %in% thisAdmin$abreviation &
    sim$table3$ecozone %in% eco, ])
  stable4 <- as.data.table(sim$table4[sim$table4$juris_id %in% thisAdmin$abreviation &
    sim$table4$ecozone %in% eco, ])
  # table5 is different since there was not have enough data to fit models for
  # all provinces. Here we are hard-coding the closest equivalent province to
  # have a complete set.
  # This first If-statement is to catch the "no-province" match


  stable5.2 <- as.data.table(sim$table5[sim$table5$juris_id %in% thisAdmin$abreviation, ])
  if (!length(unique(stable5.2$juris_id)) == length(unique(thisAdmin$abreviation))) {
    ## DANGER HARD CODED: if NFIS changes table 5, this will no longer be valid
    # juris_id: there are only 5/13 possible
    # these are the provinces available: AB BC NB NF NT
    # for the non match these would be the equivalent
    # "PE" - NB
    # "QC" - NB
    # "ON" - NB
    # "MB" - AB
    # "SK" - AB
    # "YK" - NT
    # "NU" - NT
    abreviation <- c("PE", "QC", "ON", "MB", "SK", "YK", "NU")
    t5abreviation <- c("NB", "NB", "NB", "AB", "AB", "NT", "NT")
    abreviaReplace <- data.table(abreviation, t5abreviation)
    # replace the abbreviations and select
    thisAdmin5 <- merge(abreviaReplace, thisAdmin)
    thisAdmin5[, c("abreviation", "t5abreviation") := list(t5abreviation, NULL)]
    stable5.2 <- as.data.table(sim$table5[sim$table5$juris_id %in% thisAdmin5$abreviation, ])
  }
  # This second "if-statement" is to catch is the "no-ecozone" match
  ### THIS NEEDS TO BE TESTED
  if (nrow(stable5.2) > 0) {
    stable5 <- stable5.2[ecozone %in% unique(eco), ]
  } else {
    stop(
      "There are no matches found for the parameters needed to execute the Boudewyn models.",
      "Please manually find matches for table 5."
    )
  }
  if (!length(eco) == length(unique(stable5$ecozone))) {
    # there are 9/15 ecozones
    # These are the ones in table5
    # id               name
    # 4       Taiga Plains
    # 5  Taiga Shield West
    # 6 Boreal Shield West
    # 7  Atlantic Maritime
    # 9      Boreal Plains
    # 10  Subhumid Prairies
    # 12  Boreal Cordillera
    # 13   Pacific Maritime
    # 14 Montane Cordillera

    # these are the ones that are not
    # id               name
    # 8   Mixedwood Plains  - 7  Atlantic Maritime
    # 11   Taiga Cordillera - 4 taiga plains
    # 15      Hudson Plains - 6 Boreal Shield West
    # 16  Taiga Shield East - 5  Taiga Shield West
    # 17 Boreal Shield East - 6 Boreal Shield West
    # 18  Semiarid Prairies - 10  Subhumid Prairies

    EcoBoundaryID <- c(8, 11, 15, 16, 17, 18)
    ecoNotInT5 <- c(7, 4, 6, 5, 6, 10)
    ecoReplace <- data.table(ecoNotInT5, EcoBoundaryID)
    thisAdmin5.1 <- merge(ecoReplace, thisAdmin5, by = EcoBoundaryID)
    stable5 <- as.data.table(stable5[stable5$ecozone %in% thisAdmin5.1$EcoBoundaryID, ])
  }
  if (nrow(stable5) < 1) {
    stop("There is a problem finding a parameter match in table 5.")
  }

  stable6 <- as.data.table(sim$table6[sim$table6$juris_id %in% thisAdmin$abreviation &
    sim$table6$ecozone %in% eco, ])
  stable7 <- as.data.table(sim$table7[sim$table7$juris_id %in% thisAdmin$abreviation &
    sim$table6$ecozone %in% eco, ])
  # END reducing Biomass model parameter tables -----------------------------------------------

  # Read-in user provided meta data for growth curves. This could be a complete
  # data frame with the same columns as gcMetaEg.csv OR is could be only curve
  # id and species. This format is necessary to process the curves and use the
  # resulting increments
  gcMeta <- sim$gcMeta
  ##TODO: RIA and SK don't need the same processing of this file. Need to
  ##generalize this process
   #   # help the user go from their growth curve id and leading species to the six
  #   # columns: names(gcMeta)
  #   # [1] "growth_curve_id"           "growth_curve_component_id"
  #   # [3] "species"                   "canfi_species"
  #   # [5] "genus"                     "forest_type_id"
  #   # the data frame canfi_species.csv (in userData_Defaults_spadesCBM -
  #   # https://drive.google.com/drive/folders/1OBDTp1v_3b3D3Yvg1pHLiW6A_GRklgpD?usp=sharing)
  #   # has all the possible options for canfi_species (number), genus (4 letter
  #   # code) and species (three letter code).

  ## this was for SK
  # # checking how many columns in gcMeta, if not 6, columns need to be added
  # if (!ncol(gcMeta) == 6) {
  #   gcMeta2 <- gcMeta[, .(growth_curve_id, canfi_species)]
  #   gcMeta2[, growth_curve_component_id := growth_curve_id]
  #   # check if all the species are in the canfi_species table
  #   ### THIS HAS NOT BEEN TESTED YET
  #   if (nrow(gcMeta2) == length(which(gcMeta$species %in% sim$canfi_species$name))) {
  #     spsMatch <- sim$canfi_species[
  #       , which(name %in% gcMeta2$species),
  #       .(canfi_species, genus, name, forest_type_id)
  #     ]
  #     names(spsMatch) <- c("canfi_species", "genus", "species", "forest_type_id")
  #     setkey(gcMeta2, species)
  #     setkey(spsMatch, species)
  #     gcMeta3 <- merge(gcMeta2, spsMatch) # I do not think the order of the columns matter
  #     gcMeta <- gcMeta3
  #   }
  #   ### PUT SOMETHING HERE IF THE SPECIES DONT MATCH...NOT SURE WHAT - ERROR MESSAGE?
  # }

  ## This is for the RIA fire return interval runsL using unmanged curves (VDYP)

  riaGcMeta <- gcMeta[,.(au_id, tsa, canfi_species, unmanaged_curve_id)]
  # this will be slightly different when au_id are not eqaul to
  # unmanaged_curve_id. For VDYP, those are equal.
  setnames(riaGcMeta,
           c("au_id", "tsa", "unmanaged_curve_id"),
             c("growth_curve_component_id", "TSAid","growth_curve_id"))


  # assuming gcMeta has now 6 columns, it needs a 7th: spatial_unit_id. This
  # will be used in the convertM3biom() fnct to link to the right ecozone
  # and it only needs the gc we are using in this sim.
  #spatialDt <- fread("C:/Celine/github/spadesCBM_RIA/modules/CBM_dataPrep_RIA/data/startingSpatialDT.csv")
  gcThisSim <- unique(sim$spatialDT[,.(growth_curve_component_id, spatial_unit_id, ecozones)])#as.data.table(unique(cbind(sim$spatialUnits, sim$gcids)))
  #names(gcThisSim) <- c("growth_curve_component_id","e")
  setkey(gcThisSim, growth_curve_component_id)
  setkey(riaGcMeta, growth_curve_component_id) #changed from gcMeta to riaGcMeta
  gcMeta <- merge(riaGcMeta, gcThisSim) #changed from gcMeta to riaGcMeta # adds ecozone


  browser()
  curveID <- sim$curveID
  if (!is.null(sim$level3DT)) {
    gcidsLevels <- levels(sim$level3DT$gcids)
    gcids <- factor(gcidsCreate(gcMeta[, ..curveID]), levels = gcidsLevels)
  } else {
    gcids <- factor(gcidsCreate(gcMeta[, ..curveID]))
  }


  set(gcMeta, NULL, "gcids", gcids)

  ### TODO CHECK - this in not tested
  if (!length(unique(unique(userGcM3$GrowthCurveComponentID)) ==
              length(unique(gcMeta$growth_curve_component_id)))) {
    stop("There is a missmatch in the growth curves of the userGcM3 and the gcMeta")
  }
# RIA still missing columns in gcMeta: species genus and forest_type_id
  gcMeta <- merge.data.table(gcMeta, sim$canfi_species, by = "canfi_species", all.x = TRUE)
  gcMeta[, species := NULL]
  setnames(gcMeta, "name", "species")

  ################
  warning("Modifying canfi_species 1211 ecozone to 1203")
  gcMeta[canfi_species == 1211, canfi_species := 1203]


  # START processing curves from m3/ha to tonnes of C/ha then to annual increments
  # per above ground biomass pools -------------------------------------------

  ### NEED TO MAKE SURE THE PROVIDED CURVES ARE ANNUAL
  ### if not, we need to extrapolate to make them annual

  # Matching is 1st on species, then on gcId which gives us location (admin,
  # spatial unit and ecozone)
  fullSpecies <- unique(gcMeta$species) ## RIA: change this to the canfi_sps or match??

  cumPools <- Cache(cumPoolsCreate, fullSpecies, gcMeta, userGcM3,
                             stable3, stable4, stable5, stable6, stable7, thisAdmin)

  cbmAboveGroundPoolColNames <- "totMerch|fol|other"
  colNames <- grep(cbmAboveGroundPoolColNames, colnames(cumPools), value = TRUE)

  # {
  #   # matching on species name
  #   speciesMeta <- gcMeta[species == fullSpecies[i], ]
  #   # for each species name, process one gcID at a time
  #   for (j in 1:NROW(unique(speciesMeta, on = c("growth_curve_component_id", "ecozones")))) {
  #
  #     meta <- speciesMeta[j, ]
  #     ecozone <- meta$ecozones
  #     id <- userGcM3$GrowthCurveComponentID[which(userGcM3$GrowthCurveComponentID == meta$growth_curve_component_id)][-1]
  #     ## IMPORTANT BOURDEWYN PARAMETERS FOR NOT HANDLE AGE 0 ##
  #     age <- userGcM3[GrowthCurveComponentID == meta$growth_curve_component_id, Age]
  #     age <- age[which(age>0)]
  #     # series of fncts results in curves of merch, foliage and other (SW or HW)
  #
  #     cumBiom <- as.matrix(convertM3biom(
  #       meta = meta, gCvalues = userGcM3, spsMatch = gcMeta,
  #       ecozones = thisAdmin, params3 = unique(stable3), params4 = unique(stable4),
  #       params5 = unique(stable5), params6 = unique(stable6), params7 = unique(stable7)
  #     ))
  #
  #     # going from tonnes of biomass/ha to tonnes of carbon/ha here
  #     cumBiom <- cumBiom * 0.5 ## this value is in sim$cbmData@biomassToCarbonRate
  #     # calculating the increments per year for each of the three pools (merch,
  #     # foliage and other (SW or HW))
  #     # inc <- diff(cumBiom)
  #     # CBM processes half the growth before turnover and OvermatureDecline, and
  #     # half after.
  #     # names(outInputs$allProcesses)
  #     # [1] "Disturbance"       "Growth1"           "DomTurnover"       "BioTurnover"
  #     # [5] "OvermatureDecline" "Growth2"           "DomDecay"          "SlowDecay"
  #     # [9] "SlowMixing"
  #     id_ecozone <- paste0(id, "_", ecozone)
  #     cumBiom <- cbind(id, age, cumBiom, ecozone, id_ecozone)
  #
  #     cumPools <- rbind(cumPools, cumBiom)
  #   }
  # }

  # Check models that are directly out of the Boudewyn-translation----------------------------
  # Usually, these need to be, at a minimum, smoothed out.

  # plotting the curves of the direct translation ------------------------
  # adding the zeros back in
  # cumPools <- as.data.table(cumPools)
  # cumPools[, numAge := as.numeric(age)]
  browser()
  minAgeId <- cumPools[,.(minAge = max(0, min(age) - 1)), by = "gcids"]

  fill0s <- minAgeId[,.(age = seq(from = 0, to = minAge, by = 1)), by = "gcids"]
  # might not need this
  length0s <- fill0s[,.(toMinAge = length(age)), by = "gcids"]

  # these are going to be 0s
  carbonVars <- data.table(gcids = unique(fill0s$gcids),
                                    totMerch = 0,
                                    fol = 0,
                                    other = 0 )

  fiveOf7cols <- fill0s[carbonVars, on = "gcids"]

  otherVars <- cumPools[,.(id = unique(id), ecozone = unique(ecozone)), by = "gcids"]
  add0s <- fiveOf7cols[otherVars, on = "gcids"]
  # cumPools[,numAge := NULL]

  cumPoolsRaw <- rbind(cumPools,add0s)
  set(cumPoolsRaw, NULL, "age", as.numeric(cumPoolsRaw$age))
  setorderv(cumPoolsRaw, c("gcids", "age"))

  figPath <- file.path(modulePath(sim), currentModule(sim), "figures")
  # plotting and save the plots of the raw-translation in the sim$
  if (!is.na(P(sim)$.plotInitialTime))
  sim$plotsRawCumulativeBiomass <- Cache(m3ToBiomPlots, inc = cumPoolsRaw,
                                         path = figPath,
                                         filenameBase = "rawCumBiomass_")

  # Fixing of non-smooth curves
  cumPoolsClean <- Cache(cumPoolsSmooth, cumPoolsRaw)

  # a[, totMerch := totMerchNew]
  if (!is.na(P(sim)$.plotInitialTime))
    figs <- Cache(m3ToBiomPlots, inc = cumPoolsClean,
                path = figPath,
                filenameBase = "cumPools_smoothed_postChapmanRichards")

  set(cumPoolsClean, NULL, colNames, NULL)
  colNamesNew <- grep(cbmAboveGroundPoolColNames, colnames(cumPoolsClean), value = TRUE)
  setnames(cumPoolsClean, old = colNamesNew, new = colNames)


  # Calculating Increments
  incCols <- c("incMerch", "incFol", "incOther")
  cumPoolsClean[, (incCols) := lapply(.SD, function(x) c(NA, diff(x))), .SDcols = colNames,
                by = eval("gcids")]
  colsToUse33 <- c("age", "gcids", incCols)
  if (!is.na(P(sim)$.plotInitialTime))
    rawIncPlots <- Cache(m3ToBiomPlots, inc = cumPoolsClean[, ..colsToUse33],
                       path = figPath,
                       title = "Smoothed increments merch fol other by gc id",
                       filenameBase = "Increments")




  # From: http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/

  # do.call(ggarrange, rawPlots)
  # sim$plotsRawCumulativeBiomass <- Cache(do.call,
  #   ggarrange,
  #   append(
  #     rawPlots,
  #     list(
  #       common.legend = TRUE,
  #       legend = "right",
  #       labels = names(rawPlots),
  #       font.label = list(size = 7, color = "black", face = "bold"),
  #       label.x = 0.5
  #     )
  #   )
  # )
  # # dev.new()
  # annotate_figure(rawPlots, #sim$plotsRawCumulativeBiomass,
  #                 top = text_grob("Cumulative merch fol other by gc id", face = "bold", size = 14)
  # )
  message(crayon::red("User: please inspect figures of the raw and smoothed translation of your growth curves in: ",
          figPath))

  # sim$gg
  # end plotting direct translations----------------------------------------

  # Example of smoothing: two choices presented, smoothing of the cumulative
  # carbon/ha curves or smoothing of the increment curves of carbon/ha

  # ####### Smoothing of start here#######################################################
  # ## CHOICE 1. on cumulative curves ###################################################
  # # created smooth cumulative curves using GAMs
  # sMerch <- NULL
  # sFol <- NULL
  # sOther <- NULL
  #
  # ## TODO: module not loading the package mgcv even though it is in reqdPkgs
  # #library(mgcv)
  #
  #
  # id <- unique(cumPoolsRaw$id)
  # for (val in id) {
  #   # per column
  #   # totMerch
  #   oneSet <- cumPoolsRaw[id == val, ]
  #   wts <- c(
  #     100, rep(2, (which(oneSet$totMerch == max(oneSet$totMerch)) - 1)),
  #     rep(1, (length(oneSet$age) - which(oneSet$totMerch == max(oneSet$totMerch))))
  #   )
  #   k <- 20
  #   gamMerch <- gam(oneSet$totMerch ~ s(oneSet$age, k = k), weight = wts, method = "REML")
  #   df1 <- as.data.frame(cbind(age = oneSet$age, totMerch = gamMerch$fitted.values))
  #   # User: this would be a check of the estimates from Boudewyn with the fitted GAM values
  #   # ggplot(oneSet, aes(age, totMerch)) + geom_point() +
  #   #  geom_line(data=df1,aes(color="Fitted GAM cumulative totMerch"))
  #   sMerch <- rbind(sMerch, df1)
  #   # fol
  #   # oneSet <- cbind(cumPoolsRaw[id==id[i],.(age,fol)])
  #   wts <- c(
  #     100, rep(2, (which(oneSet$fol == max(oneSet$fol)) - 1)),
  #     rep(1, (length(oneSet$age) - which(oneSet$fol == max(oneSet$fol))))
  #   )
  #   k <- 20
  #   gamMerch <- gam(oneSet$fol ~ s(oneSet$age, k = k), weight = wts, method = "REML")
  #   df2 <- as.data.frame(cbind(age = oneSet$age, fol = gamMerch$fitted.values))
  #   # User: this would be a check of the estimates from Boudewyn with the fitted GAM values
  #   # ggplot(oneSet, aes(age, fol)) + geom_point() +
  #   #   geom_line(data=df2,aes(color="Fitted GAM cumulative fol"))
  #   sFol <- rbind(sFol, df2)
  #   # other
  #   # oneSet <- cbind(cumPoolsRaw[id==id[i],.(age,other)])
  #   wts <- c(
  #     100, rep(2, (which(oneSet$other == max(oneSet$other)) - 1)),
  #     rep(1, (length(oneSet$age) - which(oneSet$other == max(oneSet$other))))
  #   )
  #   k <- 20
  #   gamMerch <- gam(oneSet$other ~ s(oneSet$age, k = k), weight = wts, method = "REML")
  #   df3 <- as.data.frame(cbind(age = oneSet$age, other = gamMerch$fitted.values))
  #   # User: this would be a check of the estimates from Boudewyn with the fitted GAM values
  #   # ggplot(oneSet, aes(age, other)) + geom_point() +
  #   #   geom_line(data=df3,aes(color="Fitted GAM cumulative other"))
  #   sOther <- rbind(sOther, df3)
  # }
  # # putting it all together
  # smoothCumPools <- as.data.table(cbind(cumPoolsRaw$id, sMerch, sFol[, 2], sOther[, 2]))
  # names(smoothCumPools) <- names(cumPoolsRaw)
  # setkey(smoothCumPools, id)
  # # half the increments are use at the begining of simulations and half later in
  # # the simulation. The order is:
  # # names(sim$allProcesses)
  # # [1] "Disturbance"       "Growth1"           "DomTurnover"
  # # [4] "BioTurnover"       "OvermatureDecline" "Growth2"
  # # [7] "DomDecay"          "SlowDecay"         "SlowMixing"
  # # here, calculating the increment, then dividing it in 2.
  # cols <- c("totMerch", "fol", "other")
  # newCols <- c("lMerch", "lFol", "lOther")
  # lagPools <- smoothCumPools[, (newCols) := lapply(.SD, data.table::shift), .SDcols = cols, by = "id"]
  # lagPools[is.na(lagPools)] <- 0
  # incCols <- c("incMerch", "incFol", "incOther")
  # incPools <- lagPools[, (incCols) := list(totMerch - lMerch, fol - lFol, other - lOther)][, .(id, age, incMerch, incFol, incOther)]
  #
  # forestType <- gcMeta[, .(id = growth_curve_component_id, forest_type_id)]
  # #       #FYI:
  # #       # cbmTables$forest_type
  # #       # id           name
  # #       # 1  1       Softwood
  # #       # 2  2      Mixedwood
  # #       # 3  3       Hardwood
  # #       # 4  9 Not Applicable
  # setkey(forestType, id)
  # incPools <- merge(incPools, forestType)
  # swCols <- c("swmerch", "swfol", "swother")
  # hwCols <- c("hwmerch", "hwfol", "hwother")
  #
  # totalIncrements <- incPools[forest_type_id == 1, (swCols) := list((incMerch), (incFol), (incOther))][forest_type_id == 3, (hwCols) := list((incMerch), (incFol), (incOther))]
  # totalIncrements[is.na(totalIncrements)] <- 0
  # outCols <- c("incMerch", "incFol", "incOther", "forest_type_id")
  # incCols <- c(swCols, hwCols)
  # totalIncrements[, (outCols) := NULL]
  # increments <- totalIncrements[, (incCols) := list(
  #   swmerch / 2, swfol / 2,
  #   swother / 2, hwmerch / 2, hwfol / 2, hwother / 2
  # )] %>% .[order(id, age), ]
  # incPlots <- m3ToBiomPlots(inc = increments)
  #
  # # From: http://www.sthda.com/english/articles/32-r-graphics-essentials/126-combine-multiple-ggplots-in-one-graph/
  # # do.call(ggarrange, rawPlots)
  # sim$checkInc <- do.call(
  #   ggarrange,
  #   append(
  #     incPlots,
  #     list(
  #       common.legend = TRUE,
  #       legend = "right",
  #       labels = names(rawPlots),
  #       font.label = list(size = 10, color = "black", face = "bold"),
  #       label.x = 0.5
  #     )
  #   )
  # )
  # # dev.new()
  # annotate_figure(sim$checkInc,
  #   top = text_grob("Halved increments for merch fol other by gc id", face = "bold", size = 14)
  # )
  # message("User: please inspect the halved increments that are used in your simulation via sim$checkInc.")

  ## CHOICE 2.on increment curves ###################################################
  # setkey(cumPoolsRaw, id)
  # # half the increments are use at the begining of simulations and half later in
  # # the simulation. The order is:
  # # names(sim$allProcesses)
  # # [1] "Disturbance"       "Growth1"           "DomTurnover"
  # # [4] "BioTurnover"       "OvermatureDecline" "Growth2"
  # # [7] "DomDecay"          "SlowDecay"         "SlowMixing"
  # # here, calculating the increment, then dividing it in 2.
  # cols <- c("totMerch", "fol", "other")
  # newCols <- c("lMerch", "lFol", "lOther")
  # lagPools <- cumPoolsRaw[, (newCols) := lapply(.SD, data.table::shift), .SDcols = cols, by = "id"]
  # lagPools[is.na(lagPools)] <- 0
  # incCols <- c("incMerch", "incFol", "incOther")
  # incPools <- lagPools[, (incCols) := list(totMerch - lMerch, fol - lFol, other - lOther)][, .(id, age, incMerch, incFol, incOther)]
  # rawIncPlots <- m3ToBiomPlots(inc = incPools)
  # plotsRawInc <- do.call(
  #   ggarrange,
  #   append(
  #     rawIncPlots,
  #     list(
  #       common.legend = TRUE,
  #       legend = "right",
  #       labels = names(rawIncPlots),
  #       font.label = list(size = 10, color = "black", face = "bold"),
  #       label.x = 0.5
  #     )
  #   )
  # )
  # # dev.new()
  # annotate_figure(plotsRawInc,
  #                 top = text_grob("RawIncrements merch fol other by gc id", face = "bold", size = 14)
  # )
  # sim$plotsRawInc <- plotsRawInc
  # ## calculate the GAMs on increments
  # # created smooth cumulative curves using GAMs
  # sIncMerch <- NULL
  # sIncFol <- NULL
  # sIncOther <- NULL
  #
  # ## TODO: module not loading the package mgcv even though it is in reqdPkgs
  # library(mgcv)
  #
  #
  # id <- unique(incPools$id)
  # for (val in id) {
  #   # per column
  #   # incMerch
  #   oneSet <- incPools[id == val, ]
  #   wts <- c(
  #     100, rep(2, (which(oneSet$incMerch == max(oneSet$incMerch)) - 1)),
  #     rep(1, (length(oneSet$age) - which(oneSet$incMerch == max(oneSet$incMerch))))
  #   )
  #   k <- 20
  #   gamMerch <- gam(oneSet$incMerch ~ s(oneSet$age, k = k), weight = wts, method = "REML")
  #   df1 <- as.data.frame(cbind(age = oneSet$age, incMerch = gamMerch$fitted.values))
  #   # User: this would be a check of the estimates from Boudewyn with the fitted GAM values
  #   # ggplot(oneSet, aes(age, incMerch)) + geom_point() +
  #   #  geom_line(data=df1,aes(color="Fitted GAM cumulative incMerch"))
  #   sIncMerch <- rbind(sIncMerch, df1)
  #   # incFol
  #   # oneSet <- cbind(cumPoolsRaw[id==id[i],.(age,incFol)])
  #   wts <- c(
  #     100, rep(2, (which(oneSet$incFol == max(oneSet$incFol)) - 1)),
  #     rep(1, (length(oneSet$age) - which(oneSet$incFol == max(oneSet$incFol))))
  #   )
  #   k <- 20
  #   gamMerch <- gam(oneSet$incFol ~ s(oneSet$age, k = k), weight = wts, method = "REML")
  #   df2 <- as.data.frame(cbind(age = oneSet$age, incFol = gamMerch$fitted.values))
  #   # User: this would be a check of the estimates from Boudewyn with the fitted GAM values
  #   # ggplot(oneSet, aes(age, incFol)) + geom_point() +
  #   #   geom_line(data=df2,aes(color="Fitted GAM cumulative incFol"))
  #   sIncFol <- rbind(sIncFol, df2)
  #   # other
  #   # oneSet <- cbind(cumPoolsRaw[id==id[i],.(age,incOther)])
  #   wts <- c(
  #     100, rep(2, (which(oneSet$incOther == max(oneSet$incOther)) - 1)),
  #     rep(1, (length(oneSet$age) - which(oneSet$incOther == max(oneSet$incOther))))
  #   )
  #   k <- 20
  #   gamMerch <- gam(oneSet$incOther ~ s(oneSet$age, k = k), weight = wts, method = "REML")
  #   df3 <- as.data.frame(cbind(age = oneSet$age, incOther = gamMerch$fitted.values))
  #   # User: this would be a check of the estimates from Boudewyn with the fitted GAM values
  #   # ggplot(oneSet, aes(age, incOther)) + geom_point() +
  #   #   geom_line(data=df3,aes(color="Fitted GAM cumulative incOther"))
  #   sIncOther <- rbind(sIncOther, df3)
  # }
  # # putting it all together
  #
  # smoothIncPools <- as.data.table(cbind(incPools$id, sIncMerch, sIncFol[, 2], sIncOther[, 2]))
  # names(smoothIncPools) <- names(incPools)
  # smoothIncPlots <- m3ToBiomPlots(inc = smoothIncPools )
  # # plot the smooth increments pools - any better?
  # plotsSmoothIncs <- do.call(
  #   ggarrange,
  #   append(
  #     smoothIncPlots,
  #     list(
  #       common.legend = TRUE,
  #       legend = "right",
  #       labels = names(smoothIncPlots),
  #       font.label = list(size = 10, color = "black", face = "bold"),
  #       label.x = 0.5
  #     )
  #   )
  # )
  # # dev.new()
  # annotate_figure(plotsSmoothIncs,
  #                 top = text_grob("RawIncrements merch fol other by gc id", face = "bold", size = 14)
  # )


  # sim$smoothIncPools <- smoothIncPools
  # sim$plotsSmoothIncs <- plotsSmoothIncs

  sim$cumPoolsClean <- cumPoolsClean # TODO: add to metadata

  colsToUseForestType <- c("growth_curve_component_id", "forest_type_id", "gcids")
  forestType <- gcMeta[, ..colsToUseForestType]
  #       #FYI:
  #       # cbmTables$forest_type
  #       # id           name
  #       # 1  1       Softwood
  #       # 2  2      Mixedwood
  #       # 3  3       Hardwood
  #       # 4  9 Not Applicable

  setkeyv(forestType, "gcids")
  cumPoolsClean <- merge(cumPoolsClean, forestType, by = "gcids")
  swCols <- c("swmerch", "swfol", "swother")
  hwCols <- c("hwmerch", "hwfol", "hwother")

  totalIncrementsSmooth <- cumPoolsClean[forest_type_id == 1, (swCols) := list((incMerch), (incFol), (incOther))]
  totalIncrementsSmooth <- totalIncrementsSmooth[forest_type_id == 3, (hwCols) := list((incMerch), (incFol), (incOther))]
  totalIncrementsSmooth[is.na(totalIncrementsSmooth)] <- 0
  outCols <- c("incMerch", "incFol", "incOther", "forest_type_id")
  incCols <- c(swCols, hwCols)
  totalIncrementsSmooth[, (outCols) := NULL]
  increments <- totalIncrementsSmooth[, (incCols) := list(
    swmerch / 2, swfol / 2,
    swother / 2, hwmerch / 2, hwfol / 2, hwother / 2
  )]
  setorderv(increments, c("gcids", "age"))
  ## TODO : clean-up this - have to remove the ecozone and id and make
  ## everything numeric to be able to hash the matrix. This will create problems
  ## when we relate the id to the id from vol curves since this id is a
  ## combination of ecozone and id with no space (to be able to make it
  ## numeric).
  incColKeep <- c("id", "age", incCols)
  set(increments, NULL, "id", as.numeric(increments[["gcids"]]))
  set(increments, NULL, setdiff(colnames(increments), incColKeep), NULL)
  setcolorder(increments, incColKeep)

  # Assertions
  if (isTRUE(P(sim)$doAssertions)) {
    # All should have same min age
    if (length(unique(increments[, min(age), by = "id"]$V1)) != 1)
      stop("All ages should start at the same age for each curveID")
    if (length(unique(increments[, max(age), by = "id"]$V1)) != 1)
      stop("All ages should end at the same age for each curveID")
  }

  sim$growth_increments <- as.matrix(increments)
  # END process growth curves -------------------------------------------------------------------------------

  sim$gcHash <- matrixHash(sim$growth_increments)
  # create a nested hash (by gcid/by age)
  ## used in SpinUp function later...
  for (item in ls(sim$gcHash)) {
    sim$gcHash[[item]] <- hash(sim$gcHash[[item]])
  }

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  # Plot(sim$object)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
Event1 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
Event2 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  # cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  cacheTags <- c(currentModule(sim), "function:.inputObjects")

  if (!suppliedElsewhere("curveID", sim)) {
    sim$curveID <- c("growth_curve_component_id", "ecozones")
  }

  # 1. tables from Boudewyn
  # these are all downloaded from the NFIS site. The NFIS however, changes the
  # tables and seems to forget parameter coumns at times. So, we added a check to
  if (!suppliedElsewhere("table3", sim)) {
    sim$table3 <- prepInputs(url = extractURL("table3"),
                             fun = "data.table::fread",
                             destinationPath = dPath,
                             overwrite = FALSE,
                             #purge = 7,
                             filename2 = "appendix2_table3.csv")
    # these are the columns needed in the functions for calculating biomass
    t3hasToHave <- c("juris_id", "ecozone", "canfi_species", "genus", "species",
                     "a", "b", "volm")
    if(!length(which(colnames(sim$table3) %in% t3hasToHave)) == length(t3hasToHave)){
      message(
        "The parameter table (appendix2_table3) does not have the expected number of columns. ",
        "This means parameters are missing. The default (older) parameter file will be used instead."
      )
      sim$table3 <- fread(file.path("modules","CBM_vol2biomass_RIA","data","appendix2_table3.csv"))
    }
  }
  if (!suppliedElsewhere("table4", sim)) {
    sim$table4 <- prepInputs(url = extractURL("table4"),
                             fun = "data.table::fread",
                             destinationPath = dPath,
                             overwrite = FALSE,
                             #purge = 7,
                             filename2 = "appendix2_table4.csv")
     t4hasToHave <- c("juris_id", "ecozone", "canfi_species", "genus", "species",
                     "a", "b", "k", "cap", "volm")
    if(!length(which(colnames(sim$table4) %in% t4hasToHave)) == length(t4hasToHave)){
      message(
        "The parameter table (appendix2_table4) does not have the expected number of columns. ",
        "This means parameters are missing. The default (older) parameter file will be used instead."
      )
    }
  }
  if (!suppliedElsewhere("table5", sim)) {
    sim$table5 <- prepInputs(url = extractURL("table5"),
                             fun = "data.table::fread",
                             destinationPath = dPath,
                             overwrite = FALSE,
                             #purge = 7,
                             filename2 = "appendix2_table5.csv")

    t5hasToHave <- c("juris_id", "ecozone", "canfi_genus", "genus", "a", "b", "k",
                     "cap", "volm")
    if(!length(which(colnames(sim$table5) %in% t5hasToHave)) == length(t5hasToHave)){
      message(
        "The parameter table (appendix2_table5) does not have the expected number of columns. ",
        "This means parameters are missing. The default (older) parameter file will be used instead."
      )
    }
  }
  if (!suppliedElsewhere("table6", sim)) {
    sim$table6 <- prepInputs(url = extractURL("table6"),
                             fun = "data.table::fread",
                             destinationPath = dPath,
                             overwrite = FALSE,
                             #purge = 7,
                             filename2 = "appendix2_table6_v2.csv")
#    sim$table6 <- fread("https://nfi.nfis.org/resources/biomass_models/appendix2_table6.csv",
#                        fill = TRUE)
    t6hasToHave <- c("juris_id", "ecozone", "canfi_species", "a1", "a2", "a3", "b1", "b2", "b3",
                     "c1", "c2", "c3" )
    if(!length(which(colnames(sim$table6) %in% t6hasToHave)) == length(t6hasToHave)){
      message(
        "The parameter table (appendix2_table6) does not have the expected number of columns. ",
        "This means parameters are missing. The default (older) parameter file will be used instead."
      )
 #     sim$table6 <- fread(file.path("modules","CBM_vol2biomass_RIA","data","appendix2_table6.csv"))
      ## TODO: use url to googleDrive (G:\userDataDefaultsCBM_SK)
    }
  }
  if (!suppliedElsewhere("table7", sim)) {
    sim$table7 <- prepInputs(url = extractURL("table7"),
                             fun = "data.table::fread",
                             destinationPath = dPath,
                             overwrite = FALSE,
                             #purge = 7,
                             filename2 = "appendix2_table7.csv")

    t7hasToHave <- c("juris_id", "ecozone", "canfi_species", "vol_min", "vol_max", "p_sw_low",
                     "p_sb_low", "p_br_low", "p_fl_low", "p_sw_high", "p_sb_high", "p_br_high",
                     "p_fl_high")
    if(!length(which(colnames(sim$table7) %in% t7hasToHave)) == length(t7hasToHave)){
      message(
        "The parameter table (appendix2_table7) does not have the expected number of columns. ",
        "This means parameters are missing. The default (older) parameter file will be used instead."
      )
    }
  }

  if (!suppliedElsewhere("gcids", sim)) {
    ## this is where the pixelGroups and their spu eco etc.
    message("No spatial information was provided for the growth curves.
            The default values (RIA NE BC) will be used to limit the number of growth curves used.")
    sim$gcids <- c(801000, 801001, 801002, 802000, 802001, 802002, 801003, 803000,
                   801008, 803002, 802008, 801006, 802006, 802003, 803003, 803001,
                   803008, 803006, 803004, 802004, 801004, 802007, 801005, 803007,
                   802005, 801007, 803005, 4002000, 4001000, 4003000, 4001002, 4001001,
                   1601005, 4003005, 4002003, 4001003, 4002005, 4003001, 4002001,
                   1601011, 4001005, 1603005, 1603007, 1601004, 1602007, 1601002,
                   4002002, 1602002, 1602004, 1601007, 1603004, 1602009, 1603009,
                   4003002, 4003003, 1603011, 1602011, 1603002, 1602005, 1601009,
                   1601001, 1601000, 1602001, 1602000, 1603000, 1601006, 1602006,
                   1603001, 1603006, 4001006, 4002006, 4003006, 1601010, 1603010,
                   2401000, 2401003, 2401005, 2401002, 2402000, 2402003, 2403003,
                   2402001, 2402004, 2403001, 2401001, 2401006, 1602010, 2401004,
                   2403004, 4002004, 4001004, 4003004, 1601003, 2403000, 1602003,
                   2402007, 1601008, 1602008, 2401007, 2402005, 2403005, 1603003,
                   1603008, 2402002, 2403006, 2402006, 2403002, 4103005, 4101000,
                   4102000, 4103003, 4103006, 4101003, 4101006, 4102003, 4102005,
                   4101005, 4101008, 4103000, 4101009, 4102007, 4101007, 4103008,
                   4102009, 4103007, 4103009, 4102002, 4103002, 4101002, 4102006,
                   2403007, 4103001, 4102008, 4102001, 4101001, 4102004, 4103004,
                   4101004)
  }

  if (!suppliedElsewhere("ecozones", sim)) {
    message("No spatial information was provided for the growth curves.
            The default values (RIA NE BC) will be used to determine which ecozones these curves are in.")
    sim$ecozones <- c(4, 12, 9, 14)
  }
  if (!suppliedElsewhere("spatialUnits", sim)) {
    message("No spatial information was provided for the growth curves.
            The default values (RIA NE BC) will be used to determine which CBM-spatial units these curves are in.")
    sim$spatialUnits <- c(38, 40, 39, 42)
  }

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

  # 2. meta info about growth and yield curves
  if (!suppliedElsewhere("gcMeta", sim)) {

    if (!suppliedElsewhere("gcMetaFile", sim)) {
      sim$gcMetaFile <- extractURL("gcMetaFile")
      }
    sim$gcMeta <- prepInputs(url = sim$gcMetaFile,
                                 fun = "data.table::fread",
                                 destinationPath = dPath,
                                 #purge = 7,
                                 filename2 = "au_table.csv")

      # message(
      #   "User has not supplied growth curves (m3 by age or the file name for the growth curves). ",
      #   "The default will be used which is for a region in Saskatchewan."
      # )
  }
  # 3. RIA: Match between TSA, spu, and eco, needed to run this independently of
  # CBM_dataPrep_RIA.r
  if (!suppliedElsewhere("TSAspuEco", sim)) {

    sim$TSAspuEco <- prepInputs(url = extractURL("TSAspuEco"),
                             fun = "data.table::fread",
                             destinationPath = dPath,
                             #purge = 7,
                             filename2 = "TSAspuEco.csv")
  }

  # 4. cbmAdmin: this is needed to match species and parameters. Boudewyn et al 2007
  # abbreviation and cbm spatial units and ecoBoundnary id is provided with the
  # adminName to avoid confusion.
  if (!suppliedElsewhere("cbmAdmin", sim)) {
    sim$cbmAdmin <-  prepInputs(url = extractURL("cbmAdmin"),
                                fun = "data.table::fread",
                                destinationPath = dPath,
                                #purge = 7,
                                filename2 = "cbmAdmin.csv")
  }

  # 5. canfi_species: for the BOudewyn parameters, the species have to be matched
  # with the ones in the Boudewyn tables. The choices HAVE to be one of these.
  # This contains three columns, canfi_species, genus and species form the
  # publication and I added (manually) one more: forest_type_id. That variable is a CBM-CFS3
  # indicator as follows:
  # cbmTables$forest_type
  # id           name
  # 1  1       Softwood
  # 2  2      Mixedwood
  # 3  3       Hardwood
  # 4  9 Not Applicable
  if (!suppliedElsewhere("canfi_species", sim)) {
    sim$canfi_species <- prepInputs(url = extractURL("canfi_species"),
                                    fun = "data.table::fread",
                                    destinationPath = dPath,
                                    #purge = 7,
                                    filename2 = "canfi_species.csv")
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
