globalVariables(c(
 'year', 'pixels'
))
#' @param rasterTomatch a template raster matching the one used to create scfmAnnualBurns
#' @param scfmAnnualBurns a data.table with an index of burned pixels
#' @param time a year in the scfmAnnualBurns table
#'
#' @return a raster layer with 1 representing a burned pixel
#'
#' @importFrom data.table data.table
#' @importFrom raster raster
rstCurrentBurnFromDT <- function(rasterToMatch, scfmAnnualBurns, time) {
  rstCurrentBurn <- raster(rasterToMatch)
  thisYearsFire <- scfmAnnualBurns[year == time]
  rstCurrentBurn[thisYearsFire$pixels]  <- 1
  return(rstCurrentBurn)
}

# scfmDT <- readRDS("outputs/scfmAnnualBurn_2015-2540.rds")
# rtm <- raster("modules/WBI_dataPrep_studyArea/data/RIA_rtm.tif")
# annualFires <- lapply(2015:2540,
#                       FUN = rstCurrentBurnFromDT,
#                       rasterToMatch = rtm,
#                       scfmAnnualBurns = scfmDT)
# annualFires <- brick(annualFires)
# writeRaster(annualFires, 'outputs/annualFires525yrs.tif', datatype = 'INT1U')
#
