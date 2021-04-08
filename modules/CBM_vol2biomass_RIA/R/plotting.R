m3ToBiomIncOnlyPlots <- function(inc=increments, id_col = "id_ecozone", nrow = 5, ncol = 5,
                                 filenameBase = "rawCumBiomass_"){
  gInc <- copy(inc)
  colsOut <- c("age","id", "ecozone")
  gInc[ ,(colsOut) := list(NULL,NULL,NULL)]
 # gInc[,c(id_col, "ageNumeric", "age")]
  gc <- data.table::melt(gInc, id.vars = c(id_col, "ageNumeric"))
  set(gc, NULL, "valueNumeric", as.numeric(gc$value))

  idSim <- unique(gc[,..id_col])[[1]]
  plots <- gc %>% # [id_ecozone %in% idSim[1:20]] %>%
    ggplot( aes(x=ageNumeric, y=valueNumeric, group=variable, color=variable)) +
    geom_line() +
    facet_wrap(facets = id_col) +
    labs(title= "Cumulative merch fol other by gc id") +
    theme(plot.title = element_text(hjust = 0.5))

  # Do first page, so that n_pages can calculate how many pages there are
  #   -- not used -- so a bit wasteful
  numPages <- ceiling(length(idSim) / (nrow * ncol))

  for (i in seq(numPages)) {
    plotsByPage <- plots + facet_wrap_paginate(facets = id_col,
                                               page = i, nrow = nrow, ncol = ncol)
    ggsave(paste0(filenameBase,i,".png"), plotsByPage)
  }
  return(plots)
}

