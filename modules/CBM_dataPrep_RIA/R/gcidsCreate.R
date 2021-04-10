gcidsCreate <- function(..., delimiter = "_") {
  do.call(paste, c(list(...)[[1]], sep= delimiter))
}
