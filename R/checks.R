# Author: Roozbeh Valavi
# contact: valavi.r@gmail.com
# Date : May 2023
# Version 0.2
# Licence GPL v3

# check for x
.check_x <- function(x, name = "x"){
  if(!methods::is(x, "sf")){
    tryCatch(
      {
        x <- sf::st_as_sf(x)
      },
      error = function(cond) {
        message(sprintf("'%s' is not convertible to an sf object!", name))
        message(sprintf("'%s' must be an sf or spatial* object.", name))
      }
    )
  }
  return(x)
}

# check for column matching colnames(x)
.check_column <- function(column, x){
  if(!is.null(column)){
    if(!column %in% colnames(x)){
      warning(sprintf("There is no column named '%s' in 'x'. Column is ignored!\n", column))
      column <- NULL
    }
  }
  return(column)
}

# column should be binary or categorical
.check_classes <- function(clen, column, th = 15){
  if(clen > th){
    warning(
      sprintf(
        paste(
          "The are too many unique values in '%s'.",
          "Use 'column' only for binary or categorical responses (ignore this if it is).\n"
        ),
        column
      )
    )
  }
}

# check for r
.check_r <- function(r, name = "r"){
  if(!methods::is(r, "SpatRaster")){
    tryCatch(
      {
        r <- terra::rast(r)
      },
      error = function(cond) {
        message(sprintf("'%s' is not convertible to a terra SpatRaster object!", name))
        message(sprintf("'%s' must be a SpatRaster, stars, Raster* object, or path to a raster file on disk.", name))
      }
    )
  }
  return(r)
}

# check for required packages
.check_pkgs <- function(pkg){
  pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
  if(length(pkgna) > 0){
    nm <- paste(pkgna, collapse = ", ")
    message("This function requires these packages: ", nm, "\nWould you like to install them now?\n1: yes\n2: no")
    user <- readline(prompt = paste0("Selection: "))
    if(tolower(user) %in% c("1", "yes", "y")){
      utils::install.packages(pkgna)
    } else{
      stop("Please install these packages for function to work: ", nm)
    }
  }
}
