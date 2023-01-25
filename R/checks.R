# Author: Roozbeh Valavi
# contact: valavi.r@gmail.com
# Date : January 2023
# Version 0.1
# Licence GPL v3

# check for x
.x_check <- function(x, name = "x"){
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
.column_check <- function(column, x){
  if(!is.null(column)){
    if(!column %in% colnames(x)){
      warning(sprintf("There is no column named '%s' in 'x'.\n", column))
      column <- NULL
    }
  }
  return(column)
}

# check for r
.r_check <- function(r, name = "r"){
  if(!methods::is(r, "SpatRaster")){
    tryCatch(
      {
        r <- terra::rast(r)
      },
      error = function(cond) {
        message(sprintf("'%s' is not convertible to a terra SpatRaster object!", name))
        message(sprintf("'%s' must be a SpatRaster, stars, Raster* object, or path to raster a file on disk.", name))
      }
    )
  }
  return(r)
}

# check for required packages
.pkg_check <- function(pkg){
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