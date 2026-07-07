# Author: Roozbeh Valavi
# contact: valavi.r@gmail.com
# Date : May 2023
# Version 0.3

# check the object is a blockCV object
.check_cv <- function(x) {
    is_cv <- inherits(x, c("cv_spatial", "cv_cluster", "cv_buffer", "cv_nndm", "cv_knndm"))
    if (!is_cv) stop("'cv' must be a blockCV cv_* object.")
}


# number of original sample points represented by a cv object (its folds index
# into these points). For leave-one-out objects the folds do not each span the
# whole data, so the distinct referenced indices are counted instead.
.cv_n_points <- function(cv){
    if(.is_loo(cv)){
        length(unique(unlist(cv$folds_list)))
    } else{
        length(unlist(cv$folds_list[[1]]))
    }
}


# the supplied sample data must match the points used to build the cv object,
# otherwise the fold indices stored in 'cv' no longer line up with 'x'
.check_x_matches_cv <- function(x, cv){
    if(nrow(x) != .cv_n_points(cv)){
        stop("Number of rows in 'x' does not match the folds in 'cv'!", call. = FALSE)
    }
    invisible(TRUE)
}


# check points fall within the raster layer
.check_within <- function(x, r) {
    bbox <- sf::st_bbox(x)
    ex <- terra::ext(r)

    y <- bbox[1] >= terra::xmin(ex) &
        bbox[3] <= terra::xmax(ex) &
        bbox[2] >= terra::ymin(ex) &
        bbox[4] <= terra::ymax(ex)

    if(!y) stop("The x's bounding box lies outside the raster extent.")
}


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

# validate presence-background data and return the row indices of the presences (1s).
# When presence_bg is FALSE the full set of row indices is returned unchanged. When
# TRUE, 'column' is required and must hold a binary 0/1 (numeric) response.
.presence_index <- function(x, column, presence_bg){
    if(!isTRUE(presence_bg)){
        return(seq_len(nrow(x)))
    }
    if(is.null(column)){
        stop("'column' must be provided for presence-background data.")
    }
    vals <- x[, column, drop = TRUE]
    unqsp <- unique(vals)
    if(!is.numeric(unqsp) || any(unqsp < 0) || any(unqsp > 1)){
        stop("Presence-background option is only for species data with 0s (backgrounds/pseudo-absences) and 1s (presences).\n",
             "The data should be numeric.\n")
    }
    which(vals == 1)
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


# check raster extent is valid
.check_ext <- function(r) {
    tryCatch(
        {
            e <- terra::ext(r)
            vals <- e[1:4]
        },
        error = function(cond) {
            stop("Failed to extract raster extent: ", conditionMessage(cond))
        }
    )

    y <- all(is.finite(vals)) &&
        (terra::xmax(e) > terra::xmin(e)) &&
        (terra::ymax(e) > terra::ymin(e))

    if (!y) stop("Invalid raster extent: values are non-finite or out of range.")
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
    # check for a valid extent
    .check_ext(r)

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

