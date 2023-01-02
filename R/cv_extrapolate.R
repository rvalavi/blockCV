#' Title
#'
#' @param cv
#' @param x
#' @param r
#' @param num_plot
#' @param mess
#' @param jitter_width
#' @param points_size
#' @param points_alpha
#' @param points_colors
#' @param progress
#'
#' @return
#' @export
#'
#' @examples
cv_extrapolate <- function(cv,
                           x,
                           r,
                           num_plot = seq_along(cv$folds_list),
                           mess = TRUE,
                           jitter_width = 0.1,
                           points_size = 2,
                           points_alpha = 0.7,
                           points_colors = NULL,
                           progress = TRUE){

  # check required packages
  pkg <- c("ggplot2", "dismo")
  .pkg_check(pkg)

  # check x is an sf object
  x <- .x_check(x)
  # change the r to terra object
  r <- .r_check(r)

  # check cv obj
  if(!class(cv) %in% c("cv_spatial", "cv_cluster", "cv_buffer")){
    stop("'cv' must be a blockCV cv_* object.")
  }

  # The iteration must be a natural number
  tryCatch(
    {
      num_plot <- abs(as.integer(num_plot))
      num_plot <- sort(num_plot)
    },
    error = function(cond) {
      message("'num_plot' must be a natural numbers.")
    }
  )

  # get the folds list
  folds_list <- cv$folds_list
  # length of the folds
  k <- length(folds_list)
  if(max(num_plot) > k){
    num_plot <- num_plot[num_plot <= k]
  }

  # extract the raster values
  points <- terra::extract(r, x, ID = FALSE)

  # to set as nrow for df; cv_buffer has only one target points unless P-BG
  n <- nrow(points)
  if(class(cv) == "cv_buffer"){
    n <- ifelse(cv$presence_background, nrow(points), 1)
  }

  m <- ncol(points)
  df <- data.frame(id = seq_len(n))

  if(progress) pb <- txtProgressBar(min = 0, max = length(num_plot), style = 3)

  # calculate MESS for testing data
  for(i in num_plot){
    df[, paste("Fold", i, sep = "")] <- NA
    train <- folds_list[[i]][[1]]
    test <- folds_list[[i]][[2]]
    mes <- sapply(1:m, function(j) dismo:::.messi3(points[test, j], points[train, j]))
    if(class(cv) == "cv_buffer"){
      mmes <- min(mes)
    } else{
      mmes <- apply(mes, 1, min, na.rm = TRUE)
    }
    df[1:length(mmes), paste("Fold", i, sep = "")] <- mmes
    if(progress) setTxtProgressBar(pb, i)
  }


  fold_names <- paste("Fold", num_plot, sep = "")
  # reshape for plotting
  mes_reshp <- stats::reshape(df,
                              direction = "long",
                              idvar = "id",
                              varying = fold_names,
                              times = fold_names,
                              v.names =  "value",
                              timevar = "folds"
  )
  # remove NAs
  mes_reshp <- mes_reshp[complete.cases(mes_reshp), ]
  if(class(cv) == "cv_buffer") mes_reshp$folds <- as.numeric(substr(mes_reshp$folds, 5, 20))
  # get the max value for color legend
  maxabs <- max(abs(mes_reshp$value))

  cols <- c("#D53E4F", "#FC8D59", "#FEE08B", "#FFFFBF", "#E6F598", "#99D594", "#3288BD")
  if(is.null(points_colors)) points_colors <- cols else points_colors

  goem_buffer <- ggplot2::geom_point(size = points_size, alpha = points_alpha)
  geom_other <- ggplot2::geom_jitter(width = jitter_width, size = points_size, alpha = points_alpha)
  geom_vio <- ggplot2::geom_violin(fill = NA)

  geom_exta <- if(class(cv) == "cv_buffer") goem_buffer else geom_other

  p1 <- ggplot2::ggplot(data = mes_reshp,
                        ggplot2::aes_string(x = "folds", y =  "value", col = "value")) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_exta +
    switch(class(cv) != "cv_buffer", geom_vio, NULL) +
    ggplot2::scale_color_gradientn(colours = points_colors,
                                   limits = c(-maxabs, maxabs),
                                   na.value = "#BEBEBE03") +
    ggplot2::labs(x = "Folds", y = "MESS Values", col = "MESS") +
    ggplot2::theme_bw()

  return(p1)
}
