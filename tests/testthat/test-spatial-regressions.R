test_that("geographic autocorrelation plots preserve short and very large ranges", {
    for(mode in 1:2){
        set.seed(if(mode == 1) 7 else 42)
        n <- if(mode == 1) 400 else 300
        lon <- runif(n, 5, if(mode == 1) 6 else 5.3)
        lat <- runif(n, 52, if(mode == 1) 53 else 52.3)
        resp <- if(mode == 1){
            sin(lon * 100) + sin(lat * 100) + rnorm(n, 0, 0.2)
        } else{
            (lon - 5) * 10 + (lat - 52) * 10 + rnorm(n, 0, 0.3)
        }
        pts <- sf::st_as_sf(data.frame(resp, lon, lat),
                            coords = c("lon", "lat"), crs = 4326)

        expect_silent(sac <- cv_spatial_autocor(
            x = pts, column = "resp", plot = FALSE, progress = FALSE
        ))
        expect_equal(sac$range, sac$variograms[[1]]$var_model[2, "range"] * 1000)
        if(mode == 1) expect_lt(sac$range, 111325 / 2)
        if(mode == 2) expect_gt(sac$range / 111325, 90 - min(lat))
        expect_true(ggplot2::is_ggplot(sac$plots))
        expect_no_error(ggplot2::ggplot_build(sac$plots))
    }
})

test_that("predefined fold counts cover every point and reject invalid IDs", {
    bbox <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 7, ymax = 1), crs = 3857)
    blocks <- sf::st_sf(fold_id = 1:7,
                        geometry = sf::st_make_grid(sf::st_as_sfc(bbox), n = c(7, 1)))
    pts <- sf::st_as_sf(data.frame(x = (0:6) + 0.5, y = 0.5),
                        coords = c("x", "y"), crs = 3857)
    make_cv <- function(...){
        cv_spatial(x = pts, user_blocks = blocks, folds_column = "fold_id",
                   selection = "predefined", hexagon = FALSE,
                   plot = FALSE, report = FALSE, progress = FALSE, ...)
    }

    expect_silent(cv <- make_cv())
    expect_equal(cv$k, 7)
    expect_equal(cv$folds_ids, 1:7)
    tested <- unlist(lapply(cv$folds_list, function(f) f[[2]]), use.names = FALSE)
    expect_equal(sort(tested), seq_len(nrow(pts)))
    expect_equal(rowSums(!cv$biomod_table), rep(1, nrow(pts)))
    expect_message(explicit <- make_cv(k = 2), "k = 7")
    expect_equal(explicit$folds_ids, cv$folds_ids)

    for(ids in list(0:6, c(1:6, 8))){
        blocks$fold_id <- ids
        expect_error(make_cv(), "consecutive integers starting at 1")
    }
})

test_that("plot labels remain correct when a fold has no test or train points", {
    pts <- sf::st_as_sf(data.frame(x = 1:3, y = 1:3),
                        coords = c("x", "y"), crs = 3857)
    for(train_only in c(TRUE, FALSE)){
        train <- if(train_only) 1:3 else integer(0)
        test <- if(train_only) integer(0) else 1:3
        cv <- list(folds_list = list(list(train, test)), folds_ids = rep(1, 3))
        x_long <- .x_to_long(pts, cv, num_plots = 1)
        expect_equal(as.character(x_long$value), rep(if(train_only) "Train" else "Test", 3))
    }
})

test_that("block-size sliders use metre limits, the long extent, and valid values", {
    skip_if_not_installed("shiny")
    slider_values <- function(app){
        html <- as.character(get("ui", envir = environment(app$httpHandler)))
        vapply(c("min", "max", "from"), function(name){
            attr <- regmatches(html, regexpr(paste0('data-', name, '="[^"]*"'), html))
            as.numeric(gsub('^[^"]*"|"$', "", attr))
        }, numeric(1))
    }
    pts <- sf::st_as_sf(data.frame(x = c(0, 1e6, 0, 1e6), y = c(0, 0, 4e5, 4e5)),
                        coords = c("x", "y"), crs = 32631)
    expect_silent(app <- cv_block_size(x = pts, min_size = 1e5, max_size = 8e5))
    expect_equal(slider_values(app), c(min = 1e5, max = 8e5, from = 1e5))
    expect_silent(app <- cv_block_size(x = pts, min_size = 1000, max_size = 10000))
    expect_equal(slider_values(app), c(min = 1000, max = 10000, from = 10000))
    expect_message(app <- cv_block_size(x = pts, max_size = 2e6), "biggest extent")
    expect_equal(slider_values(app)["max"], c(max = 1e6))

    pts <- sf::st_as_sf(data.frame(lon = c(5, 6, 5, 6), lat = c(52, 52, 53, 53)),
                        coords = c("lon", "lat"), crs = 4326)
    expect_silent(app <- cv_block_size(x = pts, min_size = 10000, max_size = 90000))
    expect_equal(slider_values(app), c(min = 10000, max = 90000, from = 11132))
    expect_silent(app <- cv_block_size(x = pts))
    expect_equal(slider_values(app), c(min = 5566, max = 111325, from = 11132))
    expect_message(
        expect_error(cv_block_size(x = pts, min_size = 2e5, max_size = 9e5),
                     "smaller than 'max_size'"),
        "biggest extent"
    )
})
