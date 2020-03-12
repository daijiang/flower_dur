if(!require(xfun)) install.packages("xfun")
xfun::pkg_attach2(c("tidyverse", "sf"))

#' Plot observation data and group by grid cells
#' 
#' @param us_map US map.
#' @param cell_size Size (meters) of the square grid cells.
#' @param dat Data frame of observed information. The first column must be longitude;
#' the second column must be latitude.
#' @param n_per_cell Minimum records per cell to be used.
#' @param show_fig Plot the figure at the same time?
#' @return A list of maps, a data frame to summarise number of records per cell,
#' and a data frame of the records that fall within cells with enough records.
#' 
plt_summary = function(us_map = readRDS("data/usa_map.rds"), 
                       cell_size = 50000, 
                       dat = tibble(long = runif(300, -110, -85),
                                    lat = runif(300, 26, 45), z = 1),
                       n_per_cell = 10,
                       show_fig = FALSE){
  # make grid over us
  grids = st_make_grid(us_map, cellsize = c(cell_size, cell_size))
  # add grid cell id
  grids = mutate(st_sf(geometry = grids), id_cells = 1:n())
  # convert lat/long to have the same crs
  dat = st_transform(st_as_sf(dat, coords = 1:2, crs = 4326), 
                         crs = st_crs(us_map)$proj4string)
  # which cell each point falls in?
  dat_cells = st_join(dat, grids)
  dat_cells_count = group_by(st_drop_geometry(dat_cells), id_cells) %>% 
    tally() %>% 
    arrange(desc(n)) %>% 
    mutate(enough_data = n >= n_per_cell)
  
  # cells with data
  cells_with_data = dplyr::filter(grids, id_cells %in% dat_cells_count$id_cells) %>% 
    left_join(dat_cells_count, by = "id_cells")
  cells_with_data = bind_cols(cells_with_data, 
                              suppressWarnings(st_centroid(cells_with_data) %>% 
                                st_transform(4326) %>% 
                                st_coordinates() %>% 
                                as.data.frame() %>% 
                                rename(long_cell = X, lat_cell = Y)))
 
  plt = ggplot() +
    geom_sf(data = us_map) +
    geom_sf(data = grids, alpha = 0, size = 0.1, color = "gray") +
    geom_sf(data = dat, size = 0.5, alpha = 0.6) +
    geom_sf(data = filter(cells_with_data, enough_data), alpha = 0, 
            size = 0.15, color = "red") +
    labs(title = paste(nrow(filter(cells_with_data, enough_data)),
                       "highlighted cells with records more than",
                       n_per_cell, 
                       "(Cell resolution:", cell_size/1000, "km by",
                       cell_size/1000, "km)",
                       collapse = " "))
  if(show_fig){
    print(plt)
  }
 
  cat(nrow(filter(cells_with_data, enough_data)), 
      "cells with records more than", n_per_cell, "\n")
  
  # records fall within cells with >= n_per_cell records
  dat_to_use = filter(dat_cells, id_cells %in%
           filter(cells_with_data, enough_data)$id_cells)
  
  list(cells_with_data = cells_with_data, dat_to_use = dat_to_use, fig = plt)
}

# usage examples ----
d = readr::read_csv("data/Claytonia virginica_inat.csv") %>% 
  dplyr::select(longitude, latitude, everything()) %>% 
  filter(flowers == 1) %>% 
  drop_na(longitude, latitude) %>% 
  rename(id_iNat = id)

cell_100k = plt_summary(cell_size = 100000, dat = d, n_per_cell = 10)
cell_100k$cells_with_data
cell_100k$dat_to_use
cell_100k$fig
