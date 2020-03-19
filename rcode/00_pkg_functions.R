if(!require(xfun)) install.packages("xfun")
if(!require(phenesse)) remotes::install_github("mbelitz/phenesse")
xfun::pkg_attach2(c("tidyverse", "sf", "lubridate", "phenesse", "pbmcapply"))

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
  if(!(grepl("^long", names(dat)[1], ignore.case = TRUE) &
              grepl("^lat", names(dat)[2], ignore.case = TRUE))){
    stop("The first two columns of dat must be longitude and latitude, respectively.")
  }
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
  # add centroid coords
  cells_with_data = bind_cols(cells_with_data, 
                              suppressWarnings(st_centroid(cells_with_data) %>% 
                                st_transform(4326) %>% 
                                st_coordinates() %>% 
                                as.data.frame() %>% 
                                rename(long_cell = X, lat_cell = Y)))
  
  # records fall within cells with >= n_per_cell records
  dat_to_use = filter(dat_cells, id_cells %in%
                        filter(cells_with_data, enough_data)$id_cells)
    
  plt_base = ggplot() +
    geom_sf(data = us_map) +
    geom_sf(data = grids, alpha = 0, size = 0.1, color = "gray") 
 
  plt = plt_base +
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

  
  list(cells_with_data = cells_with_data, grids = grids,
       dat_to_use = dat_to_use, fig = plt, fig_base = plt_base)
}

#' Take the output from plt_summary function and run phenesse on the cell of
#' interest.
#' 
#' @param plt_summary_output The named output of the plt_summary function
#' @param minimum_obs Minimum records per cell to be used.
#' @param earliest_year Earliest year to be included in analysis
#' @param latest_year Latest year to be included in analysis
#' @param flowering_cutoff Day of year to filter observations by if next year's
#' flowering maybe occurring in December
#' @param onset_perct Which percentile to use for onset?
#' @param offset_perct Which percentile to use for offset?
#' @param num_cores Number of cores to use in calculation
#' @return A dataframe of the onset, offset, and duration
#'  calculation for the cells of interest.
#'  run_phenesse(cell_100k, 100, 2018, 2018, num_cores = 4)

run_phenesse <- function(plt_summary_output, minimum_obs = 10, 
                         earliest_year = 2017, last_year = 2019, 
                         flowering_cutoff = 365, onset_perct = 0,
                         offset_perct = 1, num_cores){
  
  # make Daijiang function output into dataframe
  df <- plt_summary_output$dat_to_use %>% 
    st_drop_geometry() %>% 
    mutate(observed_year = year(observed_on)) %>% 
    mutate(observed_doy = yday(observed_on)) %>% 
    # filter data to only years of interest
    filter(observed_year >= earliest_year & observed_year <= last_year) %>% 
    filter(observed_doy <= flowering_cutoff)
  
  # count number of records for each cell x year combination
  num_of_records <- df %>% 
    group_by(observed_year, id_cells) %>% 
    summarise(count = n()) %>% ungroup() %>% 
    filter(count >= minimum_obs)
  
  # remove cell, year combinations that do not have enough records
  df_manip <- df %>% 
    group_by(observed_year, id_cells, observed_doy) %>% 
    summarise(n_obs = n()) %>% 
    ungroup() %>% 
    left_join(num_of_records,  by = c("observed_year", "id_cells")) %>% 
    filter(!is.na(count)) %>% 
    select(-count)
  
  # make list with all doy values in it for each cell x year combination
  species_cell_year <- split(df_manip, 
                             f = list(df_manip$observed_year,
                                      df_manip$id_cells),
                             drop = TRUE)
  
  # lapply functions
  setestimator <- function(x, niter = 500, perct = 0){
    tibble(est = weib_percentile(observations = x$observed_doy, 
                                 iterations = niter, percentile = perct))
  }
  
  # Estimate onseet and offset
  if(num_cores > 1){
    onset <- pbmclapply(species_cell_year, setestimator, 
                        perct = onset_perct, mc.cores = num_cores)
    offset <- pbmclapply(species_cell_year, setestimator, 
                         perct = offset_perct, mc.cores = num_cores)
  } else{
    onset <- lapply(species_cell_year, setestimator, perct = onset_perct)
    offset <- lapply(species_cell_year, setestimator, perct = offset_perct)
  }  
  
  # split outputs back to df
  onset_df = map_df(onset, ~.x, .id = "yr_cell") %>% 
    separate("yr_cell", into = c("observed_year", "id_cells"), sep = "[.]") %>% 
    mutate(id_cells = as.numeric(id_cells),
           observed_year = as.numeric(observed_year)) %>% 
    rename(onset = est)
  
  offset_df = map_df(offset, ~.x, .id = "yr_cell") %>% 
    separate("yr_cell", into = c("observed_year", "id_cells"), sep = "[.]") %>% 
    mutate(id_cells = as.numeric(id_cells),
           observed_year = as.numeric(observed_year)) %>% 
    rename(offset = est)
  
  # join estimates with original sf dataframe based on cell_ids and year
  cell_duration <- left_join(onset_df, offset_df, 
                             by = c("observed_year", "id_cells")) %>% 
    mutate(duration = offset - onset)
  
  return(cell_duration)
}

