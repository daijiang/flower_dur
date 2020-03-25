if(!require(xfun)) install.packages("xfun")
if(!require(phenesse)) remotes::install_github("mbelitz/phenesse")
xfun::pkg_attach2(c("tidyverse", "sf", "lubridate", "phenesse", "pbmcapply",
                    "raster", "corrplot", "lmerTest", "broom.mixed"))

#' Plot observation data and group by grid cells
#' 
#' @param us_map US map.
#' @param cell_size Size (meters) of the square grid cells.
#' @param dat Data frame of observed information. The first column must be longitude;
#' the second column must be latitude.
#' @param n_per_cell Minimum records per cell to be used.
#' @param days_per_cell Minimum days of data per cell to be used.
#' @param sd_cutoff The cutoff of standard deviation of the number of records of 
#' all days for each year and cell combination; if higher than this cutoff, it
#' suggests that some days have way more records than the other days; we will thin
#' records for these days.
#' @param show_fig Plot the figure at the same time?
#' @param add_lakes_map Plot lakes?
#' @return A list of maps, a data frame to summarise number of records per cell,
#' and a data frame of the records that fall within cells with enough records.
#' 
plt_summary = function(us_map = readRDS("data/usa_map.rds"), 
                       cell_size = 50000, 
                       dat = tibble(long = runif(300, -110, -85),
                                    lat = runif(300, 26, 45), z = 1),
                       n_per_cell = 10, days_per_cell = 3, sd_cutoff = 3.5,
                       show_fig = FALSE, add_lakes_map = TRUE){
  # make grid over us
  grids = st_make_grid(us_map, cellsize = c(cell_size, cell_size))
  # add grid cell id
  grids = mutate(st_sf(geometry = grids), id_cells = 1:n())
  # convert lat/long to have the same crs
  if(!(grepl("^long", names(dat)[1], ignore.case = TRUE) &
              grepl("^lat", names(dat)[2], ignore.case = TRUE))){
    stop("The first two columns of dat must be longitude and latitude, respectively.")
  }
  dat = mutate(dat, observed_year = year(observed_on))
  dat = st_transform(st_as_sf(dat, coords = 1:2, crs = 4326), 
                         crs = st_crs(us_map)$proj4string)
  # which cell each point falls in?
  dat_cells = st_join(dat, grids)
  
  # is there any cell with most obs from one day of a year? likely because of iNat city challenge
  # 2019: April 26 - 29
  # 2018: April 27 - 30
  # 2017: April 14 - 18
  # 2016: April 14 - 21
  city_challenge_days = as.Date(c("2019-04-26", "2019-04-27", "2019-04-28", "2019-04-29",
                          "2018-04-27", "2018-04-28", "2018-04-29", "2018-04-30",
                          "2017-04-14", "2017-04-15", "2017-04-16", "2017-04-17", "2017-04-18",
                          "2016-04-14", "2016-04-15", "2016-04-16", "2016-04-17", 
                          "2016-04-18", "2016-04-19", "2016-04-20", "2016-04-21"))
  
  # identify cells/year with enough records but may be peaked at specific days (based on SD)
  cell_yr_peaks = group_by(st_drop_geometry(dat_cells), id_cells, observed_year, observed_on) %>% 
    tally() %>% # count per day
    group_by(id_cells, observed_year) %>% 
    summarise(n_days = n(), 
              n_records = sum(n, na.rm = T),
              ave_records_per_day = mean(n, na.rm = T), 
              sd_record_days = sd(n, na.rm = T)) %>% 
    filter(sd_record_days > sd_cutoff, n_records >= n_per_cell, n_days >= 5) %>% 
    ungroup()
  
  # data that need to thin
  dat_cells2 = left_join(dplyr::select(cell_yr_peaks, id_cells, observed_year), 
                         dat_cells, by = c("id_cells", "observed_year")) 
  # get n record per day and check whether they are from city challenge dates
  cell_yr_peak_days = group_by(dat_cells2, id_cells, observed_year, observed_on) %>% 
    tally() %>% 
    mutate(in_city_challenge = observed_on %in% city_challenge_days) %>% 
    arrange(id_cells, observed_year, desc(n)) %>% ungroup()
  # get the max number of records that are not from city challenge
  cell_yr_peak_no_city_challenge = cell_yr_peak_days %>% 
    filter(!in_city_challenge) %>% 
    group_by(id_cells, observed_year) %>% 
    summarise(max_n_per_day_not_city_challenge = max(n)) %>% 
    ungroup()
  
  cell_yr_peaks = left_join(cell_yr_peaks, cell_yr_peak_no_city_challenge, 
                            by = c("id_cells", "observed_year"))
  # decide how many records to keep for each day
  n_to_keep = bind_rows(
    # the top 1, no matter it is from city challenge or not
    group_by(cell_yr_peak_days, id_cells, observed_year) %>% slice(1L),
    # the first top 5 that are from city challenge
    group_by(cell_yr_peak_days, id_cells, observed_year) %>% 
      slice(1L:5L) %>% filter(in_city_challenge) 
  ) %>% unique() %>% 
    ungroup() %>% 
    left_join(cell_yr_peaks, by = c("id_cells", "observed_year")) %>% 
    filter(n > max_n_per_day_not_city_challenge) %>% 
    mutate(ave_records_per_day = round(ave_records_per_day),
           # if the day with max record is the dates of city challgend
           # keep the number of records that is equal to the max number from non-city-challenge days
           # otherwise, use the average number of record per day
           n_to_keep = ifelse(in_city_challenge, 
                              max_n_per_day_not_city_challenge,
                              ave_records_per_day)) %>% 
    dplyr::select(id_cells, observed_year, observed_on, n_to_keep)
  
  # data that are good to go
  dat_cells_asis = anti_join(dat_cells, 
                         dplyr::select(n_to_keep, id_cells, observed_year, observed_on), 
                         by = c("id_cells", "observed_year", "observed_on"))
  # random thin for these days
  dat_cells_thinned = left_join(n_to_keep, dat_cells, 
                                by = c("id_cells", "observed_year", "observed_on")) %>% 
    group_by(id_cells, observed_year, observed_on) %>% 
    sample_n(unique(n_to_keep)) %>% ungroup() %>% 
    dplyr::select(-n_to_keep) %>% st_sf()
  
  dat_cells = rbind(dat_cells_asis, dat_cells_thinned) %>% unique()
  
  # count number of records per cell and year combination
  dat_cells_count = group_by(st_drop_geometry(dat_cells), id_cells, observed_year) %>% 
    summarise(n_records = n(), n_days = n_distinct(observed_on)) %>% 
    mutate(enough_data = n_records >= n_per_cell & n_days >= days_per_cell) %>% 
    ungroup()
  
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
    
  if(add_lakes_map){
    plt_base = ggplot() +
      geom_sf(data = us_map) +
      geom_sf(data = readRDS("data/lakes.rds"), color = "gray50", fill = "white") +
      geom_sf(data = grids, alpha = 0, size = 0.1, color = "gray") 
  } else {
    plt_base = ggplot() +
      geom_sf(data = us_map) +
      geom_sf(data = grids, alpha = 0, size = 0.1, color = "gray") 
  }
 
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
#' @param df The named output `dat_to_use` of the plt_summary function; after 
#' some examniation for potential outliers and decide a new base date if needed
#' to deal with flowering across calendar years.
#' @param minimum_obs Minimum records per cell to be used.
#' @param minimum_days Minimum days per cell to be used.
#' @param earliest_year Earliest year to be included in analysis
#' @param latest_year Latest year to be included in analysis
#' @param flowering_cutoff Day of year to filter observations by if next year's
#' flowering maybe occurring in December
#' @param n_item Number of iterations for phenesse.
#' @param onset_perct Which percentile to use for onset?
#' @param offset_perct Which percentile to use for offset?
#' @param num_cores Number of cores to use in calculation
#' @return A dataframe of the onset, offset, and duration
#'  calculation for the cells of interest.
#'  run_phenesse(cell_100k, 100, 2018, 2018, num_cores = 4)

run_phenesse <- function(df, minimum_obs = 10, minimum_days = 3,
                         earliest_year = 2017, last_year = 2019, 
                         flowering_cutoff = "01-01", n_item = 500,
                         onset_perct = 0, offset_perct = 1, num_cores){
  
  df = df %>% 
    mutate(base_date = as.Date(ifelse(observed_on < paste(observed_year, flowering_cutoff, sep = "-"),
                              paste(observed_year - 1, flowering_cutoff, sep = "-"),
                              paste(observed_year, flowering_cutoff, sep = "-"))),
           observed_doy = as.numeric(observed_on - base_date) + 1,
           observed_year2 = year(base_date) + 1) # new defined "year"
  
  stopifnot(df$observed_doy >= 0 & df$observed_doy < 366)
  
  df = filter(df, observed_year2 >= earliest_year & observed_year2 <= last_year)
  
  # count number of records for each cell x year combination
  num_of_records <- df %>% 
    group_by(observed_year2, id_cells) %>% 
    summarise(n_records = n(), n_days = n_distinct(observed_on)) %>% ungroup() %>% 
    filter(n_records >= minimum_obs, n_days >= minimum_days)
  
  # remove cell, year combinations that do not have enough records
  df_manip <- df %>% 
    group_by(observed_year2, id_cells, observed_doy) %>% 
    summarise(n_obs = n()) %>% 
    ungroup() %>% 
    left_join(num_of_records,  by = c("observed_year2", "id_cells")) %>% 
    filter(!is.na(n_records)) %>% 
    dplyr::select(-n_records, -n_days)
  
  # make list with all doy values in it for each cell x year combination
  species_cell_year <- split(df_manip, 
                             f = list(df_manip$observed_year2,
                                      df_manip$id_cells),
                             drop = TRUE)
  
  # lapply functions
  setestimator <- function(x, niter = n_item, perct = 0){
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
    separate("yr_cell", into = c("observed_year2", "id_cells"), sep = "[.]") %>% 
    mutate(id_cells = as.numeric(id_cells),
           observed_year2 = as.numeric(observed_year2)) %>% 
    rename(onset = est)
  
  offset_df = map_df(offset, ~.x, .id = "yr_cell") %>% 
    separate("yr_cell", into = c("observed_year2", "id_cells"), sep = "[.]") %>% 
    mutate(id_cells = as.numeric(id_cells),
           observed_year2 = as.numeric(observed_year2)) %>% 
    rename(offset = est)
  
  # join estimates with original sf dataframe based on cell_ids and year
  cell_duration <- left_join(onset_df, offset_df, 
                             by = c("observed_year2", "id_cells")) %>% 
    mutate(duration = offset - onset,
           # convert back to date
           base_date = as.Date(paste(observed_year2 - 1, flowering_cutoff, sep = "-")),
           onset_date = base_date + onset,
           offset_date = base_date + offset) %>% 
    rename(observed_year = observed_year2)
  
  if(flowering_cutoff == "01-01"){
    cell_duration$observed_year = cell_duration$observed_year - 1
  }
  
  return(cell_duration)
}

