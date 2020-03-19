#' Take the output from plt_summary function and run phenesse on the cell of
#' interest.
#' 
#' @param plt_summary_output The named output of the plt_summary function
#' @param minimum_obs Minimum records per cell to be used.
#' @param earliest_year Earliest year to be included in analysis
#' @param latest_year Latest year to be included in analysis
#' @param flowering_cutoff Day of year to filter observations by if next year's
#' flowering maybe occurring in December
#' @param num_cores Number of cores to use in calculation
#' @return A dataframe of the onset, offset, and duration
#'  calculation for the cells of interest.
#'  
#' @importFrom dplyr lubridate phenesse parallel pbmcapply

run_phenesse <- function(plt_summary_output, minimum_obs, earliest_year, last_year, 
                         flowering_cutoff,num_cores){
  # new weib_percentile function w/ Daijiang's fixes. Can get rid of these lines
  # once the new version of phenesse is on cran
  weib_percentile_dl <- function(observations, percentile, iterations = 500){
    
    # curve_intersect determines where two lines intersect
    # parameters needed are two dataframes with two columns, x and y, which could
    # be plotted.
    
    curve_intersect <- function(curve1, curve2, empirical=TRUE, domain=NULL) {
      if (!empirical & missing(domain)) {
        stop("'domain' must be provided with non-empirical curves")
      }
      
      if (!empirical & (length(domain) != 2 | !is.numeric(domain))) {
        stop("'domain' must be a two-value numeric vector, like c(0, 10)")
      }
      
      if (empirical) {
        # Approximate the functional form of both curves
        curve1_f <- stats::approxfun(curve1$x, curve1$y, rule = 2)
        curve2_f <- suppressWarnings(stats::approxfun(curve2$x, curve2$y, rule = 2))
        
        # Calculate the intersection of curve 1 and curve 2 along the x-axis
        point_x <- stats::uniroot(function(x) curve1_f(x) - curve2_f(x),
                                  c(min(curve1$x), max(curve1$x)))$root
        
        # Find where point_x is in curve 2
        point_y <- curve2_f(point_x)
      } else {
        # Calculate the intersection of curve 1 and curve 2 along the x-axis
        # within the given domain
        point_x <- stats::uniroot(function(x) curve1(x) - curve2(x), domain)$root
        
        # Find where point_x is in curve 2
        point_y <- curve2(point_x)
      }
      
      return(list(x = point_x, y = point_y))
    }
    
    # use a data frame to plot a smooth CDF from -0.001 to 1.001 given
    # our original observations
    
    create_predict_df <- function(observations){
      # previous create_cdf_ends()
      weib <- fitdistrplus::fitdist(observations, distr = "weibull",
                                    method = "mle", lower = c(0, 0))
      cdf0 <- as.numeric(weib$estimate['scale']*
                           (-log(1-0.01))^(1/weib$estimate['shape']))
      cdf100 <- as.numeric(weib$estimate['scale']*
                             (-log(1-0.99))^(1/weib$estimate['shape']))
      added_vec <- sort(append(observations, values = c(cdf0, cdf100)),
                        decreasing = FALSE)
      
      new_vec <- seq(from = min(added_vec), to = max(added_vec), by = 0.5)
      
      cdfadded <- 1 - exp(-(new_vec/weib$estimate['scale'])^weib$estimate['shape'])
      
      cdf_df <- data.frame(x = new_vec, y = cdfadded)
      ends <- data.frame(x = c(min(added_vec - 1), max(added_vec + 1)),
                         y = c(-0.001,1.001))
      cdf_df <- rbind(cdf_df, ends)
      cdf_df <- cdf_df[order(cdf_df$x, decreasing = FALSE),]
      
      return(cdf_df)
    }
    
    # calculates the theta hat value for each iteration, which
    # when averaged is used to calculate the bias value
    
    get_theta_hat_i <- function(observations, percentile){
      
      emptyvec <- vector(mode = "numeric", length = length(observations))
      df1 <- create_predict_df(observations) # move out of loop since constant
      
      for(i in seq_along(observations)){
        sim_vector <- stats::runif(n = 1, min = 0, max = 1)
        df2 <- data.frame(x = observations, y = sim_vector)
        emptyvec[i] <- curve_intersect(df1, df2)$x
      }
      
      new_vector <- sort(emptyvec, decreasing = FALSE)
      
      new_df1 <- create_predict_df(new_vector)
      
      new_df2 <- data.frame(x = observations, y = percentile)
      
      theta_hat_i <- curve_intersect(new_df1, new_df2)$x
      
      return(theta_hat_i)
    }
    
    theta_hat_df <- data.frame(x = observations, y = percentile)
    
    # calculate theta hat original value
    theta_hat <- curve_intersect(curve1 = create_predict_df(observations),
                                 curve2 = theta_hat_df)[['x']]
    # calculate bias value based off of the mean of many theta hat i calculations
    bias <- mean(replicate(n = iterations,
                           expr = get_theta_hat_i(observations = observations,
                                                  percentile = percentile)))
    # final calculation, which gives you theta bar!
    2 * theta_hat - bias
  }
  
  # make Daijiang function output into dataframe
  df <- plt_summary_output$dat_to_use
  
  df <- df %>% 
    mutate(year = year(observed_on)) %>% 
    mutate(doy = yday(observed_on))
  
  # filter data to only years of interest
  df <- df %>% 
    filter(year >= earliest_year & year <= last_year) %>% 
    filter(doy <= flowering_cutoff)
  
  # count number of records for each cell x year combination
  num_of_records <- df %>% 
    group_by(year, id_cells) %>% 
    summarise(count = n())
  
  not_enough_data <- num_of_records %>% 
    filter(count >= minimum_obs) %>% 
    as.data.frame()
  
  # remove cell, year combinations that do not have enough records
  df_manip <- as.data.frame(df) %>% 
    group_by(year, id_cells, doy) %>% 
    summarise(obs = n())
  
  enough_records <- left_join(df_manip, not_enough_data) %>% 
    filter(!is.na(count)) %>% 
    select(-geometry)
  
  # make list with all doy values in it for each cell x year combination
  species_cell_year <- split(enough_records, 
                             f = list(enough_records$year, enough_records$id_cells),
                             drop = TRUE)
  
  # lapply functions
  onsetestimator <- function(x){
    onset <- weib_percentile_dl(observations = x$doy, iterations = 500, percentile = 0)
    return(onset)
  }
  
  offsetestimator <- function(x){
    offset <- weib_percentile_dl(observations = x$doy, iterations = 500, percentile = 1.0)
    return(offset)
  }
  
  # Estimate onseet and offset
  if(num_cores > 1){
    onset <- unlist(pbmclapply(species_cell_year, FUN = onsetestimator, mc.cores = num_cores))  
    offset <- unlist(pbmclapply(species_cell_year, FUN = offsetestimator, mc.cores = num_cores))} else{
      onset <- unlist(lapply(species_cell_year, FUN = onsetestimator))  
      offset <- unlist(lapply(species_cell_year, FUN = offsetestimator)) 
    }  
  # split outputs back to df
  onset_df <- as.data.frame(split(onset, 1:1)) %>% 
    tibble::rownames_to_column(var = "rowname")
  onset_df <- onset_df %>% 
    mutate(year = substr(x = rowname, start = 0, stop = 4 )) %>% 
    mutate(id_cells = substr(x = rowname, start = 6, stop = nchar(rowname))) %>% 
    rename(onset = X1)
  onset_df$year <- as.numeric(onset_df$year)
  onset_df$id_cells <- as.numeric(onset_df$id_cells)
  offset_df <- as.data.frame(split(offset, 1:1)) %>% 
    tibble::rownames_to_column(var = "rowname")
  offset_df <- offset_df %>% 
    mutate(year = substr(x = rowname, start = 0, stop = 4 )) %>% 
    mutate(id_cells = substr(x = rowname, start = 6, stop = nchar(rowname))) %>% 
    rename(offset = X1)
  offset_df$year <- as.numeric(offset_df$year)
  offset_df$id_cells <- as.numeric(offset_df$id_cells)  
  
  # join estimates with original sf dataframe based on cell_ids and year
  cell_duration <- left_join(onset_df, offset_df) 
  cell_duration <- cell_duration %>% 
    mutate(duration = offset - onset)
  
  return(cell_duration)
}