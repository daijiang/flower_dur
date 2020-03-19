source("rcode/00_pkg_functions.R")

d = readr::read_csv("data/Claytonia virginica_inat.csv") %>% 
  dplyr::select(longitude, latitude, everything()) %>% 
  filter(flowers == 1) %>% 
  drop_na(longitude, latitude) %>% 
  rename(id_iNat = id)

cell_25k = plt_summary(cell_size = 25000, dat = d, n_per_cell = 10)


duration_claytonia_25km <- run_phenesse(plt_summary_output = cell_25k,
                                        minimum_obs = 10, earliest_year = 2017, last_year = 2019,
                                        flowering_cutoff = 300, num_cores = 30)

write.csv(x = duration_claytonia_25km, file = "phenesse_outputs/duration_claytonia_25km.csv", row.names = FALSE)

# 50 km

cell_50k = plt_summary(cell_size = 50000, dat = d, n_per_cell = 10)


duration_claytonia_50km <- run_phenesse(plt_summary_output = cell_50k,
                                        minimum_obs = 10, earliest_year = 2017, last_year = 2019,
                                        flowering_cutoff = 300, num_cores = 30)


write.csv(x = duration_claytonia_50km, file = "phenesse_outputs/duration_claytonia_50km.csv", row.names = FALSE)