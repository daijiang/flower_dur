source("Rcode/00_pkg_functions.R")

d = readr::read_csv("data/Claytonia virginica_inat.csv") %>% 
  dplyr::select(longitude, latitude, everything()) %>% 
  filter(flowers == 1) %>% 
  drop_na(longitude, latitude) %>% 
  rename(id_iNat = id)

cell_25k = plt_summary(cell_size = 25000, dat = d, n_per_cell = 10)

df_25k <- st_drop_geometry(cell_25k$dat_to_use)

message("look for potential outliers; and check whether flowering across calendar year;
          if so, when is a good cutoff for a flowering season?")

group_by(df_25k, observed_year, observed_on) %>%
  tally() %>% 
  filter(observed_year > 2014) %>% 
  ggplot(aes(x = observed_on, y = observed_year)) +
  geom_tile(aes(fill = n))


filter(df_25k, observed_year == 2017, observed_on > "2017-06-30")
# remove 2017-08-23, likely outlier
df_25k = filter(df_25k, !(observed_year == 2017 & observed_on == "2017-08-23"))


duration_claytonia_25km <- run_phenesse(df_25k, earliest_year = 2017, last_year = 2019,
                                        num_cores = 30)

write.csv(x = duration_claytonia_25km, 
          file = "phenesse_outputs/duration_claytonia_25km.csv", 
          row.names = FALSE)

# 50 km

cell_50k = plt_summary(cell_size = 50000, dat = d, n_per_cell = 10)
df_50k <- st_drop_geometry(cell_50k$dat_to_use)

message("look for potential outliers; and check whether flowering across calendar year;
          if so, when is a good cutoff for a flowering season?")

group_by(df_50k, observed_year, observed_on) %>%
  tally() %>% 
  filter(observed_year > 2014) %>% 
  ggplot(aes(x = observed_on, y = observed_year)) +
  geom_tile(aes(fill = n))


filter(df_50k, observed_year == 2017, observed_on > "2017-06-30")
# remove 2017-08-23, likely outlier
df_50k = filter(df_50k, !(observed_year == 2017 & observed_on == "2017-08-23"))


duration_claytonia_50km <- run_phenesse(df_50k,
                                        minimum_obs = 10, earliest_year = 2017, last_year = 2019,
                                        flowering_cutoff = 300, num_cores = 30)


write.csv(x = duration_claytonia_50km, file = "phenesse_outputs/duration_claytonia_50km.csv", row.names = FALSE)