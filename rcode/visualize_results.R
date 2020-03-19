source("rcode/00_pkg_functions.R")

# make function output into dataframe
cell_25k$dat_to_use
cell_25k$fig_base

est_25k <- read_csv("phenesse_outputs/duration_claytonia_25km.csv") %>% 
  select(id_cells, observed_year = year, onset, offset, dur = duration) %>% 
  left_join(cell_25k$grids, by = "id_cells") %>% 
  st_sf()

cell_25k$fig_base +
  geom_sf(data = filter(est_25k, observed_year == 2019), 
          mapping = aes(fill = onset), size = .1) +
  scale_fill_viridis_c(option = "A") +
  ggtitle("Onset")

cell_25k$fig_base +
  geom_sf(data = filter(est_25k, observed_year == 2019, offset < 250), 
          mapping = aes(fill = offset), size = .1) +
  scale_fill_viridis_c(option = "A") +
  ggtitle("Offset")

cell_25k$fig_base +
  geom_sf(data = filter(est_25k, observed_year == 2019), 
          mapping = aes(fill = dur), size = .1) +
  scale_fill_viridis_c(option = "A") +
  ggtitle("Duration")

