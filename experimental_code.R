library(here)
library(lubridate)
#install.Rtools()
library(tidyverse)
library(raster)

library(sf)
library(mapview)

##devtools::install_github("rich-iannone/splitr")

#remove.packages('splitr')
#devtools::install_github("martin-kosiik/splitr")

library(splitr)
library(devtools)
library('rich-iannone/splitr')
source('my_functions.R')


dispersion_model <-
  create_dispersion_model() %>%
  add_source(
    name = "particle",
    lat = 30, lon = 77, height = 1.5, 
    rate = 0.7, pdiam = 1.1, density = 1.5, shape_factor = 0.8,
    release_start = lubridate::ymd_hm("2014-10-01 00:00"),
    release_end = lubridate::ymd_hm("2014-10-01 00:00") + lubridate::hours(1)  
    
  ) %>%
  add_dispersion_params(
    start_time = lubridate::ymd_hm("2014-10-01 00:00"),
    end_time = lubridate::ymd_hm("2014-10-01 00:00") + lubridate::hours(48),
    direction = "forward", 
    met_type = "reanalysis",
    met_dir = here::here("met"),
    exec_dir = here::here("out"), 
  ) #%>%run_model()



dispersion_model_run <- dispersion_model %>% 
  run_model()

dispersion_model_run <- dispersion_model %>% 
  run_my_model(my_particle_num = 250, my_particle_max = 250)

class(dispersion_model)
dispersion_tbl <- dispersion_model_run %>% get_output_tbl()


species_list <- 
  dispersion_model$sources[1, ] %>% 
  dplyr::select(-c(lat, lon, height)) %>%
  as.list()


trace(hysplit_dispersion, edit = T)

library(disperseR)

.Platform$OS.type == "windows"

try_something <- hysplit_dispersion(lat = 30, lon = 76, height = 1.5,
                                    start_day = "2014-10-01", start_hour = 0, duration = 2,
                                    direction = "forward", met_type = "reanalysis", vert_motion = 0,
                                    model_height = 20000, particle_num = 100, particle_max = 100,
                                    species = species_list, disp_name = NULL, binary_path = here::here('aux_files/'), 
                                    met_dir = here::here("met"), 
                                    exec_dir = here::here("out"), clean_up = F, softrun = F)


try_something <- hysplit_dispersion(lat = 30, lon = 77, height = 1.5,
                                    start_day = "2014-10-01", start_hour = 0, duration = 4,
                                    direction = "forward", met_type = "reanalysis", vert_motion = 0,
                                    model_height = 20000, particle_num = 100, particle_max = 500,
                                    species = species_list, disp_name = NULL, binary_path = 'C:/Users/marti/Dropbox/research_projects/air_pollution/aux_files/win/', 
                                    met_dir = here::here("met"), #binary_name = 'hycs_std.exe',
                                    exec_dir = here::here("out"), clean_up = F, softrun = F)


# https://github.com/rich-iannone/splitr/blob/main/R/utils.R




set_binary_path(binary_path = 'C:/Users/marti/Dropbox/research_projects/air_pollution/aux_files/', binary_name = 'hycs_std' )

# By default, binary names should be either:
#  - hyts_std (trajectory models)
#  - hycs_std (dispersion models)

# If a user uses another binary name, the path to it should also be specified


try_something <- hysplit_dispersion_right(lat = 30, lon = 76, height = 1.5,
                                          start_day = "2014-10-01", start_hour = 0, duration = 2,
                                          direction = "forward", met_type = "reanalysis", vert_motion = 0,
                                          model_height = 20000, particle_num = 100, particle_max = 100,
                                          species = species_list, disp_name = NULL, binary_path = NULL, met_dir = here::here("met"), 
                                          exec_dir = here::here("out"), clean_up = F)

dispersion_tbl %>% dispersion_plot()


plot(x=1:3, y=1:3)


dispersion_tbl %>% 
  filter(hour == 4) %>% View()




binary_path_2 <-
  system.file(
    file.path("win", paste0('hycs_std', ".exe")),
    package = "splitr"
  )




















haryana <- read_csv("data/census_2011/haryana.csv")
names(haryana)



haryana <- st_as_sf(haryana %>% filter(district_code != "Haryana"), wkt = 'geometery_in_wkt')

haryana <- haryana %>%  st_set_crs(4326)
plot(haryana$geometery_in_wkt)






