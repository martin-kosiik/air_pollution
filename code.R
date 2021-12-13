#install.packages('splitr')
#install.packages("here")

#library(installr)
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




#the smoldering phase were reported to be coarser than those released during flaming phase
# (Ordou & Agranovski, 2019). Wardoyo et al. (2007) investigated the size distribution of smoke particles 
#in varying burning modes for some species of grass and found that the emitted smoke particles range in diameter between (30 nm and 60 nm) 
# for flaming phase and around (60 nm to 210 nm) for smoldering phase. The particles size of smoke particles released during agricultural burning
#were reported to be around 150 nm for cereal crops and 200 nm for wet fuels (e.g. montana grass) ( Zhang et al., 2011).









read_shapefile <- function(zipfile){
  temp <- tempfile()
  unzip(zipfile =  zipfile, exdir = temp)
  SHP_file <-list.files(temp, pattern = ".shp$",full.names=TRUE)
  village_map <- sf::read_sf(SHP_file)
  return(village_map)
}



hr_map_path <- 'C:/Users/marti/Dropbox/research_projects/india_water/data/india_2001_census_map/CopyOfindia-india-village-level-geospatial-socio-econ-1991-2001-hr-2001-shp.zip'

pb_map_path <- "C:\\Users\\marti\\Dropbox\\research_projects\\india_water\\data\\india_2001_census_map\\india-india-village-level-geospatial-socio-econ-1991-2001-pb-2001-shp.zip"

village_map_hr <- read_shapefile(hr_map_path)

village_map_pb <- read_shapefile(pb_map_path)



#village_map_hr <- sf::read_sf("data\\india_2001_census_map\\all_india\\all_india_village_level_map.shp")


# trims particles that are above the global max boundary value
#disp_df_trim <- disp_df[height <= 2665]


dispersion_sf <- st_as_sf(dispersion_tbl_after_day, coords = c("lon", "lat"), crs = 4326)

village_map_hr <- village_map_hr %>% 
  st_transform(crs = 4326)

#village_map_hr$geometry <- village_map_hr$geometry %>%
#  s2::s2_rebuild() %>%
 # sf::st_as_sfc()

sf::sf_use_s2(FALSE)

#(village_map_hr$particle_count <- lengths(st_intersects(village_map_hr, dispersion_sf$geometry)))



names(village_map_hr)

mapview(dispersion_sf)
mapview(village_map_hr, zcol= 'particle_count')
mapview(village_map_hr, zcol= 'part_count_per_ha')

ggplot(data = village_map_hr) +
# annotation_map_tile("osm", zoom = 8) +
  geom_sf(aes(fill = log(part_count_per_ha * 100 + 1))) +
  ggtitle(label = "Particles per ha")+
 # annotation_scale() +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "Particles per ha")

#  scale_fill_gradient2(low = "black", mid = "yellow",high = "red", midpoint = 1.5#, limits = c(0,0.15), oob = scales::squish )
#♦scale_fill_viridis_c(option = "magma")



pop_grid <- raster("data\\population_grid\\india-spatial-india-census-2011-population-grid\\india_pop\\w001001.adf") 

plot(pop_grid)

# (xmin,xmax,ymin,ymax)
b <- as(extent(74, 78, 28, 32), 'SpatialPolygons')
crs(b) <- crs(pop_grid)
pop_grid_crop <- crop(pop_grid, b)

myCol = terrain.colors(5)

plot(pop_grid_crop, breaks = c(60, 600, 6000, 60000), col = myCol)

pop_grid_crop_dis <- raster::aggregate(pop_grid_crop, fact=3, fun = sum)
plot(pop_grid_crop_dis, breaks = c(60, 600, 6000, 60000), col = myCol)

cell_size<-area(pop_grid_crop_dis, weights=FALSE, na.rm = T)
plot(cell_size)

res(pop_grid_crop_dis)
# https://gis.stackexchange.com/questions/309407/computing-number-of-points-in-a-raster-grid-cell-in-r

fires <- sf::read_sf('data\\viirs\\J1_VIIRS_C2_South_Asia_7d.shp')

fires_count <- rasterize(fires, pop_grid_crop_dis, fun='count', field = 'FRP')
plot(fires_count)


village_map_hr$area_km2 <- st_area(village_map_hr$geometry)/1000

village_map_hr <- village_map_hr %>% 
  st_transform(crs = 4326)

sf::sf_use_s2(FALSE)


get_disp_tables <- function(date_time = "2014-10-01 00:00"){
  dispersion_model <-
    create_dispersion_model() %>%
    add_source(
      name = "particle",
      lat = 30, lon = 77, height = 50, 
      rate = 5, pdiam = 1.1, density = 2, shape_factor = 0.8,
      release_start = lubridate::ymd_hm(date_time),
      release_end = lubridate::ymd_hm(date_time) + lubridate::hours(1)  
      
    ) %>%
    add_dispersion_params(
      start_time = lubridate::ymd_hm(date_time),
      end_time = lubridate::ymd_hm(date_time) + lubridate::hours(60),
      direction = "forward", 
      met_type = "reanalysis",
      met_dir = here::here("met"),
      exec_dir = here::here("out"), 
    ) %>%run_model()
  
  
  
  dispersion_tbl <- dispersion_model %>% get_output_tbl()
  
  return(dispersion_tbl)
  
}


dates_of_emission <- c('2013-10-01 00:00', '2014-10-01 00:00', '2015-10-01 00:00', '2016-10-01 00:00', '2017-10-01 00:00')

disp_tables_list <- dates_of_emission %>% 
  map(get_disp_tables)


names(disp_tables_list) <- dates_of_emission

# 59,220


disp_tables_list_after_day <- disp_tables_list %>% 
  bind_rows(.id = 'date_of_emission')


# 
# disp_tables_list_after_day <- disp_tables_list %>% 
#   map(~ . %>%   
#         group_by(particle_i) %>% 
#         mutate(reached_zero_height = cumsum(height == 0),
#               # reached_above_pbl = cumsum(height > 2665)
#                )%>%
#         ungroup() %>% 
#         mutate(reached_zero_height = (reached_zero_height > 0) *1,
#                #reached_above_pbl = (reached_above_pbl > 0) *1,
#                above_pbl = (height > 2665) * 1 )
#   )


#disp_tables_list_after_day <- disp_tables_list %>% 
#  map(~ filter(., height !=0, height < 2665))


disp_tables_list_after_day <- disp_tables_list_after_day %>% 
  group_by(date_of_emission, particle_i) %>% 
  mutate(reached_zero_height = cumsum(height == 0),
         # reached_above_pbl = cumsum(height > 2665)
  )%>%
  ungroup() %>% 
  mutate(reached_zero_height = (reached_zero_height > 0) *1,
         #reached_above_pbl = (reached_above_pbl > 0) *1,
         above_pbl = (height > 2665) * 1 )


disp_tables_list_after_day <- disp_tables_list_after_day %>% 
  filter(reached_zero_height !=1, above_pbl !=1)



get_particle_count <- function(dispersion_vec = disp_tables_list_after_day[[1]]){
  
  dispersion_sf <- st_as_sf(dispersion_vec, coords = c("lon", "lat"), crs = 4326)
  
  particle_count <- lengths(st_intersects(village_map_hr, dispersion_sf$geometry))
  
  part_count_per_ha <- particle_count/village_map_hr$AREA
  
  return(part_count_per_ha)
}

get_particle_count <- function(selected_date_of_emission = '2014-10-01 00:00'){
  
  dispersion_vec <- disp_tables_list_after_day %>% filter(date_of_emission == selected_date_of_emission)
  
  dispersion_sf <- st_as_sf(dispersion_vec, coords = c("lon", "lat"), crs = 4326)
  
  particle_count <- lengths(st_intersects(village_map_hr, dispersion_sf$geometry))
  
  part_count_per_km2 <- 1000* (particle_count/village_map_hr$area_km2)
  
  return(part_count_per_km2)
}


get_particle_count()

disp_tables_df <- dates_of_emission %>% 
  map_dfc(get_particle_count)

names(disp_tables_df) <- str_c('part_count_per_ha_', 1:5)

# tohle by mělo být try except
village_map_hr <- village_map_hr %>% 
  dplyr::select(-str_c('part_count_per_ha_', 1:5))

village_map_hr <- bind_cols(village_map_hr, disp_tables_df)
names(village_map_hr)

#village_map_hr <- village_map_hr %>% 
 # mutate_at(vars(str_c('part_count_per_ha_', 1:4)), list(mean = all_of(.)))

#village_map_hr <- village_map_hr %>% 
#  mutate(part_count_per_ha_mean = mean(all_of(str_c('part_count_per_ha_', 1:4)))) 


village_map_hr <- village_map_hr %>% 
  mutate(part_count_per_ha_mean = (part_count_per_ha_1 + part_count_per_ha_2 + part_count_per_ha_3+ part_count_per_ha_4 + part_count_per_ha_5)/5) 

names(village_map_hr)

drop_units(village_map_hr$part_count_per_ha_mean)

units(village_map_hr$part_count_per_ha_mean) <- NULL

library(viridis)


ggplot(data = village_map_hr) +
  # annotation_map_tile("osm", zoom = 8) +
  geom_sf(aes(fill = log(part_count_per_ha_mean * 100 + 1))) +
  ggtitle(label = "Particles per ha")+
  # annotation_scale() +
  scale_fill_viridis_c(option = "magma"#, limits = c(0,4.5)
                       ) +
  labs(fill = "Particles per ha")

ggsave('figures/hr_pollution_by_vilage_all_polygon_area.pdf', scale = 2)



ggplot(data = village_map_hr) +
  # annotation_map_tile("osm", zoom = 8) +
  geom_sf(aes(fill = part_count_per_ha_mean )) +
  ggtitle(label = "Particles per ha")+
  # annotation_scale() +
  scale_fill_gradient2(low = "black", mid = "yellow",high = "red", midpoint = 20, limits = c(0,40)
  )+
  #scale_fill_viridis_c(option = "magma", limits = c(0,50) ) +
  labs(fill = "Particles per ha")

ggsave('figures/hr_pollution_by_vilage_all_polygon_area_no_log.pdf', scale = 2)



ggplot(data = village_map_hr) +
  # annotation_map_tile("osm", zoom = 8) +
  geom_sf(aes(fill = part_count_per_ha_mean * 100)) +
  ggtitle(label = "Particles per ha")+
  # annotation_scale() +
  scale_fill_viridis( breaks=c(0.5,1,3,5), limits = c(0,5) ) +
  #scale_fill_viridis_c(option = "magma") +
  labs(fill = "Particles per ha")


village_map_hr <- village_map_hr %>% 
  mutate(TOT_P = as.numeric(TOT_P))

village_map_hr$TOT_P
units(village_map_hr$area_km2) <- NULL

village_map_hr %>% 
  filter(LEVEL != 'TEHSIL') %>% 
  ggplot() +
  # annotation_map_tile("osm", zoom = 8) +
  geom_sf(aes(fill = TOT_P/area_km2)) +
  ggtitle(label = "Total population (2001)")+
  # annotation_scale() +
  scale_fill_viridis_c( option = "magma", limits = c(0.001,8), trans = 'log'
  ) +
  labs(fill = "Total population (2001)")


village_map_hr %>% 
  filter(LEVEL == 'TEHSIL') %>% 
  View()

'JAGADHRI'

village_map_hr %>% 
  filter(NAME == 'JAGADHRI') %>% 
  View()


village_map_hr %>% 
  filter(LEVEL == 'TEHSIL') %>% 
  ggplot() +
  # annotation_map_tile("osm", zoom = 8) +
  geom_sf(aes(fill = TOT_P)) +
  ggtitle(label = "Total population (2001)")+
  # annotation_scale() +
  scale_fill_viridis_c(trans = "log", option = "magma"#, limits = c(0,4.5)
  ) +
  labs(fill = "Total population (2001)")




village_map_hr %>% 
  #filter(LEVEL == 'TEHSIL') %>% 
  mapview()


village_map_pb %>% 
  ggplot() +
  # annotation_map_tile("osm", zoom = 8) +
  geom_sf(aes(fill = TOT_P)) +
  ggtitle(label = "Total population (2001)")+
  # annotation_scale() +
  scale_fill_viridis_c(trans = "log", option = "magma"#, limits = c(0,4.5)
  ) +
  labs(fill = "Total population (2001)")

