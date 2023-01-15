library(here)
library(lubridate)
#install.Rtools()
library(tidyverse)
library(raster)
library(brick)
library(scales)
library(latex2exp)

setwd(here::here())
library(ncdf4)

library(mapview)
library(sf)
library(lwgeom)
library(sp)

# Change to 1 if you want to also execute the code for the approach without wind direction
compute_without_wind_direction <- 0


population <- raster("data/population_grid/ppp_2018_1km_Aggregated.tif") 

#plot(population)


b <- as(extent(70, 80, 26, 34), 'SpatialPolygons')
crs(b) <- crs(population)
#mapview(b)

population_crop <- crop(population, b)
plot(population_crop)

cell_size<-area(population_crop, weights=FALSE, na.rm = T)


pop_density <- population_crop/cell_size

#plot(pop_density)

#xy <- data.frame(x = 30, y = 75)
#xy <- data.frame(x = 75, y = 30)


# 1:ncell(r)
#pop_within_30km <- raster::extract(population_crop , y = 1:100 , buffer = 30000 , fun = sum )









###################################################
#Wind direction not taken into account

if(compute_without_wind_direction == 1){
  

x <- seq(from = 73, to = 78, by = 0.1)
y <- seq(from = 28, to = 32, by = 0.1)

## create a grid of all pairs of coordinates (as a data.frame)
xy <- expand.grid(x = x, y = y)



pop_within_30km <- raster::extract(population_crop , y = xy , buffer = 30000 , fun = sum, na.rm=T)

pop_within_30km_df <- xy
pop_within_30km_df$population <- pop_within_30km

dfr <- rasterFromXYZ(pop_within_30km_df)  #Convert first two columns as lon-lat and third as value                
plot(dfr/4)
dfr                  

plot((dfr/4) > 200000)


infant_share <- 25/1393
hist(pop_within_30km/4)

costs_of_per_km2_reduction_in_burned_area <- 2700/75 * 2.47 *100
infant_mort_increase_per_km2_burned <- 0.001

how_many_infants_we_save_per_dollar <- (dfr/4 * infant_share * infant_mort_increase_per_km2_burned)/costs_of_per_km2_reduction_in_burned_area

dollars_to_save_a_infant <- costs_of_per_km2_reduction_in_burned_area/(dfr/4 * infant_share * infant_mort_increase_per_km2_burned)

# 1 = how_many_infants_we_save_per_dollar * 

pop_within_30km_df$dollars_to_save_a_infant <- costs_of_per_km2_reduction_in_burned_area/(pop_within_30km_df$population/4 * infant_share * infant_mort_increase_per_km2_burned)

plot(how_many_infants_we_save_per_dollar)
plot(dollars_to_save_a_infant)
plot(dollars_to_save_a_infant, breaks = c(0, 100, 1000, 2000, 3000, 4000, 5000),
     col = terrain.colors(7))

ggplot() +
  geom_raster(data = pop_within_30km_df , 
              aes(x = x, y = y, 
                  fill = dollars_to_save_a_infant)) + 
  #scale_fill_viridis_c() +  
  scale_fill_viridis_c( begin = 0, end = 1, direction = -1, limits = c(0, 9000)) +  
  ggtitle("Elevation with hillshade") +
  coord_quickmap()


ggplot() +
  geom_raster(data = pop_within_30km_df , 
              aes(x = x, y = y, 
                  fill = dollars_to_save_a_infant)) + 
  #scale_fill_viridis_c() +  
  scale_fill_viridis_c( begin = 0, end = 1, direction = -1, limits = c(0, 5000)) +  
  ggtitle("Elevation with hillshade") +
  coord_quickmap()


hist(pop_within_30km_df$dollars_to_save_a_infant, xlim = c(0, 5000))

ggplot(pop_within_30km_df, aes(x = dollars_to_save_a_infant)) + geom_histogram() +   scale_x_log10()


ggplot(pop_within_30km_df, aes(x = dollars_to_save_a_infant)) + geom_histogram()  +xlim(0, 10000)

}

####################################
######################
# Wind direction taken into account ()

projected_crs <- crs(population_crop)
#Function to create circle with quadrants. Save desired projection as projected_crs
create_circle <- function(lat_x, long_y, theta_x=0, buffer_m=30000){
  #Create circle with radius buffer_m centered at (lat_x, long_y)
  circle_buffer  <- st_point(c(lat_x, long_y)) %>%  st_sfc(crs = 4326) %>% 
    st_cast("POINT")  %>% 
    st_transform(projected_crs) %>%
    st_buffer(buffer_m)
  
  #Create two orthogonal lines at origin 
  p1 <- rbind(c(lat_x,long_y - 1), c(lat_x,long_y + 1))
  p2 <- rbind(c(lat_x+1,long_y), c(lat_x-1,long_y))
  mls <- st_multilinestring(list(p1,p2))  %>% st_sfc(crs = 4326) %>% 
    st_transform(projected_crs) 
  
  #Use orthogonal lines to split circle into 4 quadrants
  x1 <- st_split(circle_buffer, mls) 
  
  #Convert origin into projected CRS
  center_in_crs  <- st_point(c(lat_x, long_y)) %>% 
    st_sfc(crs = 4326) %>%
    st_transform(projected_crs)
  
  sp_obj <- x1 %>% st_collection_extract(type="POLYGON") %>%
    #Convert to spatial to use sp functions
    as_Spatial() %>% 
    #rotate x degrees
  #  elide(rotate = theta_x + 45, center  = center_in_crs[[1]]) %>% 
    #return to sf 
    st_as_sf()
}


#install.packages('lwgeom')


#source_xy_df <- head(xy)

x <- seq(from = 73, to = 78, by = 0.1)
y <- seq(from = 28, to = 32, by = 0.1)

## create a grid of all pairs of coordinates (as a data.frame)
source_xy_df <- expand.grid(x = x, y = y)


map2(source_xy_df$x, source_xy_df$y, ~ .x + .y)

circles_list <- map2(source_xy_df$x, source_xy_df$y, ~  create_circle(lat_x = .x, long_y = .y))

pop_within_30km_quart_list <- map(circles_list, ~ raster::extract(population_crop , y = .x  , fun = sum, na.rm=T))

source_df_list <- map2(source_xy_df$x, source_xy_df$y, ~
                         data.frame(x = rep(.x, 4), y = rep(.y, 4), quadrant = c('nw', 'sw', 'se', 'ne')) )

source_df_list <- map2(source_df_list, pop_within_30km_quart_list , ~ .x %>% mutate(pop_within_30km = .y))



source_df_all <- source_df_list %>% 
  bind_rows()


# https://link.springer.com/content/pdf/10.1007/s00703-017-0512-2.pdf
# nw_share <- 0.06 + 0.16 + 0.18 + 0.03
# ne_share <- 0.005+ 0.01 + 0.03 + 0.003
# sw_share <- 0.03 + 0.05 + 0.07 + 0.005 + 0.005
# se_share <- 0.07 + 0.01 + 0.01 + 0.06+ 0.005
# 
# wd_sum <- ne_share +nw_share + sw_share + se_share
# 
# nw_share <- nw_share * (1/wd_sum)
# ne_share <- ne_share * (1/wd_sum)
# sw_share <- sw_share * (1/wd_sum)
# se_share <- se_share * (1/wd_sum)
# 
# ne_share +nw_share + sw_share + se_share

# https://www.weatheronline.in/weather/maps/city?FMM=10&FYY=2000&LMM=11&LYY=2000&WMO=42182&CONT=inin&REGION=0024&LAND=II&ART=WDR&R=0&NOREGION=0&LEVEL=162&LANG=in&MOD=tab

#nw_share <- 0.33 + 0.015 + 0.10
#ne_share <- 0.015+ 0.055
#sw_share <- 0.13 + 0.03 + 0.10
#se_share <- 0.14 + 0.03 + 0.055

#ne_share +nw_share + sw_share + se_share



wind_dir_u <- brick("data/wind_direction/uwnd.10m.mon.mean.nc")
wind_dir_v <- brick("data/wind_direction/vwnd.10m.mon.mean.nc")


all(names(wind_dir_u) == names(wind_dir_v))

layer_names <- names(wind_dir_u)
layer_months <- as.integer(substr(layer_names, 7,8))

wind_dir_u_fall <- raster::subset(wind_dir_u, which(layer_months %in% c(10, 11)  ))
wind_dir_v_fall <- raster::subset(wind_dir_v, which(layer_months  %in% c(10, 11) ))

# The meteorological convention for winds is that U component is positive for 
#a west to east flow (eastward wind) and the V component is positive for south to north flow (northward wind).

ne <- (wind_dir_u_fall > 0) * (wind_dir_v_fall >0)
nw <- (wind_dir_u_fall < 0) * (wind_dir_v_fall >0)
se <- (wind_dir_u_fall > 0) * (wind_dir_v_fall <0)
sw <- (wind_dir_u_fall < 0) * (wind_dir_v_fall <0)


ne_mean <- mean(ne)
nw_mean <- mean(nw)
se_mean <- mean(se)
sw_mean <- mean(sw)




source_df_all$ne_share <- extract(ne_mean, source_df_all %>% dplyr::select(x,y), method = 'bilinear')
source_df_all$nw_share <- extract(nw_mean, source_df_all %>% dplyr::select(x,y), method = 'bilinear')
source_df_all$se_share <- extract(se_mean, source_df_all %>% dplyr::select(x,y), method = 'bilinear')
source_df_all$sw_share <- extract(sw_mean, source_df_all %>% dplyr::select(x,y), method = 'bilinear')


#source_df_all %>% 
#  mutate(total_sum = ne_share + nw_share + se_share + sw_share) %>% 
#  View()



max(source_df_all$x)
min(source_df_all$x)
max(source_df_all$y)
min(source_df_all$y)


source_df_all %>% 
    mutate(total_sum = ne_share + nw_share + se_share + sw_share)




source_df_all <- source_df_all %>% 
  mutate(wind_dir_share = case_when(quadrant == 'nw' ~ nw_share,
                                    quadrant == 'ne' ~ ne_share,
                                    quadrant == 'sw' ~ sw_share,
                                    quadrant == 'se' ~ se_share))



source_df_all<- source_df_all %>% 
  group_by(x, y) %>% 
  summarise(weighted_pop = sum(wind_dir_share*pop_within_30km))




weighted_pop_with_30km <- rasterFromXYZ(source_df_all)  #Convert first two columns as lon-lat and third as value                
plot(weighted_pop_with_30km)
weighted_pop_with_30km                  

plot(weighted_pop_with_30km> 200000)


infant_share <- 25/1393

#costs_in_inr_per_unburned_acre <- 2700
costs_in_inr_per_unburned_acre <- 4051.3


costs_of_per_km2_reduction_in_burned_area <- costs_in_inr_per_unburned_acre/75 * 2.47 *100
infant_mort_increase_per_km2_burned <- 0.00096

how_many_infants_we_save_per_dollar <- (weighted_pop_with_30km * infant_share * infant_mort_increase_per_km2_burned)/costs_of_per_km2_reduction_in_burned_area

dollars_to_save_a_infant <- costs_of_per_km2_reduction_in_burned_area/(weighted_pop_with_30km * infant_share * infant_mort_increase_per_km2_burned)

# 1 = how_many_infants_we_save_per_dollar * 

source_df_all$dollars_to_save_a_infant <- costs_of_per_km2_reduction_in_burned_area/(source_df_all$weighted_pop * infant_share * infant_mort_increase_per_km2_burned)





km2_burned = 0.070294 * 1e6
lives_claimed <- 66000

lives_claimed/km2_burned



ggplot(source_df_all, aes(x = dollars_to_save_a_infant)) + geom_histogram() +   scale_x_log10()
ggplot(source_df_all, aes(x = dollars_to_save_a_infant)) + geom_histogram()  +xlim(0, 6000)


ggplot() +
  geom_raster(data = source_df_all , 
              aes(x = x, y = y, 
                  fill = dollars_to_save_a_infant)) + 
  #scale_fill_viridis_c() +  
  scale_fill_viridis_c( begin = 0, end = 1, direction = -1, limits = c(0, 5000)) +  
  #ggtitle("Elevation with hillshade") +
  coord_quickmap()



burned_area <- raster("data/burned_area/MCD64A1.006_Burn_Date_doy2019274_aid0001 (1).tif") 
plot(burned_area)

#crs(burned_area) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

burned_area <- projectRaster(burned_area, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )

crs(weighted_pop_with_30km) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

burned_area_cat <- (burned_area > 1)*1
plot(burned_area_cat)


resampled_weighted_pop_with_30km <- resample(weighted_pop_with_30km, burned_area_cat, method='bilinear')

plot(resampled_weighted_pop_with_30km)

resampled_weighted_pop_with_30km <- burned_area_cat * resampled_weighted_pop_with_30km

resampled_weighted_pop_with_30km[resampled_weighted_pop_with_30km <= 0] <- NA


how_many_infants_we_save_per_dollar <- (resampled_weighted_pop_with_30km * infant_share * infant_mort_increase_per_km2_burned)/costs_of_per_km2_reduction_in_burned_area

dollars_to_save_a_infant <- costs_of_per_km2_reduction_in_burned_area/(resampled_weighted_pop_with_30km * infant_share * infant_mort_increase_per_km2_burned)

#plot(dollars_to_save_a_infant)

#mapview(dollars_to_save_a_infant, layer.name = 'USD per an infant life saved')


hist(dollars_to_save_a_infant@data@values, xlab = 'USD per an infant life saved', main = 'Histogram (1 km^2 pixels are the unit of obs.)')

ggplot(data.frame(dollars_to_save_a_infant = dollars_to_save_a_infant@data@values),
       aes(x = dollars_to_save_a_infant)) + geom_histogram(col = "white") +  
  theme_minimal() + labs(x = 'USD per an infant life saved',  caption= 'pixel-level at approx. 1 km resolution') +
  theme(axis.line = element_line(size = 1), 
        # panel.grid.major.x = element_blank(), 
      #  panel.grid.minor.x = element_blank(), 
        text = element_text(size=16))
#+xlim(0, 6000)

ggsave('figures/radius_approach/hist_usd_per_life_saved_var_wind_direction.pdf')

library(ggmap)

extent(72, 80, 28, 33)

bbox <- c(left = 73, bottom = 28, right = 78, top = 33)






terrain_map <- get_stamenmap(bbox, zoom = 7, maptype = "terrain")
#map <- get_googlemap(center = c(lon= 76, lat = 30), zoom = 7, maptype = "hybrid")

#map <- get_openstreetmap(bbox, zoom = 7)


#https://stackoverflow.com/questions/48955504/how-to-overlay-a-transparent-raster-on-ggmap

# https://gis.stackexchange.com/questions/389050/how-can-you-overlay-raster-data-on-a-stamenmap-in-r


dollars_to_save_a_infant_data <- as.data.frame(rasterToPoints(dollars_to_save_a_infant))


ggmap(terrain_map) +
  geom_tile(data = dollars_to_save_a_infant_data, aes(x, y, fill = value), alpha = 0.5) +
  scale_fill_gradient(low = "green", high = "red")

ggplot() + 
  geom_raster(aes(x = x, y = y, fill = layer), data = dollars_to_save_a_infant_data)

ggplot() + 
  geom_tile(aes(x = x, y = y, fill = layer), data = dollars_to_save_a_infant_data)



ggmap(terrain_map) +
  geom_tile(aes(x = x, y = y, fill = layer), data = dollars_to_save_a_infant_data)+
  scale_fill_viridis_c(option = 'turbo', direction = -1)+ coord_quickmap()+
  labs(x = 'longitude', y = 'latitude', fill = 'USD per an infant life saved')#+  scale_fill_gradientn(limits = c(0, 4000))


ggmap(terrain_map) +
  geom_tile(aes(x = x, y = y, fill = layer), data = dollars_to_save_a_infant_data)+
  scale_fill_viridis_c(option = 'turbo', direction = -1)+ 
  labs(x = 'longitude', y = 'latitude', fill = 'USD per an infant life saved')#+  scale_fill_gradientn(limits = c(0, 4000))



ggsave('figures/radius_approach/map_usd_per_life_saved_var_wind_direction.pdf')



ggmap(terrain_map) +
  geom_tile(aes(x = x, y = y, fill = layer), data = dollars_to_save_a_infant_data)+
  scale_fill_gradient(low = "red", high = "blue")+ coord_quickmap()+
  labs(x = 'longitude', y = 'latitude', fill = 'USD per an infant life saved')
