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


population <- raster("data/population_grid/ppp_2018_1km_Aggregated.tif") 

#plot(population)


b <- as(extent(73, 78, 28, 32), 'SpatialPolygons')
crs(b) <- crs(population)
#mapview(b)

population_crop <- crop(population, b)
plot(population_crop)

cell_size<-area(population_crop, weights=FALSE, na.rm = T)


pop_density <- population_crop/cell_size

#plot(pop_density)

xy <- data.frame(x = 30, y = 75)
xy <- data.frame(x = 75, y = 30)


# 1:ncell(r)
#pop_within_30km <- raster::extract(population_crop , y = 1:100 , buffer = 30000 , fun = sum )


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



####################################
# quarters

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


# https://link.springer.com/content/pdf/10.1007/s00703-017-0512-2.pdf
nw_share <- 0.06 + 0.16 + 0.18 + 0.03
ne_share <- 0.005+ 0.01 + 0.03 + 0.003
sw_share <- 0.03 + 0.05 + 0.07 + 0.005 + 0.005
se_share <- 0.07 + 0.01 + 0.01 + 0.06+ 0.005

wd_sum <- ne_share +nw_share + sw_share + se_share

nw_share <- nw_share * (1/wd_sum)
ne_share <- ne_share * (1/wd_sum)
sw_share <- sw_share * (1/wd_sum)
se_share <- se_share * (1/wd_sum)

ne_share +nw_share + sw_share + se_share

# https://www.weatheronline.in/weather/maps/city?FMM=10&FYY=2000&LMM=11&LYY=2000&WMO=42182&CONT=inin&REGION=0024&LAND=II&ART=WDR&R=0&NOREGION=0&LEVEL=162&LANG=in&MOD=tab

#nw_share <- 0.33 + 0.015 + 0.10
#ne_share <- 0.015+ 0.055
#sw_share <- 0.13 + 0.03 + 0.10
#se_share <- 0.14 + 0.03 + 0.055

ne_share +nw_share + sw_share + se_share

source_df_all <- source_df_list %>% 
  bind_rows() %>% 
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
plot(resampled_weighted_pop_with_30km> 200000)


infant_share <- 25/1393

costs_in_inr_per_unburned_acre <- 2700
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
  ggtitle("Elevation with hillshade") +
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

plot(dollars_to_save_a_infant)

hist(dollars_to_save_a_infant@data@values)

ggplot(data.frame(dollars_to_save_a_infant = dollars_to_save_a_infant@data@values),
       aes(x = dollars_to_save_a_infant)) + geom_histogram()  +xlim(0, 6000)
