library(here)
library(lubridate)
#install.Rtools()
library(tidyverse)
library(raster)

setwd(here::here())

library(sf)


#pop_grid <- raster("data\\population_grid\\india-spatial-india-census-2011-population-grid\\india_pop\\w001001.adf") 

winter_cropped_area <- raster("data/winter_cropped_area/India_cropped-area_1km_2016.tif") 


# (xmin,xmax,ymin,ymax)
b <- as(extent(73, 78, 28, 32), 'SpatialPolygons')
crs(b) <- crs(winter_cropped_area)
winter_cropped_area <- crop(winter_cropped_area, b)

plot(winter_cropped_area)




winter_cropped_area_dis <- raster::aggregate(winter_cropped_area, fact=3, fun = mean)
plot(winter_cropped_area_dis)

winter_cropped_area_cat <- winter_cropped_area_dis > 50
plot(winter_cropped_area_cat)

cell_size<-area(winter_cropped_area_cat, weights=FALSE, na.rm = T)
plot(cell_size)

# https://gis.stackexchange.com/questions/309407/computing-number-of-points-in-a-raster-grid-cell-in-r

fires <- sf::read_sf('data/viirs/fire_archive_SV-C2_247176.shp')

fires_count <- rasterize(fires, winter_cropped_area_cat, fun='count', field = 'FRP')
plot(fires_count)
fires_count <- fires_count * winter_cropped_area_cat


fires_total_frp <- rasterize(fires, winter_cropped_area_cat, fun='sum', field = 'FRP')
fires_total_frp <- fires_total_frp * winter_cropped_area_cat


plot(fires_total_frp)

rasters_list <- list('winter_cropped_area_cont' = winter_cropped_area_dis,
                     'winter_cropped_area_cat' = winter_cropped_area_cat, 'cell_size' = cell_size,
                     'fires_count' = fires_count, 'fires_total_frp' = fires_total_frp)

raster_brick <- brick(rasters_list)

outfile <- writeRaster(raster_brick, filename='data/winter_cropped_area/raster_brick.tif', format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))


pop_exposure <- raster("data/intermed_files/pop_exposure_matrix_10_years_with_crs.tif") 
pop_exposure <- flip(pop_exposure, direction='y')

#pop_exposure_nc <- raster("data/intermed_files/pop_exposure_matrix_10_years_with_crs.nc") 


#raster::?disaggregate
plot(pop_exposure)
resampled_pop_exposure <- resample(pop_exposure, winter_cropped_area_cat, method='bilinear')
plot(resampled_pop_exposure)

resampled_pop_exposure_dot <- resampled_pop_exposure * (fires_total_frp > 5)
plot(resampled_pop_exposure_dot)


hist(resampled_pop_exposure_dot)

pop_exposure_values <- resampled_pop_exposure_dot@data@values
pop_exposure_values <- pop_exposure_values[!is.na(pop_exposure_values)]
pop_exposure_values <- pop_exposure_values[pop_exposure_values != 0]

hist(pop_exposure_values)

median(pop_exposure_values)
mean(pop_exposure_values[pop_exposure_values < median(pop_exposure_values)])
mean(pop_exposure_values[pop_exposure_values > stats::quantile(pop_exposure_values, probs=0.9)])/median(pop_exposure_values)

mean(pop_exposure_values[pop_exposure_values > stats::quantile(pop_exposure_values, probs=0.9)])/mean(pop_exposure_values)
