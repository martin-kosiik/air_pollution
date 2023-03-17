import os
import numpy as np
import pandas as pd
import datetime
import xarray as xr
import rasterio
import rioxarray
import matplotlib.pyplot as plt
#!pip install git+https://github.com/noaa-oar-arl/monetio.git
!pip install geocube


import geopandas as gpd



os.chdir(r'C:\Users\marti\Dropbox\research_projects\air_pollution')

pop_exposure_matrix = xr.open_dataarray("data/intermed_files/pop_exposure_matrix_10_years.nc")


pop_exposure_by_lonlat = pop_exposure_matrix.isel(band=0).set_index(source=['lon_source', 'lat_source'],).unstack("source")

# the layer we want LULC250K_1718


# winter cropped area https://sedac.ciesin.columbia.edu/data/set/india-india-annual-winter-cropped-area-2001-2016

#winter_cropped_area = rioxarray.open_rasterio('data/winter_cropped_area/India_cropped-area_1km_2016.tif')


winter_cropped_area_xr = rioxarray.open_rasterio('data/winter_cropped_area/India_cropped-area_1km_2016.tif')
winter_cropped_area_xr.plot()
winter_cropped_area_xr.rio.crs


winter_cropped_area_xr_mask = winter_cropped_area_xr.where(winter_cropped_area_xr > 60)
winter_cropped_area_xr_mask = winter_cropped_area_xr_mask.sel(x=slice(73, 79), y=slice(27, 24))
winter_cropped_area_xr_mask.plot()

with rasterio.open('data/winter_cropped_area/India_cropped-area_1km_2016.tif') as src:
        winter_cropped_area = src.read(1, window=crop_window)

plt.imshow(np.where(winter_cropped_area>60, 1, 0))

winter_cropped_area = rasterio.open('data/winter_cropped_area/India_cropped-area_1km_2016.tif')
np.where(winter_cropped_area.read(1)>70, 1, 0).plot()

winter_cropped_area.nodata
winter_cropped_area.nodatavals
winter_cropped_area.profile
winter_cropped_area.dataset_mask()
winter_cropped_area.colorinterp

winter_cropped_area.read(1)
plt.imshow(winter_cropped_area.read(1))
plt.show()

winter_cropped_area.read(1).shape
winter_cropped_area
winter_cropped_area.read(1)[400:3000,600:2000]

np.where(winter_cropped_area.read(1)>70, 1, 0)

from rasterio.windows import Window


crop_window = Window.from_slices((400, 1200), (400, 1200))



import rasterio
from rasterio.features import shapes
mask = None
with rasterio.Env():
    with rasterio.open('data/winter_cropped_area/India_cropped-area_1km_2016.tif') as src:
        #image = src.read(1) # first band
        image =np.where(src.read(1, window=crop_window)>60, 1, 0)
        results = (
        {'properties': {'raster_val': v}, 'geometry': s}
        for i, (s, v)
        in enumerate(
            shapes(image, mask=mask, transform=src.transform)) )

geoms = list(results)
 # first feature
print(geoms[0])


gpd_polygonized_raster  = gpd.GeoDataFrame.from_features(geoms)
# conda install -c conda-forge descartes
gpd_polygonized_raster.plot()


gpd_polygonized_raster_filter = gpd_polygonized_raster[gpd_polygonized_raster.raster_val == 1]
gpd_polygonized_raster_filter.plot()






from osgeo import ogr, gdal, osr
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import fiona
from shapely.geometry import shape
import rasterio.features

#segimg=glob.glob('Poly.tif')[0]
#src_ds = gdal.Open(segimg, GA_ReadOnly )
#srcband=src_ds.GetRasterBand(1)
#myarray=srcband.ReadAsArray()
#these lines use gdal to import an image. 'myarray' can be any numpy array

mypoly=[]
for vec in rasterio.features.shapes(myarray):
    mypoly.append(shape(vec))





viirs_fires = gpd.read_file(r'data\viirs\fire_archive_SV-C2_247176.shp')

gpd_polygonized_raster_filter = gpd_polygonized_raster_filter.set_crs(epsg=4326)
from geopandas.tools import sjoin
pointInPolys = sjoin(viirs_fires, gpd_polygonized_raster_filter, how='left')
print(pointInPolys.groupby(['index_right']).size().reset_index(name='count'))

print(pointInPolys.groupby(['index_right']).().reset_index(name='count'))


total_frp = pointInPolys.groupby(['index_right']).FRP.sum()

gpd_polygonized_raster_filter.loc[64]

gpd_polygonized_raster_filter['total_frp'] =  total_frp

gpd_polygonized_raster_filter.loc[gpd_polygonized_raster_filter['total_frp'].isna(),  'total_frp'] = 0


(gpd_polygonized_raster_filter.total_frp == 0).sum()


import geocube




fire_raster = rasterio.open('data/winter_cropped_area/raster_brick.tif')

from rasterio.plot import show

show(fire_raster.read(1, masked=True))

show(fire_raster.read(2, masked=True))

show(fire_raster.read(4, masked=True))

show(fire_raster.read(4, masked=True) > 5)




fire_raster_xr = rioxarray.open_rasterio('data/winter_cropped_area/raster_brick.tif')
fire_raster_xr.sel(band=4).plot()

fire_raster_xr_masked = fire_raster_xr.sel(band=4).where(fire_raster_xr.sel(band=4) >5 )
fire_raster_xr_masked.plot()
fire_raster_xr_masked_bin = fire_raster_xr_masked > 0
fire_raster_xr_masked_bin.plot()


fire_raster_xr_masked_bin = fire_raster_xr_masked_bin.rio.set_crs('+proj=longlat +datum=WGS84 +no_defs', inplace=False)

fire_raster_xr_masked_bin.rio.crs


pop_exposure_matrix = xr.open_dataarray("data/intermed_files/pop_exposure_matrix_10_years.nc")
#pop_exposure_matrix

pop_exposure_matrix = pop_exposure_matrix.isel(band=0).set_index(source=['lon_source', 'lat_source'],).unstack("source")
pop_exposure_matrix = pop_exposure_matrix.transpose('lat_source', 'lon_source')
pop_exposure_matrix = pop_exposure_matrix.rio.set_spatial_dims(x_dim = 'lon_source', y_dim='lat_source', inplace=False)
pop_exposure_matrix = pop_exposure_matrix.rio.set_crs('+proj=longlat +datum=WGS84 +no_defs', inplace = False) #


# pop_exposure_matrix = pop_exposure_matrix.set_index(x = 'lon_source', y='lat_source')

pop_exposure_matrix.rio.crs
fire_raster_xr_masked_bin.rio.crs

pop_exposure_matrix.plot()

rasterio.__version__

pop_exposure_matrix_match = pop_exposure_matrix.rio.reproject_match(fire_raster_xr_masked_bin)





#
