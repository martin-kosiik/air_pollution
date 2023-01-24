import os
import numpy as np
import pandas as pd
import datetime
import xarray as xr
import rasterio
import rioxarray
import matplotlib.pyplot as plt
#!pip install git+https://github.com/noaa-oar-arl/monetio.git

os.chdir('C:/Users/marti/Dropbox/research_projects/air_pollution')

region_of_interest = 'vietnam'
region_of_interest = 'india_low_res'

import monetio as mio

# hysplitfile = 'C:/HYSPLIT/working/cdumps/china/cdump_06100109'
# hysplitfile = 'C:/HYSPLIT/working/cdumps/pakistan/cdump_06100109'
# hysplitfile = 'C:/HYSPLIT/working/cdumps/india/cdump_06100109'
hysplitfile = f'hysplit/working/cdumps/{region_of_interest}/cdump_06100109'


def get_latlongrid_flat(dset, xindx, yindx):
    """
    INPUTS
    dset : xarray data set from ModelBin class
    xindx : list of integers
    yindx : list of integers
    RETURNS
    mgrid : output of numpy meshgrid function.
            Two 2d arrays of latitude, longitude.
    """
    llcrnr_lat = dset.attrs["llcrnr latitude"]
    llcrnr_lon = dset.attrs["llcrnr longitude"]
    nlat = dset.attrs["Number Lat Points"]
    nlon = dset.attrs["Number Lon Points"]
    dlat = dset.attrs["Latitude Spacing"]
    dlon = dset.attrs["Longitude Spacing"]

    lat = np.arange(llcrnr_lat, llcrnr_lat + nlat * dlat, dlat)
    lon = np.arange(llcrnr_lon, llcrnr_lon + nlon * dlon, dlon)
    print(nlat, nlon, dlat, dlon)
    print("lon shape", lon.shape)
    print("lat shape", lat.shape)
    #print(lat)
    #print(lon)
    lonlist = [lon[x - 1] for x in xindx]
    latlist = [lat[x - 1] for x in yindx]
    return (lonlist, latlist)



#my_modelbin = mio.models.hysplit.ModelBin(hysplitfile, verbose=False)


my_modelbin = mio.models.hysplit.ModelBin(hysplitfile, verbose=False)


time_period = 0

my_modelbin.dset.data_vars['0001'][time_period][0].plot(x="longitude", y="latitude")


my_modelbin.dset.data_vars['006F'][time_period][0].plot(x="longitude", y="latitude")
my_modelbin.dset.data_vars['006F'][1][0].plot(x="longitude", y="latitude")


my_modelbin.dset.attrs["Concentration Grid"]
my_modelbin.dset.attrs['llcrnr latitude']


all_source_xarray = my_modelbin.dset.to_array(dim='source')
sources_hexdec_index = list(my_modelbin.dset.data_vars.keys())

all_source_xarray['source']

int(sources_hexdec_index[2], base=16)

sources_dec_index = [int(x, base=16) for x in sources_hexdec_index]

with open(f'hysplit/working/cdumps/{region_of_interest}/CONTROL') as f:
    lines = f.readlines()



n_of_sources = int(lines[1].replace('\n', ''))

sources_coords_list = lines[2:(n_of_sources+2)]

sources_coords_list[0].split()
sources_coords_array = [x.split() for x in sources_coords_list]
sources_coords_array = np.array(sources_coords_array, dtype=np.float32)
sources_coords_array = sources_coords_array[:,0:2 ] # remove height colmun
sources_coords_df = pd.DataFrame(sources_coords_array, columns=['lat_source', 'lon_source'], index=sources_hexdec_index)

#sources_coords_df.loc['0015']



all_source_xarray = all_source_xarray.assign_coords(
        lon_source=('source', sources_coords_df.loc[all_source_xarray['source']].lon_source),
        lat_source=('source', sources_coords_df.loc[all_source_xarray['source']].lat_source)
)


all_source_xarray = all_source_xarray.isel(z=0)

all_source_xarray.shape



def read_conc_file(conc_file_path):
    hysplit_model = mio.models.hysplit.ModelBin(conc_file_path, verbose=False)
    all_source_xarray = hysplit_model.dset.to_array(dim='source')
    all_source_xarray = all_source_xarray.isel(z=0)
    all_source_xarray = all_source_xarray.assign_coords(
            lon_source=('source', sources_coords_df.loc[all_source_xarray['source']].lon_source),
            lat_source=('source', sources_coords_df.loc[all_source_xarray['source']].lat_source)
    )
    all_source_xarray = all_source_xarray.drop("longitude")
    all_source_xarray = all_source_xarray.drop("latitude")
    relative_times_array = (all_source_xarray.time.values - all_source_xarray.time.values[0]).astype('timedelta64[h]').astype(int)
    all_source_xarray = all_source_xarray.assign_coords(
            relative_time=(('time'), relative_times_array )
    ).set_index(time = ['relative_time'])
    return all_source_xarray


#year_list = ['06', '07', '08', '09', '10', '11', '12', '15', '16', '17']
#file_name_list = ['C:/HYSPLIT/working/cdumps/cdump_' + year + '10' + '01' + '00' for year in year_list]

# ten_days_length
file_path_list = os.listdir(f'hysplit/working/cdumps/{region_of_interest}/')
file_path_list = [x for x in file_path_list if x not in ['CONC.CFG', 'CONTROL', 'SETUP.CFG']]
file_path_list = [os.path.join(f'hysplit/working/cdumps/{region_of_interest}/', x) for x in file_path_list]
file_name_list = [os.path.basename(x) for x in file_path_list]
start_date_list = [datetime.datetime(int('20'+file_name[6:8]), int(file_name[8:10]), int(file_name[10:12]), hour=int(file_name[12:14])) for file_name in file_name_list]



# first is 2014
#conc_file_path_list = ['', '_2013', '_2016', '_2017']
#conc_file_path_list = ['C:/HYSPLIT/working/cdump' + x for x in conc_file_path_list]
#conc_file_path_list



conc_xarrays_list = [read_conc_file(x) for x in file_path_list]
# conc_xarrays_list_origin = conc_xarrays_list

#conc_xarrays_list = [x.chunk({'time': 'auto'}) for x in conc_xarrays_list_origin]
#conc_xarrays_list = [x.chunk({'x': 'auto', 'y': 'auto', 'source': 'auto'}) for x in conc_xarrays_list_origin]
#conc_xarrays_list = [x.chunk({'source': 'auto'}) for x in conc_xarrays_list_origin]

conc_xarrays_list = [x.chunk({'source': 'auto'}) for x in conc_xarrays_list]




import dask
dask.config.set({"array.slicing.split_large_chunks": True}) 




conc_xarray_final = xr.concat(conc_xarrays_list, pd.Index(start_date_list, name="start_date"))

del conc_xarrays_list
import gc
gc.collect()


import sys
# in GB
sys.getsizeof(conc_xarray_final)/1_000_000



 lon_vals, lat_vals = get_latlongrid_flat(conc_xarray_final, conc_xarray_final.x.values, conc_xarray_final.y.values)

 conc_xarray_final = conc_xarray_final.assign_coords(
             lon=(('x'), lon_vals ),
             lat=(('y'), lat_vals )
     ).set_index(x = ['lon'], y = ['lat'])


#conc_xarray_final = conc_xarray_final.rename({'x':'lon', 'y' :'lat'})



conc_xarray_final = conc_xarray_final.rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #
import pyproj
pyproj.show_versions()
from pyproj import Proj, transform
proj_4326 = Proj("epsg:4326")
crs_4326 = CRS("WGS84")
pyproj.datadir.get_data_dir()
from pyproj import CRS

crs=CRS('EPSG:4326')



#source_id = '0023'

#np.log(conc_xarray_final.sel(source=source_id, time=slice(0, None)).fillna(0).mean(dim='start_date').mean(dim='time')).plot()



population_raster = rioxarray.open_rasterio(r'data\population_grid\ppp_2018_1km_Aggregated.tif')

#population_raster.plot()

# look at what exact type of resampling are you using
pop_raster_matched = population_raster.rio.reproject_match(conc_xarray_final)
#pop_raster_matched.plot()

pop_raster_matched_masked = pop_raster_matched.where(pop_raster_matched >-0.1)
#pop_raster_matched_masked.plot()
#np.log(pop_raster_matched_masked).plot()


pop_raster_matched_masked_other = xr.where(pop_raster_matched >-0.1,pop_raster_matched, 0)
#conc_xarray_final
conc_xarray_final_masked = conc_xarray_final.where(pop_raster_matched >-0.1)



pop_exposure_matrix = conc_xarray_final_masked.fillna(0).mean(dim='start_date').mean(dim='time').dot(pop_raster_matched_masked.fillna(0), dims = ['x', 'y'])


pop_exposure_matrix.to_netcdf("data/intermed_files/pop_exposure_matrix_pakistan.nc")
pop_exposure_matrix.to_netcdf("data/intermed_files/pop_exposure_matrix_india.nc")


pop_exposure_matrix = xr.open_dataarray("data/intermed_files/pop_exposure_matrix_china.nc")
pop_exposure_matrix = xr.open_dataarray("data/intermed_files/pop_exposure_matrix_india.nc")


#pop_exposure_matrix.isel(band=0).set_index(source=['lon_source', 'lat_source'],).unstack("source").plot(x='lon_source', y='lat_source')

#pop_exposure_matrix.isel(band=0).plot.hist()

pop_exposure_by_lonlat = pop_exposure_matrix.isel(band=0).set_index(source=['lon_source', 'lat_source'],).unstack("source")

pop_exposure_by_lonlat = pop_exposure_by_lonlat.transpose('lat_source', 'lon_source')
pop_exposure_by_lonlat = pop_exposure_by_lonlat.rio.set_spatial_dims(x_dim = 'lon_source', y_dim='lat_source', inplace=False)
pop_exposure_by_lonlat = pop_exposure_by_lonlat.rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #


pop_exposure_by_lonlat.rio.to_raster('data/intermed_files/pop_exposure_matrix_pakistan_with_crs.tif')
pop_exposure_by_lonlat.rio.to_raster('data/intermed_files/pop_exposure_matrix_china_with_crs.tif')
pop_exposure_by_lonlat.rio.to_raster('data/intermed_files/pop_exposure_matrix_india_with_crs.tif')
pop_exposure_by_lonlat.rio.to_raster(f'data/intermed_files/pop_exposure_matrix_{region_of_interest}_deposition_with_crs.tif')



pop_exposure_by_lonlat.rio.to_raster('data/intermed_files/pop_exposure_matrix_india_72h_with_crs.tif')
pop_exposure_by_lonlat.rio.to_raster('data/intermed_files/pop_exposure_matrix_india_higher_res_with_crs.tif')



del pop_raster_matched_masked
del population_raster
import gc
gc.collect()

