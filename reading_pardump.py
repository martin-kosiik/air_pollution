import monetio
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


import monetio as mio

hysplitfile = 'C:/HYSPLIT/working/cdumps/cdump_06100100'

#hysplitfile = 'C:/HYSPLIT/working/SRM_source_conc.bin'
#one_source  = 'C:/HYSPLIT/working/conc_dump_grid/SRM_source_conc_x-86_y41.bin'




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
    llcrnr_lat = dset.attrs["Concentration Grid"]["llcrnr latitude"]
    llcrnr_lon = dset.attrs["Concentration Grid"]["llcrnr longitude"]
    nlat = dset.attrs["Concentration Grid"]["Number Lat Points"]
    nlon = dset.attrs["Concentration Grid"]["Number Lon Points"]
    dlat = dset.attrs["Concentration Grid"]["Latitude Spacing"]
    dlon = dset.attrs["Concentration Grid"]["Longitude Spacing"]

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



my_modelbin = mio.models.hysplit.ModelBin(hysplitfile, verbose=False)



dir(my_modelbin.dset.data_vars)

time_period = 0

my_modelbin.dset.data_vars['006F'][time_period][0].plot(x="longitude", y="latitude")
my_modelbin.dset.data_vars['0001'][time_period][0].plot(x="longitude", y="latitude")
my_modelbin.dset.data_vars['0079'][time_period][0].plot(x="longitude", y="latitude")



all_source_xarray = my_modelbin.dset.to_array(dim='source')
sources_hexdec_index = list(my_modelbin.dset.data_vars.keys())

all_source_xarray['source']

int(sources_hexdec_index[60], base=16)

sources_dec_index = [int(x, base=16) for x in sources_hexdec_index]

with open('config_files/CONTROL') as f:
    lines = f.readlines()



n_of_sources = int(lines[1].replace('\n', ''))

sources_coords_list = lines[2:(n_of_sources+2)]

sources_coords_list[0].split()
sources_coords_array = [x.split() for x in sources_coords_list]
sources_coords_array = np.array(sources_coords_array, dtype=np.float32)
sources_coords_array = sources_coords_array[:,0:2 ] # remove height colmun
sources_coords_df = pd.DataFrame(sources_coords_array, columns=['lat_source', 'lon_source'], index=sources_hexdec_index)

sources_coords_df.loc['0015']



all_source_xarray = all_source_xarray.assign_coords(
        lon_source=('source', sources_coords_df.loc[all_source_xarray['source']].lon_source),
        lat_source=('source', sources_coords_df.loc[all_source_xarray['source']].lat_source)
)


all_source_xarray = all_source_xarray.isel(z=0)

all_source_xarray.shape


start_date = lines[0].split()
start_date


start_year = int('20' + start_date[0])
start_month = int( start_date[1])
start_day =int( start_date[2])
start_hour =int( start_date[3])


start_date = datetime.datetime(start_year, start_month, start_day, hour=start_hour)
start_date







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


year_list = ['06', '07', '08', '09', '10', '11', '12', '15', '16', '17']
file_name_list = ['C:/HYSPLIT/working/cdumps/cdump_' + year + '10' + '01' + '00' for year in year_list]

# first is 2014
#conc_file_path_list = ['', '_2013', '_2016', '_2017']
#conc_file_path_list = ['C:/HYSPLIT/working/cdump' + x for x in conc_file_path_list]
#conc_file_path_list



conc_xarrays_list = [read_conc_file(x) for x in file_name_list]


start_year_list = [int('20'+year) for year in year_list]
start_date_list = [datetime.datetime(year, start_month, start_day, hour=start_hour) for year in start_year_list]


conc_xarray_final = xr.concat(conc_xarrays_list, pd.Index(start_date_list, name="start_date"))

del conc_xarrays_list
import gc
gc.collect()


import sys
# in GB
sys.getsizeof(conc_xarray_final)/1_000_000

#relative_times_array = (conc_xarrays_list[0].time.values - conc_xarrays_list[0].time.values[0]).astype('timedelta64[h]').astype(int)

#conc_xarrays_list[0] = conc_xarrays_list[0].assign_coords(
#        relative_time=(('time'), relative_times_array )
#).set_index(time = ['relative_time'])




#conc_xarrays_list[0].expand_dims({'new_start_time': 1}).assign_coords(new_start_time= (('new_start_time'), conc_xarrays_list[0].time.values ))




lon_vals, lat_vals = get_latlongrid_flat(conc_xarray_final, conc_xarray_final.x.values, conc_xarray_final.y.values)

conc_xarray_final = conc_xarray_final.assign_coords(
            lon=(('x'), lon_vals ),
            lat=(('y'), lat_vals )
    ).set_index(x = ['lon'], y = ['lat'])


#mgrid = get_latlongrid(conc_xarray_final, conc_xarray_final.x.values, conc_xarray_final.y.values)
#conc_xarray_final = conc_xarray_final.assign_coords(latitude=(("y", "x"), mgrid[1]))
#conc_xarray_final = conc_xarray_final.assign_coords(longitude=(("y", "x"), mgrid[0]))

# https://github.com/noaa-oar-arl/monetio/blob/master/monetio/models/hysplit.py

#years_array = conc_xarray_final.time.values.astype('datetime64[Y]').astype(int) + 1970


#start_datetime_array = [datetime.datetime(x, start_month, start_day, hour=start_hour) for x in years_array]
#start_datetime_array = np.array(start_datetime_array)

conc_xarray_final = conc_xarray_final.sel(y=slice(15, 35), x=slice(67,85))


conc_xarray_final = conc_xarray_final.rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #


conc_xarray_final









time_period = 0



#conc_xarray_final_zero = conc_xarray_final.fillna(0)


#conc_xarray_final.set_index(relative_time = ['relative_time']).sel(source='0001').sel(relative_time=0).mean(dim='start_date', skipna=True).mean(dim='time', skipna=True).plot(x="longitude", y="latitude")


#conc_xarray_final.sel(source='0001', time=12).mean(dim='start_date', skipna=True).plot(x="longitude", y="latitude")


source_id = '0023'

np.log(conc_xarray_final.sel(source=source_id, time=slice(0, None)).fillna(0).mean(dim='start_date').mean(dim='time')).plot()


conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=0).fillna(0).mean(dim='time').plot()


conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=1).fillna(0).mean(dim='time').plot()


conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=2).fillna(0).mean(dim='time').plot()

conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=3).fillna(0).mean(dim='time').plot()


import sys
# in GB
sys.getsizeof(conc_xarray_final)/1_000_000



#population_raster = rioxarray.open_rasterio( 'data\\population_grid\\india-spatial-india-census-2011-population-grid\\india_pop\\w001001.adf')
#population_raster = rioxarray.open_rasterio(r'data\population_grid\gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min.tif')
population_raster = rioxarray.open_rasterio(r'data\population_grid\ppp_2018_1km_Aggregated.tif')

#population_raster.plot()

# look at what exact type of resampling are you using
pop_raster_matched = population_raster.rio.reproject_match(conc_xarray_final)
pop_raster_matched.plot()

pop_raster_matched_masked = pop_raster_matched.where(pop_raster_matched >-0.1)
pop_raster_matched_masked.plot()
np.log(pop_raster_matched_masked).plot()

pop_raster_matched_masked_other
np.log(pop_raster_matched_masked_other).plot()

pop_raster_matched_masked_other = xr.where(pop_raster_matched >-0.1,pop_raster_matched, 0)
conc_xarray_final
conc_xarray_final_masked = conc_xarray_final.where(pop_raster_matched >-0.1)

conc_xarray_final_masked.sel(source=source_id, time=slice(None, 6)).fillna(0).mean(dim='start_date').mean(dim='time').plot()

del pop_raster_matched_masked
del conc_xarrays_list
del population_raster
import gc
gc.collect()

conc_xarray_final.dims
pop_exposure_matrix = conc_xarray_final_masked.fillna(0).mean(dim='start_date').mean(dim='time').dot(pop_raster_matched_masked.fillna(0), dims = ['x', 'y'])
#pop_exposure_matrix = conc_xarray_final.fillna(0).mean(dim='start_date').mean(dim='time').dot(pop_raster_matched_masked_other, dims = ['x', 'y'])

conc_xarray_final_masked

pop_exposure_matrix.to_netcdf("data/intermed_files/pop_exposure_matrix_10_years.nc")

pop_exposure_matrix = xr.open_dataarray("data/intermed_files/pop_exposure_matrix_10_years.nc")
pop_exposure_matrix.rio.crs

conc_xarray_final = conc_xarray_final.rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #


pop_exposure_matrix.isel(band=0).set_index(source=['lon_source', 'lat_source'],).unstack("source").plot(x='lon_source', y='lat_source')

pop_exposure_matrix.isel(band=0).plot.hist()

pop_exposure_by_lonlat = pop_exposure_matrix.isel(band=0).set_index(source=['lon_source', 'lat_source'],).unstack("source")

pop_exposure_by_lonlat.plot(x='lon_source', y='lat_source')

pop_exposure_by_lonlat = pop_exposure_by_lonlat.set_index(x = 'lon_source', y='lat_source').rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #
#pop_exposure_by_lonlat = pop_exposure_by_lonlat.rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #

pop_exposure_by_lonlat = pop_exposure_by_lonlat.set_index(x = ['lon_source'], y='lat_source')
pop_exposure_by_lonlat = pop_exposure_by_lonlat.sortby(["y", "x"]) # THIS SOLVED THE PROBLEM
pop_exposure_by_lonlat.rio.write_crs("epsg:4326", inplace=True)


pop_exposure_by_lonlat = pop_exposure_by_lonlat.transpose('lat_source', 'lon_source')
pop_exposure_by_lonlat = pop_exposure_by_lonlat.rio.set_spatial_dims(x_dim = 'lon_source', y_dim='lat_source', inplace=False)
pop_exposure_by_lonlat = pop_exposure_by_lonlat.rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #


pop_exposure_by_lonlat.rio.to_raster('data/intermed_files/pop_exposure_matrix_10_years_with_crs.tif')
pop_exposure_by_lonlat.rio.crs
pop_exposure_by_lonlat.rio.info
pop_exposure_by_lonlat.plot()


pop_exposure_by_lonlat.to_netcdf("data/intermed_files/pop_exposure_matrix_10_years_with_crs.nc")

pop_exposure_matrix = xr.open_dataarray("data/intermed_files/pop_exposure_matrix_10_years.nc")


import cartopy

import cartopy.crs as ccrs




import cartopy.io.img_tiles as cimgt

fig = plt.figure()

# create geo axes
projection = ccrs.PlateCarree()
geo_axes = plt.subplot(projection=projection)

# add open street map background
# when commenting the two following lines, the data array is plotted correctly
osm_background = cimgt.OSM()
geo_axes.add_image(osm_background, 10)
#projection.states()
# plot dataset
xr.plot.imshow(
    darray=pop_exposure_by_lonlat,
    x="lon_source",
    y="lat_source",
    ax=geo_axes,
    transform=projection,
    zorder=10
)

# show plot
plt.show()








# Now we have to do resampling of the rasters to allign the population and pollution exposure rasters
# https://pygis.io/docs/e_raster_resample.html
# we probably want to use average resampling
# https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html

lon_vals, lat_vals = get_latlongrid_flat(conc_xarray_final, conc_xarray_final.x.values, conc_xarray_final.y.values)


new_xarray = conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=3).fillna(0).mean(dim='time')
new_xarray.attrs = conc_xarray_final.attrs

lon_vals, lat_vals = get_latlongrid_flat(new_xarray, new_xarray.x.values, new_xarray.y.values)

len(new_xarray.x.values)
len(lon_vals)
len(lat_vals)
len(new_xarray.y.values)

new_xarray = new_xarray.assign_coords(
            lon=(('x'), lon_vals ),
            lat=(('y'), lat_vals )
    ).set_index(x = ['lon'], y = ['lat'])



new_xarray
lon_vals



new_xarray = new_xarray.rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #

get_latlongrid_flat()

np.log(new_xarray).plot()

new_xarray.rio.crs

new_xarray.rio.to_raster("data/planet_scope_green.tif")

!rio info data/planet_scope_green.tif





new_xarray.rio.bounds()

set_crs
rasterio.crs.CRS.from_epsg(4326)
new_xarray.rio.height
new_xarray.rio.width
new_xarray.rio._x_dim
new_xarray.rio.set_spatial_dims(x_dim='longitude', y_dim='latitude', inplace=False)

new_xarray.rio.transform()

new_xarray.rio.grid_mapping
new_xarray.encoding
rds.rio.crs
rds.rio.bounds()
rds.rio.grid_mapping
rds.rio.height
rds.rio.width
rds.rio._x_dim
set_spatial_dims()
rds.rio.transform()
rds.encoding


def print_raster(raster):
    print(
        f"shape: {raster.rio.shape}\n"
        f"resolution: {raster.rio.resolution()}\n"
        f"bounds: {raster.rio.bounds()}\n"
        f"sum: {raster.sum().item()}\n"
        f"CRS: {raster.rio.crs}\n"
    )

print_raster(rds)

print_raster(new_xarray)


rds_repr_match = rds.rio.reproject_match(new_xarray)

print_raster(rds_repr_match)

rds_repr_match.plot()

rds_repr_match_masked = rds_repr_match.where(rds_repr_match >0)

np.log(rds_repr_match_masked).plot()


new_xarray_masked = new_xarray.where(rds_repr_match >0)

new_xarray_masked.plot()

rds_repr_match_masked

new_xarray_masked.fillna(0).dot(rds_repr_match_masked.fillna(0), dims = ['x', 'y'])



crs = CRS.from_string ('+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')


rds = rioxarray.open_rasterio(
    'data\\population_grid\\india-spatial-india-census-2011-population-grid\\india_pop\\w001001.adf',
)
rds

rds.isel(x =1)
rds.rio.set_attrs({'_FillValue': 0})

rds = rds.rio.set_nodata(0)
rds.rio.nodata
rds.rio.encoded_nodata
rds.rio.write_nodata(False).isel(x=0)


conc_xarray_final.where(conc_xarray_final == np.nan, other=0)


rds[rds < 0]
rds.where(rds < 0, other=0)

dir(rds.rio)

rds.rio.estimate_utm_crs()


population_raster = rasterio.open('data\\population_grid\\india-spatial-india-census-2011-population-grid\\india_pop\\w001001.adf')

population_raster.nodata
population_raster.nodatavals
population_raster.profile
population_raster.dataset_mask()
population_raster.colorinterp


conc_xarray_final
conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=3).fillna(0).mean(dim='time').rio.crs
plot(population_raster)

rasterio.CRS.from_epsg(4326)


from matplotlib import pyplot
population_raster.read(1)
population_raster.bounds
population_raster.crs

0.008333333333*6

population_raster.res

pyplot.imshow(population_raster.read(1), cmap='viridis')
pyplot.show()

from rasterio.plot import show_hist
from rasterio.windows import Window


crop_window = Window.from_slices((0, 1500), (0, 1500))


show_hist( population_raster, bins=50, lw=0.0, stacked=False, alpha=0.3,  histtype='stepfilled', title="Histogram")

with rasterio.open('data\\population_grid\\india-spatial-india-census-2011-population-grid\\india_pop\\w001001.adf') as src:
        pop_raster = src.read(1, window=crop_window)

pop_raster

pyplot.imshow(pop_raster, cmap='viridis')

pop_raster

from rasterio.plot import show
show(population_raster)
show(pop_raster)


show_hist( pop_raster, bins=50, lw=0.0, stacked=False, alpha=0.3,  histtype='stepfilled', title="Histogram")


conc_xarray_final



#!pip install rioxarray

import rioxarray








import geopandas
my_modelbin.makegrid([300, 301], [300, 301])
my_modelbin.parse_header(hdata1 = my_modelbin)

print(mio.models.__file__)

massload =  mio.models.hysplit.hysp_massload(hxr2, threshold=0, mult=1e10)

massload =  mio.models.hysplit.hysp_heights(hxr, threshold=0, height_mult=1/1000.0, mult=1e10, mass_load=False)





from numpy import fromfile, arange

fid = open(hysplitfile_2, "rb")

fid

lines = fid.readlines()

lines










hxr = mio.models.hysplit.open_dataset(hysplitfile)
dir(mio.models)

hysplit_pardump_file = 'C:/HYSPLIT/working/PARDUMP'

par_dump_alt = mio.models.hysplit.open_dataset(hysplit_pardump_file)




#par_dump = mio.models.pardump.Pardump(hysplit_pardump_file)
#par_dump = Pardump(hysplit_pardump_file)


#read_par_dump = par_dump.read(century=1900)

import datetime

import numpy as np
import pandas as pd


# this might be interesting to explore
# https://www.ready.noaa.gov/hysplitusersguide/S350.htm
# https://github.com/noaa-oar-arl/utilhysplit/blob/master/utilhysplit/par2conc.py

# from montielio
class Pardump():
    """methods for writing and reading a pardump file.
       __init__  initializes structure of binary file.
       write   writes a pardump file.
       read    reads a pardump file. returns a dictionary.
               Keys are the date of the particle positions in YYMMDDHH.
               Values are pandas dataframe objects with the particle
               information.
    """

    def __init__(self, fname='PARINIT'):
        """
        ##initializes structures which correspond to the binary records.
        ##'p' variables coorespond to padding that fortran adds.
        """
        self.fname = fname
        self.dtfmt = "%Y%m%d%H%M"

        tp1 = '>f'  # big endian float.
        tp2 = '>i'  # big endian integer.

        # header record in fortran file.
        self.hdr_dt = np.dtype([('padding', tp2),
                                ('parnum', tp2),
                                ('pollnum', tp2),
                                ('year', tp2),
                                ('month', tp2),
                                ('day', tp2),
                                ('hour', tp2),
                                ('minute', tp2)
                                ])

        # data record in fortran file.
        self.pardt = np.dtype([('p1', tp2),
                               ('p2', tp2),
                               ('pmass', tp1),
                               ('p3', '>l'),
                               ('lat', tp1),
                               ('lon', tp1),
                               ('ht', tp1),
                               ('su', tp1),
                               ('sv', tp1),
                               ('sx', tp1),
                               ('p4', '>l'),
                               ('age', tp2),
                               ('dist', tp2),
                               ('poll', tp2),
                               ('mgrid', tp2),
                               ('sorti', tp2)])

    def write(self, numpar, pmass, lon, lat, ht, pollnum, sdate):
        """
        numpar : integer
                 number of particles
        pmass   : lists or numpy array
        lon     : list or numpy array
        lat     : list or numpy array
        ht      : list or numpy array
        pollnum : integer
                  pollutant index
        sdate   : datetime.datetime object
        """
        with open(self.fname, 'wb') as fp:
            pad1 = np.ones(numpar) * 28
            pad2 = np.ones(numpar) * 4
            pad3 = np.ones(numpar) * 17179869208
            pad4 = np.ones(numpar) * 103079215124.0
            zrs = np.zeros(numpar)
            ones = np.ones(numpar)
            sorti = np.arange(1, numpar + 1)

            print(len(lon))
            a = np.zeros((numpar, ), dtype=self.pardt)
            a['p1'] = pad1
            a['p2'] = pad2
            a['p3'] = pad3
            a['p4'] = pad4

            a['lat'] = lat
            a['lon'] = lon
            a['ht'] = ht
            a['pmass'] = pmass
            a['poll'] = pollnum

            a['age'] = zrs
            a['dist'] = zrs

            a['mgrid'] = ones
            a['sorti'] = sorti

            # print a

            hdr = np.zeros((1, ), dtype=self.hdr_dt)
            hdr['padding'] = 28
            hdr['parnum'] = numpar
            hdr['pollnum'] = 1
            hdr['year'] = sdate.year
            hdr['month'] = sdate.month
            hdr['day'] = sdate.day
            hdr['hour'] = sdate.hour
            hdr['minute'] = sdate.minute
            print(hdr)

            endrec = np.array([20], dtype='>i')

            fp.write(hdr)
            fp.write(a)
            fp.write(endrec)
            fp.write(endrec)

    # def writeascii(self, drange=[], verbose=1, century=2000, sorti=[]):
    #    read(self, drange=[], verbose=1, century=2000, sorti=[]):



