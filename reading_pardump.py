import monetio
import os
import numpy as np
import pandas as pd
import datetime
import xarray as xr


#!pip install git+https://github.com/noaa-oar-arl/monetio.git

os.chdir(r'C:\Users\marti\Dropbox\research_projects\air_pollution')


import monetio as mio

hysplitfile = 'C:/HYSPLIT/working/cdump'
hysplitfile = 'C:/HYSPLIT/working/SRM_cdump_new'

hysplitfile = 'C:/HYSPLIT/working/SRM_source_conc.bin'
one_source  = 'C:/HYSPLIT/working/conc_dump_grid/SRM_source_conc_x-86_y41.bin'



my_modelbin = mio.models.hysplit.ModelBin(hysplitfile, verbose=False)

my_modelbin.atthash['Starting Locations']
my_modelbin.atthash
my_modelbin.nonzeroconcdates
my_modelbin.dset

my_modelbin.slat
my_modelbin.sht
my_modelbin.slon
my_modelbin.llcrnr_lon
my_modelbin.nlon
my_modelbin.nlat
my_modelbin.dlat
my_modelbin.dlon
my_modelbin.levels
my_modelbin.verbose




dir(my_modelbin.dset.data_vars)

time_period = 5

my_modelbin.dset.data_vars['006F'][time_period][0].plot(x="longitude", y="latitude")


my_modelbin.dset.data_vars['0001'][time_period][0].plot(x="longitude", y="latitude")

my_modelbin.dset.data_vars['0079'][time_period][0].plot(x="longitude", y="latitude")



all_source_xarray = my_modelbin.dset.to_array(dim='source')

sources_hexdec_index = list(my_modelbin.dset.data_vars.keys())

sources_hexdec_index

int(sources_hexdec_index[60], base=16)

sources_dec_index = [int(x, base=16) for x in sources_hexdec_index]
sources_dec_index

with open('C:/HYSPLIT/working/CONTROL') as f:
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




conc_file_path_list[1]


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

    return all_source_xarray





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




# first is 2014
conc_file_path_list = ['', '_2013', '_2016', '_2017']
conc_file_path_list = ['C:/HYSPLIT/working/cdump' + x for x in conc_file_path_list]
conc_file_path_list



conc_xarrays_list = [read_conc_file(x) for x in conc_file_path_list]



conc_xarray_final.sel(start_date = datetime.datetime(2014, start_month, start_day, hour=start_hour))

conc_xarray_final = xr.concat(conc_xarrays_list, pd.Index([datetime.datetime(2014, start_month, start_day, hour=start_hour),
                                                            datetime.datetime(2013, start_month, start_day, hour=start_hour),
                                                            datetime.datetime(2016, start_month, start_day, hour=start_hour),
                                                            datetime.datetime(2017, start_month, start_day, hour=start_hour)], name="start_date"))


relative_times_array = (conc_xarrays_list[0].time.values - conc_xarrays_list[0].time.values[0]).astype('timedelta64[h]').astype(int)


conc_xarrays_list[0] = conc_xarrays_list[0].assign_coords(
        relative_time=(('time'), relative_times_array )
).set_index(time = ['relative_time'])




conc_xarrays_list[0].expand_dims({'new_start_time': 1}).assign_coords(new_start_time= (('new_start_time'), conc_xarrays_list[0].time.values ))
xr.concat()

conc_xarrays_list[0].assign_coords({'new_start_time': conc_xarrays_list[0].time.values[0]})



mgrid = get_latlongrid(conc_xarray_final, conc_xarray_final.x.values, conc_xarray_final.y.values)
#newhxr = newhxr.drop("longitude")
#newhxr = newhxr.drop("latitude")
conc_xarray_final = conc_xarray_final.assign_coords(latitude=(("y", "x"), mgrid[1]))
conc_xarray_final = conc_xarray_final.assign_coords(longitude=(("y", "x"), mgrid[0]))

# https://github.com/noaa-oar-arl/monetio/blob/master/monetio/models/hysplit.py
(conc_xarray_final.time[1] - conc_xarray_final.start_date[1]).values.astype('timedelta64[h]')

conc_xarray_final.time[1].values.year

years_array = conc_xarray_final.time.values.astype('datetime64[Y]').astype(int) + 1970


start_datetime_array = [datetime.datetime(x, start_month, start_day, hour=start_hour) for x in years_array]
start_datetime_array = np.array(start_datetime_array)


(conc_xarray_final.time - conc_xarray_final.start_date).astype('timedelta64[h]')

conc_xarray_final.sel(start_date = '2016-10-01T00:00:00.000000000', time='2013-10-01T00:00:00.000000000')


(conc_xarray_final.time - start_datetime_array.astype('datetime64')).values.astype('timedelta64[h]').astype(int)
conc_xarray_final = conc_xarray_final.assign_coords(
        relative_time=(('time'), (conc_xarray_final.time - start_datetime_array.astype('datetime64')).values.astype('timedelta64[h]').astype(int))
)




conc_xarray_final = conc_xarray_final.assign_coords(
        relative_time=(('time'), (conc_xarray_final['time'] - conc_xarray_final['start_date']).astype('timedelta64[h]').astype(int) )
)

conc_xarray_final.time - np.repeat(conc_xarray_final.start_date.values, 6)

conc_xarray_final.start_date.values

conc_xarray_final = conc_xarray_final.assign_coords(
        {'relative_time':conc_xarray_final.time - conc_xarray_final.start_date}
)

np.repeat(conc_xarray_final.start_date.values, 6)



time_period = 0

conc_xarray_final.assign_coords({"relative_time": (((da.lon + 180) % 360) - 180)})
conc_xarray_final.sel(source='0001')

conc_xarray_final.sel(relative_time=43200000000000).sel(source='0001').mean(dims='').plot(x="longitude", y="latitude")


conc_xarray_final_zero = conc_xarray_final.fillna(0)

conc_xarray_final.where(conc_xarray_final.relative_time==0, drop=True)

conc_xarray_final.set_index(relative_time = ['relative_time']).sel(source='0001').sel(relative_time=0).mean(dim='start_date', skipna=True).mean(dim='time', skipna=True).plot(x="longitude", y="latitude")


conc_xarray_final.sel(source='0001', time=12).mean(dim='start_date', skipna=True).plot(x="longitude", y="latitude")


source_id = '0072'

conc_xarray_final.sel(source=source_id, time=slice(1, None)).fillna(0).mean(dim='start_date').mean(dim='time').plot(x="longitude", y="latitude")


conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=0).fillna(0).mean(dim='time').plot(x="longitude", y="latitude")


conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=1).fillna(0).mean(dim='time').plot(x="longitude", y="latitude")


conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=2).fillna(0).mean(dim='time').plot(x="longitude", y="latitude")

conc_xarray_final.sel(source=source_id, time=slice(1, None)).isel(start_date=3).fillna(0).mean(dim='time').plot(x="longitude", y="latitude")




# Now we have to do resampling of the rasters to allign the population and pollution exposure rasters
# https://pygis.io/docs/e_raster_resample.html
# we probably want to use average resampling
# https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html











conc_xarray_final.sel(source='0072', time=60).fillna(0).mean(dim='start_date').plot(x="longitude", y="latitude")


conc_xarray_final.sel(source='0001').mean(dim='time', skipna=True).mean(dim='start_date', skipna=True).plot(x="longitude", y="latitude")

conc_xarray_final.sel(source='0079').mean(dim='time', skipna=True).mean(dim='start_date', skipna=True).plot(x="longitude", y="latitude")



conc_xarray_final_zero.sel(source='0001').mean(dim='time').mean(dim='start_date').plot(x="longitude", y="latitude")



my_modelbin.dset.data_vars['0001'][time_period][0].plot(x="longitude", y="latitude")

my_modelbin.dset.data_vars['0079'][time_period][0].plot(x="longitude", y="latitude")





conc_xarray_final








def get_latlongrid(dset, xindx, yindx):
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
    print(lat)
    print(lon)
    lonlist = [lon[x - 1] for x in xindx]
    latlist = [lat[x - 1] for x in yindx]
    mgrid = np.meshgrid(lonlist, latlist)
    return mgrid







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


        recs = define_struct()
        rec1, rec2, rec3, rec4a = recs[0], recs[1], recs[2], recs[3]
        rec4b, rec5a, rec5b, rec5c = recs[4], recs[5], recs[6], recs[7]
        rec6, rec8a, rec8b, rec8c = recs[8], recs[9], recs[10], recs[11]
        # rec7 = rec6
        # start_loc in rec1 tell how many rec there are.
        tempzeroconcdates = []
        # Reads header data. This consists of records 1-5.
        hdata1 = fromfile(fid, dtype=rec1, count=1)
        nstartloc = parse_header(hdata1)

        hdata2 = fromfile(fid, dtype=rec2, count=nstartloc)
        century = parse_hdata2(hdata2, nstartloc, 1900)
        recs = self.define_struct()
        rec1, rec2, rec3, rec4a = recs[0], recs[1], recs[2], recs[3]
        rec4b, rec5a, rec5b, rec5c = recs[4], recs[5], recs[6], recs[7]
        rec6, rec8a, rec8b, rec8c = recs[8], recs[9], recs[10], recs[11]
        # rec7 = rec6
        # start_loc in rec1 tell how many rec there are.
        tempzeroconcdates = []
        # Reads header data. This consists of records 1-5.
        hdata1 = fromfile(fid, dtype=rec1, count=1)
        nstartloc = self.parse_header(hdata1)

        hdata2 = fromfile(fid, dtype=rec2, count=nstartloc)
        century = self.parse_hdata2(hdata2, nstartloc, century)












    def define_struct():
        """Each record in the fortran binary begins and ends with 4 bytes which
        specify the length of the record. These bytes are called pad below.
        They are not used here, but are thrown out. The following block defines
        a numpy dtype object for each record in the binary file. """
        from numpy import dtype

        real4 = ">f"
        int4 = ">i"
        int2 = ">i2"
        char4 = ">a4"

        rec1 = dtype(
            [
                ("pad1", int4),
                ("model_id", char4),  # meteorological model id
                ("met_year", int4),  # meteorological model starting time
                ("met_month", int4),
                ("met_day", int4),
                ("met_hr", int4),
                ("met_fhr", int4),  # forecast hour
                ("start_loc", int4),  # number of starting locations
                ("conc_pack", int4),  # concentration packing flag (0=no, 1=yes)
                ("pad2", int4),
            ]
        )

        # start_loc in rec1 tell how many rec there are.
        rec2 = dtype(
            [
                ("pad1", int4),
                ("r_year", int4),  # release starting time
                ("r_month", int4),
                ("r_day", int4),
                ("r_hr", int4),
                ("s_lat", real4),  # Release location
                ("s_lon", real4),
                ("s_ht", real4),
                ("r_min", int4),  # release startime time (minutes)
                ("pad2", int4),
            ]
        )

        rec3 = dtype(
            [
                ("pad1", int4),
                ("nlat", int4),
                ("nlon", int4),
                ("dlat", real4),
                ("dlon", real4),
                ("llcrnr_lat", real4),
                ("llcrnr_lon", real4),
                ("pad2", int4),
            ]
        )

        rec4a = dtype(
            [
                ("pad1", int4),
                ("nlev", int4),  # number of vertical levels in concentration grid
            ]
        )

        rec4b = dtype([("levht", int4)])  # height of each level (meters above ground)

        rec5a = dtype(
            [
                ("pad1", int4),
                ("pad2", int4),
                ("pollnum", int4),  # number of different pollutants
            ]
        )

        rec5b = dtype([("pname", char4)])  # identification string for each pollutant

        rec5c = dtype([("pad2", int4)])

        rec6 = dtype(
            [
                ("pad1", int4),
                ("oyear", int4),  # sample start time.
                ("omonth", int4),
                ("oday", int4),
                ("ohr", int4),
                ("omin", int4),
                ("oforecast", int4),
                ("pad3", int4),
            ]
        )

        # rec7 has same form as rec6.            #sample stop time.

        # record 8 is pollutant type identification string, output level.

        rec8a = dtype(
            [
                ("pad1", int4),
                ("poll", char4),  # pollutant identification string
                ("lev", int4),
                ("ne", int4),  # number of elements
            ]
        )

        rec8b = dtype(
            [
                ("indx", int2),  # longitude index
                ("jndx", int2),  # latitude index
                ("conc", real4),
            ]
        )

        rec8c = dtype([("pad2", int4)])
        recs = (
            rec1,
            rec2,
            rec3,
            rec4a,
            rec4b,
            rec5a,
            rec5b,
            rec5c,
            rec6,
            rec8a,
            rec8b,
            rec8c,
        )
        return recs

    def parse_header(hdata1):
        """
        hdata1 : dtype
        Returns
        nstartloc : int
           number of starting locations in file.
        """
        if len(hdata1["start_loc"]) != 1:
            print(
                "WARNING in ModelBin _readfile - number of starting locations "
                "incorrect"
            )
            print(hdata1["start_loc"])
        # in python 3 np.fromfile reads the record into a list even if it is
        # just one number.
        # so if the length of this record is greater than one something is
        # wrong.
        nstartloc = hdata1["start_loc"][0]
        #self.atthash["Meteorological Model ID"] = hdata1["model_id"][0].decode("UTF-8")
        #self.atthash["Number Start Locations"] = nstartloc
        return nstartloc

    def parse_hdata2(hdata2, nstartloc, century):

        # Loop through starting locations
        for nnn in range(0, nstartloc):
            # create list of starting latitudes, longitudes and heights.
            slat.append(hdata2["s_lat"][nnn])
            slon.append(hdata2["s_lon"][nnn])
            self.sht.append(hdata2["s_ht"][nnn])
            self.atthash["Starting Locations"].append(
                (hdata2["s_lat"][nnn], hdata2["s_lon"][nnn])
            )

            # try to guess century if century not given
            if century is None:
                if hdata2["r_year"][0] < 50:
                    century = 2000
                else:
                    century = 1900
                print(
                    "WARNING: Guessing Century for HYSPLIT concentration file", century
                )
            # add sourcedate which is datetime.datetime object
            sourcedate = datetime.datetime(
                century + hdata2["r_year"][nnn],
                hdata2["r_month"][nnn],
                hdata2["r_day"][nnn],
                hdata2["r_hr"][nnn],
                hdata2["r_min"][nnn],
            )
            self.sourcedate.append(sourcedate)
            self.atthash["Source Date"].append(sourcedate)
            return century





hxr = mio.models.hysplit.open_dataset(hysplitfile)
dir(mio.models)

hysplit_pardump_file = 'C:/HYSPLIT/working/PARDUMP'

par_dump_alt = mio.models.hysplit.open_dataset(hysplit_pardump_file)




par_dump = mio.models.pardump.Pardump(hysplit_pardump_file)
par_dump = Pardump(hysplit_pardump_file)


read_par_dump = par_dump.read(century=1900)

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

    def read(self, drange=None, verbose=1, century=2000, sorti=None):
        """
        daterange should be a list of two datetime.datetime objects
        indicating the beginning
        and ending date of the particle positions of interest.
        Returns
           pframehash :  dictionary.
        The key is the date of the particle positions in YYMMDDHH.
        The value is a pandas dataframe object with the particle information.
        sorti is a list of sort indices. If sorti not None then will
        only return particles
        with those sort indices.
        nsort keeps track of which particle it is throughout the time.
        Could use this to keep track of initial height and time of release.
        """

        imax = 100
        # returns a dictionary of pandas dataframes. Date valid is the key.
        pframe_hash = {}
        with open(self.fname, 'rb') as fp:
            i = 0
            testf = True
            while testf:
                hdata = np.fromfile(fp, dtype=self.hdr_dt, count=1)
                if verbose:
                    print('Record Header ', hdata)
                if not hdata:
                    print('Done reading ', self.fname)
                    break
                if hdata['year'] < 1000:
                    year = hdata['year'] + century
                else:
                    year = hdata['year']
                pdate = datetime.datetime(
                    year,
                    hdata['month'],
                    hdata['day'],
                    hdata['hour'],
                    hdata['minute'])
                # if drange==[]:
                #   drange = [pdate, pdate]
                parnum = hdata['parnum']
                data = np.fromfile(fp, dtype=self.pardt, count=parnum[0])
                # n = parnum - 1
                # padding at end of each record
                np.fromfile(fp, dtype='>i', count=1)
                if verbose:
                    print('Date ', pdate, ' **** ', drange)

                testdate = False
                if not drange:
                    testdate = True
                elif pdate >= drange[0] and pdate <= drange[1]:
                    testdate = True

               # Only store data if it is in the daterange specified.
                if testdate:
                    print('Adding data ', hdata, pdate)
                    # otherwise get endian error message when create dataframe.
                    ndata = data.byteswap().newbyteorder()
                    par_frame = pd.DataFrame.from_records(
                        ndata)  # create data frame
                    # drop the fields which were padding
                    par_frame.drop(['p1', 'p2', 'p3', 'p4'],
                                   inplace=True,
                                   axis=1)
                    par_frame.drop(['su', 'sv', 'sx', 'mgrid'],
                                   inplace=True,
                                   axis=1)  # drop other fields
                    # drop where the lat field is 0. because
                    par_frame = par_frame.loc[par_frame['lat'] != 0]
                    # in pardump file particles which have not been
                    # released yet

                    if sorti:
                        # returns only particles with
                        par_frame = par_frame.loc[par_frame['sorti'].isin(
                            sorti)]
                        # sort index in list sorti
                    par_frame['date'] = pdate
                    # par_frame.sort('ht', inplace=True)  # sort by height
                    par_frame = pd.concat(
                        [par_frame], keys=[
                            self.fname])  # add a filename key
                    # create dictionary key for output.
                    datekey = pdate.strftime(self.dtfmt)
                    # Add value to dictionary.
                    pframe_hash[datekey] = par_frame

                # Assume data is written sequentially by date.
                i += 1

                if drange:
                    if pdate > drange[1]:
                        testf = False
                        if verbose:
                            print("Past date. Closing file.", drange[1], pdate)
                    # elif  drange[0] < pdate:
                    #   testf=False
                    #   if verbose:
                    #      print "Before date. Closing file"
                if i > imax:
                    print('Read pardump. Limited to 100 iterations. Stopping')
                    testf = False
        return pframe_hash




import monetio.models.pardump as pardump

par_dump_2 = pardump.open_dataset(hysplit_pardump_file, century=1900)


read_par_dump = par_dump.read(century=1900)

# https://hysplitbbs.arl.noaa.gov/viewtopic.php?f=28&t=1245
# https://www.ready.noaa.gov/hysplitusersguide/S350.htm

class VolcPar:
    def __init__(self, fdir="./", fname="PARDUMP.A"):
        fname = fname
        self.tname = os.path.join(fdir, fname)
        self.strfmt = "%Y%m%d%H%M"
        self.start = datetime.datetime.now()
        self.delta = datetime.timedelta(hours=1)
        self.delt = 5
        self.ymax = 9
        self.pdict = {}

    def read_pardump(self, drange=None, century=2000):
        self.df = pardump.open_dataset(fname=self.tname, drange=drange, century=century)
        # pd = pardump.Pardump(fname=self.tname)
        # self.df = pd.read(century=century)

    # def key2time(self):
    #    datelist=[]
    #    for key in self.pdict.keys:
    #        datelist.append(datetime.datetime.strptime(key, self.strfmt))
    #    #self.datetlist = datelist
    #    return datelist

    def getbytime(self,time):
        # returns a dataframe for that time.
        dstr = time.strftime(self.strfmt)
        return self.pdict[dstr]

    def findsource(self, sorti):
        done = False
        iii = 0
        while not done:
            d1 = self.start + iii * self.delta
            print("find source", d1)
            df1 = self.getbytime(d1)
            print("Ages", df1["age"].unique())
            # keep only the new particles.
            df1 = df1[df1["age"] == self.delt]

            # if no particles released then
            if df1.empty and iii != 0:
                done = True
                print("empty", iii)
                continue
            df1 = df1[df1["sorti"].isin(sorti)]
            if iii == 0:
                dfsource = df1
            else:
                dfsource = pd.concat([dfsource, df1])
            iii += 1
        return dfsource

    def plotsource(self, dfhash):
        x = []
        y = []
        for key in dfhash.keys():
            sorti = dfhash[key]["sorti"]
            dfsource = self.findsource(sorti)
            x.extend(dfsource["date"])
            y.extend(dfsource["ht"])
        x2 = time2int(x)
        # put time into minutes since start
        x2 = np.array(x2) / 60.0
        # put height into km
        y = np.array(y) / 1000.0
        xbins = np.arange(0, 200, 5)
        ybins = np.arange(0, 4 * 9) / 4.0
        cb = plt.hist2d(x2, y, bins=[xbins, ybins])
        plt.colorbar(cb[3])
        # sns.heatmap(x,y)
        return x2, y

def time2int(timelist):
    newlist = []
    tmin = np.min(timelist)
    for ttt in timelist:
        val = ttt - tmin
        newlist.append(val.seconds)
    return newlist


def average_mfitlist(mfitlist, dd=None, dh=None, buf=None, lat=None, lon=None, ht=None):
    """
    mfitlist : list of MassFit objects
    returns xarray DataArray
    """
    logger.debug("Running average_mfitlist in par2conc")
    concra = combine_mfitlist(mfitlist, dd, dh, buf, lat, lon, ht)
    concra = concra.mean(dim="time")
    return concra


def combine_mfitlist(
    mfitlist, dd=None, dh=None, buf=None, lat=None, lon=None, ht=None,
):
    """
    mfitlist : list of MassFit objects.
    finds concentrations from each fit and combines into
    one xarray along dimension called 'time'. Although in
    some cases that dimension may represent something other than time.
    e.g. ensemble member number.
    returns xarray DataArray
    """
    logger.debug("Running combine_mfitlist in par2conc")
    iii = 0
    concra = xr.DataArray(None)
    templist = []
    # conclist = []
    minlat = 90
    minlon = 180
    maxlat = -90
    maxlon = -180

    # fit all time periods and put arrays in a list.
    # keep track of range of latitude and longitude so
    # lat-lon grid can be created later.
    for mfit in mfitlist:
        latra, lonra, htra = mfit.get_grid(dd, dh, buf, lat, lon, ht)
        conc = mfit.get_conc2(dd=dd, dh=dh, latra=latra, lonra=lonra, htra=htra)
        minlat = np.min([minlat, np.min(latra)])
        minlon = np.min([minlon, np.min(lonra)])
        maxlat = np.max([maxlat, np.max(latra)])
        maxlon = np.max([maxlon, np.max(lonra)])
        templist.append(conc)
    # for xr.align to work properly, the coordinates
    # need to be integers.
    # This list comprehension changes the lat-lon coordinates to ints.
    # by applying the reindex function to all xarrays in templist.
    # re-index all the arrays to the largest grid.
    nlat = np.abs(np.ceil((maxlat - minlat) / dd)) + 1
    nlon = np.abs(np.ceil((maxlon - minlon) / dd)) + 1
    conclist = [reindex(x, minlat, minlon, nlat, nlon, dd, dd) for x in templist]
    # for conc in templist:
    # nlat = np.abs(np.ceil((maxlat - minlat) / dd)) + 1
    # nlon = np.abs(np.ceil((maxlon - minlon) / dd)) + 1
    #    conc2 = reindex(conc, minlat, minlon, nlat, nlon, dd, dd)
    #    conclist.append(conc2)
    #    print('MFIT', conc2)
    # create large xarray to align to.
    iii = 0
    for conc in conclist:
        if iii == 0:
            xnew = conc.copy()
        else:
            a, xnew = xr.align(conc, xnew, join="outer")
        iii += 1

    # align all arrays to largest one.
    iii = 0
    templist = []
    for temp in conclist:
        aaa, bbb = xr.align(temp, xnew, join="outer")
        aaa = aaa.fillna(0)
        aaa.expand_dims("time")
        aaa["time"] = iii
        templist.append(aaa)
        iii += 1

    # concatenate the aligned arrays.
    concra = xr.concat(templist, "time")

    # fill nans with 0
    concra = concra.fillna(0)

    # add time dimesion if not there.
    # if 'time' not in concra.dims:
    #    concra = concra.expand_dims('time')

    # add the lat lon coordinates back in
    concra = concra.drop("latitude")
    concra = concra.drop("longitude")

    # np.clip was added because sometimes due to the floating point
    # arithmetic arange returns an array with an extra number.
    # clip replaces any numbers larger than the max with the max value.
    # then remove duplicate values at the end.
    latra = np.clip(np.arange(minlat, maxlat + dd, dd), None, maxlat)
    lonra = np.clip(np.arange(minlon, maxlon + dd, dd), None, maxlon)
    if latra[-1] == latra[-2]:
        latra = latra[0:-1]
    if lonra[-1] == lonra[-2]:
        lonra = lonra[0:-1]
    mgrid = np.meshgrid(lonra, latra)

    concra = concra.assign_coords(longitude=(("y", "x"), mgrid[0]))
    concra = concra.assign_coords(latitude=(("y", "x"), mgrid[1]))

    return concra


def height_correction(zaprime, zter, msl=True, zmdl=25000):
    """
    zaprime : float : height from pardump
    zter : terrain height
    zmdl : model top
    za : actual height.
    """

    za = zaprime * (zmdl - zter) / float(zmdl)
    if msl:
        za += zter
    return za


def process_under(under):
    """
    under : xarray
    helper function for shift_underground and reflect_underground.
    """
    lastz = under.z.values[-1]
    under = under.sum(dim="z")
    under = under.assign_coords(z=lastz)
    under = under.expand_dims("z")
    return under


def shift_underground(dra):
    """
    dra : xarray
    Takes all mass that is underground and puts it in the first level.
    """
    import math

    iii = 0
    height = 1e-12
    for val in dra.z.values:
        if height < val:
            break
        if math.isclose(height, val):
            break
        iii += 1
    under = dra.isel(z=np.arange(0, iii + 1))
    under = process_under(under)
    above = dra.isel(z=np.arange(iii + 1, len(dra.z.values), 1))
    new = xr.concat([under, above], dim="z")
    return new


def reflect_underground(dra):
    """
    dra : xarray
    reflects mass that is in first 3 levels underground
    onto first three levels above ground.
    """
    import math

    iii = 0
    height = 0
    for val in dra.z.values:
        if height < val:
            break
        if math.isclose(height, val):
            break
        iii += 1
    print("INDEX", iii)
    under1 = dra.isel(z=[iii - 1, iii])
    print("under1", under1)
    under1 = process_under(under1)
    under = under1
    jjj = iii + 1
    if iii - 2 > 0:
        under2 = dra.isel(z=[iii - 2, iii + 1])
        print("under2", under2)
        under2 = process_under(under2)
        jjj = iii + 2
        under = xr.concat([under, under2], dim="z")
    if iii - 3 > 0:
        under3 = dra.isel(z=[iii - 3, iii + 2])
        print("under3", under3)
        under3 = process_under(under3)
        jjj = iii + 3
        under = xr.concat([under, under3], dim="z")

    # under = under1
    print("under", under)
    above = dra.isel(z=np.arange(jjj, len(dra.z.values), 1))
    print("above", above.z.values)

    # lastz = under.z.values[-1]
    # under = under.sum(dim='z')

    # under = under.assign_coords(z=lastz)
    # under = under.expand_dims('z')

    new = xr.concat([under, above], dim="z")

    return new


def getvolume(mean, cov, x, y, z, dh, dz, verbose=False):
    """
    mean and cov of gaussian
    rv is multivariate_normal object
    x, y, z : float. center position.
    dh : float. half width/ length
    dz : float. half height.
    Returns
    volume under portion of 3d Gaussian.
    """
    rv = multivariate_normal(mean, cov)

    aa1 = rv.cdf([x + dh, y + dh, z + dz])
    aa5 = rv.cdf([x - dh, y + dh, z - dz])
    aa8 = rv.cdf([x - dh, y + dh, z + dz])
    aa4 = rv.cdf([x + dh, y + dh, z - dz])

    aa6 = rv.cdf([x + dh, y - dh, z + dz])
    aa2 = rv.cdf([x - dh, y - dh, z - dz])
    aa3 = rv.cdf([x + dh, y - dh, z - dz])
    aa7 = rv.cdf([x - dh, y - dh, z + dz])

    v1 = (aa1 + aa5) - (aa8 + aa4)
    v2 = (aa6 + aa2) - (aa3 + aa7)

    if verbose:
        print(aa1, aa5)
    if verbose:
        print(aa8, aa4)
    if verbose:
        print(aa6, aa2)
    if verbose:
        print(aa7, aa3)
    volume = v1 - v2
    if verbose:
        print("volume, V1, V2", volume, "=", v1, "-", v2)
    return volume