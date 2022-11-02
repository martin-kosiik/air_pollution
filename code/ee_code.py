import ee
import os

os.chdir('C:/Users/marti/Dropbox/research_projects/air_pollution')

# Trigger the authentication flow.
ee.Authenticate()

# Initialize the library.
ee.Initialize()

# Print metadata for a DEM dataset.
print(ee.Image('USGS/SRTMGL1_003').getInfo())


year = '2016'
dataset = ee.ImageCollection('MODIS/061/MCD64A1').filter(ee.Filter.date(year+ '-09-01',str((int(year)+1))+ '-08-25'))
burnedArea = dataset.select('BurnDate')


imageVisParam4 = {"opacity":1,"bands":["BurnDate"],"min":0.0004509199352469295,"max":0.0016086751129478216,"palette":["3730a5","ff0000"]}

burnedAreaVis = {min: 30.0, max: 341.0, 'palette': ['4e0400', '951003', 'c61503', 'ff1901'],}
burnedAreaVis = {min: 30.0, max: 341.0, 'palette': ['4e0400', '951003', 'c61503', 'ff1901'],}
croplandVis = {min: 0, max: 1, 'palette': ['brown', 'green'],}


pop_exposure_orig = ee.Image('projects/ee-martinkosiik/assets/pop_exposure_matrix_india_72h_with_crs')
mcd12q1 = ee.ImageCollection("MODIS/006/MCD12Q1")
pop_exposure = ee.Image("projects/ee-martinkosiik/assets/pop_exposure_matrix_india_96h_with_crs")
punjab_har_table = ee.FeatureCollection("projects/ee-martinkosiik/assets/punjab_haryana_village")

burnedArea = dataset.select('BurnDate')
pop_exp_b1 = pop_exposure.select('b1')


modis_projection = burnedArea.first().projection()
pop_exp_b1_reproject = pop_exp_b1.resample('bicubic').reproject(modis_projection)

mcd12q1Yr = ee.Image(mcd12q1.filter(ee.Filter.calendarRange(int(year),int(year),'year')).first()).select('LC_Type2').eq(12)
mcd12q1Yr = mcd12q1Yr.updateMask(mcd12q1Yr.neq(0))


burnedAreaBinary = burnedArea.max().gt(1.0)
burnedAreaBinary = burnedAreaBinary.multiply(mcd12q1Yr)

burnedAreaM2 = burnedAreaBinary.multiply(ee.Image.pixelArea())
mcd12q1YrM2 = mcd12q1Yr.multiply(ee.Image.pixelArea())


FinalImpact = burnedAreaBinary.multiply(pop_exp_b1_reproject)
FinalImpact



Map = geemap.Map()
Map = geemap.Map(center=[30, 76], zoom=7)
url = 'https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}&hl=en'
Map.add_tile_layer(url, name='Google Maps', attribution='Google')

Map.addLayer(mcd12q1Yr, croplandVis, 'Cropland', opacity=0.68)
Map.addLayer(burnedAreaBinary, burnedAreaVis, 'Burned Area', opacity=0.75)
#Map.addLayer(punjab_har_table, {'color': 'blue'}, 'table')
#Map.addLayer(FinalImpact, imageVisParam4, 'Final impact')

Map.add_legend(legend_title="Burned area", legend_dict={'Cropland': '008000', 'Burned area': 'ff1901'})

Map

#Map.to_image(filename=png_file)
