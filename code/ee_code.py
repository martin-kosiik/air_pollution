import ee
import os
#import geemap

import folium
import eefolium

import io
#from PIL import Image
import branca.colormap as cmp


def add_categorical_legend(folium_map, title, colors, labels):
    import folium
    if len(colors) != len(labels):
        raise ValueError("colors and labels must have the same length.")

    color_by_label = dict(zip(labels, colors))
    
    legend_categories = ""     
    for label, color in color_by_label.items():
        legend_categories += f"<li><span style='background:{color}'></span>{label}</li>"
        
    legend_html = f"""
    <div id='maplegend' class='maplegend'>
      <div class='legend-title'>{title}</div>
      <div class='legend-scale'>
        <ul class='legend-labels'>
        {legend_categories}
        </ul>
      </div>
    </div>
    """
    script = f"""
        <script type="text/javascript">
        var oneTimeExecution = (function() {{
                    var executed = false;
                    return function() {{
                        if (!executed) {{
                             var checkExist = setInterval(function() {{
                                       if ((document.getElementsByClassName('leaflet-top leaflet-right').length) || (!executed)) {{
                                          document.getElementsByClassName('leaflet-top leaflet-right')[0].style.display = "flex"
                                          document.getElementsByClassName('leaflet-top leaflet-right')[0].style.flexDirection = "column"
                                          document.getElementsByClassName('leaflet-top leaflet-right')[0].innerHTML += `{legend_html}`;
                                          clearInterval(checkExist);
                                          executed = true;
                                       }}
                                    }}, 100);
                        }}
                    }};
                }})();
        oneTimeExecution()
        </script>
      """
   

    css = """

    <style type='text/css'>
      .maplegend {
        z-index:9999;
        float:right;
        background-color: rgba(255, 255, 255, 1);
        border-radius: 5px;
        border: 2px solid #bbb;
        padding: 10px;
        font-size:12px;
        positon: relative;
      }
      .maplegend .legend-title {
        text-align: left;
        margin-bottom: 5px;
        font-weight: bold;
        font-size: 90%;
        }
      .maplegend .legend-scale ul {
        margin: 0;
        margin-bottom: 5px;
        padding: 0;
        float: left;
        list-style: none;
        }
      .maplegend .legend-scale ul li {
        font-size: 80%;
        list-style: none;
        margin-left: 0;
        line-height: 18px;
        margin-bottom: 2px;
        }
      .maplegend ul.legend-labels li span {
        display: block;
        float: left;
        height: 16px;
        width: 30px;
        margin-right: 5px;
        margin-left: 0;
        border: 0px solid #ccc;
        }
      .maplegend .legend-source {
        font-size: 80%;
        color: #777;
        clear: both;
        }
      .maplegend a {
        color: #777;
        }
    </style>
    """

    folium_map.get_root().header.add_child(folium.Element(script + css))

    return folium_map


os.chdir('C:/Users/marti/Dropbox/research_projects/air_pollution')

import sys
 
# adding Folder_2/subfolder to the system path
sys.path.insert(0, '/code')
 
# importing the hello
from ee_utilities import add_categorical_legend
from ee_utilities import hello
import ee_utilities


# calling hello function
hello()

# Trigger the authentication flow.
ee.Authenticate()

# Initialize the library.
ee.Initialize()

# Print metadata for a DEM dataset.


m = folium.Map(location=[45.5236, -122.6750])
m

import folium
map = folium.Map([51., 12.], zoom_start=6,control_scale=True)
#folium.GeoJson(data).add_to(map)
map.save('figures/map.html')


img_data = m._to_png(5)
img = Image.open(io.BytesIO(img_data))
img.save('figures/image.png')

map._to_png() 



import eefolium



year = '2016'
dataset = ee.ImageCollection('MODIS/061/MCD64A1').filter(ee.Filter.date(year+ '-09-01',str((int(year)+1))+ '-08-25'))
burnedArea = dataset.select('BurnDate')

AreaOfInterest = ee.Geometry.Rectangle([73.5,27,78,32],'EPSG:4326', False)

0.0014647410716861486

imageVisParam2 = {"opacity":1,"bands":["b1"],"max":0.001464119297452271,"min":0.0001823857455747202,"palette":["92a58d","0e2c80"]}
imageVisParam4 = {"opacity":1,"bands":["BurnDate"],"max":0.001464119297452271,"min":0.00040931659168563783,"palette":["3730a5","ff0000"]}

burnedAreaVis = {min: 30.0, max: 341.0, 'palette': ['4e0400', '951003', 'c61503', 'ff1901'],}
burnedAreaVis = {min: 30.0, max: 341.0, 'palette': ['4e0400', '951003', 'c61503', 'ff1901'],}
croplandVis = {min: 0, max: 1, 'palette': ['brown', 'green'],}


#pop_exposure_orig = ee.Image('projects/ee-martinkosiik/assets/pop_exposure_matrix_india_72h_with_crs')
mcd12q1 = ee.ImageCollection("MODIS/006/MCD12Q1")
#pop_exposure = ee.Image("projects/ee-martinkosiik/assets/pop_exposure_matrix_india_96h_with_crs")
pop_exposure = ee.Image("projects/ee-martinkosiik/assets/pop_exposure_matrix_india_higher_res_with_crs")

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


#################################
# 
Map = eefolium.Map()
Map = eefolium.Map(center=[29.5, 76], zoom=7)
url = 'https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}&hl=en'
Map.add_tile_layer(url, name='Google Maps', attribution='Google')

Map.addLayer(mcd12q1Yr, croplandVis, 'Cropland', opacity=0.68)
Map.addLayer(burnedAreaBinary, burnedAreaVis, 'Burned Area', opacity=0.75)
Map.addLayer(AreaOfInterest, {'pallete': 'black'}, 'Output region', opacity=0.55)

#palette = ['blue', 'purple', 'cyan', 'green','red']
#Map.create_colorbar(width=250, height=30, palette=palette, vertical=False,add_labels=True, font_size=20, labels=[-40, 35])
add_categorical_legend(Map, '',
                             colors = ['black', '#008000','#ff1901'],
                           labels = ['Area of interest', 'Cropland', 'Burned area'])



img_data = Map._to_png(5)
img = Image.open(io.BytesIO(img_data))
img.save('figures/burned_area_folium.png')






# Final impact
Map = eefolium.Map()
Map = eefolium.Map(center=[30, 76], zoom=7.2)
url = 'https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}&hl=en'
Map.add_tile_layer(url, name='Google Maps', attribution='Google')

#Map.addLayer(punjab_har_table, {'color': 'blue'}, 'table')
Map.addLayer(FinalImpact, imageVisParam4, 'Final impact')

imageVisParam4['palette'] = ['#3730a5', '#ff0000']
linear = cmp.LinearColormap(
    imageVisParam4['palette'],
    vmin=imageVisParam4['min'], vmax=imageVisParam4['max'],
    caption='Source impact (alpha)' #Caption for Color scale or Legend
)
linear.add_to(Map)   #adds colorscale or legend



img_data = Map._to_png(5)
img = Image.open(io.BytesIO(img_data))
img.save('figures/final_impact_folium.png')


# source impact no interpolation
Map = eefolium.Map()
Map = eefolium.Map(center=[30, 76], zoom=7.2)
Map.add_tile_layer(url, name='Google Maps', attribution='Google')

Map.addLayer(pop_exp_b1, {'bands': 'b1', 'palette': imageVisParam4['palette'], 
                          'min': imageVisParam2['min'], 'max': imageVisParam2['max']}, 'Final impact', opacity=0.92)


linear.add_to(Map)   #adds colorscale or legend



img_data = Map._to_png(5)
img = Image.open(io.BytesIO(img_data))
img.save('figures/source_impact_no_interpolation_folium.png')


# source impact with interpolation
Map = eefolium.Map()
Map = eefolium.Map(center=[30, 76], zoom=7.2)
Map.add_tile_layer(url, name='Google Maps', attribution='Google')

Map.addLayer(pop_exp_b1_reproject, {'bands': 'b1', 'palette': imageVisParam4['palette'], 
                          'min': imageVisParam2['min'], 'max': imageVisParam2['max']}, 'Final impact', opacity=0.92)


linear.add_to(Map)   #adds colorscale or legend



img_data = Map._to_png(5)
img = Image.open(io.BytesIO(img_data))
img.save('figures/source_impact_with_interpolation_folium.png')











#Map.to_image(filename='figures/burned_area_folium.png')




linear = cmp.LinearColormap(
    ['yellow', 'green', 'purple'],
    vmin=3, vmax=10,
    caption='Color Scale for Map' #Caption for Color scale or Legend
)
linear
linear.add_to(Map)   #adds colorscale or legend



Map.add_legend(legend_title="Burned area", legend_dict={'Cropland': '008000', 'Burned area': 'ff1901'})

Map

Map.to_image(filename='figures/map.png')





def crop_images(file_name='burned_area_folium.png'):

    im = Image.open(r"figures/" + file_name)
 
    # Size of the image in pixels (size of original image)
    # (This is not mandatory)
    width, height = im.size
 
# Setting the points for cropped image
    left = round(0.18 * width)
    top = 0 
    right = width
    bottom = round(0.95 * height)
 
# Cropped image of above dimension
# (It will not change original image)
    im1 = im.crop((left, top, right, bottom))
 
# Shows the image in image viewer
    #im1.show()

    im1.save('figures/cropped/' + file_name)



crop_images()

figures_to_crop = ['burned_area_folium.png', 'final_impact_folium.png', 
                    'source_impact_no_interpolation_folium.png', 
                    'source_impact_with_interpolation_folium.png']


[crop_images(file_name=image_path) for image_path in figures_to_crop]

# Setting the points for cropped image
left = 5
top = height / 4
right = 164
bottom = 3 * height / 4
