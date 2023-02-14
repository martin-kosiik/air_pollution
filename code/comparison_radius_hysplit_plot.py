import os
import numpy as np
import pandas as pd
import xarray as xr
import rasterio
import rioxarray
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf


os.chdir('C:/Users/marti/Dropbox/research_projects/air_pollution')

#region_of_interest = 'vietnam'
#region_of_interest = 'india_low_res'
region_of_interest = 'india'


population_raster = rioxarray.open_rasterio(r'data\population_grid\ppp_2018_1km_Aggregated.tif')


hysplit_pop_exposure = rioxarray.open_rasterio(f'data/intermed_files/pop_exposure_matrix_{region_of_interest}_deposition_with_crs.tif')

radius_pop_exposure = rioxarray.open_rasterio(f'data/intermed_files/wind_direction_weighted_pop_within_30km.tif')

infant_share =  26/1210
infant_mort_increase_per_km2_burned = 1.060/1000
infant_mort_increase_per_ha_burned = infant_mort_increase_per_km2_burned* 0.01
radius_pop_exposure = radius_pop_exposure*infant_share*infant_mort_increase_per_ha_burned 

pm25_increase = 0.45
infat_mortality_increase = infant_mort_increase_per_km2_burned/pm25_increase
#Crop yield in kg/ha
cy = 3774
# Residue-to-crop weight ratio (unitless)
rc = 1.5
# Dry matter fraction of the crop
f_dm = 0.85
# Combustion completeness (fraction of the dry matter burned)
f_cc = 0.67
# Emission factor (g of species per kg of dry matter)
ef = 5.1

# phi in mg of species per ha of crop
phi_mg_per_ha = cy * rc * f_dm * f_cc * ef

# adjust the concentration (since we the HYSPLIT simulations ran only for 4 days)
# but we are interested in the effect on monthly PM2.5 concentrations
conc_adj = 4/30


#(hysplit_pop_exposure*infant_share*phi_mg_per_ha*infat_mortality_increase*conc_adj*10).plot()

hysplit_pop_exposure = hysplit_pop_exposure*infant_share*phi_mg_per_ha*infat_mortality_increase*conc_adj

hysplit_pop_exposure.plot()
(radius_pop_exposure*infant_share*infant_mort_increase_per_ha_burned).plot()
radius_pop_exposure.plot()
population_raster_cropped = population_raster.rio.clip_box(
    minx=73.5,
    miny=27,
    maxx=78,
    maxy=32,
    crs="EPSG:4326",
)


ccrs.epsg(4326)
ccrs.epsg(26910)

hysplit_pop_exposure_matched = hysplit_pop_exposure.rio.reproject_match(population_raster_cropped, 
                                                           resampling=rasterio.enums.Resampling.cubic)


radius_pop_exposure_matched = radius_pop_exposure.rio.reproject_match(population_raster_cropped, 
                                                           resampling=rasterio.enums.Resampling.cubic)




burned_area = rioxarray.open_rasterio("data/burned_area/MCD64A1.006_Burn_Date_doy2019274_aid0001 (1).tif") 
#burned_area = burned_area.rio.set_crs(rasterio.crs.CRS.from_epsg(4326), inplace = False) #
burned_area = (burned_area > 0)*1
burned_area.rio.crs
burned_area.rio.nodata
burned_area.rio.nodata

burned_area_reproject = burned_area.rio.reproject(rasterio.crs.CRS.from_epsg(4326))



burned_area_matched = burned_area_reproject.rio.reproject_match(population_raster_cropped, 
                                                           resampling=rasterio.enums.Resampling.cubic)
burned_area_matched.plot(ylim=[28,32])


burned_area_matched = (burned_area_matched > 0)*1


burned_area.plot()
burned_area_reproject.plot()


xlim_end = 0.3
bins =np.linspace(start = 0, stop = xlim_end, num = 50)
hysplit_pop_exposure_matched.plot.hist(bins = bins, xlim=[0, xlim_end], alpha=0.55)
radius_pop_exposure_matched.plot.hist(bins = bins,xlim=[0, xlim_end], alpha=0.55)
plt.legend(['HYSPLIT average dispersion approach', 'Reduced-form approach'], loc="upper right")
plt.xlabel('Expected infant deaths per hectare of burned land')
plt.title('')
plt.savefig("figures/radius_approach/hist_comparison_exp_inf_deaths.png")



ax = plt.axes(projection = ccrs.Mercator())
ax.add_feature(cf.COASTLINE)
ax.add_feature(cf.BORDERS)
p = radius_pop_exposure_matched.plot(alpha=0.5, transform=ccrs.Mercator())

plt.title('Expected infant deaths per hectare of burned land')
plt.ylabel('latitude')
plt.xlabel('longitude')



#[73, 78.5, 26.75, 32.25]

ax = plt.axes(projection = ccrs.PlateCarree())
ax.set_extent([73.5, 78, 27.0, 32], crs=ccrs.PlateCarree())
ax.coastlines(resolution='50m')
ax.add_feature(cf.BORDERS, linewidth=2.5)
ax.add_feature(cf.STATES, linewidth=0.5)
gl = ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,
              linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
ax.plot(77.216721, 28.644800, 'bo', markersize=2, transform=ccrs.PlateCarree())
ax.text(77, 28.3, 'New Delhi', transform=ccrs.PlateCarree())
radius_pop_exposure_matched.plot(alpha=1, transform=ccrs.PlateCarree())
plt.title('Expected infant deaths per hectare of burned land')
plt.savefig("figures/radius_approach/radius_pop_exposure_exp_inf_deaths.png")


ax = plt.axes(projection = ccrs.PlateCarree())
ax.set_extent([73.5, 78, 27.0, 32], crs=ccrs.PlateCarree())
ax.coastlines(resolution='50m')
ax.add_feature(cf.BORDERS, linewidth=2.5)
ax.add_feature(cf.STATES, linewidth=0.5)
gl = ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,
              linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
gl.top_labels = False
gl.right_labels = False
ax.plot(77.216721, 28.644800, 'bo', markersize=2, transform=ccrs.PlateCarree())
ax.text(77, 28.3, 'New Delhi', transform=ccrs.PlateCarree())
hysplit_pop_exposure_matched.plot(alpha=1, transform=ccrs.PlateCarree())
plt.title('Expected infant deaths per hectare of burned land')
plt.savefig("figures/radius_approach/hysplit_pop_exposure_exp_inf_deaths.png")





hysplit_pop_exposure_matched_masked = hysplit_pop_exposure_matched.where(burned_area_matched >0)
radius_pop_exposure_matched_masked = radius_pop_exposure_matched.where(burned_area_matched >0)

hysplit_pop_exposure_matched_masked.plot()
radius_pop_exposure_matched_masked.plot()


xlim_end = 0.3
bins =np.linspace(start = 0, stop = xlim_end, num = 50)
hysplit_pop_exposure_matched_masked.plot.hist(bins = bins, xlim=[0, xlim_end], alpha=0.55)
radius_pop_exposure_matched_masked.plot.hist(bins = bins,xlim=[0, xlim_end], alpha=0.55)
plt.legend(['HYSPLIT average dispersion approach', 'Reduced-form approach'], loc="upper right")
plt.xlabel('Expected infant deaths per hectare of burned land')
plt.title('')
plt.savefig("figures/radius_approach/hist_comparison_exp_inf_deaths_filter.png")
