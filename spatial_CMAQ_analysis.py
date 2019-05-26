#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:51:30 2019
@author: ehsan (ehsanm@dri.edu , ehsan.mosadegh@gmail.com)
purpose: spatial plot for CMAQ
"""

# import libraries
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap , cm
from osgeo import gdal, gdal_array, osr


# define functions that calculate concentrations of pollutants
def function_co ( domain_rows , domain_cols , cmaq_data , lay ):  # the order of argumenrs is important when input.
	
	data_mesh = np.empty( shape=( domain_rows , domain_cols ) )

	# start CMAQ algorithm
	for row in range(0,domain_rows,1):

		print('--------------------------------------')
		print('   new row starts   ')
		print('--------------------------------------')

		for col in range(0, domain_cols,1):

			aconc_24hr_cell_list = []

			for tstep in range(0,24,1):

				print('--------------------------------------')

				print('loop for row=%s col=%s time-step=%s' %(row,col,tstep))

				hrly_aconc = cmaq_data[tstep][lay][row][col]

				aconc_24hr_cell_list.append(hrly_aconc)


			aconc_24hr_cell_array = np.array( aconc_24hr_cell_list)

			cell_mean = aconc_24hr_cell_array.mean()

			data_mesh[row][col] = cell_mean

			del aconc_24hr_cell_list

	return data_mesh


# some settings
cmaq_pol = 'CO'
lay = 0
domain_cols = 250
domain_rows = 265

# import input files
cmaq_file = '/storage/ehsanm/USFS_CA_WRF_1km/plots/CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen1_mpi_standard_20160901.nc'
mcip_file = '/storage/ehsanm/USFS_CA_WRF_1km/plots/GRIDDOT2D_161001'

mcip_input = Dataset( mcip_file )
cmaq_input = Dataset( cmaq_file )

# get some info
print('-> MCIP file dimensions: %s' %str( mcip_input.variables['LATD'].dimensions ) )
print('-> shape of each dimension: %s' %( str(mcip_input.variables['LATD'].shape ) ))

# extract lat and lon parameteres
lat_mesh = np.array( mcip_input.variables['LATD'][0][0][:][:] ) # select only rosws and cols for the 1st timestep and layer
lon_mesh = np.array( mcip_input.variables['LOND'][0][0][:][:] )

#data_mesh = np.random.rand(265,250)*10
cmaq_data = cmaq_input.variables[cmaq_pol]

# functions for each pollutant - the output will be data_mesh array
data_mesh = function_co( domain_rows , domain_cols , cmaq_data , lay )





# data_mesh = np.empty( shape=(domain_rows , domain_cols) )

# # start CMAQ algorithm
# for row in range(0,domain_rows,1):

# 	print('--------------------------------------')
# 	print('   new row starts   ')
# 	print('--------------------------------------')

# 	for col in range(0, domain_cols,1):

# 		aconc_24hr_cell_list = []

# 		for tstep in range(0,24,1):

# 			print('--------------------------------------')

# 			print('loop for row=%s col=%s time-step=%s' %(row,col,tstep))

# 			hrly_aconc = cmaq_data[tstep][lay][row][col]

# 			aconc_24hr_cell_list.append(hrly_aconc)


# 		aconc_24hr_cell_array = np.array( aconc_24hr_cell_list)

# 		cell_mean = aconc_24hr_cell_array.mean()

# 		data_mesh[row][col] = cell_mean

# 		del aconc_24hr_cell_list



# create raster file
xmin,ymin,xmax,ymax = [lon_mesh.min(),lat_mesh.min(),lon_mesh.max(),lat_mesh.max()]

nrows,ncols = np.shape(data_mesh)

xres = (xmax-xmin)/float(ncols)
yres = (ymax-ymin)/float(nrows)
geotransform=(xmin,xres,0,ymax,0, -yres)
# That's (top left x, w-e pixel resolution, rotation (0 if North is up),
#         top left y, rotation (0 if North is up), n-s pixel resolution)
# I don't know why rotation is in twice???

output_raster = gdal.GetDriverByName('GTiff').Create('myraster.tif' , ncols , nrows , 1 , gdal.GDT_Float32)  # Open the file

output_raster.SetGeoTransform(geotransform)  # Specify its coordinates

srs = osr.SpatialReference()                 # Establish its coordinate encoding

srs.ImportFromEPSG(4326)                     # This one specifies WGS84 lat long.

output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system
                                                   # to the file
output_raster.GetRasterBand(1).WriteArray(data_mesh)   # Writes my array to the raster

output_raster.FlushCache()




# plot dots from grid coordinates of the dots
print('-> plotting the data...')

#plt.plot( lon_mesh , lat_mesh , marker='.' , color='b' , linestyle= 'none' )

# plot the domain/region borders
xcent =-120.806 # degrees
ycent =40.000 # degrees

NROWS = 265*1000 # meters
NCOLS = 250*1000 # meters

llcornerx=-117500 # meters
llcornery=-265500 # meters

urcornerx=132500 # meters
urcornery=-500 # meters

# create a Basemap class/model instance for a specific projection
#basemap_instance = Basemap(projection='lcc' , lat_0=ycent , lon_0=xcent , height=NROWS , width=NCOLS , resolution='i') # , area_thresh=0.1) # latlon=True for when x and y are not in map proj. coordinates
basemap_instance = Basemap(projection='lcc' ,
                           llcrnrx=llcornerx , llcrnry=llcornery , urcrnrx=urcornerx , urcrnry=urcornery ,
                           lat_0=ycent , lon_0=xcent , height=NROWS , width=NCOLS ,
                           resolution='f')

basemap_instance.bluemarble()

x_mesh, y_mesh = basemap_instance(lon_mesh , lat_mesh) # order: x , y, transforms from degree to meter for LCC

basemap_instance.drawmapboundary(color='k' ) #, fill_color='aqua')
basemap_instance.drawcoastlines(color = '0.15')
basemap_instance.drawcounties()
basemap_instance.drawstates()

#basemap_instance.fillcontinents(lake_color='aqua')

image1 = basemap_instance.pcolormesh(x_mesh , y_mesh , data_mesh , cmap=plt.cm.OrRd , shading='flat')
#im2 = basemap_instance.pcolormesh(lon_mesh , lat_mesh , data_mesh , cmap=plt.cm.jet , shading='flat')

cb = basemap_instance.colorbar(image1 , 'bottom' , label='CO concentration [ppmV]')

#cs = basemap_instance.contourf(lon_mesh , lat_mesh , data_mesh)
#cbar = basemap_instance.colorbar(cs, location='bottom')

fig_dir = '/storage/ehsanm/USFS_CA_WRF_1km/plots/CMAQ_analysis/'
fig_name = 'co_spatial_1.png'

out_fig = fig_dir+fig_name
plt.savefig(out_fig)

#plt.show() # opens a window to show the results - after saving

print('-> output figure is stored at:')
print(out_fig)

plt.close()

mcip_input.close()
cmaq_input.close()
