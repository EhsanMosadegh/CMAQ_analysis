#!/usr/bin/env python3
#====================================================
#import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
#from mpl_toolkits.basemap import Basemap
#from osgeo import gdal, gdal_array, osr , ogr
from osgeo import gdal
import ogr, os, osr
#import rasterio
#from rasterio.transform import from_origin
#import time

#====================================================

home_dir = '/Volumes/Ehsanm_DRI/cmaq_usfs/'   # '/' at the end
mcip_dir = home_dir   # '/' at the end
raster_dir = home_dir+'/raster_dir/'

aconc_file = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen1_mpi_standard_20161110.nc'
acon_file_dir = '/Volumes/Ehsanm_DRI/cmaq_usfs/scen_1/'

mcip_file = 'GRIDDOT2D_161001'

cmaq_pol = 'CO'
domain_rows = 265
domain_cols = 250
lay=0

test_with_original_data = 'off' # on OR off
test_with_sample_data = 'on'

#====================================================
print('-> opening MCIP file...')
mcip_in = mcip_dir + mcip_file
mcip_open = Dataset( mcip_in , 'r')

aconc_open = Dataset ( acon_file_dir+aconc_file , 'r')
print('-> opening ACONC file...')

### get some info
print('-> MCIP file info:')
print('-> MCIP file dimensions: %s' %str( mcip_open.variables['LATD'].dimensions ) )
print('-> shape of each dimension: %s' %( str(mcip_open.variables['LATD'].shape ) ))

### extract lat and lon parameteres
lat_mesh = np.array( mcip_open.variables['LATD'][ 0 , 0 , : , : ] ) # select only rosws and cols for the 1st timestep and layer = [ tstep=0 , lay=0]
lon_mesh = np.array( mcip_open.variables['LOND'][ 0 , 0 , : , : ] )

# get the min of each array as the origin of array
array_origin_lat = lat_mesh.min()
array_origin_lon = lon_mesh.min()

print('-> closing MCIP file...')
mcip_open.close()
print(" ")

#====================================================
if ( test_with_original_data == 'on' ):

	print('-> processing for single plot ...')
	### create an empty tensor for each cell and day as container of daily 24-hr t-step concentrations
	daily_tensor_scen = np.empty ( shape= ( 24 , domain_rows , domain_cols ))  # when assignin gby index, z-dim should be 1.

	#print(f'-> shape of daily tseries array={daily_tensor_scen.shape }')
	#print('-> traversing each cell and extract pollutants ...')

	### traverse each cell in the C-storing style for each day: row and then col
	for row in range( 0 , domain_rows , 1 ) :

		for col in range( 0 , domain_cols , 1 ) :

			#if ( processing_pol == 'co' ) :

			print( f'-> extracting cell for single POL at row= {row} and col={col} ... ' )


			cell_24hr_series_list = []
			# extract all 24 t-step
			cell_24hr_series_list = aconc_open.variables[ cmaq_pol ][ : , lay , row , col ]
			# change daily list to daily np.Array
			cell_24hr_series_array = np.array( cell_24hr_series_list )

			#cell_24hr_tseries_for_singlePol = function_cell_24hr_timeSeries_singlePOL( aconc_open_scen , cmaq_pol , lay , row , col )
			#print(f'--> cell tseries is= {cell_24hr_tseries_for_singlePol}')
			daily_tensor_scen [:,row,col]  = cell_24hr_series_array



	daily_2d_array_scen = np.ndarray( shape= ( domain_rows , domain_cols ) )

	for row in range( 0 , daily_tensor_scen.shape[1] , 1 ) :

		for col in range( 0 , daily_tensor_scen.shape[2] , 1 ) :

			print(f'-> processing the daily-tensor @ row= {row} and col= {col} ')

			cell_z_axis = daily_tensor_scen [ : , row , col ]

			#print(f'-> size of z-axis is= { cell_z_axis.shape } ')

			### take average of each z-axis
			cell_z_axis_mean = cell_z_axis.mean()
			### asign the cell mean to 2D array
			daily_2d_array_scen [ row , col ] = cell_z_axis_mean

	output_array = daily_2d_array_scen



#====================================================
if ( test_with_sample_data == 'on') :

	output_array = np.array([[ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
	                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
	                      [ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
	                      [ 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
	                      [ 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
	                      [ 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1],
	                      [ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1],
	                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
	                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
	                      [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]])


raster_name = 'raster_'+cmaq_pol+'_singleday.tif'
path = raster_dir + raster_name

print( f'-> raster name= {raster_name}')
print( f'-> raster path= {path}')

# center of domain
xcent =-120.806 # degrees
ycent =40.000 # degrees
# upper_lat = 40.450
# lower_lat = 36.450

no_of_bands = 1
datatype = gdal.GDT_Float32
#epsg_code = 4326 # output coord-ref-sys

raster_origin = ( array_origin_lon , array_origin_lat ) # unit? metere or degree? it can be either meter, or degrees --> ImportFromProj4(+units=m)

pixelWidth = 1000 # meters
pixelHeight = 1000 # meters

Xorig = raster_origin[0] # -117500 #
Yorig =  raster_origin[1] # -265500

rows = output_array.shape[0]
cols = output_array.shape[1]

geotransform = ( Xorig , pixelWidth , 0 , Yorig , 0 , pixelHeight ) # units? meter or degrees?

# get the class of coordinate reference system
crs = osr.SpatialReference()
#crs.ImportFromEPSG( epsg_code )
# ??? how define/set projection parameters for lcc for my dataset???
crs.ImportFromProj4( '+proj=lcc +lat_0=40.000  +lon_0=-120.806  +lat_1=40.450  +lat_2=36.450  +units=m , +datum=WGS84 +no_defs' )

# Initialize driver & create file
driver = gdal.GetDriverByName('GTiff')
# create output raster/matrix dataseet to put data into it
out_raster = driver.Create( path, cols, rows, no_of_bands , datatype ) # (path, cols, rows, bands, dtype-> GDAL data type arg)

out_raster.SetGeoTransform(geotransform)  # Specify its coordinates
# ??? how set projection for lcc???
out_raster.SetProjection( crs.ExportToWkt() )

print( f'-> writing the raster ...')
out_raster.GetRasterBand( no_of_bands ).WriteArray( output_array )  # write my array to the raster
out_raster.FlushCache()

# Once we're done, close the dataset properly
out_raster = None
print(" ")
