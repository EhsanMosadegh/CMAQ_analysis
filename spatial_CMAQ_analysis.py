#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#====================================================================================================
# import libraries

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
#from osgeo import gdal, gdal_array, osr , ogr
from osgeo import gdal
import ogr, os, osr
#import rasterio
#from rasterio.transform import from_origin
import time

#====================================================================================================
# run time setting

### get the starting time
start = time.time()

### run time settings
cmaq_file_month = '10'		 # 07, 08, 09, 10, 11
sim_month = 'oct'  			# jul, aug, sep, oct, nov

cmaq_file_year = '2016'
days_to_run_in_month = 1
scenario = '4' 			# 1-5, baseline
mcip_date_tag = '161001'

cmaq_pol = 'CO'  # for plot title 'CO' ro 'PM2.5'
pol_unit = 'ppmV'		#'[ppmV]' or [ug/m^3]
max_conc_threshold = 0.2  # for Basemap plot

### spatial plot
spatial_plotting = 'yes' # yes or no
processing_pol = 'co' 		# 'co' or 'pm2.5'
processing_method = 'single_plot' 	# 'single_plot' or 'diff_plot'
produce_raster = 'yes' 	# 'yes' OR 'no'

### time-series plot
timeseries_plotting = 'no' # yes or not
favorite_row = 5
favorite_col = 5

platform = 'Mac'  # 'Mac' or 'cluster'


# ### Basemap plot setting
# center of domain
xcent =-120.806 # degrees
ycent =40.000 # degrees
# # domain size
pixel_size = 1000 # meters
NROWS = 265*pixel_size # meters
NCOLS = 250*pixel_size # meters
# # lower-left corner
llx = -117500 # meters
lly = -265500 # meters
ulx = llx
uly = lly + NROWS
# # upper-right corner
# urcornerx=132500 # meters
# urcornery=-500 # meters

### Basemap plot setting to zoom
### center of domain
xcent_zoom =-120.0324 # degrees
ycent_zoom =39.09 # degrees
### domain size
NROWS_zoom = 100000#265*1000 # meters
NCOLS_zoom = 100000#250*1000 # meters

### domain settings
lay = 0
domain_cols = 250
domain_rows = 265

# xorig = -117500.000 # degree?
# yorig = -265500.000 # degree?

#raster_origin = ( xorig , yorig )
pixelWidth = NROWS
pixelHeight = NCOLS
newRasterfn = 'co_test_raster.tif'

print( f'-> CMAQ year= {cmaq_file_year}')
print( f'-> CMAQ month of analysis= {cmaq_file_month}')
print( f'-> LANDIS scenarios= {scenario}')
print( f'-> processing pollutant= {processing_pol}')
print( f'-> processing method= {processing_method}')
print( f'-> number of days to run= {days_to_run_in_month}')
print( f'-> platform is= {platform}')
print( f'-> spatial plotting= {spatial_plotting}')
print( f'-> produce raster= {produce_raster}')
print( f'-> time-series plotting= {timeseries_plotting}')
print(" ")

#====================================================================================================
# main
#====================================================================================================

def main() :

	### input directory setting
	if ( platform == 'Mac' ) :

		# input_dir = '/Volumes/USFSdata/no_oct/'   # '/' at the end
		# mcip_dir = '/Volumes/USFSdata/'   # '/' at the end
		# fig_dir = '/Volumes/USFSdata/no_oct/cmaq_figs/'  # '/' at the end
		# raster_dir =

		input_dir = '/Volumes/Ehsanm_DRI/cmaq_usfs/'   # '/' at the end
		mcip_dir = '/Volumes/Ehsanm_DRI/cmaq_usfs/'   # '/' at the end
		fig_dir = '/Volumes/Ehsanm_DRI/cmaq_usfs/cmaq_figs/'  # '/' at the end
		raster_dir = '/Volumes/Ehsanm_DRI/cmaq_usfs/raster_dir/'

	elif ( platform == 'cluster' ) :

		input_dir = '/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/cmaq_data/'
		mcip_dir = '/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/'
		fig_dir = '/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/cmaq_data/cmaq_figs'
		raster_dir = '/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/cmaq_data/raster_dir/'

	else:

		print( '-> ERROR: specify running platform ' )
		print('-> exiting ...')
		raise SystemExit()

	### set input pathes
	input_path_scen = input_dir + 'scen_' + scenario + '/' #+ sim_month + '/'
	input_path_base = input_dir + 'scen_baseline' + '/' #+ sim_month + '/'

	print('-> CMAQ input directory is:')
	print(input_path_scen)
	print(input_path_base)

	print('-> MCIP input directory is:')
	print(mcip_dir)
	print(" ")

# def function_3D_mesh_maker ( days_to_run_in_month , domain_rows , domain_cols , cmaq_file_month , scenario , input_path_scen , input_path_base ) :
# 	" processes the files and days and returns two 3D meshes"

	print('-> month of analysis is=' , cmaq_file_month )

	# ### define the monthly array mesh
	# if ( processing_pol == 'co' ) :

	#monthly_tseries_tensor_from_scen = np.ndarray( shape=( 24 , domain_rows , domain_cols ) )
	monthly_tseries_tensor_from_scen = np.empty( shape=( 0 , domain_rows , domain_cols ) ) # use zero for concatenate method
	monthly_tseries_tensor_from_base = np.empty( shape=( 0 , domain_rows , domain_cols ) ) # zero means there is no cell in z-dir

	# elif ( processing_pol == 'pm2.5' ) :

	# 	monthly_tseries_tensor_from_scen = np.empty( shape=( days_to_run_in_month , domain_rows , domain_cols ) )
	# 	monthly_tseries_tensor_from_base = np.empty( shape=( days_to_run_in_month , domain_rows , domain_cols ) )

	# else:
	# 	pass

	## create a day list for a month to create file-date-tag, use an argument-unpacking operator * to unpack the list
	day_list = [*range( 1 , days_to_run_in_month+1 , 1)] # don't forget the [] around range function to create the list

	# # to run for specific day
	# day_list = [21]  # use the favorite day
	#====================================================================================================
	### traverse the list for each day

	for day_of_the_month in day_list :

		print('-----------------------------')
		print( " ")
		print('-> no. of days to process are:')
		print(day_list)
		print( f'-> processing for day= {day_of_the_month}' )
		print(" ")

		# prepare the day flags
		if day_of_the_month <= 9 :
			# if jday is less than 10, add zero before it
			day_count = '0'+str(day_of_the_month)

		else:
			# if jday is bigger than 9, use it.
			day_count = str(day_of_the_month)

		### define file tags
		file_date_tag = cmaq_file_year + cmaq_file_month + day_count

		#====================================================================================================
		### define input files

		aconc_scen = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+scenario+'_mpi_standard_'+file_date_tag+'.nc'
		pmdiag_scen = 'CCTM_PMDIAG_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+scenario+'_mpi_standard_'+file_date_tag+'.nc'

		aconc_base = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonpt_baseline_AgBioNonpt_mpi_standard_'+file_date_tag+'.nc'
		pmdiag_base = 'CCTM_PMDIAG_v52_CA_WRF_1km_griddedAgBioNonpt_baseline_AgBioNonpt_mpi_standard_'+file_date_tag+'.nc'

		#====================================================================================================
		### define netcdf files based on each processing method and pollutant

		if ( processing_pol == 'co') : # (processing_pol == 'no2')  # we need only 1 file: "aconc"
			# for single scenario plot
			if ( processing_method == 'single_plot' ) :

				# define input files
				print('-> setting path for single POL and single_plot...')
				aconc_input_scen = input_path_scen + aconc_scen
				aconc_open_scen = Dataset( aconc_input_scen , 'r' )
				print('-> opening/reading CMAQ files:')
				print( aconc_input_scen )
				print( " ")

			# for difference between 2 scenarios
			elif ( processing_method == 'diff_plot' ) :  # open 2 netcdf files: "aconc_scen" and "aconc_baseline"
				# define input files
				print('-> setting path for single POL and diff_plot...')
				aconc_input_scen = input_path_scen + aconc_scen
				aconc_input_base = input_path_base + aconc_base
				aconc_open_scen = Dataset( aconc_input_scen , 'r' )
				aconc_open_base = Dataset( aconc_input_base , 'r' )
				print('-> opening/reading CMAQ files for single POL:')
				print( aconc_input_scen )
				print( aconc_input_base )
				print( " ")

			else:

				print( '-> WARNING: define/check single processing POL and method first! ')
				print('-> exiting ...')
				raise SystemExit()

		elif ( processing_pol == 'pm2.5' ) :   # we need 2 files: "aconc" and "pmdiag" files

			if ( processing_method == 'single_plot' ) :

				# define input files
				aconc_input_scen = input_path_scen + aconc_scen
				pmdiag_input_scen = input_path_scen + pmdiag_scen

				print('-> opening/reading CMAQ files for pm2.5:')
				print( aconc_input_scen )
				print( pmdiag_input_scen )
				print( " ")

				# define netcdf file
				aconc_open_scen = Dataset( aconc_input_scen , 'r' )
				pmdiag_open_scen = Dataset( pmdiag_input_scen , 'r' )

			elif ( processing_method == 'diff_plot' ) :

				# define input files
				aconc_input_scen = input_path_scen + aconc_scen
				pmdiag_input_scen = input_path_scen + pmdiag_scen
				aconc_input_base = input_path_base + aconc_base
				pmdiag_input_base = input_path_base + pmdiag_base

				print('-> opening/reading CMAQ files for pm2.5:')
				print(aconc_input_scen)
				print(pmdiag_input_scen)
				print(aconc_input_base)
				print(pmdiag_input_base)

				# define netcdf files
				aconc_open_scen = Dataset( aconc_input_scen , 'r' )
				pmdiag_open_scen = Dataset( pmdiag_input_scen , 'r' )
				aconc_open_base = Dataset( aconc_input_base , 'r' )
				pmdiag_open_base = Dataset( pmdiag_input_base , 'r' )

			else:
				print('-> ERROR: select single or diff methid! ')
				print('-> exiting ...')
				raise SystemExit()

		else:

			print( '-> ERROR: define/check POL first! ')
			print('-> exiting ...')
			raise SystemExit()

		#====================================================================================================
		# process netcdf files for each cell with a specific function
		# single plot:

		if ( processing_method == 'single_plot' ) :

			### create an empty tensor for each cell and day as container of daily 24-hr t-step concentrations
			daily_tensor_scen = np.empty ( shape= ( 24 , domain_rows , domain_cols ))  # when assignin gby index, z-dim should be 1.

			#print(f'-> shape of daily tseries array={daily_tensor_scen.shape }')
			#print('-> traversing each cell and extract pollutants ...')

			### traverse each cell in the C-storing style for each day: row and then col
			for row in range( 0 , domain_rows , 1 ) :

				for col in range( 0 , domain_cols , 1 ) :

					if ( processing_pol == 'co' ) :

						#print( f'-> extracting cell for single POL - singlePlot - at row= {row} and col={col} ... ' )

						cell_24hr_tseries_for_singlePol = function_cell_24hr_timeSeries_singlePOL( aconc_open_scen , cmaq_pol , lay , row , col )
						#print(f'--> cell tseries is= {cell_24hr_tseries_for_singlePol}')
						daily_tensor_scen [:,row,col]  = cell_24hr_tseries_for_singlePol

					elif ( processing_pol == 'pm2.5') :

						#print( f'-> extracting cell for pm2.5 at row= {row} and col={col} ... ' )

						cell_24hr_tseries_for_pm25 = function_pm25_daily_cell_tseries( aconc_open_scen , pmdiag_open_scen , lay , row , col )

						daily_tensor_scen [:,row,col] = cell_24hr_tseries_for_pm25 # fill tensor for all cells in domain

					else:

						print( '-> WARNING: define/check single POL or pm2.5 settings and processing method first! ')
						print('-> exiting ...')
						raise SystemExit()

			if ( produce_raster == 'yes' ) :
				print( " " )
				print( f'-> producing raster file ...')

				daily_2d_array_scen = function_3Dto2D( domain_rows , domain_cols , daily_tensor_scen )

				array2raster( raster_dir , processing_pol , file_date_tag , daily_2d_array_scen )
				# print( f'-> raster file = { output_raster }')
				# print( " " )

			### now we concatenate daily timeseries mesh to monthly tensor
			monthly_tseries_tensor_from_scen = np.concatenate( ( monthly_tseries_tensor_from_scen ,  daily_tensor_scen ) , axis=0 )

			### after all days are extracted, add the daily frame to monthly frame
			### now we concatenate daily timeseries mesh to monthly tensor
			monthly_tseries_tensor_from_scen = np.concatenate( ( monthly_tseries_tensor_from_scen , daily_tensor_scen ) , axis=0 )

		#====================================================================================================
		# diff plot

		elif ( processing_method == 'diff_plot' ) :

			### create an empty tensor for each cell and day as container of daily 24-hr t-step concentrations
			daily_tensor_scen = np.empty ( shape=( 24 , domain_rows , domain_cols ) )
			daily_tensor_base = np.empty ( shape=( 24 , domain_rows , domain_cols ) )


			print('-> traversing each cell and extract pollutants ...')
			### traverse each cell in the C-storing style for each day: row and then col
			for row in range( 0 , domain_rows , 1 ) :

				for col in range( 0 , domain_cols , 1 ) :

					if ( processing_pol == 'co' ) :

						#print( f'-> extracting cell for single POL - diff - at row= {row} and col={col} ... ' )

						# we calculate cell means for each scenario
						cell_24hr_timeSeries_array_scen = function_cell_24hr_timeSeries_singlePOL( aconc_open_scen , cmaq_pol , lay , row , col )
						cell_24hr_timeSeries_array_base = function_cell_24hr_timeSeries_singlePOL( aconc_open_base , cmaq_pol , lay , row , col )

					elif ( processing_pol == 'pm2.5' ) :

						# we calculate cell means
						cell_24hr_timeSeries_array_scen = function_pm25_daily_cell_tseries( aconc_open_scen , pmdiag_open_scen , lay , row , col )
						cell_24hr_timeSeries_array_base = function_pm25_daily_cell_tseries( aconc_open_base , pmdiag_open_base , lay , row , col )

					else:
						pass

					### fill daily tensors cells with 24-hr time-series
					daily_tensor_scen [: , row , col] = cell_24hr_timeSeries_array_scen
					daily_tensor_base [: , row , col] = cell_24hr_timeSeries_array_base

			### now we add/concatenate each daily tensor to monthly tensor
			#print( f'-> add/pin each cell 24-hr t-series to monthly_tseries_tensor_from_scen at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )

			monthly_tseries_tensor_from_scen = np.concatenate( ( monthly_tseries_tensor_from_scen ,  daily_tensor_scen ) , axis=0 )

			#print( f'-> add/pin each cell 24hr t-series to monthly_tseries_tensor_from_base at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )

			monthly_tseries_tensor_from_base = np.concatenate( ( monthly_tseries_tensor_from_base ,  daily_tensor_base ) , axis=0 )

		else:
			print('-> ERROR: processing method NOT defined!')

		#====================================================================================================
		### closing nc files for each day

		if ( processing_pol == 'co' ) :

			if ( processing_method == 'single_plot' ) :

				print('-> closing ACONC file ...')
				aconc_open_scen.close()

			elif ( processing_method == 'diff_plot' ) :

				print('-> closing ACONC scen and baseline files ...')
				aconc_open_scen.close()
				aconc_open_base.close()

			else:
				pass

		elif ( processing_pol == 'pm25' ) :

			if ( processing_method == 'single_plot' ) :

				print('-> closing ACONC and PMDIAG files ...')
				aconc_open_scen.close()
				pmdiag_open_scen.close()

			elif ( processing_method == 'diff_plot' ) :

				print('-> closing ACONC and PMDIAG files for both scen and baseline ...')
				aconc_open_scen.close()
				pmdiag_open_scen.close()
				aconc_open_base.close()
				pmdiag_open_base.close()

		else:
			pass

	#====================================================================================================
	### change 3D to 2D array to make monthly_mean_mesh_2d: used for spatial plots

	if ( spatial_plotting == 'yes' ) :

		monthly_tseries_tensor_from_scen_intermed = monthly_tseries_tensor_from_scen  # intermediate file is used in spatial plotting

		if ( processing_method == 'single_plot' ) :

			### look at the 3D tensors for scenario
			print('-----------------------------')
			print('-> 3D data mesh info:')
			print( f'-> LANDIS scenario 3D monthly tensor: number of dimensions= { monthly_tseries_tensor_from_scen_intermed.ndim}' )
			print( f'-> LANDIS scenario 3D monthly tensor: shape of data-mesh= { monthly_tseries_tensor_from_scen_intermed.shape}' )
			print('-----------------------------')

			### we only have one 3D tensor

			print('-> change mesh 3D to 2D for single plot ...')

			monthly_mean_mesh_2d = function_3Dto2D( domain_rows , domain_cols , monthly_tseries_tensor_from_scen_intermed )

		elif ( processing_method == 'diff_plot' ) :

			### look at the 3D meshes
			print('-----------------------------')
			print('-> check 3D data mesh info:')
			print( f'-> LANDIS scenario 3D monthly mean mesh: number of dimensions= {monthly_tseries_tensor_from_scen_intermed.ndim}' )
			print( f'-> LANDIS scenario 3D monthly mean mesh: shape of data-mesh= {monthly_tseries_tensor_from_scen_intermed.shape}' )
			print( f'-> baseline 3D monthly mesh: number of dimensions= {monthly_tseries_tensor_from_base.ndim}' )
			print( f'-> baseline 3D monthly mesh: shape of data-mesh= {monthly_tseries_tensor_from_base.shape}' )
			print('-----------------------------')

			### we have 2 tensors
			print('-> change mesh 3D-to-2D for diff-plot for mesh-3D-LANDIS scenario ...')
			monthly_mean_2d_mesh_scen = function_3Dto2D( domain_rows , domain_cols , monthly_tseries_tensor_from_scen_intermed )

			print('-> change mesh 3D-to-2D for diff-plot for mesh-3D-baseline scenario ...')
			monthly_mean_2d_mesh_base = function_3Dto2D( domain_rows , domain_cols , monthly_tseries_tensor_from_base )

			# now subtract 2 meshes to get the diff mesh for spatial plotting
			monthly_mean_mesh_2d = monthly_mean_2d_mesh_scen - monthly_mean_2d_mesh_base

		else:
			print('-> ERROR: check processing method for 3Dto2D ...')
			print('-> exiting ...')
			raise SystemExit()

	#====================================================================================================
	### open MCIP file to get lon-lat of domain

	print('-> opening MCIP file...')
	mcip_file = 'GRIDDOT2D_'+mcip_date_tag
	mcip_in = mcip_dir + mcip_file
	mcip_open = Dataset( mcip_in )

	### get some info
	print('-> MCIP file info:')
	print('-> MCIP file dimensions: %s' %str( mcip_open.variables['LATD'].dimensions ) )
	print('-> shape of each dimension: %s' %( str(mcip_open.variables['LATD'].shape ) ))

	### extract lat and lon parameteres
	lat_mesh = np.array( mcip_open.variables['LATD'][ 0 , 0 , : , : ] ) # select only rosws and cols for the 1st timestep and layer = [ tstep=0 , lay=0]
	lon_mesh = np.array( mcip_open.variables['LOND'][ 0 , 0 , : , : ] )
	print('-> closing MCIP file...')
	mcip_open.close()
	print(" ")

	#====================================================================================================
	### time-series plotting
	#====================================================================================================

	if ( timeseries_plotting == 'yes') :

		print('-> time-series plotting is yes ...')

		if( processing_method != 'single_plot') :

			print('-> WARNING: single plot is not set for time-series ...')

		else :

			row_ = favorite_row
			col_ = favorite_col
			conc_timeseries_list = monthly_tseries_tensor_from_scen [ : , row_ , col_ ]

			x_ = [ *range(1, ((days_to_run_in_month*24)+1) , 1) ]  # x bar is no. of hours in aconc files
			y_ = conc_timeseries_list		#.tolist  # y-bar is hourly concentrations/timeseries

			print(f'-> size of x_ = {len(x_)}')
			print(f'-> x_ = {x_}' )
			print(f'-> size of y_ = {len(y_)}')
			print(f'-> y_ = {y_}' )

			plt.plot( x_ , y_ )
			#plt.show()


			plot_name = cmaq_pol+'_timeseries'+'_scen_'+scenario+'_'+cmaq_file_year+'-'+cmaq_file_month+'_summed_'+str(days_to_run_in_month)+'_days'+'.png'
			plot_dir = fig_dir
			saved_plot = fig_dir+plot_name

			plt.savefig( saved_plot , dpi=1200 , format='png' )
			plt.close()


	#====================================================================================================
	### Spatial plotting
	#====================================================================================================
	# use Basemap library and make spatial plots

	if ( spatial_plotting == 'yes') :

		print('-> spatial plotting ...')

		### plot dots from grid coordinates of the dots
		#plt.plot( lon_mesh , lat_mesh , marker='.' , color='b' , linestyle= 'none' )

		# ### create a Basemap class/model instance for a specific projection
		# # basemap_instance = Basemap(projection='lcc' , lat_0=ycent , lon_0=xcent , height=NROWS , width=NCOLS , resolution='i') # , area_thresh=0.1) # latlon=True for when x and y are not in map proj. coordinates
		# theMap = Basemap(projection='lcc' ,
		# 													 llcrnrx=llx , llcrnry=lly , urcrnrx=urcornerx , urcrnry=urcornery ,
		# 													 lat_0=ycent , lon_0=xcent , height=NROWS , width=NCOLS ,
		# 													 resolution='f' , area_thresh=0.5)


		### create Basemap model instance from its class, it is a map that color mesh sits on it.
		theMap_zoomed = Basemap(projection='lcc' , lat_0=ycent_zoom , lon_0=xcent_zoom , height=NROWS_zoom , width=NCOLS_zoom , resolution='f' , area_thresh=0.5)

		x_mesh, y_mesh = theMap_zoomed(lon_mesh , lat_mesh) # order: x , y; Basemap model transforms lon/lat from degree to meter for LCC projection map

		#basemap_instance.fillcontinents(lake_color='aqua')

		#my_levels = [ 0.02 , 0.05 ]
		#my_colors = ( 'g' , 'b' , 'r' )

		### create a color mesh image from basemap model instance, the color mesh is constant, cos it is plotted from lon/lat values
		colorMesh = theMap_zoomed.pcolormesh( x_mesh , y_mesh , monthly_mean_mesh_2d , cmap=plt.cm.OrRd , shading='flat' , vmin=0.0 , vmax=max_conc_threshold ) #levels=my_levels , colors=my_colors

		#im2 = basemap_instance.pcolormesh(lon_mesh , lat_mesh , data_mesh , cmap=plt.cm.jet , shading='flat')

		theMap_zoomed.drawmapboundary(color='k' ) #, fill_color='aqua')
		theMap_zoomed.drawcoastlines(color = '0.15')
		theMap_zoomed.drawstates()
		theMap_zoomed.drawcounties(linewidth=0.5 , color='k')

		#theMap_zoomed.bluemarble()

		### create colorbar
		colorbar = theMap_zoomed.colorbar( colorMesh , 'bottom' , label= f'{cmaq_pol} mean concentration {pol_unit}' )
		#cs = basemap_instance.contourf(lon_mesh , lat_mesh , data_mesh)
		#colorbar = basemap_instance.colorbar(cs, location='bottom')
		#plt.subplot( figsize=(10,10) )
		if ( processing_method == 'single_plot' ) :

			plt.title(f' {cmaq_pol} monthly mean concentrations for {sim_month}, {cmaq_file_year} - LANDIS scenario {scenario}' , fontsize=7 )

		elif ( processing_method == 'diff_plot' ) :

			plt.title(f' {cmaq_pol} monthly mean concentration difference between LANDIS scenario-{scenario} and baseline for {sim_month}, {cmaq_file_year} ' , fontsize=7 )

		print(" ")

		#====================================================================================================
		# save the plots
		#====================================================================================================

		### path for saving plots
		print('-> fig directory is:')
		print(fig_dir)

		### plot name
		if ( processing_method == 'single_plot' ) :

			fig_name = cmaq_pol + '_monthlyMean' + '_scen_' + scenario + '_' + cmaq_file_year+'-'+cmaq_file_month + '_summed_' + str(days_to_run_in_month) + '_days' + '.png'

		elif ( processing_method == 'diff_plot' ) :

			fig_name = cmaq_pol + '_monthlyMean' + '_scen_' + scenario + '_difference_from_baseline_' + cmaq_file_year+'-'+cmaq_file_month + '_summed_' + str(days_to_run_in_month) + '_days' + '.png'

		else:
			pass

		### plot full path
		out_fig = fig_dir + fig_name
		print('-> output figure is stored at:')
		print(out_fig)
		### export the figure
		plt.savefig( out_fig , dpi=1200 , format='png')
		### opens a window to show the results - after savefig
		#plt.show()
		### close the plot
		plt.close()

	end = time.time()
	print( f'-> run time of main function= { (( end - start ) / 60 ) :.2f} min' )  # f-string

#====================================================================================================
### function to change 3D to 2D array

def function_3Dto2D ( domain_rows , domain_cols , mesh_3d  ) :
	" returns monthly mean of each cell, changes 3D mesh of daily mean conc to a 2D mesh of monthly mean conc"

	### define a 2d array
	mesh_2d = np.ndarray( shape= ( domain_rows , domain_cols ) )

	for row in range( 0 , mesh_3d.shape[1] , 1 ) :

		for col in range( 0 , mesh_3d.shape[2] , 1 ) :

			#print(f'-> processing mesh-3D-monthly mean @ row= {row} and col= {col} ')

			cell_z_axis = mesh_3d [ : , row , col ]

			#print(f'-> size of z-axis is= { cell_z_axis.shape } ')

			### take average of each z-axis
			cell_z_axis_mean = cell_z_axis.mean()
			### asign the cell mean to 2D array
			mesh_2d [ row ][ col ] = cell_z_axis_mean

	print(" ")
	print( f'-> shape of 2D array output of function: 3Dto2D= { mesh_2d.shape }' )
	print(" ")
	### function returns a 2D array to be used for plotting
	return mesh_2d

#====================================================================================================
# function to produce raster image

def array2raster( raster_dir , processing_pol , file_date_tag , output_array ) :

	raster_name = 'raster_'+processing_pol+'_'+file_date_tag+'.tiff'
	path = raster_dir + raster_name

	print( f'-> raster name= {raster_name}')
	print( f'-> raster path= {path}')

	no_of_bands = 1
	datatype = gdal.GDT_Float32
	#epsg_code = 4326 # output coord-ref-sys

	raster_origin = ( llx , lly ) # unit? metere or degree? it can be either meter, or degrees --> ImportFromProj4(+units=m)

	pixelWidth = 1000 # meters
	pixelHeight = 1000 # meters

	Xorig = raster_origin[0]
	Yorig = raster_origin[1]

	rows = output_array.shape[0]
	cols = output_array.shape[1]

	geotransform = ( Xorig , pixelWidth , 0 , Yorig , 0 , pixelHeight ) # units? either meter or degrees

	# get the class of coordinate system
	srs = osr.SpatialReference()
	#srs.ImportFromEPSG( epsg_code )
	srs.ImportFromProj4( f'+proj=lcc +lat_0={ycent} +lon_0={xcent} +units=m ' )

	# Initialize driver & create file
	driver = gdal.GetDriverByName('GTiff')
	# create output raster/matrix dataseet to put data into it
	out_raster = driver.Create( path, cols, rows, no_of_bands , datatype ) # (path, cols, rows, bands, dtype-> GDAL data type arg)

	out_raster.SetGeoTransform(geotransform)  # Specify its coordinates
	out_raster.SetProjection( srs.ExportToWkt() )

	print( f'-> writing the raster ...')
	out_raster.GetRasterBand( no_of_bands ).WriteArray( output_array )  # write my array to the raster
	out_raster.FlushCache()

	# Once we're done, close the dataset properly
	out_raster = None

	del out_raster

# def array2raster( input_array , NCOLS , NROWS , fig_dir ) :
# 	" array > raster"

# 	raster_file_name = 'co_raster_test.tif'
# 	path_to_raster_file = fig_dir
# 	raster_full_path = path_to_raster_file+raster_file_name

# 	print( f'-> creating raster file= {raster_full_path}')

# 	transformation_matrix = from_origin( ulx , uly , pixel_size , pixel_size )  # coordinates of uppper-left corner

# 	new_dataset = rasterio.open ( raster_full_path , 'w' , driver='GTiff' , height=NROWS , width=NCOLS , count=1 , dtype=str(input_array.dtype) , crs= '+proj=lcc' , transform=transformation_matrix )  # a dataset to store our grid

# 	print( f'-> writing raster file ...')
# 	new_dataset.write( input_array , 1 )
# 	new_dataset.close()

#====================================================================================================
### function to

def function_cell_24hr_timeSeries_singlePOL ( aconc_open , cmaq_pol , lay , row , col ):  # the order of argumenrs is important when input.
	" returns 24-hr time series of singlePOL"

	cell_24hr_series_list = []
	# extract all 24 t-step
	cell_24hr_series_list = aconc_open.variables[ cmaq_pol ][ : , lay , row , col ]
	# change daily list to daily np.Array
	cell_24hr_series_array = np.array( cell_24hr_series_list )
	# # get the mean of each cell
	# daily_cell_mean_for_singlePOL = cell_24hr_series_array.mean()

	# delete daily list
	del cell_24hr_series_list

	# function returns mean of the pollutant for each cell
	return cell_24hr_series_array

#====================================================================================================
### function to

def function_pm25_daily_cell_tseries ( aconc_open , pmdiag_open , lay , row , col ) : # arg are the variables that are defined insdie this function
	" returns daily timeseries for pm2.5 for each cell"

	print( f'-> processing row= {row} and col= {col}' )

	# loop inside 24 time-steps and extract pm concentrations
	# extract PM2.5 species from input files
	#print('-> extracting several species from CMAQ files for pm2.5 ...')
	# species from aconc [1]
	AH3OPI = aconc_open.variables['AH3OPI'][:,lay,row,col]
	#print( f'-> shape of AORGCJ= {AH3OPI.shape }')
	#print(f'-> type of AH3= { type(AH3OPI) }')
	AH3OPI = np.array(AH3OPI)
	# print(f'-> AH3OPI: {AH3OPI}')
	# print(f'-> AH3OPI mean= {AH3OPI_mean}')

	AH3OPJ = aconc_open.variables['AH3OPJ'][:,lay,row,col]
	AH3OPJ = np.array(AH3OPJ)

	AH3OPK = aconc_open.variables['AH3OPK'][:,lay,row,col]
	AH3OPK = np.array(AH3OPK)

	ACLI = aconc_open.variables['ACLI'][:,lay,row,col]
	ACLI = np.array(ACLI)

	ACLJ = aconc_open.variables['ACLJ'][:,lay,row,col]
	ACLJ = np.array(ACLJ)

	ACLK = aconc_open.variables['ACLK'][:,lay,row,col]
	ACLK = np.array(ACLK)

	AECI = aconc_open.variables['AECI'][:,lay,row,col]
	AECI = np.array(AECI)

	AECJ = aconc_open.variables['AECJ'][:,lay,row,col]
	AECJ = np.array(AECJ)

	ANAI = aconc_open.variables['ANAI'][:,lay,row,col]
	ANAI = np.array(ANAI)

	ANAJ = aconc_open.variables['ANAJ'][:,lay,row,col]
	ANAJ = np.array(ANAJ)

	AMGJ = aconc_open.variables['AMGJ'][:,lay,row,col]
	AMGJ = np.array(AMGJ)

	AKJ = aconc_open.variables['AKJ'][:,lay,row,col]
	AKJ = np.array(AKJ)

	ACAJ = aconc_open.variables['ACAJ'][:,lay,row,col]
	ACAJ = np.array(ACAJ)

	ANH4I = aconc_open.variables['ANH4I'][:,lay,row,col]
	ANH4I = np.array(ANH4I)

	ANH4J = aconc_open.variables['ANH4J'][:,lay,row,col]
	ANH4J = np.array(ANH4J)

	ANO3I = aconc_open.variables['ANO3I'][:,lay,row,col]
	ANO3I = np.array(ANO3I)

	ANO3J = aconc_open.variables['ANO3J'][:,lay,row,col]
	ANO3J = np.array(ANO3J)

	ASOIL = aconc_open.variables['ASOIL'][:,lay,row,col]
	ASOIL = np.array(ASOIL)

	ASO4I = aconc_open.variables['ASO4I'][:,lay,row,col]
	ASO4I = np.array(ASO4I)

	ASO4J = aconc_open.variables['ASO4J'][:,lay,row,col]
	ASO4J = np.array(ASO4J)

	ALVPO1I = aconc_open.variables['ALVPO1I'][:,lay,row,col]
	ALVPO1I = np.array(ALVPO1I)

	ASVPO1I = aconc_open.variables['ASVPO1I'][:,lay,row,col]
	ASVPO1I = np.array(ASVPO1I)

	ASVPO2I = aconc_open.variables['ASVPO2I'][:,lay,row,col]
	ASVPO2I = np.array(ASVPO2I)

	ALVOO1I = aconc_open.variables['ALVOO1I'][:,lay,row,col]
	ALVOO1I = np.array(ALVOO1I)

	ALVOO2I = aconc_open.variables['ALVOO2I'][:,lay,row,col]
	ALVOO2I = np.array(ALVOO2I)

	ASVOO1I = aconc_open.variables['ASVOO1I'][:,lay,row,col]
	ASVOO1I = np.array(ASVOO1I)

	ASVOO2I = aconc_open.variables['ASVOO2I'][:,lay,row,col]
	ASVOO2I = np.array(ASVOO2I)

	ALVPO1J = aconc_open.variables['ALVPO1J'][:,lay,row,col]
	ALVPO1J = np.array(ALVPO1J)

	ASVPO1J = aconc_open.variables['ASVPO1J'][:,lay,row,col]
	ASVPO1J = np.array(ASVPO1J)

	ASVPO2J = aconc_open.variables['ASVPO2J'][:,lay,row,col]
	ASVPO2J = np.array(ASVPO2J)

	ASVPO3J = aconc_open.variables['ASVPO3J'][:,lay,row,col]
	ASVPO3J = np.array(ASVPO3J)

	AIVPO1J = aconc_open.variables['AIVPO1J'][:,lay,row,col]
	AIVPO1J = np.array(AIVPO1J)

	AXYL1J = aconc_open.variables['AXYL1J'][:,lay,row,col]
	AXYL1J = np.array(AXYL1J)

	AXYL2J = aconc_open.variables['AXYL2J'][:,lay,row,col]
	AXYL2J = np.array(AXYL2J)

	AXYL3J = aconc_open.variables['AXYL3J'][:,lay,row,col]
	AXYL3J = np.array(AXYL3J)

	ATOL1J = aconc_open.variables['ATOL1J'][:,lay,row,col]
	ATOL1J = np.array(ATOL1J)

	ATOL2J = aconc_open.variables['ATOL2J'][:,lay,row,col]
	ATOL2J = np.array(ATOL2J)

	ATOL3J = aconc_open.variables['ATOL3J'][:,lay,row,col]
	ATOL3J = np.array(ATOL3J)

	ABNZ1J = aconc_open.variables['ABNZ1J'][:,lay,row,col]
	ABNZ1J = np.array(ABNZ1J)

	ABNZ2J = aconc_open.variables['ABNZ2J'][:,lay,row,col]
	ABNZ2J = np.array(ABNZ2J)

	ABNZ3J = aconc_open.variables['ABNZ3J'][:,lay,row,col]
	ABNZ3J = np.array(ABNZ3J)

	AISO1J = aconc_open.variables['AISO1J'][:,lay,row,col]
	AISO1J = np.array(AISO1J)

	AISO2J = aconc_open.variables['AISO2J'][:,lay,row,col]
	AISO2J = np.array(AISO2J)

	AISO3J = aconc_open.variables['AISO3J'][:,lay,row,col]
	AISO3J = np.array(AISO3J)

	ATRP1J = aconc_open.variables['ATRP1J'][:,lay,row,col]
	ATRP1J = np.array(ATRP1J)

	ATRP2J = aconc_open.variables['ATRP2J'][:,lay,row,col]
	ATRP2J = np.array(ATRP2J)

	ASQTJ = aconc_open.variables['ASQTJ'][:,lay,row,col]
	ASQTJ = np.array(ASQTJ)

	AALK1J = aconc_open.variables['AALK1J'][:,lay,row,col]
	AALK1J = np.array(AALK1J)

	AALK2J = aconc_open.variables['AALK2J'][:,lay,row,col]
	AALK2J = np.array(AALK2J)

	AORGCJ = aconc_open.variables['AORGCJ'][:,lay,row,col]
	AORGCJ = np.array(AORGCJ)

	AOLGBJ = aconc_open.variables['AOLGBJ'][:,lay,row,col]
	AOLGBJ = np.array(AOLGBJ)

	AOLGAJ = aconc_open.variables['AOLGAJ'][:,lay,row,col]
	AOLGAJ = np.array(AOLGAJ)

	APAH1J = aconc_open.variables['APAH1J'][:,lay,row,col]
	APAH1J = np.array(APAH1J)

	APAH2J = aconc_open.variables['APAH2J'][:,lay,row,col]
	APAH2J = np.array(APAH2J)

	APAH3J = aconc_open.variables['APAH3J'][:,lay,row,col]
	APAH3J = np.array(APAH3J)

	ALVOO1J = aconc_open.variables['ALVOO1J'][:,lay,row,col]
	ALVOO1J = np.array(ALVOO1J)

	ALVOO2J = aconc_open.variables['ALVOO2J'][:,lay,row,col]
	ALVOO2J = np.array(ALVOO2J)

	ASVOO1J = aconc_open.variables['ASVOO1J'][:,lay,row,col]
	ASVOO1J = np.array(ASVOO1J)

	ASVOO2J = aconc_open.variables['ASVOO2J'][:,lay,row,col]
	ASVOO2J = np.array(ASVOO2J)

	ASVOO3J = aconc_open.variables['ASVOO3J'][:,lay,row,col]
	ASVOO3J = np.array(ASVOO3J)

	APCSOJ = aconc_open.variables['APCSOJ'][:,lay,row,col]
	APCSOJ = np.array(APCSOJ)

	AALJ = aconc_open.variables['AALJ'][:,lay,row,col]
	AALJ = np.array(AALJ)

	ASIJ = aconc_open.variables['ASIJ'][:,lay,row,col]
	ASIJ = np.array(ASIJ)

	AFEJ = aconc_open.variables['AFEJ'][:,lay,row,col]
	AFEJ = np.array(AFEJ)

	ATIJ = aconc_open.variables['ATIJ'][:,lay,row,col]
	ATIJ = np.array(ATIJ)

	AOTHRI = aconc_open.variables['AOTHRI'][:,lay,row,col]
	AOTHRI = np.array(AOTHRI)

	AOTHRJ = aconc_open.variables['AOTHRJ'][:,lay,row,col]
	AOTHRJ = np.array(AOTHRJ)

	ACORS = aconc_open.variables['ACORS'][:,lay,row,col]
	ACORS = np.array(ACORS)

	ASEACAT = aconc_open.variables['ASEACAT'][:,lay,row,col]
	ASEACAT = np.array(ASEACAT)

	ASO4K = aconc_open.variables['ASO4K'][:,lay,row,col]
	ASO4K = np.array(ASO4K)

	ANO3K = aconc_open.variables['ANO3K'][:,lay,row,col]
	ANO3K = np.array(ANO3K)

	ANH4K = aconc_open.variables['ANH4K'][:,lay,row,col]
	ANH4K = np.array(ANH4K)

	AMNJ = aconc_open.variables['AMNJ'][:,lay,row,col]
	AMNJ = np.array(AMNJ)

	# species from pmdiag [3]
	PM25AT = pmdiag_open.variables['PM25AT'][:,lay,row,col]
	PM25AT = np.array(PM25AT)

	PM25AC = pmdiag_open.variables['PM25AC'][:,lay,row,col]
	PM25AC = np.array(PM25AC)

	PM25CO = pmdiag_open.variables['PM25CO'][:,lay,row,col]
	PM25CO = np.array(PM25CO)

	 # species calculated inside SpecDef file [0]
	 # perform arithmetic operations on arrays
	ANAK = 0.8373*ASEACAT + 0.0626*ASOIL + 0.0023*ACORS
	AMGK = 0.0997*ASEACAT + 0.0170*ASOIL + 0.0032*ACORS
	AKK = 0.0310*ASEACAT + 0.0242*ASOIL + 0.0176*ACORS
	ACAK = 0.0320*ASEACAT + 0.0838*ASOIL + 0.0562*ACORS

	APOCI = ALVPO1I /1.39 + ASVPO1I /1.32 + ASVPO2I /1.26
	ASOCI = ALVOO1I /2.27 + ALVOO2I /2.06 + ASVOO1I /1.88 + ASVOO2I /1.73
	APOCJ = ALVPO1J /1.39 + ASVPO1J /1.32 + ASVPO2J /1.26 +ASVPO3J /1.21 + AIVPO1J /1.17
	ASOCJ = AXYL1J /2.42  + AXYL2J /1.93  + AXYL3J /2.30 \
	+ATOL1J /2.26  + ATOL2J /1.82  + ATOL3J /2.70 \
	+ABNZ1J /2.68  + ABNZ2J /2.23  + ABNZ3J /3.00 \
	+AISO1J /2.20  + AISO2J /2.23  + AISO3J /2.80 \
	+ATRP1J /1.84  + ATRP2J /1.83  + ASQTJ /1.52  \
	+AALK1J /1.56  + AALK2J /1.42                   \
	+AORGCJ /2.00  + AOLGBJ /2.10  + AOLGAJ /2.50 \
	+APAH1J /1.63  + APAH2J /1.49  + APAH3J /1.77 \
	+ALVOO1J /2.27 + ALVOO2J /2.06 + ASVOO1J /1.88\
	+ASVOO2J /1.73 + ASVOO3J /1.60 + APCSOJ  /2.00

	APOMI = ALVPO1I  + ASVPO1I  + ASVPO2I
	ASOMI = ALVOO1I  + ALVOO2I  + ASVOO1I  + ASVOO2I
	APOMJ = ALVPO1J  + ASVPO1J  + ASVPO2J  +ASVPO3J  + AIVPO1J
	ASOMJ = AXYL1J   + AXYL2J   + AXYL3J   + ATOL1J  \
	+ATOL2J   + ATOL3J   + ABNZ1J   + ABNZ2J  \
	+ABNZ3J   + AISO1J   + AISO2J   + AISO3J  \
	+ATRP1J   + ATRP2J   + ASQTJ    + AALK1J  \
	+AALK2J   + APAH1J   + APAH2J   + APAH3J  \
	+AORGCJ   + AOLGBJ   + AOLGAJ               \
	+ALVOO1J  + ALVOO2J  + ASVOO1J  + ASVOO2J \
	+ASVOO3J  + APCSOJ

	AOCI = APOCI + ASOCI
	AOCJ = APOCJ + ASOCJ
	AOMI = APOMI    + ASOMI
	AOMJ = APOMJ    + ASOMJ
	ASOILJ = 2.20*AALJ + 2.49*ASIJ + 1.63*ACAJ + 2.42*AFEJ + 1.94*ATIJ
	ATOTI = ASO4I + ANO3I + ANH4I + ANAI + ACLI + AECI + AOMI + AOTHRI
	ATOTJ = ASO4J + ANO3J + ANH4J + ANAJ + ACLJ + AECJ + AOMJ + AOTHRJ + AFEJ + ASIJ + ATIJ + ACAJ + AMGJ + AMNJ + AALJ + AKJ
	ATOTK = ASOIL + ACORS + ASEACAT + ACLK + ASO4K + ANO3K + ANH4K
# !! PM2.5 species computed using modeled size distribution,
# reference: https://github.com/USEPA/CMAQ/blob/5.2/CCTM/src/MECHS/cb6r3_ae6_aq/SpecDef_cb6r3_ae6_aq.txt
	PM25_HP      = (AH3OPI * PM25AT + AH3OPJ * PM25AC + AH3OPK * PM25CO) * (1.0/19.0)
	PM25_CL      = ACLI * PM25AT + ACLJ * PM25AC + ACLK * PM25CO
	PM25_EC      = AECI * PM25AT + AECJ * PM25AC
	PM25_NA      = ANAI * PM25AT + ANAJ * PM25AC + ANAK * PM25CO
	PM25_MG      = AMGJ * PM25AC + AMGK * PM25CO
	PM25_K       = AKJ * PM25AC + AKK * PM25CO
	PM25_CA      = ACAJ * PM25AC + ACAK * PM25CO
	PM25_NH4     = ANH4I * PM25AT + ANH4J * PM25AC + ANH4K * PM25CO
	PM25_NO3     = ANO3I * PM25AT + ANO3J * PM25AC + ANO3K * PM25CO
	PM25_OC      = AOCI * PM25AT + AOCJ * PM25AC
	PM25_OM      = AOMI * PM25AT + AOMJ * PM25AC
	PM25_SOIL    = ASOILJ * PM25AC + ASOIL * PM25CO
	PM25_SO4     = ASO4I * PM25AT + ASO4J * PM25AC + ASO4K * PM25CO
	PM25_TOT     = ATOTI * PM25AT + ATOTJ * PM25AC + ATOTK * PM25CO
	PM25_UNSPEC1 = PM25_TOT - ( PM25_CL + PM25_EC + PM25_NA + PM25_NH4 + PM25_NO3 + PM25_OC + PM25_SOIL + PM25_SO4 )

	# now sum all species to get hourly PM2.5 concentratiosn
	pm25_cell_daily_tseries = PM25_HP + PM25_CL + PM25_EC + PM25_NA + PM25_MG + PM25_K + PM25_CA + \
					PM25_NH4 + PM25_NO3 + PM25_OC + PM25_OM + PM25_SOIL + PM25_SO4 + PM25_TOT + PM25_UNSPEC1

	# function returns the mean of pm2.5 for each cell
	return pm25_cell_daily_tseries

#====================================================================================================
# start running main()

if __name__ == '__main__' :  # __name__ is special variable by python interpreter; order of func is not importent anymore

	main()
