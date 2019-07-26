#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################################
# author: Ehsan Mosadegh
# usage: to analyze CMAQ output files spatially
# date: June 1, 2019
# email: ehsan.mosadegh@gmail.com, ehsanm@dri.edu
# NOTES: 
# in all print statements attach equal sign "=" to its preceding word for future grep searching.
############################################################

#====================================================================================================
# import libraries

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
#from osgeo import gdal, gdal_array, osr , ogr
from osgeo import gdal
import ogr, os, osr
from os import environ
#import rasterio
#from rasterio.transform import from_origin
import time

#====================================================================================================
# main

def main() :

	#====================================================================================================
	# run time setting

	### get the starting time
	start = time.time()

	cctm_process = environ.get('CCTM_PROCESS')  #'dep' 				# 'atm' or 'dep' 
	dep_type= environ.get('DEP_TYPE') #'DRYDEP'

	### run time settings
	cmaq_file_month= environ.get('CMAQ_MONTH_NUMBER')																			#  07, 08, 	09,  10,  11
	sim_month= environ.get('CMAQ_MONTH_STRING')   																				# Jul, Aug, Sep, Oct, Nov
	cmaq_file_year= '2016'
	mcip_date_tag= '161001'

	scenario= environ.get('LANDIS_SCENARIO') 	 																						# 1-5, baseline
	days_to_run_in_month= environ.get('DAYS_IN_MONTH_TO_RUN') 
	cmaq_pol= environ.get('CMAQ_POL')														# for plot title 'CO','PM2.5','NH3','O3','HNO3','NO2','SO2'
	processing_pollutant= environ.get('PROCESSING_POLLUTANT') # 'pm2.5' OR 'single_pollutant'== nh3,o3,no2,no,co
	pol_unit= environ.get('POL_UNIT') 												#	'ppmV' or 'ug/m^3'
	include_pmdiag_file= 'yes' 					 											# 'yes' OR 'no'

	### spatial plot
	spatial_plotting= environ.get('SPATIAL_PLOTTING_KEY')		 	# yes or no
	plot_method= environ.get('PLOT_METHOD')										# 'single_plot' or 'diff_plot'
	colorbar_method= environ.get('COLOR_METHOD')							# 'zero_to_max' , 'min_to_max' , 'minus_abs_max_to_max'
	minus_abs_max_diffPlot= environ.get('MINUS_ABS_MAX_DIFF')
	abs_max_diffPlot= environ.get('ABS_MAX_DIFF')
	# my_vmin_for_singlePlot= -0.4
	# my_vmax_for_singlePlot= 0.4
	produce_raster= environ.get('PRODUCE_RASTER')							# 'yes' OR 'no'

	### set mapping parameters for spatial plotting
	mapping= 'no' # 'yes' OR 'no'
	lower_bound_mapping_conc= 0.0
	upper_bound_mapping_conc= 0.120

	### time-series plot
	timeseries_plotting= environ.get('TIMESERIES_PLOTTING') 	# yes or not

	platform= 'cluster'  # 'Mac' or 'cluster'
	storage= '10T' # 'personal' OR '10T'
	dpi_scale=300

	# run time setting
	#====================================================================================================
	
	#====================================================================================================
	# Basemap plot setting
	
	### center of domain
	xcent =-120.806 # degrees
	ycent =40.000 # degrees
	
	### domain size
	pixel_size = 1000 # meters
	NROWS = 265*pixel_size # meters
	NCOLS = 250*pixel_size # meters
	
	### lower-left corner
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
	NROWS_zoom = 100000		#265*1000 # meters
	NCOLS_zoom = 100000		#250*1000 # meters

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

	# Basemap plot setting
	#====================================================================================================

	#====================================================================================================
	# print the setting

	print( f'-> scenario= {scenario}')
	print( f'-> CMAQ month= {cmaq_file_month}')
	print( f'-> number of days to run= {days_to_run_in_month}')
	print( f'-> processing pollutant= {cmaq_pol}')
	print( f'-> CCTM process= {cctm_process}')
	if ( cctm_process == 'dep' ):
		print( f'-> deposition type is= {dep_type}')
	print( f'-> pollutat unit= {pol_unit}')
	print( f'-> CMAQ year= {cmaq_file_year}')
	print( f'-> platform is= {platform}')
	print( f'-> timeSeries plotting= {timeseries_plotting}')
	print( f'-> produce raster= {produce_raster}')
	print( f'-> spatial plotting= {spatial_plotting}')
	if ( spatial_plotting == 'yes'):
		print( f'-> processing method for spatial plot= {plot_method}')
		if (plot_method=='single_plot'):
			#print( f'-> for single plot: vmin= {my_vmin_for_singlePlot} and vmax= {my_vmax_for_singlePlot} ')
			print('-> NOTE: for single plot we plot for min/max value of the region')
		if (plot_method=='diff_plot'):
			print( f'-> NOTE: colorBar method for spatial diff plot= {colorbar_method}')
			if ( colorbar_method == 'minus_abs_max_to_max' ):
				print( f'-> for diff plot: minus absolute Max. values= {minus_abs_max_diffPlot}')
				print( f'-> for diff plot: plus absolute Max. values= {abs_max_diffPlot}')				
	print( f'-> mapping spatial data= {mapping}')

	print(" ")

	# print the setting
	#====================================================================================================

	#====================================================================================================
	# input directory setting
	
	if ( platform == 'Mac' ) :

		if ( storage == '10T') :

			home_dir = '/Volumes/USFSdata/'   # '/' at the end
			mcip_dir = home_dir+'mcip_files/'   # '/' at the end
			fig_dir = home_dir+'cmaq_figs/'  # '/' at the end
			cmaq_output_dir = home_dir+'cmaq_output/'
			raster_dir = home_dir+'raster_dir/'

		if ( storage == 'personal') :

			home_dir = '/Volumes/Ehsanm_DRI/cmaq_usfs/'   # '/' at the end
			mcip_dir = home_dir+'mcip_files/'   # '/' at the end
			fig_dir = home_dir+'cmaq_figs/'  # '/' at the end
			cmaq_output_dir = home_dir+'cmaq_output/'
			raster_dir = home_dir+'raster_dir/'

	elif ( platform == 'cluster' ) :

		home_dir = '/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/'
		mcip_dir = home_dir+'mcip_files/'
		fig_dir = home_dir+'cmaq_figs/'
		cmaq_output_dir = home_dir+'cmaq_output/'
		raster_dir = home_dir+'raster_dir/'

	else:

		print( '-> ERROR: specify running platform ' )
		print('-> exiting ...')
		raise SystemExit()

	### set input pathes
	input_path_scen = cmaq_output_dir + 'scen_' + scenario + '/'
	input_path_base = cmaq_output_dir + 'scen_baseline' + '/' 

	print('-> CMAQ input directory is:')
	print(input_path_scen)
	print(input_path_base)
	print(" ")
	print('-> MCIP input directory is:')
	print(mcip_dir)
	print(" ")

	# input directory setting
	#====================================================================================================

	#====================================================================================================
	# open MCIP file to get lon-lat of domain

	print('-> opening MCIP files ...')

	grid_dot_file_name='GRIDDOT2D_'+mcip_date_tag
	grid_cross_file_name='GRIDCRO2D_'+mcip_date_tag

	grid_dot_file=mcip_dir + grid_dot_file_name
	grid_cross_file=mcip_dir + grid_cross_file_name

	print(f'-> GRID_DOT file is= {grid_dot_file}')
	print(f'-> GRID_CRO file is= {grid_cross_file}')

	grid_dot_open= Dataset(grid_dot_file , 'r')
	grid_cross_open=Dataset( grid_cross_file , 'r')

	### get some info
	print('-> MCIP file info:')

	print('-----------------------------------------------------------')
	print(f'-> info for file {grid_dot_file_name} is=')
	print('-> dimension is= %s' %( str(grid_dot_open.variables['LOND'].dimensions ))) # f-string does not work for this syntax
	print('-> shape is= %s' %str(grid_dot_open.variables['LOND'].shape )) # and this
	print('-----------------------------------------------------------')
	print(f'-> info for file {grid_cross_file_name} is=')
	print('-> dimension is= %s' %( str(grid_cross_open.variables['LON'].dimensions ))) # f-string does not work for this syntax
	print('-> shape is= %s' %str(grid_cross_open.variables['LON'].shape )) # and this
	print('-----------------------------------------------------------')

	### extract lat and lon arrays from DOT file
	lon_dot_array=np.array( grid_dot_open.variables['LOND'][0,0,:,:] ) # select only rosws and cols for the 1st timestep and layer = [ tstep=0 , lay=0]
	lat_dot_array=np.array( grid_dot_open.variables['LATD'][0,0,:,:] )

	### extract lat and lon arrays from CROSS file
	lon_cross_array=np.array( grid_cross_open.variables['LON'][0,0,:,:] )
	lat_cross_array=np.array( grid_cross_open.variables['LAT'][0,0,:,:] )

	# get the min of each array as the origin of array
	array_origin_lat = lat_dot_array.min()
	array_origin_lon = lon_dot_array.min()

	print('-> closing MCIP files ...')
	grid_dot_open.close()
	grid_cross_open.close()
	print(" ")

	# open MCIP file to get lon-lat of domain
	#====================================================================================================

	#====================================================================================================
# def function_3D_mesh_maker ( days_to_run_in_month , domain_rows , domain_cols , cmaq_file_month , scenario , input_path_scen , input_path_base ) :
# 	" processes the files and days and returns two 3D meshes"

	print('-> month of analysis is=' , cmaq_file_month )

	# ### define the monthly array mesh

	#monthly_tseries_tensor_atm_scen = np.ndarray( shape=( 24 , domain_rows , domain_cols ) )
	monthly_tseries_tensor_atm_scen = np.empty( shape=( 0 , domain_rows , domain_cols ) ) # use zero for concatenate method
	monthly_tseries_tensor_atm_base = np.empty( shape=( 0 , domain_rows , domain_cols ) ) # zero means there is no cell in z-dir
	monthly_tseries_tensor_dep_scen = np.empty( shape=( 0 , domain_rows , domain_cols ) ) # zero means there is no cell in z-dir
	monthly_tseries_tensor_dep_base = np.empty( shape=( 0 , domain_rows , domain_cols ) ) # zero means there is no cell in z-dir

	## create a day list for a month to create file-date-tag, use an argument-unpacking operator * to unpack the list
	days_to_run_in_month= int(days_to_run_in_month) # to change type from str as an env.var. to int
	day_list = [*range( 1 , days_to_run_in_month+1 , 1)] # don't forget the [] around range function to create the list

	### to run for specific day
	# day_list = [21]  # use the favorite day
	#====================================================================================================

	#====================================================================================================
	# traverse the list for each day

	for day_of_the_month in day_list :

		print('-----------------------------------------------------------')
		print( " ")
		print('-> no. of days to process are:')
		print(day_list)
		print(" ")
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
		
		#====================================================================================================
		# define input files

		aconc_scen = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+scenario+'_mpi_standard_'+file_date_tag+'.nc'
		pmdiag_scen = 'CCTM_PMDIAG_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+scenario+'_mpi_standard_'+file_date_tag+'.nc'
		dep_scen = 'CCTM_'+dep_type+'_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+scenario+'_mpi_standard_'+file_date_tag+'.nc'

		aconc_base = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonpt_baseline_AgBioNonpt_mpi_standard_'+file_date_tag+'.nc'
		pmdiag_base = 'CCTM_PMDIAG_v52_CA_WRF_1km_griddedAgBioNonpt_baseline_AgBioNonpt_mpi_standard_'+file_date_tag+'.nc'
		dep_base = 'CCTM_'+dep_type+'_v52_CA_WRF_1km_griddedAgBioNonpt_baseline_AgBioNonpt_mpi_standard_'+file_date_tag+'.nc'

		# define input files
		#====================================================================================================

		#====================================================================================================
		# define input files based on each processing method, plot, and pollutant

		if ( processing_pollutant == 'single_pollutant') : # we need only 1 file: "aconc"
			if ( plot_method == 'single_plot' ) :

				if ( cctm_process == 'atm') :

					print('-> opening/reading CMAQ files:')
					print( f'-> setting path for {cmaq_pol} for {cctm_process} and single_plot...')
					aconc_input_scen = input_path_scen + aconc_scen
					aconc_open_scen = Dataset( aconc_input_scen , 'r' )
					print( aconc_input_scen )

				if ( cctm_process == 'dep') :

					print( f'-> setting path for {cmaq_pol} for {cctm_process} and single_plot...')
					dep_input_scen = input_path_scen + dep_scen
					dep_open_scen = Dataset( dep_input_scen , 'r')
					print( dep_input_scen )

				print(" ")

			# for difference between 2 scenarios
			elif ( plot_method == 'diff_plot' ) :  # open 2 netcdf files: "aconc_scen" and "aconc_baseline"
				# define input files
				if ( cctm_process == 'atm' ) :

					print( f'-> setting path for {cmaq_pol} {cctm_process} and diff_plot...')

					aconc_input_scen = input_path_scen + aconc_scen
					aconc_input_base = input_path_base + aconc_base

					aconc_open_scen = Dataset( aconc_input_scen , 'r' )
					aconc_open_base = Dataset( aconc_input_base , 'r' )

					print( f'-> opening/reading CMAQ files for {cmaq_pol} {cctm_process} ')

					print( aconc_input_scen )
					print( aconc_input_base )

				if ( cctm_process == 'dep' ) :

					print( f'-> setting path for {cmaq_pol} {cctm_process} and diff_plot...')

					dep_input_scen = input_path_scen + dep_scen
					dep_input_base = input_path_base + dep_base

					dep_open_scen = Dataset( dep_input_scen , 'r' )
					dep_open_base = Dataset( dep_input_base , 'r' )

					print( f'-> opening/reading CMAQ files for {cmaq_pol} {cctm_process} ')

					print( dep_input_scen )
					print( dep_input_base )

				print(" ")

			else:

				print( '-> WARNING: define/check single processing POL and method first! ')
				print('-> exiting ...')
				raise SystemExit()

		elif ( processing_pollutant == 'pm2.5' ) :   # we need 2 files: "aconc" and "pmdiag" files

			if ( plot_method == 'single_plot' ) :

				print('-> setting path for PM2.5 and single_plot...')
				# define input files
				aconc_input_scen = input_path_scen + aconc_scen
				pmdiag_input_scen = input_path_scen + pmdiag_scen

				print('-> opening/reading CMAQ files for pm2.5:')
				print( aconc_input_scen )
				print( pmdiag_input_scen )
				print(" ")

				# define netcdf file
				aconc_open_scen = Dataset( aconc_input_scen , 'r' )
				pmdiag_open_scen = Dataset( pmdiag_input_scen , 'r' )

			elif ( plot_method == 'diff_plot' ) :

				print('-> setting path for PM2.5 and diff_plot...')
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
				print(" ")

			else:
				print('-> ERROR: select single or diff method! ')
				print('-> exiting ...')
				raise SystemExit()

		else:

			print( '-> ERROR: define/check POL first! ')
			print('-> exiting ...')
			raise SystemExit()

		# define input files based on each processing method and pollutant
		#====================================================================================================

		#====================================================================================================
		# process input files for each cell with a specific function

		# single plot

		if ( plot_method == 'single_plot' ) :

			print('-> processing for single plot ...')
			### create an empty tensor for each cell and day as container of daily 24-hr t-step concentrations
			daily_tseries_tensor_atm_scen = np.empty ( shape= ( 24 , domain_rows , domain_cols ))  # when assignin gby index, z-dim should be 1.
			daily_tseries_tensor_dep_scen = np.empty ( shape= ( 24 , domain_rows , domain_cols ))  # when assignin gby index, z-dim should be 1.

			print(f'-> shape of initial daily tseries tensor - atm scen - before filling = {daily_tseries_tensor_atm_scen.shape }')
			print(f'-> shape of initial daily tseries tensor - dep scen - before filling = {daily_tseries_tensor_dep_scen.shape }')
			# print(f'-> chk: initial empty daily tseries tensor- atm scen =')
			# print(daily_tseries_tensor_atm_scen[:1,:1])
			# print(f'-> chk: initial empty daily tseries tensor- dep scen =')
			# print(daily_tseries_tensor_dep_scen[:10])

			print(f'-> traversing each cell to extract species for {plot_method} ...')
			### traverse each cell in the C-storing style for each day: row and then col
			for row in range( 0 , domain_rows , 1 ) :
				for col in range( 0 , domain_cols , 1 ) :

					if ( processing_pollutant == 'single_pollutant' ) :

						#print( f'-> extracting cell for {cmaq_pol} - singlePlot - at row= {row} and col={col} ... ' )
						if ( cctm_process == 'atm' ) :

							cell_24hr_tseries_for_singlePol = function_cell_24hr_timeSeries_singlePOL( aconc_open_scen , cmaq_pol , lay , row , col )
							daily_tseries_tensor_atm_scen [:,row,col]  = cell_24hr_tseries_for_singlePol

						elif ( cctm_process == 'dep') :

							#print(f'-> we are going to use time-series function for dep ... ')
							cell_24hr_tseries_for_singlePol_dep = function_cell_24hr_timeSeries_singlePOL( dep_open_scen , cmaq_pol , lay , row , col )
							daily_tseries_tensor_dep_scen [:,row,col]  = cell_24hr_tseries_for_singlePol_dep

						else:
							print( '-> WARNING: define/check process type first! ')
							print('-> exiting ...')
							raise SystemExit()

					elif ( processing_pollutant == 'pm2.5') :

						print( f'-> inside loop: extracting cell for pm2.5 at row= {row} and col={col} ... ' )

						cell_24hr_tseries_for_pm25 = function_pm25_daily_cell_tseries( include_pmdiag_file , aconc_open_scen , pmdiag_open_scen , lay , row , col )
						daily_tseries_tensor_atm_scen [:,row,col] = cell_24hr_tseries_for_pm25 # fill tensor for all cells in domain
						# for dep files
						#
					else:

						print( '-> WARNING: define/check single POL or pm2.5 settings and processing method first! ')
						print('-> exiting ...')
						raise SystemExit()

			if ( produce_raster == 'yes' ) :

				print( " " )
				print( f'-> producing raster file ...')
				print( f'-> note: blah blah')
				### array to raster and write out a raster 
				daily_2d_array_scen_and_mapped = function_3Dto2D( domain_rows , domain_cols , daily_tseries_tensor_atm_scen )
				daily_2d_array_scen = daily_2d_array_scen_and_mapped[0]
				array2raster( raster_dir , cmaq_pol , file_date_tag , daily_2d_array_scen , array_origin_lon , array_origin_lat )
				# print( f'-> raster file = { output_raster }')
				# print( " " )


			### after all days are extracted, now we concatenate the daily timeseries tensor to monthly timeseries tensor
			if ( cctm_process == 'atm') :

				print( f'-> note: shape of initial daily tseries tensor - atm scen - after filling is = {daily_tseries_tensor_atm_scen.shape}')
				print( f'-> note: shape of monthly tseries tensor before filling is = {monthly_tseries_tensor_atm_scen.shape}')
				monthly_tseries_tensor_atm_scen = np.concatenate( ( monthly_tseries_tensor_atm_scen ,  daily_tseries_tensor_atm_scen ) , axis=0 )

			if( cctm_process == 'dep') :

				print( f'-> note: shape of initial daily tseries tensor - dep scen - after filling is = {daily_tseries_tensor_dep_scen.shape}')
				# print( f'-> daily tseries tensor after filling is = ')
				# print( daily_tseries_tensor_dep_scen )
				print( f'-> note: shape of monthly tseries tensor before filling is = {monthly_tseries_tensor_dep_scen.shape}')

				monthly_tseries_tensor_dep_scen = np.concatenate( ( monthly_tseries_tensor_dep_scen ,  daily_tseries_tensor_dep_scen ) , axis=0 )

		# diff plot

		elif ( plot_method == 'diff_plot' ) :

			print('-> processing for diff plot ...')
			### create an empty tensor for each day as container of daily 24-hr t-step concentrations and then fill each cell up for each cell
			daily_tseries_tensor_atm_scen = np.empty ( shape=( 24 , domain_rows , domain_cols ) )
			daily_tseries_tensor_atm_base = np.empty ( shape=( 24 , domain_rows , domain_cols ) )
			daily_tseries_tensor_dep_scen = np.empty ( shape=( 24 , domain_rows , domain_cols ) )
			daily_tseries_tensor_dep_base = np.empty ( shape=( 24 , domain_rows , domain_cols ) )

			print('-> traversing each cell to extract species ...')
			### traverse each cell in the C-storing style for each day: row and then col
			for row in range( 0 , domain_rows , 1 ) :
				for col in range( 0 , domain_cols , 1 ) :

					if ( processing_pollutant == 'single_pollutant' ) :

						#print( f'-> extracting cell for {cmaq_pol} - diff - at row= {row} and col={col} ... ' )
						if ( cctm_process == 'atm' ) :

							# we extract tseries from each cell
							cell_24hr_timeSeries_array_scen = function_cell_24hr_timeSeries_singlePOL( aconc_open_scen , cmaq_pol , lay , row , col )
							cell_24hr_timeSeries_array_base = function_cell_24hr_timeSeries_singlePOL( aconc_open_base , cmaq_pol , lay , row , col )
							
							### fill daily tensors cells with daily time-series
							daily_tseries_tensor_atm_scen [: , row , col] = cell_24hr_timeSeries_array_scen
							daily_tseries_tensor_atm_base [: , row , col] = cell_24hr_timeSeries_array_base


						if ( cctm_process == 'dep' ) :

							# we extract tseries from each cell
							cell_24hr_timeSeries_array_scen = function_cell_24hr_timeSeries_singlePOL( dep_open_scen , cmaq_pol , lay , row , col )
							cell_24hr_timeSeries_array_base = function_cell_24hr_timeSeries_singlePOL( dep_open_base , cmaq_pol , lay , row , col )
							
							### fill daily tensors cells with daily time-series
							daily_tseries_tensor_dep_scen [: , row , col] = cell_24hr_timeSeries_array_scen
							daily_tseries_tensor_dep_base [: , row , col] = cell_24hr_timeSeries_array_base


					elif ( processing_pollutant == 'pm2.5' ) :

						print(" ")
						print('-----------------------------------------------------------')
						if ( cctm_process == 'atm' ) :

							### we extract timeseries 
							#print( f'-> inside loop: extracting cell for PM2.5 - diff - Landis - at row= {row} and col= {col} ... ' )
							cell_24hr_timeSeries_array_scen = function_pm25_daily_cell_tseries( include_pmdiag_file , aconc_open_scen , pmdiag_open_scen , lay , row , col )

							#print( f'-> inside loop: extracting cell for PM2.5 - diff - baseline - at row= {row} and col= {col} ... ' )
							cell_24hr_timeSeries_array_base = function_pm25_daily_cell_tseries( include_pmdiag_file , aconc_open_base , pmdiag_open_base , lay , row , col )

							### fill daily tensors cells with 24-hr time-series
							daily_tseries_tensor_atm_scen [: , row , col] = cell_24hr_timeSeries_array_scen
							daily_tseries_tensor_atm_base [: , row , col] = cell_24hr_timeSeries_array_base


						if ( cctm_process == 'dep' ) :

							print(f'-> NOTE: we do not need to calculate pm2.5 dep on the lake for now.\
								in future if we need to work on this selction, it takes a lot of time, needs a new function maybe ???')


						print('-----------------------------------------------------------')
						print(" ")

					else:
						pass



			### after all days are extracted, now we concatenate the daily timeseries tensor to monthly timeseries tensor
			if ( cctm_process == 'atm' ) :

				#print( f'-> add/pin each cell 24-hr t-series to monthly_tseries_tensor_atm_scen at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )
				monthly_tseries_tensor_atm_scen = np.concatenate( ( monthly_tseries_tensor_atm_scen ,  daily_tseries_tensor_atm_scen ) , axis=0 ) # along axis=0 of the tensor
				#print( f'-> add/pin each cell 24hr t-series to monthly_tseries_tensor_atm_base at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )
				monthly_tseries_tensor_atm_base = np.concatenate( ( monthly_tseries_tensor_atm_base ,  daily_tseries_tensor_atm_base ) , axis=0 ) # along axis=0 of the tensor

			if ( cctm_process == 'dep' ) :

				#print( f'-> add/pin each cell 24-hr t-series to monthly_tseries_tensor_atm_scen at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )
				monthly_tseries_tensor_dep_scen = np.concatenate( ( monthly_tseries_tensor_dep_scen ,  daily_tseries_tensor_dep_scen ) , axis=0 ) # along axis=0 of the tensor
				#print( f'-> add/pin each cell 24hr t-series to monthly_tseries_tensor_atm_base at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )
				monthly_tseries_tensor_dep_base = np.concatenate( ( monthly_tseries_tensor_dep_base ,  daily_tseries_tensor_dep_base ) , axis=0 ) # along axis=0 of the tensor

		else:
			print('-> ERROR: processing method NOT defined!')

		# process input files for each cell with a specific function
		#====================================================================================================

		#====================================================================================================
		# closing nc files for each day

		if ( processing_pollutant == 'single_pollutant' ) :
			if ( plot_method == 'single_plot' ) :
				if ( cctm_process == 'atm' ) :

					print('-> closing ACONC file ...')
					aconc_open_scen.close()

				if ( cctm_process == 'dep' ) :

					print('-> closing DEP file ...')
					dep_open_scen.close()

			elif ( plot_method == 'diff_plot' ) :
				if ( cctm_process == 'atm' ) :

					print('-> closing ACONC scen & baseline files ...')
					aconc_open_scen.close()
					aconc_open_base.close()

				if ( cctm_process == 'dep' ) :

					print('-> closing DEP scen & baseline files ...')
					dep_open_scen.close()
					dep_open_base.close()

			else:
				pass

		elif ( processing_pollutant == 'pm25' ) :

			if ( plot_method == 'single_plot' ) :

				print('-> closing ACONC and PMDIAG files ...')
				aconc_open_scen.close()
				pmdiag_open_scen.close()

			elif ( plot_method == 'diff_plot' ) :

				print('-> closing ACONC and PMDIAG files for both scen and baseline ...')
				aconc_open_scen.close()
				pmdiag_open_scen.close()
				aconc_open_base.close()
				pmdiag_open_base.close()

		else:
			pass

	# closing nc files for each day
	#====================================================================================================

	#====================================================================================================
	# change 3D to 2D array to make monthly_mean_2d_mesh: used for spatial plots; the final product for spatial plotting is: monthly_mean_2d_mesh

	# NOTE: we only need one scen and ase tensor, we direct different scen files to one intermediate file and work on that

	### region dictionary, starting point == origin == ll point of the bounding box, and then ul, ur, lr, ll
	regions_dict={

	'SouthTahoe':			[-120.05 , 38.88 , 8 , 8 		],
	'NorthTahoe':			[-120.05 , 39.217 , 8 , 8 	],
	'LakeTahoeBasin':	[-120.30 , 38.87 , 50 , 38 	],
	'Ndep':						[-120.152 , 38.93 , 37 , 20 ]  # row, col = range_in_raw and col 
	}


	# intermed file is directed to be used in spatial plotting; and the original tensor is used for time-series plotting
	# monthly_tseries_tensor_scen_intermed = monthly_tseries_tensor_atm_scen
	# monthly_tseries_tensor_scen_intermed_drydep = monthly_tseries_tensor_dep_scen
	if ( cctm_process == 'atm' ) :
		monthly_tseries_tensor_scen_intermediate = monthly_tseries_tensor_atm_scen  # change atm scen to intermediate 

		if ( plot_method == 'diff_plot') :
			monthly_tseries_tensor_baseline_intermediate = monthly_tseries_tensor_atm_base

	if ( cctm_process == 'dep' ) :
		monthly_tseries_tensor_scen_intermediate = monthly_tseries_tensor_dep_scen
		
		if ( plot_method == 'diff_plot') :
			monthly_tseries_tensor_baseline_intermediate = monthly_tseries_tensor_dep_base


	### before starting plot merhod, we define list for min/max
	min_list_for_all_regions = []
	max_list_for_all_regions = []


	if ( plot_method == 'single_plot' ) :


		if ( cctm_process == 'atm' ) :
			# for o3 and N and so2, change ppm to ppb for plotting
			if any([ cmaq_pol=='O3',cmaq_pol=='NH3',cmaq_pol=='HNO3',cmaq_pol=='NO2',cmaq_pol=='SO2' ]) :

				print( f'-> changing ppm to ppb in place for scenario pol={cmaq_pol} ...' )
				pol_unit = 'ppb'

				monthly_tseries_tensor_scen_intermediate = monthly_tseries_tensor_scen_intermediate * 1000

		if ( cctm_process == 'dep' ) :

				print( f'-> for {cmaq_pol} in {cctm_process} for single_plot the unit will be updated from "{pol_unit}" to "kg/hectare" ' )
				pol_unit = 'kg/hectare'


		### look at the 3D tensors for scenario
		print('-----------------------------------------------------------')
		print('-> 3D data mesh info:')
		print( f'-> monthly tensor from LANDIS scenario: dimensions= { monthly_tseries_tensor_scen_intermediate.ndim} and shape of data-mesh= { monthly_tseries_tensor_scen_intermediate.shape}' )

		#print( f'-> monthly tensor for LANDIS drydep scenario: dimensions= { monthly_tseries_tensor_scen_intermed_drydep.ndim} and shape of data-mesh= { monthly_tseries_tensor_scen_intermed_drydep.shape}' )
		print('-----------------------------------------------------------')

		### we only have one 3D tensor

		print('-> changing monthly 3D mesh of time-series to 2D for single plotting ...')

		monthly_mean_2d_mesh_and_mapped = function_3Dto2D( domain_rows , domain_cols , monthly_tseries_tensor_scen_intermediate , mapping )

		if ( mapping == 'yes' ) :
			print('-> mapping abs concentrations to mapped domain ...')

			monthly_mean_2d_mesh = monthly_mean_2d_mesh_and_mapped[1]

		else :
			print('-> using the absolute concentrations from 2D mesh array...')

			monthly_mean_2d_mesh = monthly_mean_2d_mesh_and_mapped[0]

		print( f'-> quick look at final product of 3D-to-2D for single plotting: monthly mean 2D mesh=')
		print(monthly_mean_2d_mesh)

		print( " ")
		print( '-> monthly mean mesh statistics:')

		print( f'-> shape= {monthly_mean_2d_mesh.shape } and dimension= {monthly_mean_2d_mesh.ndim }')

		regions_list = [ *regions_dict.keys() ]
		print(f'-> region list is= {regions_list}')

		for region in regions_list:

			ll_lon = regions_dict[region][0]
			ll_lat = regions_dict[region][1]
			range_in_raw = regions_dict[region][2]
			range_in_col = regions_dict[region][3]

			
			tuple_of_stats = function_stats_of_desired_region(monthly_mean_2d_mesh , timeseries_plotting , lon_cross_array , lat_cross_array , ll_lon , ll_lat , range_in_raw , range_in_col , region )

			print('-----------------------------------------------------------')
			print( f'-> stats: regionNameSingleMesh= 			{region} ')
			print( f'-> stats: minSingleMesh{region}= 		{ tuple_of_stats[0]	}')
			print( f'-> stats: meanSingleMesh{region}= 		{ tuple_of_stats[1]	}')
			print( f'-> stats: medianSingleMesh{region}= 	{ tuple_of_stats[2]	}')
			print( f'-> stats: stdSingleMesh{region}= 		{ tuple_of_stats[3]	}')
			print( f'-> stats: maxSingleMesh{region}= 		{ tuple_of_stats[4]	}')
			# print( f'-> row no. of max value{region}= { round(row_of_max_cell,6) }')
			# print( f'-> col no. of max value{region}= { round(col_of_max_cell,6) }')
			print('-----------------------------------------------------------')

			min_list_for_all_regions.append(tuple_of_stats[0])
			max_list_for_all_regions.append(tuple_of_stats[4])

		#print( f'-> applying the mapp function...')

		### to here

	elif ( plot_method == 'diff_plot' ) :

		if ( cctm_process == 'atm' ) :
			if any([cmaq_pol == 'O3' , cmaq_pol == 'NH3' , cmaq_pol == 'HNO3' , cmaq_pol == 'NO2' , cmaq_pol == 'SO2']) :

				print( f'-> changing ppm to ppb in place for scenario pol={cmaq_pol} ...' )
				pol_unit = 'ppb'

			monthly_tseries_tensor_scen_intermediate = 	monthly_tseries_tensor_scen_intermediate * 1000
			monthly_tseries_tensor_baseline_intermediate = 	monthly_tseries_tensor_baseline_intermediate * 1000

		if ( cctm_process == 'dep' ) :

			print( f'-> for {cmaq_pol} in {cctm_process} for diff_plot the unit will be updated from "{pol_unit}" to "kg/hectare" ' )
			pol_unit = 'kg/hectare'


		### look at the 3D meshes
		print('---------------------------------------------------------------')
		print('-> check 3D data mesh info:')
		print( f'-> LANDIS scenario, monthly mean tensor: dimensions= {monthly_tseries_tensor_scen_intermediate.ndim} and shape of tensor= {monthly_tseries_tensor_scen_intermediate.shape} ' )
		print( f'-> baseline, monthly mean tensor: dimensions= {monthly_tseries_tensor_baseline_intermediate.ndim} and shape of tensor= {monthly_tseries_tensor_baseline_intermediate.shape} ' )
		print('---------------------------------------------------------------')

		### we have 2 tensors
		print( f'-> first, changing monthly tensor of time-series to 2D mesh for diff-plotting for 3D LANDIS scenario ({scenario}) mesh ...')
		monthly_mean_2d_mesh_scen_and_mapped = function_3Dto2D( domain_rows , domain_cols , monthly_tseries_tensor_scen_intermediate , mapping )
		monthly_mean_2d_mesh_scen = monthly_mean_2d_mesh_scen_and_mapped[0]

		print('-> now, changing monthly tensor of time-series to 2D mesh for diff-plotting for 3D "baseline" scenario mesh ...')
		monthly_mean_2d_mesh_base_and_mapped = function_3Dto2D( domain_rows , domain_cols , monthly_tseries_tensor_baseline_intermediate , mapping )
		monthly_mean_2d_mesh_base = monthly_mean_2d_mesh_base_and_mapped[0]

		print ('-> now subtract 2 meshes to get the diff mesh for spatial plotting ')

		monthly_mean_2d_mesh = monthly_mean_2d_mesh_scen - monthly_mean_2d_mesh_base # if + then over-estimate; if - then underestimate
		print( " ")
		print( f'-> quick look at final product of 3D-to-2D for - diff - plotting: monthly mean 2D mesh=')

		print( f'-> monthly mesh scenario= ')
		print( monthly_mean_2d_mesh_scen )
		print( " ")
		print( f'-> monthly mesh baseline= ')
		print( monthly_mean_2d_mesh_base )
		print( " ")
		print( f'-> monthly mean diff-mesh= ')
		print( monthly_mean_2d_mesh )
		print( " ")
		print( '-> diff matrix statistics:')

		print( f'-> shape= {monthly_mean_2d_mesh.shape } and dimension= {monthly_mean_2d_mesh.ndim }')

		regions_list = [ *regions_dict.keys() ]
		print(f'-> region list is = {regions_list}')

		for region in regions_list:

			ll_lon = regions_dict[region][0]
			ll_lat = regions_dict[region][1]
			range_in_raw = regions_dict[region][2]
			range_in_col = regions_dict[region][3]


			tuple_of_stats = function_stats_of_desired_region(monthly_mean_2d_mesh , timeseries_plotting , lon_cross_array , lat_cross_array , ll_lon , ll_lat , range_in_raw , range_in_col , region )

			print('-----------------------------------------------------------')
			print( f'-> stats: regionNameDiffMesh= 			{region} ')
			print( f'-> stats: minDiffMesh{region}= 		{ tuple_of_stats[0] }')
			print( f'-> stats: meanDiffMesh{region}= 		{ tuple_of_stats[1] }')
			print( f'-> stats: medianDiffMesh{region}= 	{ tuple_of_stats[2]}')
			print( f'-> stats: stdDiffMesh{region}= 		{ tuple_of_stats[3]	}')
			print( f'-> stats: maxDiffMesh{region}= 		{ tuple_of_stats[4]	}')
			# print( f'-> row no. of max value{region}= { round(row_of_max_cell,6) }')
			# print( f'-> col no. of max value{region}= { round(col_of_max_cell,6) }')
			print('-----------------------------------------------------------')

			min_list_for_all_regions.append(tuple_of_stats[0])
			max_list_for_all_regions.append(tuple_of_stats[4])

	else:
		print('-> ERROR: check processing method for 3Dto2D ...')
		print('-> exiting ...')
		raise SystemExit()


	min_of_single_region = min(min_list_for_all_regions)
	max_of_single_region = max(max_list_for_all_regions)
	min_of_diff_region = min(min_list_for_all_regions)
	max_of_diff_region = max(max_list_for_all_regions)


	# change 3D to 2D array to make monthly_mean_2d_mesh: used for spatial plots; the final product for spatial plotting is: monthly_mean_2d_mesh
	#====================================================================================================
	# time-series plotting

	if ( timeseries_plotting == 'yes') :
		print('-> time-series plotting is YES, so we plot time-series...')

		row_ = row_of_max_cell
		col_ = col_of_max_cell
		conc_timeseries_list = monthly_tseries_tensor_atm_scen [ : , row_ , col_ ]

		if ( cmaq_pol == 'O3' ) :

			conc_timeseries_list = conc_timeseries_list*1000 # change to ppb

		x_ = [ *range(1, ((days_to_run_in_month*24)+1) , 1) ]  # x bar is no. of hours in aconc files
		y_ = conc_timeseries_list		#.tolist  # y-bar is hourly concentrations/timeseries

		print(f'-> size of x_ axis list= {len(x_)}')
		#print(f'-> x_ = {x_}' )
		print(f'-> size of y_ axis list= {len(y_)}')
		#print(f'-> y_ = {y_}' )

		plt.plot( x_ , y_ )
		#plt.show()
		xticks_position = [ i*24 for i in range(0 , days_to_run_in_month+1 , 1) ]
		xticks_labels = [ f'{i+1}' for i in range (0 , days_to_run_in_month , 1)]
		plt.xticks( xticks_position , xticks_labels , fontsize=6 )
		plt.yticks(fontsize=6)

		plt.xlabel(f' days in {sim_month}')
		plt.ylabel( f' {cmaq_pol} concentration')
		plt.title( f'time-series of {cmaq_pol} for {sim_month}, 2016')
		plt.grid(True)

		plot_name = cmaq_pol+'_timeseries'+'_scen_'+scenario+'_'+cmaq_file_year+'-'+cmaq_file_month+'_summed_'+str(days_to_run_in_month)+'_days'+'.png'
		saved_plot = fig_dir+plot_name

		plt.savefig( saved_plot , dpi=dpi_scale , format='png' )
		plt.close()
		print('-> time-series plot saved at:')
		print(saved_plot)
		print(" ")

	# time-series plotting
	#====================================================================================================

	#====================================================================================================
	# Spatial plotting with Basemap

	# use Basemap library and make spatial plots

	if ( spatial_plotting == 'yes') :

		print(f'-> spatial plotting for CCTM process= {cctm_process}...')

		### plot dots from grid coordinates of the dots
		#plt.plot( lon_mesh , lat_mesh , marker='.' , color='b' , linestyle= 'none' )

		# ### create a Basemap class/model instance for a specific projection
		# # basemap_instance = Basemap(projection='lcc' , lat_0=ycent , lon_0=xcent , height=NROWS , width=NCOLS , resolution='i') # , area_thresh=0.1) # latlon=True for when x and y are not in map proj. coordinates
		# theMap = Basemap(projection='lcc' ,
		# 													 llcrnrx=llx , llcrnry=lly , urcrnrx=urcornerx , urcrnry=urcornery ,
		# 													 lat_0=ycent , lon_0=xcent , height=NROWS , width=NCOLS ,
		# 													 resolution='f' , area_thresh=0.5)


		### create Basemap model instance from its class, it is a map that color mesh sits on it.
		theMap = Basemap(projection='lcc' , lat_0=ycent_zoom , lon_0=xcent_zoom , height=NROWS_zoom , width=NCOLS_zoom , resolution='f' , area_thresh=0.5)

		x_mesh, y_mesh = theMap(lon_dot_array , lat_dot_array) # order: x , y; Basemap model transforms lon/lat from degree to meter for LCC projection map; or just use latlon=True inpcolormesh function and use lat-lon values

		#basemap_instance.fillcontinents(lake_color='aqua')

		#my_levels = [ 0.02 , 0.05 ]
		#my_colors = ( 'g' , 'b' , 'r' )

		# set the color map
		if ( plot_method == 'single_plot' ) :

			color_mapping_function = 'Reds'

		if ( plot_method == 'diff_plot' ) :

			color_mapping_function = 'RdBu_r'

		### create a color mesh image from basemap model instance, the color mesh is constant, cos it is plotted from lon/lat values
		print(" ")
		print( '-> making the colormesh ...')
		print( f'-> shape of x_mesh= {x_mesh.shape }')
		print( f'-> shape of y_mesh= {y_mesh.shape }')
		print( f'-> shape of monthly_mean_diff_mesh= {monthly_mean_2d_mesh.shape }')
		print(" ")
		# define the image first
		colorImage = theMap.pcolormesh( x_mesh , y_mesh , monthly_mean_2d_mesh , cmap=color_mapping_function , shading='flat' )# , vmin=-5e-5 , vmax=5e-5 )

		lon_of_text , lat_of_text = theMap( xcent_zoom , ycent_zoom )
		plt.text( lon_of_text-45000 , lat_of_text+45000 , 'scen %s %s' %(scenario , sim_month) )

		#im2 = basemap_instance.pcolormesh(lon_mesh , lat_mesh , data_mesh , cmap=plt.cm.jet , shading='flat')

		# then, set the color limit
		if ( plot_method == 'single_plot' ) :
			print( f'-> plot method= {plot_method}')
			
			if ( mapping == 'yes' ) :
				plt.clim( lower_bound_mapping_conc , upper_bound_mapping_conc )

			### you can choose either of these options
			plt.clim( vmin=min_of_single_region ,  vmax=max_of_single_region )
			#plt.clim( vmin=my_vmin_for_singlePlot , vmax=my_vmax_for_singlePlot )


		if ( plot_method == 'diff_plot' ) :
			print( f'-> plot method= {plot_method}')
			
			if ( colorbar_method == 'zero_to_max' ) :
				
				print(f'-> colorbar method= {colorbar_method}')
				plt.clim( 0.0 ,  max_of_diff_region )
				print( f'-> plot the image for vmin= {0.0} and vmax= {max_of_diff_region}')

			if ( colorbar_method == 'min_to_max' ) :

				print(f'-> colorbar method= {colorbar_method}')
				plt.clim( min_of_diff_region ,  max_of_diff_region )
				print( f'-> plot the image for vmin= {min_of_diff_region} and vmax= {max_of_diff_region}')

			if ( colorbar_method == 'minus_abs_max_to_max' ) :

				print(f'-> colorbar method= {colorbar_method}')
				plt.clim( minus_abs_max_diffPlot ,  abs_max_diffPlot )
				print( f'-> plot the image for vmin= {minus_abs_max_diffPlot} and vmax= {abs_max_diffPlot}')
		print(" ")
			# plt.clim( my_vmin_for_singlePlot , vmax_mine )
			# print( f'-> plot the image for vmin={vmin_mine} and vmax={vmax_mine}')

		# set the map features
		theMap.drawmapboundary(color='k' ) #, fill_color='aqua')
		theMap.drawcoastlines(color = '0.15')
		theMap.drawstates()
		theMap.drawcounties(linewidth=0.5 , color='k' )

		#theMap.bluemarble()

		### create colorbar/legend in a seperate obj to play with it later
		colorbar = theMap.colorbar( colorImage , 'bottom' , size='4%' )

		colorbar_label = f' {cmaq_pol} concentration [{pol_unit}] ' 
		colorbar.set_label( label=colorbar_label , size=6 )

		colorbar.ax.tick_params( labelsize= 5 ) # provide access to the usual axis methods including tick formatting
		#cs = basemap_instance.contourf(lon_mesh , lat_mesh , data_mesh)
		#colorbar = basemap_instance.colorbar(cs, location='bottom')
		#plt.subplot( figsize=(10,10) )
		if ( plot_method == 'single_plot' ) :
			if ( cctm_process == 'atm' ) :
				plt.title(f' {cctm_process} {cmaq_pol} monthly mean concentrations for {sim_month}, {cmaq_file_year} - LANDIS scenario {scenario}' , fontsize=6 )

			if ( cctm_process == 'dep' ) :
				plt.title(f' {cmaq_pol} {dep_type} monthly mean concentrations for {sim_month}, {cmaq_file_year} - LANDIS scenario {scenario}' , fontsize=6 )

		elif ( plot_method == 'diff_plot' ) :
			if ( cctm_process == 'atm' ) :
				plt.title(f' {cctm_process} {cmaq_pol} monthly mean concentration for {sim_month}, scenario{scenario}_minus_baseline' , fontsize=6 )

			if ( cctm_process == 'dep' ) :
				plt.title(f' {cmaq_pol} {dep_type} monthly mean concentration for {sim_month}, scenario{scenario}_minus_baseline' , fontsize=6 )


		print(" ")

		# Spatial plotting with Basemap
		#====================================================================================================

		#====================================================================================================
		# save the plots

		### path for saving plots
		# print('-> fig directory is:')
		# print(fig_dir)

		### plot name
		if ( plot_method == 'single_plot' ) :
			if( cctm_process == 'atm' ) :

				fig_name = cctm_process+'_conc_'+cmaq_pol+'_monthlyMean_scen_'+scenario+'_singlePlot_month_'+cmaq_file_month+'_summed_'+str(days_to_run_in_month)+'_days'+'_dpi_'+str(dpi_scale)+'.png'
			
			if ( cctm_process == 'dep') :

				fig_name = dep_type+'_'+cmaq_pol+'_monthlyMean_scen_'+scenario+'_singlePlot_month_'+cmaq_file_month+'_summed_'+str(days_to_run_in_month)+'_days'+'_dpi_'+str(dpi_scale)+'.png'

		elif ( plot_method == 'diff_plot' ) :
			if( cctm_process == 'atm' ) :

				fig_name = cctm_process+'_conc_'+cmaq_pol+'_monthlyMean_scen_'+scenario+'_diff_from_baseline_month_'+cmaq_file_month+'_summed_'+str(days_to_run_in_month)+'_days_colorbarMethod_'+colorbar_method+'_dpi_'+str(dpi_scale)+'.png'

			if ( cctm_process == 'dep' ) :

				fig_name = dep_type+'_'+cmaq_pol+'_monthlyMean_scen_'+scenario+'_diff_from_baseline_month_'+cmaq_file_month+'_summed_'+str(days_to_run_in_month)+'_days_colorbarMethod_'+colorbar_method+'_dpi_'+str(dpi_scale)+'.png'
		else:
			pass

		### full path of the plot
		out_fig = fig_dir + fig_name
		print('-> output figure is stored at:')
		print(out_fig)
		### export the figure
		plt.savefig( out_fig , dpi=dpi_scale , format='png')
		### opens a window to show the results - after savefig
		#plt.show()
		### close the plot
		plt.close()

		# save the plots
		#====================================================================================================
		

	print('------------------------------------')
	end = time.time()
	print( f'-> run time of main function= { (( end - start ) / 60 ) :.2f} min' )  # f-string
	print('*** MAIN COMPLETED SUCCESSFULLY! ***')
	print('------------------------------------')

# main
#====================================================================================================


#====================================================================================================
# other functions
#====================================================================================================


#====================================================================================================
# function to map concentrations from abs value to map space

def mapping_function( abs_conc , lower_bound_mapping_conc , upper_bound_mapping_conc ) :

	mapped_conc = abs_conc * ( (abs_conc - lower_bound_mapping_conc) / (upper_bound_mapping_conc - lower_bound_mapping_conc) )

	return mapped_conc

# function to map concentrations from abs value to map space
#====================================================================================================

#====================================================================================================
# function to calculate statistics of a 2D array region

def function_stats_of_desired_region( monthly_mean_2d_mesh , timeseries_plotting , lon_cross_array , lat_cross_array , ll_lon_of_region , ll_lat_of_region , range_in_row , range_in_col , region_name ) :

	print(" ")
	print( f'-> getting the statistics of a region/mesh == {region_name} ')

	list_from_region = []

	( mesh_row , mesh_col ) = monthly_mean_2d_mesh.shape

	print( f'-> input mesh has total number of rows= {mesh_row} and total number of col= {mesh_col} ')

	### NOITE: totasl diff array is computed from CROSS points (= center of the cells)
	lon_diff_arr = lon_cross_array - ll_lon_of_region
	lat_diff_arr = lat_cross_array - ll_lat_of_region

	#print(f'-> type of {lon_diff_arr} is= {type(lon_diff_arr)} ')

	### create the abs total array of abs min distances
	total_diff_arr = np.abs( lon_diff_arr ) + np.abs( lat_diff_arr )
	### find the row/col with min distance
	tuple_of_row_col_of_cell = np.argwhere( total_diff_arr == np.min(total_diff_arr) )[0]

	print(f'-> row/col of lower-left cell is= {tuple_of_row_col_of_cell} ')
	row_of_cell_with_point = tuple_of_row_col_of_cell[0]
	col_of_cell_with_marker = tuple_of_row_col_of_cell[1]

	# print(f'-> row of the cell contaning the point is= {row_of_cell_with_point}')
	# print(f'-> column of the cell contaning the point is= {col_of_cell_with_marker}')
	print(f'-> now plot the mesh from the origin == starting cell that includes the point == (row,col) = {row_of_cell_with_point , col_of_cell_with_marker} ')
	print(f'-> range in row is= {range_in_row} ')
	print(f'-> range in col is= {range_in_col} ')

	for cell_row in range( row_of_cell_with_point , row_of_cell_with_point + range_in_row , 1 ) :
		for cell_col in range( col_of_cell_with_marker , col_of_cell_with_marker + range_in_col , 1 ) :

			#print( f'-> loop for row= {cell_row} and col= {cell_col}')

			### extract the desired region from lower-left (ll) point + number of cells to add
			cell_val = monthly_mean_2d_mesh[ cell_row , cell_col ]
			list_from_region.append( cell_val )


	print( f'-> size of list_from_region= { len(list_from_region) } ')

	min_of_region= np.min( list_from_region )
	mean_of_region= np.mean( list_from_region )
	median_of_region= np.median( list_from_region )
	std_of_region= np.std( list_from_region )
	max_of_region= np.max( list_from_region )

	# loop again to find the row-col of max cell
	if ( timeseries_plotting == 'yes' ) :

		for row in range(mesh_row) :
			for col in range(mesh_col) :

				cell_val = monthly_mean_2d_mesh[ row , col ]

				if ( max_of_region == cell_val ) :

					row_of_max = row
					col_of_max = col
					break

		return min_of_region , mean_of_region , median_of_region , std_of_region , max_of_region , row_of_max , col_of_max

	if ( timeseries_plotting == 'no' ) :
		return min_of_region , mean_of_region , median_of_region , std_of_region , max_of_region


# function to calculate statistics of a 2D array region
#====================================================================================================

#====================================================================================================
# function to change 3D to 2D array

def function_3Dto2D ( domain_rows , domain_cols , monthly_tseries_tensor , mapping ) :
	" returns monthly mean of each cell, changes 3D mesh of daily mean conc to a 2D mesh of monthly mean conc"

	### define a 2d array
	mesh_2d = np.ndarray( shape= ( domain_rows , domain_cols ) )
	mesh_2d_mapped = np.ndarray( shape= ( domain_rows , domain_cols ) )

	for row in range( 0 , monthly_tseries_tensor.shape[1] , 1 ) :

		for col in range( 0 , monthly_tseries_tensor.shape[2] , 1 ) :

			#print(f'-> processing mesh-3D-monthly mean @ row= {row} and col= {col} ')

			cell_z_axis = monthly_tseries_tensor [ : , row , col ]

			# print( f'-> z axis = {cell_z_axis}' )
			#print(f'-> size of z-axis is= { cell_z_axis.shape }' )

			### take average of each z-axis
			cell_z_axis_mean = cell_z_axis.mean()

			# print( f'-> mean= {cell_z_axis_mean} for {cmaq_pol} @ row= {row} and col= {col}')
			if ( mapping == 'yes' ) :
				#print('-> apply mapping function...')

				conc_mapped = mapping_function( cell_z_axis_mean , lower_bound_mapping_conc , upper_bound_mapping_conc )

				mesh_2d_mapped [ row , col ] = conc_mapped

			### asign the cell mean to 2D array
			mesh_2d [ row , col ] = cell_z_axis_mean

	print(" ")
	print( f'-> shape of monthly 2D array, output of function_3Dto2D = { mesh_2d.shape }' )
	print(" ")
	### function returns a 2D array to be used for plotting
	return mesh_2d , mesh_2d_mapped

# function to change 3D to 2D array
#====================================================================================================

#====================================================================================================
# function to produce raster image

def array2raster( raster_dir , cmaq_pol , file_date_tag , output_array , array_origin_lon , array_origin_lat ) :

	raster_name = 'raster_'+cmaq_pol+'_'+file_date_tag+'.tif'
	path = raster_dir + raster_name

	print( f'-> raster name= {raster_name}')
	print( f'-> raster path= {path}')

# center of domain
	xcent =-120.806 # degrees
	ycent =40.000 # degrees
	upper_lat = 40.450
	lower_lat = 36.450

	no_of_bands = 1
	datatype = gdal.GDT_Float32
	#epsg_code = 4326 # output coord-ref-sys

	raster_origin = ( array_origin_lon , array_origin_lat ) # unit? metere or degree? it can be either meter, or degrees --> ImportFromProj4(+units=m)

	pixelWidth = 1000 # meters
	pixelHeight = 1000 # meters

	Xorig = raster_origin[0]
	Yorig = raster_origin[1]

	rows = output_array.shape[0]
	cols = output_array.shape[1]

	geotransform = ( Xorig , pixelWidth , 0 , Yorig , 0 , pixelHeight ) # units? either meter or degrees

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

#====================================================================================================
# function to

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

#====================================================================================================
# function to

def function_pm25_daily_cell_tseries ( include_pmdiag_file , aconc_open , pmdiag_open , lay , row , col ) : # arg are the variables that are defined insdie this function
	" returns daily timeseries for pm2.5 for each cell"
	print(" ")
	# loop inside 24 time-steps and extract pm concentrations
	# extract PM2.5 species from input files
	print('-> inside pm2.5 daily function: extracting several species from CMAQ files for pm2.5 ...')
	print( f'-> processing row= {row} and col= {col}' )

	# species from aconc [1]
	AH3OPI = aconc_open.variables['AH3OPI'][:,lay,row,col]
	AH3OPI = np.array(AH3OPI)

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

	if ( include_pmdiag_file == 'yes' ) :

		print( '-> PMDIAG file is included in PM2.5 calculations ...' )

		### species from pmdiag [3]
		PM25AT = pmdiag_open.variables['PM25AT'][:,lay,row,col]
		PM25AT = np.array(PM25AT)
		#print( f'-> PM25AT time-series= {PM25AT}')
		#print( f'-> mean of PM25AT= {PM25AT.mean()}')

		PM25AC = pmdiag_open.variables['PM25AC'][:,lay,row,col]
		PM25AC = np.array(PM25AC)
		#print( f'-> PM25AC time-series= {PM25AC}')
		#print( f'-> mean of PM25AC= {PM25AC.mean()}')

		PM25CO = pmdiag_open.variables['PM25CO'][:,lay,row,col]
		PM25CO = np.array(PM25CO)
		#print( f'-> PM25CO time-series= {PM25CO}')
		#print( f'-> mean of PM25CO= {PM25CO.mean()}')

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

	else:

		print( '-> PMDIAG file is _NOT_ included in PM2.5 calculations ...' )
		print( f'-> PM2.5 concentrations will be over-estimated!')

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
		PM25_HP      = (AH3OPI  + AH3OPJ  + AH3OPK ) * (1.0/19.0)
		PM25_CL      = ACLI  + ACLJ  + ACLK
		PM25_EC      = AECI  + AECJ
		PM25_NA      = ANAI  + ANAJ  + ANAK
		PM25_MG      = AMGJ  + AMGK
		PM25_K       = AKJ  + AKK
		PM25_CA      = ACAJ  + ACAK
		PM25_NH4     = ANH4I  + ANH4J  + ANH4K
		PM25_NO3     = ANO3I  + ANO3J  + ANO3K
		PM25_OC      = AOCI  + AOCJ
		PM25_OM      = AOMI  + AOMJ
		PM25_SOIL    = ASOILJ + ASOIL
		PM25_SO4     = ASO4I  + ASO4J  + ASO4K
		PM25_TOT     = ATOTI  + ATOTJ  + ATOTK
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
