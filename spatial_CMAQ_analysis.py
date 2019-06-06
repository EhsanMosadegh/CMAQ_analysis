#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################################
# import libraries

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap , cm
#from osgeo import gdal, gdal_array, osr , ogr
from osgeo import gdal
import time

###################################################################################

###################################################################################
# define functions that calculate concentrations of pollutants

def function_3D_mesh_maker ( days_to_run_in_month , domain_rows , domain_cols , cmaq_file_month , scenario , input_path_scen , input_path_base ) :
	" processes the files and days and returns two 3D meshes"

	print('-> month of analysis is=' , cmaq_file_month)

	### define the monthly array mesh
	if ( processing_pol == 'co' ) :

		#monthly_3D_tensor_scen = np.ndarray( shape=( 24 , domain_rows , domain_cols ) )
		monthly_3D_tensor_scen = np.empty( shape=( 0 , domain_rows , domain_cols ) )
		monthly_3D_tensor_base = np.empty( shape=( 0 , domain_rows , domain_cols ) )

	elif ( processing_pol == 'pm2.5' ) :

		monthly_3D_tensor_scen = np.empty( shape=( days_to_run_in_month , domain_rows , domain_cols ) )
		monthly_3D_tensor_base = np.empty( shape=( days_to_run_in_month , domain_rows , domain_cols ) )

	else:
		pass

	# monthly_3D_tensor_scen = np.zeros( shape=( 31 , domain_rows , domain_cols ) )
	# monthly_3D_tensor_base = np.zeros( shape=( 31 , domain_rows , domain_cols ) )
	# print(f'-> shape of zero array= {monthly_3D_tensor_scen.shape}')

	## create a day list for a month to create file-date-tag, use an argument-unpacking operator * to unpack the list
	day_list = [*range( 1 , days_to_run_in_month+1 , 1)] # don't forget the [] around range function to create the list

	# # to run for specific day
	# day_list = [21]  # use the favorite day
	############################################################################################
	### traverse the list for each day

	for day_of_the_month in day_list :

		print( " ")
		print('-> we are processing for following days:')
		print(day_list)
		print( f'-> processing for day= {day_of_the_month}' )

		# prepare the day flags
		if day_of_the_month <= 9 :
			# if jday is less than 10, add zero before it
			day_count = '0'+str(day_of_the_month)

		else:
			# if jday is bigger than 9, use it.
			day_count = str(day_of_the_month)

		### define file tags
		file_date_tag = cmaq_file_year + cmaq_file_month + day_count

		############################################################################################
		### define input files

		aconc_scen = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+scenario+'_mpi_standard_'+file_date_tag+'.nc'
		pmdiag_scen = 'CCTM_PMDIAG_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+scenario+'_mpi_standard_'+file_date_tag+'.nc'

		aconc_base = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonpt_baseline_AgBioNonpt_mpi_standard_'+file_date_tag+'.nc'
		pmdiag_base = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonpt_baseline_AgBioNonpt_mpi_standard_'+file_date_tag+'.nc'

		############################################################################################
		### open netcdf files based on each processing method and pollutant

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

				# open netcdf file
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

				# open netcdf files
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

		############################################################################################
		### process the nc opened files for each cell with a specific function

		if ( processing_method == 'single_plot' ) :

			### create an empty tensor for each cell and day as container of daily 24-hr t-step concentrations 

			#daily_3d_mesh_scen = np.ndarray ( shape=( 24 , domain_rows , domain_cols ) )
			daily_3d_mesh_scen = np.empty ( shape=( 24 , domain_rows , domain_cols ) )

			#print(f'-> shape of daily tseries array={daily_3d_mesh_scen.shape }')
			#print('-> traversing each cell and extract pollutants ...')
			
			### traverse each cell in the C-storing style for each day: row and then col
			for row in range( 0 , domain_rows , 1 ) :

				for col in range( 0 , domain_cols , 1 ) :

					if ( processing_pol == 'co' ) :

						print( f'-> extracting cell for single POL - single - at row= {row} and col={col} ... ' )

						cell_24hr_timeSeries_array = function_cell_24hr_timeSeries_singlePOL( aconc_open_scen , cmaq_pol , lay , row , col )
						#print(f'--> cell tseries is= {cell_24hr_timeSeries_array}')

					elif ( processing_pol == 'pm2.5') :

						print( f'-> extracting cell for pm2.5 at row= {row} and col={col} ... ' )

						daily_cell_mean = function_daily_cell_mean_pm25( aconc_open_scen , pmdiag_open_scen , lay , row , col )

					else:

						print( '-> WARNING: define/check single POL or pm2.5 settings and processing method first! ')
						print('-> exiting ...')
						raise SystemExit()

					### fill daily data-mesh with each cell data, based on the order: z, x, y == layer, row, col in mesh_3d_monthly array

					#print( f'-> add/pin each cell mean value to monthly_3D_tensor_scen at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )
					
					daily_3d_mesh_scen [:,row,col]  = cell_24hr_timeSeries_array

					#print(f'--> shape of daily 24hr tseries is= {daily_3d_mesh_scen.shape}')

			# ### after each day is extracted, add the daily frame to monthly frame
			# if ( day_of_the_month == '01' ) :
			# 	### replace the array for 1st day
			# 	monthly_3D_tensor_scen = daily_3d_mesh_scen # because the montly mesh is not empy!

			# else:
			# 	# daily mesh is filled with 1st time step and now we concatenate the new mesh to previous one
			monthly_3D_tensor_scen = np.concatenate( ( monthly_3D_tensor_scen , daily_3d_mesh_scen ) , axis=0 )


		elif ( processing_method == 'diff_plot' ) :

			### create an empty tensor for each cell and day as container of daily 24-hr t-step concentrations 
			daily_3d_mesh_scen = np.empty ( shape=( 24 , domain_rows , domain_cols ) )
			daily_3d_mesh_base = np.empty ( shape=( 24 , domain_rows , domain_cols ) )


			print('-> traversing each cell and extract pollutants ...')
			### traverse each cell in the C-storing style for each day: row and then col
			for row in range( 0 , domain_rows , 1 ):

				for col in range( 0 , domain_cols , 1 ):

					if ( processing_pol == 'co' ) :

						print( f'-> extracting cell for single POL - diff - at row= {row} and col={col} ... ' )
						
						# we calculate cell means for each scenario
						cell_24hr_timeSeries_array_scen = function_cell_24hr_timeSeries_singlePOL( aconc_open_scen , cmaq_pol , lay , row , col )
						cell_24hr_timeSeries_array_base = function_cell_24hr_timeSeries_singlePOL( aconc_open_base , cmaq_pol , lay , row , col )

						### fill daily tensors cells with 24-hr time-series
						daily_3d_mesh_scen [:,row,col] = cell_24hr_timeSeries_array_scen
						daily_3d_mesh_base [:,row,col] = cell_24hr_timeSeries_array_base

			### now we fill 2 meshes
			print( f'-> add/pin each cell 24-hr t-series to monthly_3D_tensor_scen at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )

			monthly_3D_tensor_scen = np.concatenate( ( monthly_3D_tensor_scen ,  daily_3d_mesh_scen ) , axis=0 )
			#monthly_3D_tensor_scen [ day_of_the_month-1 ][ row ][ col ] = cell_24hr_timeSeries_array_scen

			print( f'-> add/pin each cell 24hr t-series to monthly_3D_tensor_base at sheet(=day-1)= {day_of_the_month-1} , row= {row} , col= {col}' )
			
			monthly_3D_tensor_base = np.concatenate( ( monthly_3D_tensor_base ,  daily_3d_mesh_base ) , axis=0 )

						#monthly_3D_tensor_base [ day_of_the_month-1 ][ row ][ col ] = cell_24hr_timeSeries_array_base

					# elif ( processing_pol == 'pm25' ) :

					# 	# we calculate cell means
					# 	daily_cell_mean_scen = function_daily_cell_mean_pm25( aconc_open_scen , pmdiag_open_scen , lay , row , col )
					# 	daily_cell_mean_base = function_daily_cell_mean_pm25( aconc_open_base , pmdiag_open_base , lay , row , col )

					# 	# now we fill the 3D mesh with cell means
					# 	monthly_3D_tensor_scen [ day_of_the_month-1 ][ row ][ col ] = daily_cell_mean_scen
					# 	monthly_3D_tensor_base [ day_of_the_month-1 ][ row ][ col ] = daily_cell_mean_base




					# else:
					# 	pass
		else:
			print('-> ERROR: processing method NOT defined!')
		
		############################################################################################
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
	
	############################################################################################
	### return two 3D meshs
	if ( processing_method == 'single_plot' ) :

		return monthly_3D_tensor_scen

	else:
	 	return monthly_3D_tensor_scen , monthly_3D_tensor_base

	# # this function retuns data-mesh-2D for plotting
	# return monthly_mean_mesh_2d

############################################################################################
### function to change 3D to 2D array

def function_3Dto2D ( domain_rows , domain_cols , mesh_3d_monthly  ) :
	" returns monthly mean of each cell, changes 3D mesh of daily mean conc to a 2D mesh of monthly mean conc"

	### define a 2d array
	cell_monthly_mean_2d_mesh = np.ndarray( shape= ( domain_rows , domain_cols ) )

	for row in range( 0 , mesh_3d_monthly.shape[1] , 1 ) :

		for col in range( 0 , mesh_3d_monthly.shape[2] , 1 ) :

			#print(f'-> processing mesh-3D-monthly mean @ row= {row} and col= {col} ')

			cell_z_axis = mesh_3d_monthly [ : , row , col ]

			#print(f'-> size of z-axis is= { cell_z_axis.shape } ')

			### take average of each z-axis
			cell_z_axis_mean = cell_z_axis.mean()
			### asign the cell mean to 2D array
			cell_monthly_mean_2d_mesh [ row ][ col ] = cell_z_axis_mean
	print(" ")
	print( f'-> shape of 2D array output of function: 3Dto2D= { cell_monthly_mean_2d_mesh.shape }' )
	print(" ")
	### function returns a 2D array to be used for plotting
	return cell_monthly_mean_2d_mesh

############################################################################################
### function

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

############################################################################################
### function

def function_daily_cell_mean_pm25 ( aconc_open , pmdiag_open , lay , row , col ) : # arg are the variables that are defined insdie this function
	" returns daily mean of pm2.5 for each cell"

	print( f'-> processing row= {row} and col= {col}' )

	# loop inside 24 time-steps and extract pm cons
	# extract PM2.5 species from input files
	print('-> extracting several species from CMAQ files for pm2.5 ...')
	# species from aconc [1]
	AH3OPI = aconc_open.variables['AH3OPI'][:,lay,row,col]
	AH3OPI = np.array(AH3OPI).mean()

	AH3OPJ = aconc_open.variables['AH3OPJ'][:,lay,row,col]
	AH3OPJ = np.array(AH3OPJ).mean()

	AH3OPK = aconc_open.variables['AH3OPK'][:,lay,row,col]
	AH3OPK = np.array(AH3OPK).mean()

	ACLI = aconc_open.variables['ACLI'][:,lay,row,col]
	ACLI = np.array(ACLI).mean()

	ACLJ = aconc_open.variables['ACLJ'][:,lay,row,col]
	ACLJ = np.array(ACLJ).mean()

	ACLK = aconc_open.variables['ACLK'][:,lay,row,col]
	ACLK = np.array(ACLK).mean()

	AECI = aconc_open.variables['AECI'][:,lay,row,col]
	AECI = np.array(AECI).mean()

	AECJ = aconc_open.variables['AECJ'][:,lay,row,col]
	AECJ = np.array(AECJ).mean()

	ANAI = aconc_open.variables['ANAI'][:,lay,row,col]
	ANAI = np.array(ANAI).mean()

	ANAJ = aconc_open.variables['ANAJ'][:,lay,row,col]
	ANAJ = np.array(ANAJ).mean()

	AMGJ = aconc_open.variables['AMGJ'][:,lay,row,col]
	AMGJ = np.array(AMGJ).mean()

	AKJ = aconc_open.variables['AKJ'][:,lay,row,col]
	AKJ = np.array(AKJ).mean()

	ACAJ = aconc_open.variables['ACAJ'][:,lay,row,col]
	ACAJ = np.array(ACAJ).mean()

	ANH4I = aconc_open.variables['ANH4I'][:,lay,row,col]
	ANH4I = np.array(ANH4I).mean()

	ANH4J = aconc_open.variables['ANH4J'][:,lay,row,col]
	ANH4J = np.array(ANH4J).mean()

	ANO3I = aconc_open.variables['ANO3I'][:,lay,row,col]
	ANO3I = np.array(ANO3I).mean()

	ANO3J = aconc_open.variables['ANO3J'][:,lay,row,col]
	ANO3J = np.array(ANO3J).mean()

	ASOIL = aconc_open.variables['ASOIL'][:,lay,row,col]
	ASOIL = np.array(ASOIL).mean()

	ASO4I = aconc_open.variables['ASO4I'][:,lay,row,col]
	ASO4I = np.array(ASO4I).mean()

	ASO4J = aconc_open.variables['ASO4J'][:,lay,row,col]
	ASO4J = np.array(ASO4J).mean()

	ALVPO1I = aconc_open.variables['ALVPO1I'][:,lay,row,col]
	ALVPO1I = np.array(ALVPO1I).mean()

	ASVPO1I = aconc_open.variables['ASVPO1I'][:,lay,row,col]
	ASVPO1I = np.array(ASVPO1I).mean()

	ASVPO2I = aconc_open.variables['ASVPO2I'][:,lay,row,col]
	ASVPO2I = np.array(ASVPO2I).mean()

	ALVOO1I = aconc_open.variables['ALVOO1I'][:,lay,row,col]
	ALVOO1I = np.array(ALVOO1I).mean()

	ALVOO2I = aconc_open.variables['ALVOO2I'][:,lay,row,col]
	ALVOO2I = np.array(ALVOO2I).mean()

	ASVOO1I = aconc_open.variables['ASVOO1I'][:,lay,row,col]
	ASVOO1I = np.array(ASVOO1I).mean()

	ASVOO2I = aconc_open.variables['ASVOO2I'][:,lay,row,col]
	ASVOO2I = np.array(ASVOO2I).mean()

	ALVPO1J = aconc_open.variables['ALVPO1J'][:,lay,row,col]
	ALVPO1J = np.array(ALVPO1J).mean()

	ASVPO1J = aconc_open.variables['ASVPO1J'][:,lay,row,col]
	ASVPO1J = np.array(ASVPO1J).mean()

	ASVPO2J = aconc_open.variables['ASVPO2J'][:,lay,row,col]
	ASVPO2J = np.array(ASVPO2J).mean()

	ASVPO3J = aconc_open.variables['ASVPO3J'][:,lay,row,col]
	ASVPO3J = np.array(ASVPO3J).mean()

	AIVPO1J = aconc_open.variables['AIVPO1J'][:,lay,row,col]
	AIVPO1J = np.array(AIVPO1J).mean()

	AXYL1J = aconc_open.variables['AXYL1J'][:,lay,row,col]
	AXYL1J = np.array(AXYL1J).mean()

	AXYL2J = aconc_open.variables['AXYL2J'][:,lay,row,col]
	AXYL2J = np.array(AXYL2J).mean()

	AXYL3J = aconc_open.variables['AXYL3J'][:,lay,row,col]
	AXYL3J = np.array(AXYL3J).mean()

	ATOL1J = aconc_open.variables['ATOL1J'][:,lay,row,col]
	ATOL1J = np.array(ATOL1J).mean()

	ATOL2J = aconc_open.variables['ATOL2J'][:,lay,row,col]
	ATOL2J = np.array(ATOL2J).mean()

	ATOL3J = aconc_open.variables['ATOL3J'][:,lay,row,col]
	ATOL3J = np.array(ATOL3J).mean()

	ABNZ1J = aconc_open.variables['ABNZ1J'][:,lay,row,col]
	ABNZ1J = np.array(ABNZ1J).mean()

	ABNZ2J = aconc_open.variables['ABNZ2J'][:,lay,row,col]
	ABNZ2J = np.array(ABNZ2J).mean()

	ABNZ3J = aconc_open.variables['ABNZ3J'][:,lay,row,col]
	ABNZ3J = np.array(ABNZ3J).mean()

	AISO1J = aconc_open.variables['AISO1J'][:,lay,row,col]
	AISO1J = np.array(AISO1J).mean()

	AISO2J = aconc_open.variables['AISO2J'][:,lay,row,col]
	AISO2J = np.array(AISO2J).mean()

	AISO3J = aconc_open.variables['AISO3J'][:,lay,row,col]
	AISO3J = np.array(AISO3J).mean()

	ATRP1J = aconc_open.variables['ATRP1J'][:,lay,row,col]
	ATRP1J = np.array(ATRP1J).mean()

	ATRP2J = aconc_open.variables['ATRP2J'][:,lay,row,col]
	ATRP2J = np.array(ATRP2J).mean()

	ASQTJ = aconc_open.variables['ASQTJ'][:,lay,row,col]
	ASQTJ = np.array(ASQTJ).mean()

	AALK1J = aconc_open.variables['AALK1J'][:,lay,row,col]
	AALK1J = np.array(AALK1J).mean()

	AALK2J = aconc_open.variables['AALK2J'][:,lay,row,col]
	AALK2J = np.array(AALK2J).mean()

	AORGCJ = aconc_open.variables['AORGCJ'][:,lay,row,col]
	AORGCJ = np.array(AORGCJ).mean()

	AOLGBJ = aconc_open.variables['AOLGBJ'][:,lay,row,col]
	AOLGBJ = np.array(AOLGBJ).mean()

	AOLGAJ = aconc_open.variables['AOLGAJ'][:,lay,row,col]
	AOLGAJ = np.array(AOLGAJ).mean()

	APAH1J = aconc_open.variables['APAH1J'][:,lay,row,col]
	APAH1J = np.array(APAH1J).mean()

	APAH2J = aconc_open.variables['APAH2J'][:,lay,row,col]
	APAH2J = np.array(APAH2J).mean()

	APAH3J = aconc_open.variables['APAH3J'][:,lay,row,col]
	APAH3J = np.array(APAH3J).mean()

	ALVOO1J = aconc_open.variables['ALVOO1J'][:,lay,row,col]
	ALVOO1J = np.array(ALVOO1J).mean()

	ALVOO2J = aconc_open.variables['ALVOO2J'][:,lay,row,col]
	ALVOO2J = np.array(ALVOO2J).mean()

	ASVOO1J = aconc_open.variables['ASVOO1J'][:,lay,row,col]
	ASVOO1J = np.array(ASVOO1J).mean()

	ASVOO2J = aconc_open.variables['ASVOO2J'][:,lay,row,col]
	ASVOO2J = np.array(ASVOO2J).mean()

	ASVOO3J = aconc_open.variables['ASVOO3J'][:,lay,row,col]
	ASVOO3J = np.array(ASVOO3J).mean()

	APCSOJ = aconc_open.variables['APCSOJ'][:,lay,row,col]
	APCSOJ = np.array(APCSOJ).mean()

	AALJ = aconc_open.variables['AALJ'][:,lay,row,col]
	AALJ = np.array(AALJ).mean()

	ASIJ = aconc_open.variables['ASIJ'][:,lay,row,col]
	ASIJ = np.array(ASIJ).mean()

	AFEJ = aconc_open.variables['AFEJ'][:,lay,row,col]
	AFEJ = np.array(AFEJ).mean()

	ATIJ = aconc_open.variables['ATIJ'][:,lay,row,col]
	ATIJ = np.array(ATIJ).mean()

	AOTHRI = aconc_open.variables['AOTHRI'][:,lay,row,col]
	AOTHRI = np.array(AOTHRI).mean()

	AOTHRJ = aconc_open.variables['AOTHRJ'][:,lay,row,col]
	AOTHRJ = np.array(AOTHRJ).mean()

	ACORS = aconc_open.variables['ACORS'][:,lay,row,col]
	ACORS = np.array(ACORS).mean()

	ASEACAT = aconc_open.variables['ASEACAT'][:,lay,row,col]
	ASEACAT = np.array(ASEACAT).mean()

	ASO4K = aconc_open.variables['ASO4K'][:,lay,row,col]
	ASO4K = np.array(ASO4K).mean()

	ANO3K = aconc_open.variables['ANO3K'][:,lay,row,col]
	ANO3K = np.array(ANO3K).mean()

	ANH4K = aconc_open.variables['ANH4K'][:,lay,row,col]
	ANH4K = np.array(ANH4K).mean()

	AMNJ = aconc_open.variables['AMNJ'][:,lay,row,col]
	AMNJ = np.array(AMNJ).mean()

	# species from pmdiag [3]
	PM25AT = pmdiag_open.variables['PM25AT'][:,lay,row,col]
	PM25AT = np.array(PM25AT).mean()

	PM25AC = pmdiag_open.variables['PM25AC'][:,lay,row,col]
	PM25AC = np.array(PM25AC).mean()

	PM25CO = pmdiag_open.variables['PM25CO'][:,lay,row,col]
	PM25CO = np.array(PM25CO).mean()

	 # species calculated inside SpecDef file [0]
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
	PM25_HP      = (AH3OPI * PM25AT + AH3OPJ * PM25AC + AH3OPK * PM25CO) * 1.0/19.0
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
	PM25_UNSPEC1 = PM25_TOT - (PM25_CL + PM25_EC + PM25_NA + PM25_NH4 + PM25_NO3 + PM25_OC + PM25_SOIL + PM25_SO4 )

	# now sum all species to get hourly PM2.5 concentratiosn
	cell_pm25_daily_mean = PM25_HP + PM25_CL + PM25_EC + PM25_NA + PM25_MG + PM25_K + PM25_CA + \
					PM25_NH4 + PM25_NO3 + PM25_OC + PM25_OM + PM25_SOIL + PM25_SO4 + PM25_TOT + PM25_UNSPEC1

	# function returns the mean of pm2.5 for each cell
	return cell_pm25_daily_mean

###################################################################################
### run time setting
###################################################################################

### get the starting time
start = time.time()

### run time settings
cmaq_file_year = '2016'
cmaq_file_month = '10'
sim_month = 'oct'
days_to_run_in_month = 2
scenario = '4' # 1-5, baseline
mcip_date_tag = '161001'

cmaq_pol = 'CO'  # for plot title
pol_unit = '[ppmV]'
max_conc_threshold = 0.3  # for Basemap plot

processing_pol = 'co' 		# 'co' or 'pm2.5'
processing_method = 'single_plot' 	# 'single_plot' or 'diff_plot'
platform = 'Mac'  # 'Mac' or 'cluster'

print(f'-> CMAQ year= {cmaq_file_year}')
print(f'-> CMAQ month of analysis= {cmaq_file_month}')
print(f'-> LANDIS scenarios= {scenario}')
print(f'-> processing pollutant= {processing_pol}')
print(f'-> processing method= {processing_method}')
print(f'-> number of days to run= {days_to_run_in_month}')
print(f'-> platform is= {platform}')
print(" ")


### domain settings
lay = 0
domain_cols = 250
domain_rows = 265


# ### Basemap plot setting
# # center of domain
# xcent =-120.806 # degrees
# ycent =40.000 # degrees
# # domain size
# NROWS = 265*1000 # meters
# NCOLS = 250*1000 # meters
# # lower-left corner
# llcornerx=-117500 # meters
# llcornery=-265500 # meters
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


### input directory setting
if ( platform == 'Mac' ) :

	input_dir = '/Users/ehsan/Documents/Python_projects/CMAQ_analysis/cmaq_inputs/'
	mcip_dir = '/Users/ehsan/Documents/Python_projects/CMAQ_analysis/cmaq_inputs/'
	fig_dir = '/Users/ehsan/Documents/Python_projects/CMAQ_analysis/cmaq_figs/'

elif ( platform == 'cluster' ) :

	input_dir = '/storage/ehsanm/USFS_CA_WRF_1km/plots/cmaq_usfs_data/'
	mcip_dir = '/storage/ehsanm/USFS_CA_WRF_1km/plots/'
	fig_dir = '/storage/ehsanm/USFS_CA_WRF_1km/plots/CMAQ_analysis/cmaq_figs/'

else:

	print( '-> ERROR: specify running platform ' )
	print('-> exiting ...')
	raise SystemExit()

### set input pathes
input_path_scen = input_dir + 'scen_' + scenario + '/' + sim_month + '/'
input_path_base = input_dir + 'scen_baseline' + '/' + sim_month + '/'

print('-> CMAQ input directory is:')
print(input_path_scen)
print(input_path_base)

print('-> MCIP input directory is:')
print(mcip_dir)
print(" ")

###################################################################################
### main
###################################################################################
# processing CMAQ files

### extract necessary data from CMAQ for each mesh and calculate data_mesh_3d
print('-> calculating monthly tensor ...')

if ( processing_method == 'single_plot' ) :
#monthly_mean_mesh_2d = function_3D_mesh_maker( days_to_run_in_month , domain_rows , domain_cols , cmaq_file_month , scenario , input_path_scen , input_path_base )
	monthly_3D_tensor_scen = function_3D_mesh_maker( days_to_run_in_month , domain_rows , domain_cols , cmaq_file_month , scenario , input_path_scen , input_path_base )

else:
	monthly_3D_tensor_scen , monthly_3D_tensor_base = function_3D_mesh_maker( days_to_run_in_month , domain_rows , domain_cols , cmaq_file_month , scenario , input_path_scen , input_path_base )

############################################################################################
### change 3D to 2D array

if ( processing_method == 'single_plot' ) :

	### look at the 3D tensors for scenario
	print('-----------------------------')
	print('-> 3D data mesh info:')
	print( f'-> LANDIS scenario 3D monthly tensor: number of dimensions= { monthly_3D_tensor_scen.ndim}' )
	print( f'-> LANDIS scenario 3D monthly tensor: shape of data-mesh= { monthly_3D_tensor_scen.shape}' )
	print('-----------------------------')

	### we pnly have one 3D tensor

	print('-> change mesh 3D to 2D for single plot ...')

	monthly_mean_mesh_2d = function_3Dto2D( domain_rows , domain_cols , monthly_3D_tensor_scen )

elif ( processing_method == 'diff_plot' ) :

	### look at the 3D meshes
	print('-----------------------------')
	print('-> check 3D data mesh info:')
	print( f'-> LANDIS scenario 3D monthly mean mesh: number of dimensions= {monthly_3D_tensor_scen.ndim}' )
	print( f'-> LANDIS scenario 3D monthly mean mesh: shape of data-mesh= {monthly_3D_tensor_scen.shape}' )
	print( f'-> baseline 3D monthly mesh: number of dimensions= {monthly_3D_tensor_base.ndim}' )
	print( f'-> baseline 3D monthly mesh: shape of data-mesh= {monthly_3D_tensor_base.shape}' )
	print('-----------------------------')

	### we have 2 3D meshes
	print('-> change mesh 3D-to-2D for diff-plot for mesh-3D-LANDIS ...')
	monthly_mean_2d_mesh_scen = function_3Dto2D( domain_rows , domain_cols , monthly_3D_tensor_scen )

	print('-> change mesh 3D-to-2D for diff-plot for mesh-3D-baseline ...')
	monthly_mean_2d_mesh_base = function_3Dto2D( domain_rows , domain_cols , monthly_3D_tensor_base )

	# now subtract 2 meshes to get the diff mesh
	monthly_mean_mesh_2d = monthly_mean_2d_mesh_scen - monthly_mean_2d_mesh_base

else:
	print('-> ERROR: ')
	print('-> exiting ...')
	raise SystemExit()

###################################################################################
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

###################################################################################
### plotting
###################################################################################
# use Basemap library and make spatial plots
print('-> plotting the data...')

### plot dots from grid coordinates of the dots
#plt.plot( lon_mesh , lat_mesh , marker='.' , color='b' , linestyle= 'none' )

# ### create a Basemap class/model instance for a specific projection
# # basemap_instance = Basemap(projection='lcc' , lat_0=ycent , lon_0=xcent , height=NROWS , width=NCOLS , resolution='i') # , area_thresh=0.1) # latlon=True for when x and y are not in map proj. coordinates
# theMap = Basemap(projection='lcc' ,
# 													 llcrnrx=llcornerx , llcrnry=llcornery , urcrnrx=urcornerx , urcrnry=urcornery ,
# 													 lat_0=ycent , lon_0=xcent , height=NROWS , width=NCOLS ,
# 													 resolution='f' , area_thresh=0.5)


### create Basemap model instance from its class, it is a map that color mesh sits on it.
theMap_zoomed = Basemap(projection='lcc' , lat_0=ycent_zoom , lon_0=xcent_zoom , height=NROWS_zoom , width=NCOLS_zoom , resolution='f' , area_thresh=0.5)

theMap_zoomed.bluemarble()
x_mesh, y_mesh = theMap_zoomed(lon_mesh , lat_mesh) # order: x , y; Basemap model transforms lon/lat from degree to meter for LCC projection map
theMap_zoomed.drawmapboundary(color='k' ) #, fill_color='aqua')
theMap_zoomed.drawcoastlines(color = '0.15')
theMap_zoomed.drawcounties(linewidth=0.5 , color='k')
theMap_zoomed.drawstates()
#basemap_instance.fillcontinents(lake_color='aqua')

#my_levels = [ 0.02 , 0.05 ]
#my_colors = ( 'g' , 'b' , 'r' )
### create a color mesh image from basemap model instance, the color mesh is constant, cos it is plotted from lon/lat values
colorMesh = theMap_zoomed.pcolormesh( x_mesh , y_mesh , monthly_mean_mesh_2d , cmap=plt.cm.OrRd , shading='flat' , vmin=0.0 , vmax=max_conc_threshold ) #levels=my_levels , colors=my_colors
#im2 = basemap_instance.pcolormesh(lon_mesh , lat_mesh , data_mesh , cmap=plt.cm.jet , shading='flat')

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

###################################################################################

###################################################################################
# save the plots

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
print( f'-> run time= { (( end - start ) / 60 ) :.2f} min' )  # f-string

###################################################################################
