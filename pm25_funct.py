#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np

# each function reads its own files inside the function and writes out the data-mesh array

def function_pm25_monthly_mean ( days_in_month , domain_rows , domain_cols ) :

	print('-> month of analysis is=' , cmaq_file_month)

	pm25_mesh = np.ndarray( shape=(days_in_month , domain_rows , domain_cols) )
	# create a list in a range, use an argument-unpacking operator * to unpack the list
	day_list = [*range( 1 , days_in_month+1 , 1)] # don't forget the [] around range function to create the list

	for day_of_the_month in day_list :

		print('-> we are analyzing the follwoing days:')
		print(day_list)
		print( '-> calculating for day= %s' %day_of_the_month )
		# prepare the count days
		if day_of_the_month <= 9 :
			# if jday is less than 10, add zero before it
			day_count = '0'+str(day_of_the_month)

		else:
			# if jday is bigger than 9, use it.
			day_count = str(day_of_the_month)

		file_date_tag = '2016'+cmaq_file_month+day_count
		# setting the input files
		aconc_file_name = 'CCTM_ACONC_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+Landis_scenario+'_mpi_standard_'+file_date_tag+'.nc'

		pmdiag_file_name = 'CCTM_PMDIAG_v52_CA_WRF_1km_griddedAgBioNonptPtfire_scen'+Landis_scenario+'_mpi_standard_'+file_date_tag+'.nc'

		print('-> reading CMAQ files:')
		print( aconc_file_name )
		print( pmdiag_file_name )
		# set input directory
		input_dir = '/storage/ehsanm/USFS_CA_WRF_1km/plots/'
		# define input files
		aconc_input = input_dir + aconc_file_name

		pmdiag_input = input_dir + pmdiag_file_name
		# read in cmaq and pmdiag input files
		aconc_in = Dataset( aconc_input )
		pmdiag_in = Dataset( pmdiag_input )

		# read each cell in the C-storing style: row and then col
		for row in range( 0 , domain_rows , 1 ):

			print('-> reading row=' , row)

			for col in range( 0 , domain_cols , 1 ):

				print('-> reading col=' , col)
				# create a list for each row-coll cell inside the loop
				collected_24hr_pm_conc_list = []
				# loop inside 24 time-steps and extract pm cons
				for tstep in range( 0, 24, 1 ): # check the time steps inside each file: is 24 there?
					# extract PM2.5 species from input files
					print('-> extracting species from CMAQ files...')
					# species from aconc [1]
					AH3OPI = aconc_in.variables['AH3OPI'][tstep][lay][row][col]
					AH3OPJ = aconc_in.variables['AH3OPJ'][tstep][lay][row][col]
					AH3OPK = aconc_in.variables['AH3OPK'][tstep][lay][row][col]
					ACLI = aconc_in.variables['ACLI'][tstep][lay][row][col]
					ACLJ = aconc_in.variables['ACLJ'][tstep][lay][row][col]
					ACLK = aconc_in.variables['ACLK'][tstep][lay][row][col]
					AECI = aconc_in.variables['AECI'][tstep][lay][row][col]
					AECJ = aconc_in.variables['AECJ'][tstep][lay][row][col]
					ANAI = aconc_in.variables['ANAI'][tstep][lay][row][col]
					ANAJ = aconc_in.variables['ANAJ'][tstep][lay][row][col]
					AMGJ = aconc_in.variables['AMGJ'][tstep][lay][row][col]
					AKJ = aconc_in.variables['AKJ'][tstep][lay][row][col]
					ACAJ = aconc_in.variables['ACAJ'][tstep][lay][row][col]
					ANH4I = aconc_in.variables['ANH4I'][tstep][lay][row][col]
					ANH4J = aconc_in.variables['ANH4J'][tstep][lay][row][col]
					ANO3I = aconc_in.variables['ANO3I'][tstep][lay][row][col]
					ANO3J = aconc_in.variables['ANO3J'][tstep][lay][row][col]
					ASOIL = aconc_in.variables['ASOIL'][tstep][lay][row][col]
					ASO4I = aconc_in.variables['ASO4I'][tstep][lay][row][col]
					ASO4J = aconc_in.variables['ASO4J'][tstep][lay][row][col]
					ALVPO1I = aconc_in.variables['ALVPO1I'][tstep][lay][row][col]
					ASVPO1I = aconc_in.variables['ASVPO1I'][tstep][lay][row][col]
					ASVPO2I = aconc_in.variables['ASVPO2I'][tstep][lay][row][col]
					ALVOO1I = aconc_in.variables['ALVOO1I'][tstep][lay][row][col]
					ALVOO2I = aconc_in.variables['ALVOO2I'][tstep][lay][row][col]
					ASVOO1I = aconc_in.variables['ASVOO1I'][tstep][lay][row][col]
					ASVOO2I = aconc_in.variables['ASVOO2I'][tstep][lay][row][col]
					ALVPO1J = aconc_in.variables['ALVPO1J'][tstep][lay][row][col]
					ASVPO1J = aconc_in.variables['ASVPO1J'][tstep][lay][row][col]
					ASVPO2J = aconc_in.variables['ASVPO2J'][tstep][lay][row][col]
					ASVPO3J = aconc_in.variables['ASVPO3J'][tstep][lay][row][col]
					AIVPO1J = aconc_in.variables['AIVPO1J'][tstep][lay][row][col]
					AXYL1J = aconc_in.variables['AXYL1J'][tstep][lay][row][col]
					AXYL2J = aconc_in.variables['AXYL2J'][tstep][lay][row][col]
					AXYL3J = aconc_in.variables['AXYL3J'][tstep][lay][row][col]
					ATOL1J = aconc_in.variables['ATOL1J'][tstep][lay][row][col]
					ATOL2J = aconc_in.variables['ATOL2J'][tstep][lay][row][col]
					ATOL3J = aconc_in.variables['ATOL3J'][tstep][lay][row][col]
					ABNZ1J = aconc_in.variables['ABNZ1J'][tstep][lay][row][col]
					ABNZ2J = aconc_in.variables['ABNZ2J'][tstep][lay][row][col]
					ABNZ3J = aconc_in.variables['ABNZ3J'][tstep][lay][row][col]
					AISO1J = aconc_in.variables['AISO1J'][tstep][lay][row][col]
					AISO2J = aconc_in.variables['AISO2J'][tstep][lay][row][col]
					AISO3J = aconc_in.variables['AISO3J'][tstep][lay][row][col]
					ATRP1J = aconc_in.variables['ATRP1J'][tstep][lay][row][col]
					ATRP2J = aconc_in.variables['ATRP2J'][tstep][lay][row][col]
					ASQTJ = aconc_in.variables['ASQTJ'][tstep][lay][row][col]
					AALK1J = aconc_in.variables['AALK1J'][tstep][lay][row][col]
					AALK2J = aconc_in.variables['AALK2J'][tstep][lay][row][col]
					AORGCJ = aconc_in.variables['AORGCJ'][tstep][lay][row][col]
					AOLGBJ = aconc_in.variables['AOLGBJ'][tstep][lay][row][col]
					AOLGAJ = aconc_in.variables['AOLGAJ'][tstep][lay][row][col]
					APAH1J = aconc_in.variables['APAH1J'][tstep][lay][row][col]
					APAH2J = aconc_in.variables['APAH2J'][tstep][lay][row][col]
					APAH3J = aconc_in.variables['APAH3J'][tstep][lay][row][col]
					ALVOO1J = aconc_in.variables['ALVOO1J'][tstep][lay][row][col]
					ALVOO2J = aconc_in.variables['ALVOO2J'][tstep][lay][row][col]
					ASVOO1J = aconc_in.variables['ASVOO1J'][tstep][lay][row][col]
					ASVOO2J = aconc_in.variables['ASVOO2J'][tstep][lay][row][col]
					ASVOO3J = aconc_in.variables['ASVOO3J'][tstep][lay][row][col]
					APCSOJ = aconc_in.variables['APCSOJ'][tstep][lay][row][col]
					AALJ = aconc_in.variables['AALJ'][tstep][lay][row][col]
					ASIJ = aconc_in.variables['ASIJ'][tstep][lay][row][col]
					AFEJ = aconc_in.variables['AFEJ'][tstep][lay][row][col]
					ATIJ = aconc_in.variables['ATIJ'][tstep][lay][row][col]
					AOTHRI = aconc_in.variables['AOTHRI'][tstep][lay][row][col]
					AOTHRJ = aconc_in.variables['AOTHRJ'][tstep][lay][row][col]
					ACORS = aconc_in.variables['ACORS'][tstep][lay][row][col]
					ASEACAT = aconc_in.variables['ASEACAT'][tstep][lay][row][col]
					ASO4K = aconc_in.variables['ASO4K'][tstep][lay][row][col]
					ANO3K = aconc_in.variables['ANO3K'][tstep][lay][row][col]
					ANH4K = aconc_in.variables['ANH4K'][tstep][lay][row][col]
					AMNJ = aconc_in.variables['AMNJ'][tstep][lay][row][col]

					# species from pmdiag [3]
					PM25AT = pmdiag_in.variables['PM25AT'][tstep][lay][row][col]
					PM25AC = pmdiag_in.variables['PM25AC'][tstep][lay][row][col]
					PM25CO = pmdiag_in.variables['PM25CO'][tstep][lay][row][col]

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
					PM25_HP     = (AH3OPI * PM25AT + AH3OPJ * PM25AC + AH3OPK * PM25CO) * 1.0/19.0
					PM25_CL     = ACLI * PM25AT + ACLJ * PM25AC + ACLK * PM25CO
					PM25_EC     =    AECI * PM25AT + AECJ * PM25AC
					PM25_NA       =    ANAI * PM25AT + ANAJ * PM25AC + ANAK * PM25CO
					PM25_MG      = AMGJ * PM25AC + AMGK * PM25CO
					PM25_K       =     AKJ * PM25AC + AKK * PM25CO
					PM25_CA       =         ACAJ * PM25AC + ACAK * PM25CO
					PM25_NH4      =    ANH4I * PM25AT + ANH4J * PM25AC + ANH4K * PM25CO
					PM25_NO3    =    ANO3I * PM25AT + ANO3J * PM25AC + ANO3K * PM25CO
					PM25_OC      =   AOCI * PM25AT + AOCJ * PM25AC
					PM25_OM       =    AOMI * PM25AT + AOMJ * PM25AC
					PM25_SOIL       =    ASOILJ * PM25AC + ASOIL * PM25CO
					PM25_SO4       =    ASO4I * PM25AT + ASO4J * PM25AC + ASO4K * PM25CO
					PM25_TOT       =     ATOTI * PM25AT + ATOTJ * PM25AC + ATOTK * PM25CO
					PM25_UNSPEC1   =    PM25_TOT - (PM25_CL + PM25_EC + PM25_NA + PM25_NH4 + PM25_NO3 + PM25_OC + PM25_SOIL + PM25_SO4 )

					# now sum all species to get hourly PM2.5 concentratiosn
					hrly_pm25_conc = PM25_HP + PM25_CL + PM25_EC + PM25_NA + PM25_MG + PM25_K + PM25_CA + PM25_NH4 + PM25_NO3 + PM25_OC + PM25_OM + PM25_SOIL + PM25_SO4 + PM25_TOT + PM25_UNSPEC1
					# append sum of species to houtly list
					collected_24hr_pm_conc_list.append( hrly_pm25_conc )

				# make an array from list when 24 time-step is finished and we get out of the t-step loop
				total_24hr_pm_conc_array = np.array( collected_24hr_pm_conc_list )
				# get the mean of all 24-hour-summed species
				cell_daily_mean_pm25 = total_24hr_pm_conc_array.mean()
				# fill the data-mesh with data, based on the order: z, x, y
				print('-> pin the data at day= %s , row= %s , col= %s' %( day_of_the_month-1 , row , col ) )
				# add daily mean pm to each layer, row, col of pm25 mesh array
				pm25_mesh [ day_of_the_month-1 ][ row ][ col ] = cell_daily_mean_pm25
				# remove the 24-hr list for each row-col cell, it will be created for the next cell
				del collected_24hr_pm_conc_list
		# close nc file
		aconc_in.close()
		pmdiag_in.close()
	# function returns the data-mesh to use in plotting
	return pm25_mesh


cmaq_file_month = '09'
days_in_month = 2
Landis_scenario = '1'
#cmaq_pol = 'CO'
lay = 0
domain_cols = 250
domain_rows = 265

# use the function
data_mesh = function_pm25_monthly_mean( days_in_month , domain_rows , domain_cols )
print('-> number of dimensions=' , data_mesh.ndim )
print('-> shape of data-mesh=' , data_mesh.shape )







