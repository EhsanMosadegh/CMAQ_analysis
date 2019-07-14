#!/bin/env python3
###############################################################
# NOTE: in this script change the following params for each run:
# pollutant
# stats_property
###############################################################
import numpy as np
import os
#--------------------------------------------------------

pollutant='CO'
stats_property='max'
region='LakeTahoeBasin' # 'LakeTahoeBasin' or 'NorthTahoe' or 'SouthTahoe'

#--------------------------------------------------------

print('------------------------------------------------')
print('-> now we are inside python...')
print( f'-> we are at= {os.getcwd() }')
print('------------------------------------------------')
print( f'-> chk: pollutant is= {pollutant}')
print( f'-> chk: statisticsl property is= {stats_property}')
print( f'-> chk: region is= {region} ')
print('------------------------------------------------')
#--------------------------------------------------------

file_name=stats_property+'DiffMesh'+region+'_list_total_for_'+pollutant+'_in_'+region+'.txt'

#home_dir='/Users/ehsan/Documents/Python_projects/CMAQ_analysis'
#home_dir='/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/CMAQ_analysis'

#file_path=home_dir+'/logs/'

#log_file=file_path+file_name

print(f'-> looking inside log file= {file_name}')
print( " ")
with open( file_name , 'r') as input_file :

	input_list=[]

	for line in input_file:
		print( f'-> line to append is= {line}' )
		input_list.append( line )

	# change the data type from string to float!
	input_list=np.array(input_list).astype(np.float)

	if (stats_property=='min') :
		print('----------------------------------------')
		print( f'-> chk: abs min for ({pollutant}) in region ({region}) is= {np.min(input_list)}' )

	elif (stats_property=='max') :
		print('----------------------------------------')
		print( f'-> chk: abs max for ({pollutant}) in region ({region}) is= {np.max(input_list)}' )

	else:
		print('----------------------------------------')
		print('-> ERROR: check the stats_property for min/max')
