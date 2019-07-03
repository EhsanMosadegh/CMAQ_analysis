#!/bin/env python3
import numpy as np
import os
#--------------------------------------------------------

pollutant='co'
stats_property='min'

#--------------------------------------------------------

print('------------------------------------------------')
print('-> now we are inside python...')
print( f'-> we are at= {os.getcwd() }')
print('------------------------------------------------')
print( f'-> chk: pollutant is= {pollutant}')
print( f'-> chk: statisticsl property is= {stats_property}')
print('------------------------------------------------')
#--------------------------------------------------------

file_name=stats_property+'DiffMesh_list_total_for_'+pollutant+'.txt'

#home_dir='/Users/ehsan/Documents/Python_projects/CMAQ_analysis'
home_dir='/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/CMAQ_analysis'

file_path=home_dir+'/logs/'

log_file=file_path+file_name

print(f'-> log file is= {log_file}')
print( " ")
with open( log_file , 'r') as input_file :

	input_list=[]

	for line in input_file:
		print( f'-> line to append is= {line}' )
		input_list.append( line )

	# change the data type from string to float!
	input_list=np.array(input_list).astype(np.float)

	if (stats_property=='min') :

		print( f'-> chk: min of the list= {np.min(input_list)}' )
#		print( np.min(input_list) )

	elif (stats_property=='max') :

		print( f'-> chk: max of the list= {np.max(input_list)}' )
#		print( np.max(input_list) )

	else:

		print('-> ERROR: check the stats_property for min/max')
