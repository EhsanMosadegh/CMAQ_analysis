#!/bin/env python3

file='/Users/ehsan/Documents/Python_projects/CMAQ_analysis/logs/minDiffMesh_list_total_for_co.txt'
# input_file=open(file , 'r')

with open( file , 'r') as input_file :

	input_list=[]

	for line in input_file:
		print( line )
		input_list.append( line )

	print( '-> min of the list= ')
	print( min(input_list))
