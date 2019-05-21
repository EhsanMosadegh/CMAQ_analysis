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
from mpl_toolkits.basemap import Basemap as bm

# import input files
#cmaq_file =
mcip_file = '/Users/ehsan/Documents/Python_projects/USFS_fire/inputs/cmaq_inputs/GRIDDOT2D_161001'

mcip_input = Dataset( mcip_file )
#cmaq_input = Dataset( cmaq_input )

# some info
print('-> MCIP file dimensions: %s' %str( mcip_input['LATD'].dimensions ) )
print('-> shape of each dimension: %s' %( str(mcip_input['LATD'].shape ) ))

# extract lat and lon parameteres
lat_list = np.array( mcip_input['LATD'][0][0][:][:] ) # select only rosws and cols for the 1st timestep and layer
lon_list = np.array( mcip_input['LOND'][0][0][:][:])
data = np.random.rand(265,250)

# plot dots from grid coordinates of the dots
#plt.plot( lon_list , lat_list , marker='.' , color='b' , linestyle= 'none' )

# plot the domain/region borders
xcent = -120.805999755859 # degrees
ycent = 40.0 # degrees
NROWS = 266*1000 # meters
NCOLS = 251*1000 # meters

bmap = bm( projection='lcc', lat_0= ycent , lon_0= xcent , height= NROWS , width= NCOLS )
bmap.drawmapboundary(color='k' ) #, fill_color='aqua')
bmap.drawcoastlines()#color = '0.15')
bm.pcolormesh(lat_list , lon_list , data )
plt.show()