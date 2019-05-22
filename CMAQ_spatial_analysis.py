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

# import input files
#cmaq_file =
mcip_file = '/Users/ehsan/Documents/Python_projects/USFS_fire/inputs/cmaq_inputs/GRIDDOT2D_161001'

mcip_input = Dataset( mcip_file )
#cmaq_input = data_meshset( cmaq_input )

# some info
print('-> MCIP file dimensions: %s' %str( mcip_input['LATD'].dimensions ) )
print('-> shape of each dimension: %s' %( str(mcip_input['LATD'].shape ) ))
print('-> plotting the data...')

# extract lat and lon parameteres
lat_mesh = np.array( mcip_input['LATD'][0][0][:][:] ) # select only rosws and cols for the 1st timestep and layer
lon_mesh = np.array( mcip_input['LOND'][0][0][:][:])
data_mesh = np.random.rand(265,250)*10

# plot dots from grid coordinates of the dots
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
                           resolution='c')

basemap_instance.bluemarble()

x_mesh, y_mesh = basemap_instance(lon_mesh , lat_mesh) # order: x , y, transforms from degree to meter for LCC

basemap_instance.drawmapboundary(color='k' ) #, fill_color='aqua')
basemap_instance.drawcoastlines(color = '0.15')
basemap_instance.drawcounties()
basemap_instance.drawstates()

#basemap_instance.fillcontinents(lake_color='aqua')

im1 = basemap_instance.pcolormesh(x_mesh , y_mesh , data_mesh , cmap=plt.cm.Reds_r , shading='flat')
#im2 = basemap_instance.pcolormesh(lon_mesh , lat_mesh , data_mesh , cmap=plt.cm.jet , shading='flat')

cb = basemap_instance.colorbar(im1 , 'bottom')

#cs = basemap_instance.contourf(lon_mesh , lat_mesh , data_mesh)
#cbar = basemap_instance.colorbar(cs, location='bottom')


plt.show()