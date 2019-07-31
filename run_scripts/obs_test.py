#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

obs_station_name = 'TahoeCity'

cmaq_file_month = '07'
number_of_days = 2

home_dir = '/Volumes/USFSdata/' 
obs_dir = home_dir+'obs/'

obs_file = 'ozone_'+obs_station_name+'_month_'+cmaq_file_month+'.csv'
obs_input = obs_dir + obs_file

print(f'-> obs input file is= {obs_input} ')

input_df = pd.read_csv( obs_input , sep=',' , header=0 )

#x_ = input_df['date']
y_ = input_df['value'][0:number_of_days*24]

plt.plot( y_ , label='obs o3' , color='red')

plt.legend(title='obs legend')
plt.show()
