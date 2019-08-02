#!/usr/bin/env python3

import cv2
import os

#==========================================================

def get_the_day(list_elem):
	"gets each element and returns the first 2 characters"
	return list_elem[0:2]

#==========================================================
scenario='5'
month='nov'
pollutant='co'

image_folder = '/Users/ehsan/Documents/Python_projects/CMAQ_analysis/animations/daily_plots/scen_'+scenario+'/'+month+'/'
#image_folder = '/storage/ehsanm/~/daily_plots/'+'scen_'+scenario+'/'+month

print(f'-> scenario is= {scenario}')
print(f'-> month is= {month}')
print(f'-> pollutant is= {pollutant}')
print(f'-> image folder is: {image_folder}')

video_name = 'cmaq_'+pollutant+'_daily_plots_scen_'+scenario+'_'+month+'.avi' # only avi

#==========================================================

print('-> making the video...')
# sorting based on fig number??????

image_list_unsorted = [img for img in os.listdir(image_folder) if img.endswith(".png")]
#print('-> original list=')
#print(image_list_unsorted)

image_list_sorted = sorted(image_list_unsorted , key=get_the_day)  # sorts the list based on first 2 days on the nametag
#print('-> sorted list=')
#print(image_list_sorted)

#==========================================================

frame = cv2.imread(os.path.join(image_folder, image_list_sorted[0]))
height, width, layers = frame.shape
fps=0.5
video = cv2.VideoWriter(video_name, 0, fps, (width,height))

print('-> writing the video from the following plots...')

for image in image_list_sorted:

		#print('-> adding the following plot to the video:')
		print(image)

		video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()