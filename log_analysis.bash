#!/bin/bash -f

log_dir='/Users/ehsan/Documents/Python_projects/CMAQ_analysis/log_dir/'
#log_dir = '/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/CMAQ_analysis/logs/'

current_dir=$(pwd)
echo '-> we are currently at='
echo $current_dir

echo '-> log directory=' 
echo $log_dir

echo '-> change to log dir...'
cd $log_dir

work_dir=$(pwd)
echo '-> we are at work directory=' 
echo $work_dir

echo '-> list of files at current dir='
ls -la .

echo '-> list of log files in work directory='
log_list=$(ls *.txt)


echo $log_list



#foreach logfile ( $log_list ) 

#	echo '-> working on log file=' $logfile

#end


#set vmin_list = 'grep -irnH 'vmin' log*'

#foreach log_line ( vmin_list )

#	echo '-> line is=' $ log_line