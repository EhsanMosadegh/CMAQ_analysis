#!/bin/csh -f

set root_dir = '/Users/ehsan/Documents/Python_projects/CMAQ_analysis/log_test/'

cd root_dir

echo '-> root directory=' $root_dir

echo '-> list of log files in root directory='

ls log*

set log_list = ( ls log* )

foreach logfile ( log_list ) 

	echo '-> log file is='
	echo $logfile

#end
