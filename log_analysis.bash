#!/bin/bash -f

#log_dir='/Users/ehsan/Documents/Python_projects/CMAQ_analysis/log_dir'
log_dir = '/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/CMAQ_analysis/logs/'

stats_property=min
stats_pattern='minDiffMesh'
log_file_pattern=log.o3*

mon=10
scenario=1
output_file_name=$stats_pattern'_scen_'$scenario'_month_'$mon.txt

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

#echo '-> log files in work directory='

#log_files=$(grep -irnH 'vmin' .)
echo '-> get the list of log files in the log directory and write to file'
ls $log_file_pattern > log_list.txt

echo '-> get the number of lines in log list'
wc -l log_list.txt

# check and remove the file is it exist before
echo '-> output file will be removed first:' $output_file_name
rm ./$output_file_name

echo '-> loop and read in log_list'
while read log_list_line
do 
	echo '-> log file is=' $log_list_line
	echo '-> capturing stats_pattern=' $stats_pattern

	thePattern=$(grep -o $stats_pattern.* $log_list_line)
	echo '-> extracting pattern=' $thePattern

	grep -o $stats_pattern.* $log_list_line | cut -f2 -d' ' >> $output_file_name

done < log_list.txt

echo '-> size of the' $stats_pattern 'list is=' 
cat $output_file_name | py -l 'print(len(l))'

echo '-> absolute' $stats_property 'value from the file:' $output_file_name 'is='
cat $output_file_name | py -l 'min(l)'


#--- method 1 to read each line

# echo '-> open file and do for loop to show the lines ...'
#??? how open the file first?

# for word in $log_files; do

# 	echo the new word= $word
# done

#--- method 2 to read each line

# filename=$work_dir
# n=1
# echo '-> read file content line-by-line'

# while read eachLogFile
# do 
# 	echo '-> line no.=' $n
# 	echo '-> the line is=' $eachLogFile
# 	new_stats_pattern=grep 'vmin=' $eachLogFile
# 	echo $new_stats_pattern
# 	#n=$((n+1))
# 	let n=$n+1

# done < logfile.txt