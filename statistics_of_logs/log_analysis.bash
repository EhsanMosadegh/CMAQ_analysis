#!/bin/bash -f
#------------------------------------------------------------------------
# run-time setting
# NOTE: this script writes its output to <logs/> directory
# to run: bash "script_name" > POL_logs/log.POL.minOrMax.txt
# to check the end result either open the log file and check the end of the log file for abs min/max values, 
# or do-> grep 'chk' on_the_logfile_of_bash_script
# I defined 'chk' as a keyword to check importent parameters and end value
#------------------------------------------------------------------------

pollutant='SO2'
stats_property='max'

stats_pattern=$stats_property'DiffMesh'
log_file_pattern=log.${pollutant}.scen*

#---------------------------------

#home_dir='/Users/ehsan/Documents/Python_projects/CMAQ_analysis'
home_dir='/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/CMAQ_analysis/'
log_dir=${home_dir}'logs/min_max_logs/'
statistics_dir=${home_dir}'statistics_of_logs/'

#---------------------------------

echo '-------------------------------'
echo '-> chk: pollutant is=' ${pollutant}
echo '-> chk: statistical property is=' ${stats_property}
echo '-> stats patter is=' ${stats_pattern}
echo '-------------------------------'
echo '  '

#---------------------------------

output_file_name=$stats_pattern'_list_total_for_'$pollutant.txt

#---------------------------------

current_dir=$(pwd)
echo '-> we are currently at='
echo $current_dir

echo '-> log directory=' 
echo $log_dir

echo '-> we change dir to <logs/>(work) dir...'
cd $log_dir

work_dir=$(pwd)
echo '-> now we are at=' 
echo $work_dir

echo '-> we list of files at current dir='
ls -la .

# check and remove the file is it exist before
echo '-> NOTE: old output files will be removed first:' 
echo $output_file_name 
echo $'log_list_for_'${pollutant}_${stats_property}.txt

rm ./$output_file_name
rm ./log_list_for_${pollutant}_${stats_property}.txt

#echo '-> log files in work directory='

#log_files=$(grep -irnH 'vmin' .)
echo '-> now we get the list of *log* files in the log directory and write to file'
ls $log_file_pattern > log_list_for_${pollutant}_${stats_property}.txt

echo '-> number of lines in log list=' 
wc -l log_list_for_${pollutant}_${stats_property}.txt 

echo '-> loop and read in log_list... '
while read log_list_line
do 
	echo ' '
	echo '-> opening log file=' $log_list_line
	echo '-> capturing stats_pattern=' $stats_pattern

	thePattern=$(grep -o $stats_pattern.* $log_list_line)
	echo '-> extracting pattern=' $thePattern

	grep -o $stats_pattern.* $log_list_line | cut -f2 -d' ' >> $output_file_name

done < log_list_for_${pollutant}_${stats_property}.txt
echo ' '
echo '-> size of the' $stats_pattern 'output list is=' 
cat $output_file_name | py -l 'print(len(l))'
echo '------------------------------------------------'
echo '-> now we do the arithmetic operations in python since shell 
					does _not_ understand arithmetic operations specially with scientific notation!!!'
echo ' '
echo '-> we are currently at=' $(pwd)
python ${statistics_dir}/get_min_max.py
echo '------------------------------------------------'

#if [ $stats_property == 'min' ]; then
	
#	echo '-> absolute' $stats_property 'value from the file:' $output_file_name ', is='
#	cat $output_file_name | py -l 'min(l)'

#elif [ $stats_property == 'max' ]; then

#        echo '-> absolute' $stats_property 'value from the file:' $output_file_name ', is='
#        cat $output_file_name | py -l 'max(l)'

#else

#	echo '-> WARNING: set the stats_property correctly...'

#fi

#echo '-----------------------------------------------------------------'

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
