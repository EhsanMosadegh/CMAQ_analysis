#!/bin/bash -f
#------------------------------------------------------------------------
# run-time setting
# NOTE: this script writes its output to <logs/> directory

# to run: 
# set the first 3 variables in both this script and python scripts
# check the directory paths 
# bash "script_name" > log_file_name.txt (POL_logs/log.POL.minOrMax.txt)

# to check the end result either open the log file and check the end of the log file for abs min/max values, 
# or do-> grep 'chk' on_the_logfile_of_bash_script
# I defined 'chk' as a keyword to check importent parameters and end value
#------------------------------------------------------------------------

pollutant='NO'
stats_property='max'
region='Ndep'  # 'LakeTahoeBasin' or 'NorthTahoe' or 'SouthTahoe' or 'Ndep'

#stats_pattern=$stats_property'DiffMesh'
stats_pattern=$stats_property'DiffMesh'$region

log_file_pattern=log.${pollutant}.scen*

#---------------------------------

#home_dir='/Users/ehsan/Documents/Python_projects/CMAQ_analysis/'
home_dir='/storage/ehsanm/USFS_CA_WRF_1km_project/data_analysis/CMAQ_analysis/' # do not add logs/ at the end 
log_dir_name='drydep/'
log_dir_full_path=${home_dir}'logs/'${log_dir_name}  # logs/ dir is added here
statistics_dir=${home_dir}'statistics_of_logs/'

#---------------------------------

echo '-------------------------------'
echo '-> chk: pollutant is=' ${pollutant}
echo '-> chk: statistical property is=' ${stats_property}
echo '-> chk: stats patter is=' ${stats_pattern}
echo '-> chk: region is=' ${region}
echo '-------------------------------'
echo '  '

#---------------------------------

#output_file_name=$stats_pattern'_list_total_for_'$pollutant'.txt'
output_file_name=$stats_pattern'_list_total_for_'$pollutant'_in_'${region}.txt
echo '-> output log file name is =' ${output_file_name}
echo ' ' 
#---------------------------------
echo '-> home dir=' ${home_dir}
echo '-> log dir=' ${log_dir_name}
echo '-> log dir full path=' ${log_dir_full_path}

current_dir=$(pwd)
echo '-> we are currently at=' $current_dir

echo '-> log directory=' $log_dir_name

echo '-> we change dir to <logs/>(work) dir...' 
cd $log_dir_full_path

log_dir=$(pwd)
echo '-> now we are at=' $log_dir

echo '-> list of files at current (hopefully be logs/) dir='
ls -la .

# check and remove the file is it exist before
echo '-> NOTE: following old output files will be removed first:' 
echo $output_file_name 
echo $'log_list_for_'${pollutant}_${stats_property}.txt

rm ./$output_file_name
rm ./log_list_for_${pollutant}_${stats_property}.txt

echo ' '

echo '-> now we get the list of *log* files in the log directory and write to file'
ls $log_file_pattern > log_list_for_${pollutant}_${stats_property}.txt
echo ' '
echo '-> number of lines in log list=' 
wc -l log_list_for_${pollutant}_${stats_property}.txt 
echo ' '
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
echo '-> size of the' $stats_pattern 'output log file with no. of captured patterns is=' 
cat $output_file_name | py -l 'print(len(l))'
echo '------------------------------------------------'
echo '-> now we do the arithmetic operations in python since shell \
	does _not_ understand arithmetic operations specially with scientific notation!!!'

echo ' '
echo '-> we are currently at=' $(pwd)
echo '-> deploying python script ...'
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
