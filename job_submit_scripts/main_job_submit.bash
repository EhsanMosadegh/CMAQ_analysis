#!/bin/bash -f

for scen_no in 1 2 3 4 5
do
	for month_name in 'jul' 'aug' 'sep' 'oct' 'nov'
	do
		qsub ./job_co_scen${scen_no}_${month_name}.csh
		#echo '-> job script= job_co_scen'${scen_no}'_'${month_name}'.csh'
	done
done