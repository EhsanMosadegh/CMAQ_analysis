#!/bin/csh -f

echo '-> set environmental variables first'
#setenv 




echo '-> looping for scen and month...'

foreach scen_no (1 2 3 4 5)

	foreach month_name ('jul' 'aug' 'sep' 'oct' 'nov')

		qsub ./job_co_scen${scen_no}_${month_name}.csh
		#echo '-> job script= job_co_scen'${scen_no}'_'${month_name}'.csh'
	end
end
