#!/bin/csh -f

echo '-> set environmental variables first'

setenv CCTM_PROCESS 'dep'  				# 'atm' or 'dep'
setenv DEP_TYPE 'DRYDEP'				# we only support DRYDEP for depossition
setenv PROCESSING_POLLUTANT 'single_pollutant'  	# 'pm2.5' OR 'single_pollutant'== nh3,o3,no2,no,co - either for 'atm' or 'dep'
setenv POL_UNIT 'kg/hectare' 				# 'ppmV' or 'ug/m^3' 	or for depossition- 'kg/hectare'
setenv SPATIAL_PLOTTING_KEY 'yes'			# yes or no
setenv PLOT_METHOD 	'diff_plot'			# 'single_plot' or 'diff_plot'
setenv COLOR_METHOD 'min_to_max'		# 'zero_to_max' , 'min_to_max' , 'minus_abs_max_to_max'
setenv PRODUCE_RASTER 'no'				# 'yes' OR 'no'
setenv TIMESERIES_PLOTTING 'no' 			# yes or not

# if ==  minus_abs_max_to_max
setenv MINUS_ABS_MAX_DIFF '-8.306249300485834e-07'
setenv ABS_MAX_DIFF '8.306249300485834e-07'
echo '------------------------------------------------------'
echo '-> cctm process=' ${CCTM_PROCESS}
if ( $CCTM_PROCESS == 'dep' ) then
	echo '-> depossition type=' ${DEP_TYPE}
endif
echo '-> processing pollutant=' ${PROCESSING_POLLUTANT}
echo '-> pollutant unit=' ${POL_UNIT}
echo '-> spatial plotting=' ${SPATIAL_PLOTTING_KEY}
echo '-> plotting mehtod=' ${PLOT_METHOD}
echo '-> colorbar method=' ${COLOR_METHOD}
echo '-> produce raster=' ${PRODUCE_RASTER}
echo '-> time-series plotting=' ${TIMESERIES_PLOTTING}
echo '-> minus-abs max values for diff plot=' ${MINUS_ABS_MAX_DIFF}
echo '-> abs max values for diff plot=' ${ABS_MAX_DIFF}
echo '------------------------------------------------------'

#foreach species ('NH3')		#	( 'CO' 'O3' 'NH3' 'HNO3' 'NO2' 'SO2' ) # 'PM2.5' ) 	# atm
foreach species ( 'HNO3' )       # ( 'HONO' 'N2O5'  'NO3' 'NO'  'NO2' 'HNO3' 'NH3' )  # dep      
 

	setenv CMAQ_POL ${species}                                    # for plot title 'CO','PM2.5','NH3','O3','HNO3','NO2','SO2'

		foreach scen_no ( 1 2 3 4 5 )  # ( 1 2 3 4 5 )

			setenv LANDIS_SCENARIO ${scen_no}

			foreach month_name ( 'jul' 'aug' 'sep' 'oct' 'nov' ) # ( 'jul' 'aug' 'sep' 'oct' 'nov' ) 
				
				if ($month_name == 'jul') then

					setenv CMAQ_MONTH_STRING 'Jul'
					setenv CMAQ_MONTH_NUMBER '07'
					setenv DAYS_IN_MONTH_TO_RUN '31'  #31

					setenv JOB_NAME 'S'${scen_no}${CMAQ_MONTH_STRING}

				else if ($month_name == 'aug') then

					setenv CMAQ_MONTH_STRING 'Aug'
					setenv CMAQ_MONTH_NUMBER '08'
					setenv DAYS_IN_MONTH_TO_RUN '31' # 31

					setenv JOB_NAME 'S'${scen_no}${CMAQ_MONTH_STRING}

				else if ($month_name == 'sep') then

					setenv CMAQ_MONTH_STRING 'Sep'
					setenv CMAQ_MONTH_NUMBER '09'
					setenv DAYS_IN_MONTH_TO_RUN '30' #30

                                        setenv JOB_NAME 'S'${scen_no}${CMAQ_MONTH_STRING}

				else if ($month_name == 'oct') then

					setenv CMAQ_MONTH_STRING 'Oct'
					setenv CMAQ_MONTH_NUMBER '10'
					setenv DAYS_IN_MONTH_TO_RUN '31'	#'31'

                                        setenv JOB_NAME 'S'${scen_no}${CMAQ_MONTH_STRING}

				else if ($month_name == 'nov') then

					setenv CMAQ_MONTH_STRING 'Nov'
					setenv CMAQ_MONTH_NUMBER '11'
					setenv DAYS_IN_MONTH_TO_RUN '29' #'29'

                                        setenv JOB_NAME 'S'${scen_no}${CMAQ_MONTH_STRING}

				else

					echo '-> ERROR: month setting wrong!'

				endif

				echo '-> job name=' ${JOB_NAME}

				qsub ./job_scen${LANDIS_SCENARIO}_${CMAQ_MONTH_STRING}.csh
				echo '-> job script= job_scen'${scen_no}'_'${month_name}'.csh' '-->' ${CMAQ_POL}
		
			end
		end
end
