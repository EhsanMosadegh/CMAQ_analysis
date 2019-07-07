#!/bin/csh -f

echo '-> set environmental variables first'

setenv PROCESSING_POLLUTANT 'pm2.5'			# 'pm2.5' OR 'single_pollutant'== nh3,o3,no2,no,co
setenv POL_UNIT 'ug/m^3'				# 'ppmV' or 'ug/m^3'
setenv SPATIAL_PLOTTING_KEY 'yes'			# yes or no
setenv PLOT_METHOD 	'diff_plot'			# 'single_plot' or 'diff_plot'
setenv COLOR_METHOD 'min_to_max'	  		# 'zero_to_max' , 'min_to_max' , 'minus_abs_max_to_max'
setenv PRODUCE_RASTER 'no'				# 'yes' OR 'no'
setenv TIMESERIES_PLOTTING 'yes' 			# yes or not

# if ==  minus_abs_max_to_max
setenv MINUS_ABS_MAX_DIFF ''
setenv ABS_MAX_DIFF ''

echo '-> processing pollutant=' ${PROCESSING_POLLUTANT}
echo '-> pollutant unit=' ${POL_UNIT}
echo '-> spatial plotting=' ${SPATIAL_PLOTTING_KEY}
echo '-> plotting mehtod=' ${PLOT_METHOD}
echo '-> colorbar method=' ${COLOR_METHOD}
echo '-> produce raster=' ${PRODUCE_RASTER}
echo '-> time-series plotting=' ${TIMESERIES_PLOTTING}
echo '-> minus-abs max values for diff plot=' ${MINUS_ABS_MAX_DIFF}
echo '-> abs max values for diff plot=' ${ABS_MAX_DIFF}


foreach species ('PM2.5')  # ('CO' 'O3' 'NH3' 'HNO3' 'NO2' 'SO2' 'PM2.5')

	setenv CMAQ_POL ${species}                                    # for plot title 'CO','PM2.5','NH3','O3','HNO3','NO2','SO2'

		foreach scen_no (1 2)

			foreach month_name ('jul' 'aug' 'sep' 'oct')
			
				qsub ./job_scen${scen_no}_${month_name}.csh
				echo '-> job script= job_scen'${scen_no}'_'${month_name}'.csh' '-->' ${CMAQ_POL}
		
			end
		end
end