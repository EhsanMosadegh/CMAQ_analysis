
##########################################################################################
# author: Ehsan Mosadegh
# usage: to plot observation data from AQ surface stations
# date: May 3, 2018
# email: ehsanm@dri.edu
# notes: change dtype of some columns from float to object (string) to be able to compare to str. dtype was float before.
#
##########################################################################################


import os
from   netCDF4 import Dataset
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
from   matplotlib.dates import drange, DateFormatter
import time

##########################################################################################
start_time = time.time()
#----------------------------------------------|
#           set input parameters               |
#----------------------------------------------|
pol_index = 3            # set for a pollutant |
stn_index = 3            # set for a station   |
VAR_unit =        'Micrograms.m$^{-3}$' # based on pollutant    O3-CO->'PPM'  other PPB
#----------------------------------------------|
cmaq_layer = 1          # set for a sfc layer  |
utc_conversion = 7  # for August; keep in UTC  |
#----------------------------------------------|

# define pollutant data

pol_index_ref = [0      ,1                ,2      ,3      ,4                 ,5               ]
pol_name =      ['Ozone','Carbon_monoxide','PM10' ,'PM2.5','Nitrogen_dioxide','Sulfor_dioxide']
cmaq_VAR =      ['O3'   ,'CO'             ,'PM10' ,'PM2.5','NO2'             ,'SO2'           ] # check with CMAQ
obs_PolCode =   ['44201','42101'          ,'81102','88101','42602'           ,'42401'         ]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define station data

stn_index_ref =[0                 ,1         ,2           ,3           ,4           ,5              ,6              ,7           ,8           ,9               ,10                 ]
stn_name =     ['Incline_Village' ,'Carson'  ,'Reno3'     ,'Sparks'    ,'EchoSummit','LemmonValley' ,'SouthReno'    ,'TahoeCity' ,'Toll'      ,'SouthLakeTahoe','SE_Gardnerville'  ]
stn_lat =      [39.250409         ,39.1447   ,39.525083   ,39.540917   ,38.81161    ,39.645264      ,39.469219      ,39.166017   ,39.399837   ,38.944979       ,38.897557          ]
stn_lon =      [-119.956738       ,-119.7661 ,-119.807717 ,-119.746761 ,-120.033084 ,-119.840025    ,-119.775354    ,-120.148833 ,-119.739606 ,-119.970609     ,-119.732507        ]
StateCode =    ['32'              ,'32'      ,'32'        ,'32'        ,'06'        ,'32'           ,'32'           ,'06'        ,'32'        ,'06'            ,'32'               ]
CountyCode =   ['031'             ,'510'     ,'031'       ,'031'       ,'017'       ,'031'          ,	'031'          ,'061'       ,'031'       ,'017'           ,'005'              ]
site_no =      ['2002'            ,'0020'    ,'0016'      ,'1005'      ,'0012'      ,'2009'         ,'0020'         ,'1004'      ,'0025'      ,'0011'          ,'0007'             ]
#obs_file_name = stn_name[stn_index]+'_'+StateCode[stn_index]+'_'+CountyCode[stn_index]+'_'+site_no[stn_index]+'.txt'

print('****************************************************************************************************')
print('-> doing calculaiton for: \
pollutant: "%s"; station name: "%s"; station ID: "%s"'\
%( pol_name[pol_index] , stn_name[stn_index] , (StateCode[stn_index]+'_'+CountyCode[stn_index]+'_'+site_no[stn_index])))
print('****************************************************************************************************')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define CMAQ input files

cmaq_conc_file_list = [	'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140810.nc',\
                       	'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140811.nc',\
                    	'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140812.nc',\
                    	'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140813.nc',\
                    	'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140814.nc',\
                    	'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140815.nc']

cmaq_pmdiag_file_list = ['CCTM_PMDIAG_v5.2_cb6_intel_Tahoe_WRF_1km_20140810.nc',\
                    'CCTM_PMDIAG_v5.2_cb6_intel_Tahoe_WRF_1km_20140811.nc',\
                    'CCTM_PMDIAG_v5.2_cb6_intel_Tahoe_WRF_1km_20140812.nc',\
                    'CCTM_PMDIAG_v5.2_cb6_intel_Tahoe_WRF_1km_20140813.nc',\
                    'CCTM_PMDIAG_v5.2_cb6_intel_Tahoe_WRF_1km_20140814.nc',\
                    'CCTM_PMDIAG_v5.2_cb6_intel_Tahoe_WRF_1km_20140815.nc']


#cmaq_apmdiag_file_list = ['CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140810.nc',\
#                    'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140811.nc']#,\
#                    'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140812.nc',\
#                    'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140813.nc',\
#                    'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140814.nc',\
#                    'CCTM_CONC_v5.2_cb6_intel_Tahoe_WRF_1km_20140815.nc']


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set pathes based on work directory:

work_dir = '/Users/ehsan/Documents/my_Python_codes/cmaq_sanity'
script_dir = work_dir+'/scripts'
output_dir = work_dir+'/outputs'
figure_dir = work_dir+'/figures'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


print ('-> work directory is: "%s"' %work_dir)
print ('-> script directory is: "%s"' %script_dir)
print('-> output directory is %s' %output_dir)
print ('-> changing directory to work directory ...')
os.chdir(work_dir)

# set input parameters
##########################################################################################


##########################################################################################
# read time-invariant VAR (lat/lon) from MCIP file

mcip_file_name = 'GRIDCRO2D_LakeTahoe'
mcip_file_path = '/Volumes/Ehsanm_DRI/MCIP_4_3'
mcip_file_full_path = os.path.join(mcip_file_path,mcip_file_name)

if (os.path.isfile(mcip_file_full_path)!=True):
    print('-> ERROR: MCIP input file does NOT exist at path: "%s"' %mcip_file_path)
else:
    print('-> MCIP input file exists!')

    # read in MCIP file
    mcip_read = Dataset(mcip_file_full_path,'r')

    # get lat/lon of center of CMAQ cells from MCIP file
    lat_mcip = np.array(mcip_read['LAT'])
    lon_mcip = np.array(mcip_read['LON'])

    # find the CMAQ cell that our station is inside it
    total_diff = np.abs(lon_mcip - stn_lon[stn_index]) + np.abs(lat_mcip - stn_lat[stn_index])
    np.shape(total_diff)
    [row_index,col_index] = np.argwhere(total_diff == np.min(total_diff))[0,2:4]
    print('-> x-index of station cell is: %s'%row_index)
    print('-> y-index of station cell is: %s'%col_index)
    mcip_read.close()

# read time-invariant VAR (lat/lon) from MCIP file
##########################################################################################


##########################################################################################
# read CMAQ files

# set CMAQ directory pathes
cmaq_input_dir = '/Volumes/Ehsanm_DRI/CMAQ_files/cmaq_baseline'
cmaq_conc_file_name = cmaq_conc_file_list[0]
cmaq_CONC_file_full_path = os.path.join(cmaq_input_dir,cmaq_conc_file_name)
print('-> CMAQ input dir. full path is: %s' %cmaq_input_dir)


print('-> check to see if CMAQ CONC file exists...')
if (os.path.isfile(cmaq_CONC_file_full_path)!=True):
    print('-> ERROR: CMAQ CONC input file does NOT exist at path: "%s"' %cmaq_input_dir)
else:
    print('-> NOTE: CMAQ CONC input file exists!')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # total date-time list for the modeling period with CMAQ daily/hourly date-time
    dt_CMAQ_CONC_LocalTime_TotalList = []
    dt_CMAQ_PMDIAG_LocalTime_TotalList = []

    # total list for VAR, extracted from CMAQ daily files
    CMAQ_VAR_TimeSeries_4WholePeriod_List = []

    PM10AT_3_list = []
    PM10AC_3_list = []
    PM10CO_3_list = []

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('-> there are (%s) CONC files inside CMAQ input directory.' %(len(cmaq_conc_file_list)))
    print('-> the following CONC files are in CMAQ input directory:')

    for cmaq_conc_file in cmaq_conc_file_list:
        print(cmaq_conc_file)

    print('-> read/loop over every CMAQ daily input files:')
    for cmaq_conc_file_no in range(len(cmaq_conc_file_list)):

        cmaq_conc_file_name = cmaq_conc_file_list[cmaq_conc_file_no]
        cmaq_CONC_file_full_path = os.path.join(cmaq_input_dir,cmaq_conc_file_name)


        print('**************************************************************************')
        print('-> opening CMAQ input file no.:(%s)' %(cmaq_conc_file_no+1))

        # open and read each cmaq CONC file seperately; should use the full file path in Dateset;
        read_from_daily_CMAQ_CONC_file = Dataset(cmaq_CONC_file_full_path,'r')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # loop for CONC date-times: TFLAGs

        cmaq_conc_dt_daily_list = []
        print('-> date-time loop over CMAQ file no. %s...' %(cmaq_conc_file_no+1))
        cmaq_conc_tflag_i = read_from_daily_CMAQ_CONC_file.variables['TFLAG']
        # extract date-time for one day/loop
        cmaq_conc_datetimes_i = cmaq_conc_tflag_i[:,0,:] # for only one VAR is enough='0'; only for TSTEP and DATE-TIME

        # loop for 24-hour
        for i in range(len(cmaq_conc_datetimes_i)):
            # calculate date:
            date_conc_i = cmaq_conc_datetimes_i[i,0]
            date_conc_i_obj = dt.datetime.strptime(str(date_conc_i),'%Y%j').date() # to get the date

            # calculate time:
            time_conc_i = cmaq_conc_datetimes_i[i,1]
            if time_conc_i < 100000:
                if time_conc_i == 0 :
                    time_conc_i_str = str(time_conc_i).rjust(6,'0') # final width=6, pad with '0''
                    #time_conc_i_str = time_conc_i_str+'0000'
                else:
                    time_conc_i_str = str(time_conc_i).rjust(6,'0') # final width=6, pad with '0''
            elif time_conc_i >= 100000:
                time_conc_i_str = str(time_conc_i)

            # make dt obj from str
            time_conc_i_obj = dt.datetime.strptime(time_conc_i_str,'%H%M%S').time() # just get the time

            # combine date-obj and time-obj together
            dt_conc_obj_i = dt.datetime.combine(date_conc_i_obj,time_conc_i_obj)

            # add each dt_obj_hr to a 24hr list
            cmaq_conc_dt_daily_list.append(dt_conc_obj_i - dt.timedelta(hours=utc_conversion))

            # now delete the 1st hr=0, so the list starts from 1-24
            cmaq_conc_dt_daily_list_new = np.delete(cmaq_conc_dt_daily_list , 0)

        # append each dt_conc_obj_i to the total list
        # now appending CMAQ daily date-time object to its total list...')
        dt_CMAQ_CONC_LocalTime_TotalList.extend(cmaq_conc_dt_daily_list_new)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # loop to extract CMAQ VAR
        # read CMAQ file and extract VARs
        VAR = cmaq_VAR[pol_index]

        if VAR == 'PM2.5' or VAR=='PM10':  # for 2 or more conditions put in tuple
            print('-> NOTE: VAR is "%s", so we start extracting all species for "%s"...\
            and will delete the 1st row=hr=0 from each CONC file to make 24 hrs in each file,\
            so each CONC species will have same no. of time-steps (1-0) similar to PMDIAG species' %(VAR,VAR))

            # extract PM species from CMAQ_CONC_* files ...')
            # delete 0 from all files to start from 1-24
            ANAI = read_from_daily_CMAQ_CONC_file['ANAI'][:,cmaq_layer,row_index,col_index]
            ANAI = np.delete(ANAI , 0)

            ACLI = read_from_daily_CMAQ_CONC_file['ACLI'][:,cmaq_layer,row_index,col_index]
            ACLI = np.delete(ACLI , 0)

            AOTHRI = read_from_daily_CMAQ_CONC_file['AOTHRI'][:,cmaq_layer,row_index,col_index]
            AOTHRI = np.delete(AOTHRI , 0)

            ANAJ = read_from_daily_CMAQ_CONC_file['ANAJ'][:,cmaq_layer,row_index,col_index]
            ANAJ = np.delete(ANAJ , 0)

            ACLJ = read_from_daily_CMAQ_CONC_file['ACLJ'][:,cmaq_layer,row_index,col_index]
            ACLJ = np.delete(ACLJ , 0)

            AOTHRJ = read_from_daily_CMAQ_CONC_file['AOTHRJ'][:,cmaq_layer,row_index,col_index]
            AOTHRJ = np.delete(AOTHRJ , 0)

            AFEJ = read_from_daily_CMAQ_CONC_file['AFEJ'][:,cmaq_layer,row_index,col_index]
            AFEJ = np.delete(AFEJ , 0)

            ASIJ = read_from_daily_CMAQ_CONC_file['ASIJ'][:,cmaq_layer,row_index,col_index]
            ASIJ = np.delete(ASIJ , 0)

            ATIJ = read_from_daily_CMAQ_CONC_file['ATIJ'][:,cmaq_layer,row_index,col_index]
            ATIJ = np.delete(ATIJ , 0)

            ACAJ = read_from_daily_CMAQ_CONC_file['ACAJ'][:,cmaq_layer,row_index,col_index]
            ACAJ = np.delete(ACAJ , 0)

            AMGJ = read_from_daily_CMAQ_CONC_file['AMGJ'][:,cmaq_layer,row_index,col_index]
            AMGJ = np.delete(AMGJ , 0)

            AMNJ = read_from_daily_CMAQ_CONC_file['AMNJ'][:,cmaq_layer,row_index,col_index]
            AMNJ = np.delete(AMNJ , 0)

            AKJ = read_from_daily_CMAQ_CONC_file['AKJ'][:,cmaq_layer,row_index,col_index]
            AKJ = np.delete(AKJ , 0)

            AALJ = read_from_daily_CMAQ_CONC_file['AALJ'][:,cmaq_layer,row_index,col_index]
            AALJ = np.delete(AALJ , 0)

            ASEACAT = read_from_daily_CMAQ_CONC_file['ASEACAT'][:,cmaq_layer,row_index,col_index]
            ASEACAT = np.delete(ASEACAT , 0)

            ACORS = read_from_daily_CMAQ_CONC_file['ACORS'][:,cmaq_layer,row_index,col_index]
            ACORS = np.delete(ACORS , 0)

            ASOIL = read_from_daily_CMAQ_CONC_file['ASOIL'][:,cmaq_layer,row_index,col_index]
            ASOIL = np.delete(ASOIL , 0)

            ACLK = read_from_daily_CMAQ_CONC_file['ACLK'][:,cmaq_layer,row_index,col_index]
            ACLK = np.delete(ACLK , 0)

            ASO4K = read_from_daily_CMAQ_CONC_file['ASO4K'][:,cmaq_layer,row_index,col_index]
            ASO4K = np.delete(ASO4K , 0)

            ANO3K = read_from_daily_CMAQ_CONC_file['ANO3K'][:,cmaq_layer,row_index,col_index]
            ANO3K = np.delete(ANO3K , 0)

            ANH4K = read_from_daily_CMAQ_CONC_file['ANH4K'][:,cmaq_layer,row_index,col_index]
            ANH4K = np.delete(ANH4K , 0)

            ALVPO1I = read_from_daily_CMAQ_CONC_file['ALVPO1I'][:,cmaq_layer,row_index,col_index]
            ALVPO1I = np.delete(ALVPO1I , 0)

            ASVPO1I = read_from_daily_CMAQ_CONC_file['ASVPO1I'][:,cmaq_layer,row_index,col_index]
            ASVPO1I = np.delete(ASVPO1I , 0)

            ASVPO2I = read_from_daily_CMAQ_CONC_file['ASVPO2I'][:,cmaq_layer,row_index,col_index]
            ASVPO2I = np.delete(ASVPO2I , 0)

            ALVOO1I = read_from_daily_CMAQ_CONC_file['ALVOO1I'][:,cmaq_layer,row_index,col_index]
            ALVOO1I = np.delete(ALVOO1I , 0)

            ALVOO2I = read_from_daily_CMAQ_CONC_file['ALVOO2I'][:,cmaq_layer,row_index,col_index]
            ALVOO2I = np.delete(ALVOO2I , 0)

            ASVOO2I = read_from_daily_CMAQ_CONC_file['ASVOO2I'][:,cmaq_layer,row_index,col_index]
            ASVOO2I = np.delete(ASVOO2I , 0)

            ASVOO1I = read_from_daily_CMAQ_CONC_file['ASVOO1I'][:,cmaq_layer,row_index,col_index]
            ASVOO1I = np.delete(ASVOO1I , 0)

            ASO4J = read_from_daily_CMAQ_CONC_file['ASO4J'][:,cmaq_layer,row_index,col_index] # PM2.5 sulphate aerosol  PM2.5_SO4=ASO4J+ASO4I
            ASO4J = np.delete(ASO4J , 0)

            ASO4I = read_from_daily_CMAQ_CONC_file['ASO4I'][:,cmaq_layer,row_index,col_index] #
            ASO4I = np.delete(ASO4I , 0)

            ANH4J = read_from_daily_CMAQ_CONC_file['ANH4J'][:,cmaq_layer,row_index,col_index] # PM2.5 ammonium aerosol  PM2.5_NH4=ANH4J+ANH4I
            ANH4J = np.delete(ANH4J , 0)

            ANH4I = read_from_daily_CMAQ_CONC_file['ANH4I'][:,cmaq_layer,row_index,col_index] #
            ANH4I = np.delete(ANH4I , 0)

            ANO3J = read_from_daily_CMAQ_CONC_file['ANO3J'][:,cmaq_layer,row_index,col_index] # PM2.5 nitrate aerosol   PM2.5_NO3=ANO3J+ANO3I
            ANO3J = np.delete(ANO3J , 0)

            ANO3I = read_from_daily_CMAQ_CONC_file['ANO3I'][:,cmaq_layer,row_index,col_index] #
            ANO3I = np.delete(ANO3I , 0)

            AECJ = read_from_daily_CMAQ_CONC_file['AECJ'][:,cmaq_layer,row_index,col_index] # PM2.5 elemental carbon aerosol  PM2.5_EC=AECJ+AECI
            AECJ = np.delete(AECJ , 0)

            AECI = read_from_daily_CMAQ_CONC_file['AECI'][:,cmaq_layer,row_index,col_index] #
            AECI = np.delete(AECI , 0)

            ALVPO1J = read_from_daily_CMAQ_CONC_file['ALVPO1J'][:,cmaq_layer,row_index,col_index] #
            ALVPO1J = np.delete(ALVPO1J , 0)

            ASVPO1J = read_from_daily_CMAQ_CONC_file['ASVPO1J'][:,cmaq_layer,row_index,col_index] #
            ASVPO1J = np.delete(ASVPO1J , 0)

            ASVPO2J = read_from_daily_CMAQ_CONC_file['ASVPO2J'][:,cmaq_layer,row_index,col_index] #
            ASVPO2J = np.delete(ASVPO2J , 0)

            ASVPO3J = read_from_daily_CMAQ_CONC_file['ASVPO3J'][:,cmaq_layer,row_index,col_index] #
            ASVPO3J = np.delete(ASVPO3J , 0)

            AIVPO1J = read_from_daily_CMAQ_CONC_file['AIVPO1J'][:,cmaq_layer,row_index,col_index] #
            AIVPO1J = np.delete(AIVPO1J , 0)

            AXYL1J = read_from_daily_CMAQ_CONC_file['AXYL1J'][:,cmaq_layer,row_index,col_index] #
            AXYL1J = np.delete(AXYL1J , 0)

            AXYL2J = read_from_daily_CMAQ_CONC_file['AXYL2J'][:,cmaq_layer,row_index,col_index] #
            AXYL2J = np.delete(AXYL2J , 0)

            AXYL3J = read_from_daily_CMAQ_CONC_file['AXYL3J'][:,cmaq_layer,row_index,col_index] #
            AXYL3J = np.delete(AXYL3J , 0)

            ATOL1J = read_from_daily_CMAQ_CONC_file['ATOL1J'][:,cmaq_layer,row_index,col_index] #
            ATOL1J = np.delete(ATOL1J , 0)

            ATOL2J = read_from_daily_CMAQ_CONC_file['ATOL2J'][:,cmaq_layer,row_index,col_index] #
            ATOL2J = np.delete(ATOL2J , 0)

            ATOL3J = read_from_daily_CMAQ_CONC_file['ATOL3J'][:,cmaq_layer,row_index,col_index] #
            ATOL3J = np.delete(ATOL3J , 0)

            ABNZ1J = read_from_daily_CMAQ_CONC_file['ABNZ1J'][:,cmaq_layer,row_index,col_index] #
            ABNZ1J = np.delete(ABNZ1J , 0)

            ABNZ2J = read_from_daily_CMAQ_CONC_file['ABNZ2J'][:,cmaq_layer,row_index,col_index] #
            ABNZ2J = np.delete(ABNZ2J , 0)

            ABNZ3J = read_from_daily_CMAQ_CONC_file['ABNZ3J'][:,cmaq_layer,row_index,col_index] #
            ABNZ3J = np.delete(ABNZ3J , 0)

            AISO1J = read_from_daily_CMAQ_CONC_file['AISO1J'][:,cmaq_layer,row_index,col_index] #
            AISO1J = np.delete(AISO1J , 0)

            AISO2J = read_from_daily_CMAQ_CONC_file['AISO2J'][:,cmaq_layer,row_index,col_index] #
            AISO2J = np.delete(AISO2J , 0)

            ATRP1J = read_from_daily_CMAQ_CONC_file['ATRP1J'][:,cmaq_layer,row_index,col_index] #
            ATRP1J = np.delete(ATRP1J , 0)

            AISO3J = read_from_daily_CMAQ_CONC_file['AISO3J'][:,cmaq_layer,row_index,col_index] #
            AISO3J = np.delete(AISO3J , 0)

            ASQTJ = read_from_daily_CMAQ_CONC_file['ASQTJ'][:,cmaq_layer,row_index,col_index] #
            ASQTJ = np.delete(ASQTJ , 0)

            ATRP2J = read_from_daily_CMAQ_CONC_file['ATRP2J'][:,cmaq_layer,row_index,col_index] #
            ATRP2J = np.delete(ATRP2J , 0)

            AALK1J = read_from_daily_CMAQ_CONC_file['AALK1J'][:,cmaq_layer,row_index,col_index] #
            AALK1J = np.delete(AALK1J , 0)

            AALK2J = read_from_daily_CMAQ_CONC_file['AALK2J'][:,cmaq_layer,row_index,col_index] #
            AALK2J = np.delete(AALK2J , 0)

            APAH1J = read_from_daily_CMAQ_CONC_file['APAH1J'][:,cmaq_layer,row_index,col_index] #
            APAH1J = np.delete(APAH1J , 0)

            APAH2J = read_from_daily_CMAQ_CONC_file['APAH2J'][:,cmaq_layer,row_index,col_index] #
            APAH2J = np.delete(APAH2J , 0)

            APAH3J = read_from_daily_CMAQ_CONC_file['APAH3J'][:,cmaq_layer,row_index,col_index] #
            APAH3J = np.delete(APAH3J , 0)

            AOLGBJ = read_from_daily_CMAQ_CONC_file['AOLGBJ'][:,cmaq_layer,row_index,col_index] #
            AOLGBJ = np.delete(AOLGBJ , 0)

            AORGCJ = read_from_daily_CMAQ_CONC_file['AORGCJ'][:,cmaq_layer,row_index,col_index] #
            AORGCJ = np.delete(AORGCJ , 0)

            AOLGAJ = read_from_daily_CMAQ_CONC_file['AOLGAJ'][:,cmaq_layer,row_index,col_index] #
            AOLGAJ = np.delete(AOLGAJ , 0)

            ALVOO1J = read_from_daily_CMAQ_CONC_file['ALVOO1J'][:,cmaq_layer,row_index,col_index] #
            ALVOO1J = np.delete(ALVOO1J , 0)

            ALVOO2J = read_from_daily_CMAQ_CONC_file['ALVOO2J'][:,cmaq_layer,row_index,col_index] #
            ALVOO2J = np.delete(ALVOO2J , 0)

            ASVOO1J = read_from_daily_CMAQ_CONC_file['ASVOO1J'][:,cmaq_layer,row_index,col_index] #
            ASVOO1J = np.delete(ASVOO1J , 0)

            ASVOO2J = read_from_daily_CMAQ_CONC_file['ASVOO2J'][:,cmaq_layer,row_index,col_index] #
            ASVOO2J = np.delete(ASVOO2J , 0)

            APCSOJ = read_from_daily_CMAQ_CONC_file['APCSOJ'][:,cmaq_layer,row_index,col_index] #
            APCSOJ = np.delete(APCSOJ , 0)

            ASVOO3J = read_from_daily_CMAQ_CONC_file['ASVOO3J'][:,cmaq_layer,row_index,col_index] #
            ASVOO3J = np.delete(ASVOO3J , 0)




#            CONC_file_species_list = [AKJ]# , AMNJ , AMGJ , ACAJ , ATIJ , ASIJ , AFEJ , AOTHRJ , ACLJ ,\
#                                      ANAJ , AOTHRI , ACLI , ANAI , ASVOO3J , APCSOJ , ASO4J , ASO4I,\
#                                      ANH4J , ANH4I , ANO3J , ANO3I , AECJ , AECI , AALK1J , AALK2J ,\
#                                      AXYL1J , AXYL2J , AXYL3J , ATOL1J , ATOL2J , ATOL3J , ABNZ1J ,\
#                                      ABNZ2J , ABNZ3J , AOLGAJ , AORGCJ , ATRP1J , ATRP2J , AISO1J , \
#                                      AISO2J , AISO3J , ASQTJ , AOLGBJ , ASVOO2J , ASVOO1J , ASVOO1I ,\
#                                      ALVOO2J , ALVOO1J , APAH3J , APAH2J , ASVOO2I , ALVOO2I , ALVOO1I ,\
#                                      APAH1J , AIVPO1J , ASVPO3J , ASVPO2J , ASVPO1J , ALVPO1J , ASVPO2I ,\
#                                      ASVPO1I , ALVPO1I , ACLK , ASO4K , ANO3K , ANH4K , AALJ , ASEACAT ,\
#                                      ACORS , ASOIL]

#            # delete the first row of CONC species to be in shape with PMDIAG files; 25:0-0 == 24:1-0
#            for ispecies in CONC_file_species_list:
#                print('-> delete the first row for %s' %ispecies)
#                ispecies = np.delete(ispecies , 0)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # for PM10 we should read PMDIAG file
            if VAR == 'PM10':

                print('-> now opening/reading PMDIAG ncfile...')

                # use the same loop (CONC files) for PM10
                cmaq_pmdiag_file_name = cmaq_pmdiag_file_list[cmaq_conc_file_no]
                cmaq_pmdiag_file_full_path = os.path.join(cmaq_input_dir,cmaq_pmdiag_file_name)

                print('-> check to see if CMAQ PMDIAG file exists...')
                if (os.path.isfile(cmaq_pmdiag_file_full_path)!=True):
                    print('-> ERROR: CMAQ PMDIAG input file does NOT exist at path: "%s"' %cmaq_input_dir)
                else:
                    print('-> NOTE: CMAQ PMDIAG input file exists!')
                    print('-> file name is: "%s"'%cmaq_pmdiag_file_name)


                    # open and read each cmaq PMDIAG file seperately; should use the full file path in Dateset;
                    cmaq_read_pmdiag_file_i = Dataset(cmaq_pmdiag_file_full_path,'r')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    # loop for PMDIAG date-times: TFLAGs

                    print('-> date-time loop over PMDIAG file no. %s...' %(cmaq_conc_file_no+1))
                    cmaq_pmdiag_tflag_i = cmaq_read_pmdiag_file_i.variables['TFLAG']
                    # extract date-time for one day/loop
                    cmaq_pmdiag_datetimes_i = cmaq_pmdiag_tflag_i[:,0,:] # for only one VAR is enough='0'; only for TSTEP and DATE-TIME

                    # loop for 24-hour
                    for hr in range(len(cmaq_pmdiag_datetimes_i)):
                        # calculate date:
                        date_pmdiag_i = cmaq_pmdiag_datetimes_i[hr,0]
                        date_pmdiag_i_obj = dt.datetime.strptime(str(date_pmdiag_i),'%Y%j').date() # to get the date
                        #print('date is %s'%date_pmdiag_i)

                        # calculate time:
                        time_pmdiag_i = cmaq_pmdiag_datetimes_i[hr,1]
                        if time_pmdiag_i < 100000:
                            if time_pmdiag_i == 0 :
                                time_pmdiag_i_str = str(time_pmdiag_i).rjust(6,'0') # final width=6, pad with '0''
                                #time_pmdiag_i_str = time_pmdiag_i_str+'0000'
                            else:
                                time_pmdiag_i_str = str(time_pmdiag_i).rjust(6,'0') # final width=6, pad with '0''
                        elif time_pmdiag_i >= 100000:
                            time_pmdiag_i_str = str(time_pmdiag_i)


                        #print('time is %s'%time_pmdiag_i)
                        # make dt obj from str
                        time_pmdiag_i_obj = dt.datetime.strptime(time_pmdiag_i_str,'%H%M%S').time() # just get the time

                        # combine date-obj and time-obj together
                        dt_pmdiag_obj_hr = dt.datetime.combine(date_pmdiag_i_obj,time_pmdiag_i_obj)
                        # append each dt_pmdiag_obj_hr to the total list

                        # now appending PMDIAG daily date-time object to its total list...')
                        dt_CMAQ_PMDIAG_LocalTime_TotalList.append(dt_pmdiag_obj_hr - dt.timedelta(hours=utc_conversion))
#                        mpdiag_list_stack = [ , ]
#                        np.hstack(dt_pmdiag_obj_hr - dt.timedelta(hours=utc_conversion))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    print('-> now extracting PM10 species...')
                    PM10AT_3 = cmaq_read_pmdiag_file_i['PM10AT'][:,cmaq_layer,row_index,col_index]
                    PM10AC_3 = cmaq_read_pmdiag_file_i['PM10AC'][:,cmaq_layer,row_index,col_index]
                    PM10CO_3 = cmaq_read_pmdiag_file_i['PM10CO'][:,cmaq_layer,row_index,col_index]

#                    PM10AT_3_list.extend( PM10AT_3 )
#                    PM10AC_3_list.extend( PM10AC_3 )
#                    PM10CO_3_list.extend( PM10CO_3 )


                   # PM10 definition from EPA (!! PM10.0 and Coarse-Sized Species; "_0" means VAR is defined in the file
                    print('-> now calculate PM10-species concentrations for each day from CONC and PMDIAG species...')

                    APOMI_0 = ALVPO1I + ASVPO1I + ASVPO2I
                    ASOMI_0 = ALVOO1I + ALVOO2I + ASVOO1I + ASVOO2I

                    AOMI_0 = APOMI_0 + ASOMI_0



                    APOMJ_0 = ALVPO1J + ASVPO1J + ASVPO2J +ASVPO3J + AIVPO1J
                    ASOMJ_0 = AXYL1J  + AXYL2J  + AXYL3J  + ATOL1J    \
                               +ATOL2J  + ATOL3J  + ABNZ1J  + ABNZ2J  \
                               +ABNZ3J  + AISO1J  + AISO2J  + AISO3J  \
                               +ATRP1J  + ATRP2J  + ASQTJ   + AALK1J  \
                               +AALK2J  + APAH1J  + APAH2J  + APAH3J  \
                               +AORGCJ  + AOLGBJ  + AOLGAJ            \
                               +ALVOO1J + ALVOO2J + ASVOO1J + ASVOO2J \
                               +ASVOO3J + APCSOJ

                    AOMJ_0 = APOMJ_0 + ASOMJ_0



                    ATOTI_0 = ASO4I + ANO3I + ANH4I + ANAI + ACLI + AECI + AOMI_0 + AOTHRI
                    ATOTJ_0 = ASO4J + ANO3J + ANH4J + ANAJ + ACLJ + AECJ + AOMJ_0 + AOTHRJ \
                            + AFEJ + ASIJ + ATIJ + ACAJ + AMGJ + AMNJ + AALJ + AKJ
                    ATOTK_0 = ASOIL + ACORS + ASEACAT + ACLK + ASO4K + ANO3K + ANH4K

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    print('-> now calculating PM10 daily timeseries...') # check to see if we have to sum columns or this operation does that for us.
                    PM10_daily_timeseries = (ATOTI_0 * PM10AT_3) + (ATOTJ_0 * PM10AC_3) + (ATOTK_0 * PM10CO_3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    VAR_daily_TimeSeries_extract = PM10_daily_timeseries
                    CMAQ_VAR_TimeSeries_4WholePeriod_List.extend(VAR_daily_TimeSeries_extract)  #???

                    print('-> closing PMDIAG ncfile ...')
                    cmaq_read_pmdiag_file_i.close()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # for PM2.5
            else:
                # for PM2.5; new version from EPA')
                   #pm25_ot_0 = ATOTI_0 * PM25AT_3 + ATOTJ_0 * PM25AC_3 + ATOTK_0 * PM25CO_3

                # old method from John Mejia

                print('-> calculating PM2.5 species...')
                AALK1J = read_from_daily_CMAQ_CONC_file['AALK1J'][:,cmaq_layer,row_index,col_index] # PM2.5 anthropogenic SOA PM2.5_SOAA=AORGAJ+AORGAI or PM2.5_SOAA=AALK1J+AALK2J+AXYL1J+AXYL2J+AXYL3J+ATOL1J+ATOL2J+ATOL3J+ABNZ1J+ABNZ2J+ABNZ3J+AOLGAJ+AORGCJ/2.
                AALK1J = AALK1J*0.5
                AALK1J = np.delete(AALK1J , 0)

                AALK2J = read_from_daily_CMAQ_CONC_file['AALK2J'][:,cmaq_layer,row_index,col_index] #
                AALK2J = AALK2J*0.5
                AALK2J = np.delete(AALK2J , 0)

                AXYL1J = read_from_daily_CMAQ_CONC_file['AXYL1J'][:,cmaq_layer,row_index,col_index] #
                AXYL1J = AXYL1J*0.5
                AXYL1J = np.delete(AXYL1J , 0)

                AXYL2J = read_from_daily_CMAQ_CONC_file['AXYL2J'][:,cmaq_layer,row_index,col_index] #
                AXYL2J = AXYL2J*0.5
                AXYL2J = np.delete(AXYL2J , 0)

                AXYL3J = read_from_daily_CMAQ_CONC_file['AXYL3J'][:,cmaq_layer,row_index,col_index] #
                AXYL3J = AXYL3J*0.5
                AXYL3J = np.delete(AXYL3J , 0)

                ATOL1J = read_from_daily_CMAQ_CONC_file['ATOL1J'][:,cmaq_layer,row_index,col_index] #
                ATOL1J = ATOL1J*0.5
                ATOL1J = np.delete(ATOL1J , 0)

                ATOL2J = read_from_daily_CMAQ_CONC_file['ATOL2J'][:,cmaq_layer,row_index,col_index] #
                ATOL2J = ATOL2J*0.5
                ATOL2J = np.delete(ATOL2J , 0)

                ATOL3J = read_from_daily_CMAQ_CONC_file['ATOL3J'][:,cmaq_layer,row_index,col_index] #
                ATOL3J = ATOL3J*0.5
                ATOL3J = np.delete(ATOL3J , 0)

                ABNZ1J = read_from_daily_CMAQ_CONC_file['ABNZ1J'][:,cmaq_layer,row_index,col_index] #
                ABNZ1J = ABNZ1J*0.5
                ABNZ1J = np.delete(ABNZ1J , 0)

                ABNZ2J = read_from_daily_CMAQ_CONC_file['ABNZ2J'][:,cmaq_layer,row_index,col_index] #
                ABNZ2J = ABNZ2J*0.5
                ABNZ2J = np.delete(ABNZ2J , 0)

                ABNZ3J = read_from_daily_CMAQ_CONC_file['ABNZ3J'][:,cmaq_layer,row_index,col_index] #
                ABNZ3J = ABNZ3J*0.5
                ABNZ3J = np.delete(ABNZ3J , 0)

                AOLGAJ = read_from_daily_CMAQ_CONC_file['AOLGAJ'][:,cmaq_layer,row_index,col_index] #
                AOLGAJ = AOLGAJ*0.5
                AOLGAJ = np.delete(AOLGAJ , 0)

                AORGCJ = read_from_daily_CMAQ_CONC_file['AORGCJ'][:,cmaq_layer,row_index,col_index] #
                AORGCJ = AORGCJ*0.5
                AORGCJ = np.delete(AORGCJ , 0)

                ATRP1J = read_from_daily_CMAQ_CONC_file['ATRP1J'][:,cmaq_layer,row_index,col_index] # PM2.5 biogenic SOA  PM2.5_SOAB=AORGBJ+AORGBI or PM2.5_SOAB=ATRP1J+ATRP2J+AISO1J+AISO2J+AISO3J+ASQTJ+AOLGBJ+AORGCJ/2.
                ATRP1J = ATRP1J*0.5
                ATRP1J = np.delete(ATRP1J , 0)

                ATRP2J = read_from_daily_CMAQ_CONC_file['ATRP2J'][:,cmaq_layer,row_index,col_index] #
                ATRP2J = ATRP2J*0.5
                ATRP2J = np.delete(ATRP2J , 0)

                AISO1J = read_from_daily_CMAQ_CONC_file['AISO1J'][:,cmaq_layer,row_index,col_index] #
                AISO1J = AISO1J*0.5
                AISO1J = np.delete(AISO1J , 0)

                AISO2J = read_from_daily_CMAQ_CONC_file['AISO2J'][:,cmaq_layer,row_index,col_index] #
                AISO2J = AISO2J*0.5
                AISO2J = np.delete(AISO2J , 0)

                AISO3J = read_from_daily_CMAQ_CONC_file['AISO3J'][:,cmaq_layer,row_index,col_index] #
                AISO3J = AISO3J*0.5
                AISO3J = np.delete(AISO3J , 0)

                ASQTJ = read_from_daily_CMAQ_CONC_file['ASQTJ'][:,cmaq_layer,row_index,col_index] #
                ASQTJ = ASQTJ*0.5
                ASQTJ = np.delete(ASQTJ , 0)

                AOLGBJ = read_from_daily_CMAQ_CONC_file['AOLGBJ'][:,cmaq_layer,row_index,col_index] #
                AOLGBJ = AOLGBJ*0.5
                AOLGBJ = np.delete(AOLGBJ , 0)

                AORGCJ = read_from_daily_CMAQ_CONC_file['AORGCJ'][:,cmaq_layer,row_index,col_index] #
                AORGCJ = AORGCJ*0.5
                AORGCJ = np.delete(AORGCJ , 0)




                PM25_daily_timeseries = ASO4J + ASO4I + ANH4J + ANH4I + ANO3J + ANO3I + AECJ    \
                + AECI + AALK1J + AALK2J + AXYL1J + AXYL2J + AXYL3J + ATOL1J + ATOL2J + ATOL3J  \
                + ABNZ1J + ABNZ2J + ABNZ3J + AOLGAJ + AORGCJ + ATRP1J + ATRP2J + AISO1J + AISO2J\
                + AISO3J + ASQTJ + AOLGBJ + AORGCJ

                VAR_daily_TimeSeries_extract = PM25_daily_timeseries
                CMAQ_VAR_TimeSeries_4WholePeriod_List.extend(VAR_daily_TimeSeries_extract)


                # define species lists that join/stand together

#                PM2_5_species_list = [ASO4J,ASO4I,ANH4J,ANH4I,ANO3J,ANO3I,AECJ,AECI,AALK1J,AALK2J,AXYL1J,AXYL2J,AXYL3J,ATOL1J,ATOL2J,ATOL3J,ABNZ1J,ABNZ2J,ABNZ3J,AOLGAJ,AORGCJ,ATRP1J,ATRP2J,AISO1J,AISO2J,AISO3J,ASQTJ,AOLGBJ,AORGCJ]
#
#                # vertically stack all species lists to sum them vertically for each hour, cos hourly PM2.5=sum of all species at each hour
#                PM2_5_species_array = np.vstack(PM2_5_species_list )#  vstack (vertically stacks lists) vs. hstack (series lists after each other)
#
#                # sum vertically over each column of total array, which is each hour
#                PM25_daily_timeseries = PM2_5_species_array.sum(axis=0) # axis=0 => vertically
#
#                # extend the total list, extend adds elements
#                CMAQ_VAR_TimeSeries_4WholePeriod_List.extend(PM25_daily_timeseries)
#                read_from_daily_CMAQ_CONC_file.close()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#        elif VAR in ('NO2','SO2','O3','CO'):
        elif VAR=='NO2' or VAR=='SO2' or VAR=='O3' or VAR=='CO':

            print('-> now extract VAR=%s from CMAQ station cell for file no. %s' %(VAR,(cmaq_conc_file_no+1)))
            print(read_from_daily_CMAQ_CONC_file[VAR])

            print('-> dimension of VAR=%s is:' %VAR)
            print(read_from_daily_CMAQ_CONC_file[VAR].dimensions)
            print('-> shape of VAR=%s is:' %VAR)
            print(read_from_daily_CMAQ_CONC_file[VAR].shape)

            #ozone = cmaq_read_conc_file_i['O3']  # ozone is still an object in this line;
            #ozone = cmaq_read_conc_file_i['O3'][:]  # Read the entire array into memory (this is what "[:]" is doing).

            # now read CMAQ VAR
            VAR_daily_TimeSeries_extract_25hr = read_from_daily_CMAQ_CONC_file[VAR][:,cmaq_layer,row_index,col_index] # slice a time series for a cell

            # now delete the 1st hr=0, so the list starts from 1-24
            VAR_daily_TimeSeries_extract = np.delete(VAR_daily_TimeSeries_extract_25hr , 0)

            if VAR in ('NO2','SO2'):
                VAR_daily_TimeSeries_extract = VAR_daily_TimeSeries_extract*1000 # change ppmv from CMAQ to ppb for NO2, SO2; O3 is in PPM

            # now extend/update daily list to total period list=CMAQ_VAR_TimeSeries_4WholePeriod_List
            CMAQ_VAR_TimeSeries_4WholePeriod_List = np.concatenate([CMAQ_VAR_TimeSeries_4WholePeriod_List , VAR_daily_TimeSeries_extract],axis=0)

        else: # O3 and CO are in PPM
            print('-> ERROR: VAR= %s NOT inside default VAR list' %VAR)

        print('-> closing CMAQ input file no.:(%s)' %(cmaq_conc_file_no+1))
        print('**************************************************************************')
        read_from_daily_CMAQ_CONC_file.close()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('-> CHECKING CMAQ files to see if they start from 0-24 or 1-24...')
    print('-> CMAQ CONC date-time file starts from:   %s and has %s elements!' %( dt_CMAQ_CONC_LocalTime_TotalList[0] , len(dt_CMAQ_CONC_LocalTime_TotalList) ))

    if VAR == 'PM10':
        print('-> CMAQ PMDIAG file, date-time starts from: %s and has %s elements!' %( dt_CMAQ_PMDIAG_LocalTime_TotalList[0] , len(dt_CMAQ_PMDIAG_LocalTime_TotalList) ))

    print('-> CMAQ VAR="%s" time-series list has ................................ %s elements!' %(VAR,len(CMAQ_VAR_TimeSeries_4WholePeriod_List)))
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # make CMAQ DF for plotting durnal cycle
    cmaq_df = pd.DataFrame({'cmaq_datetimes':dt_CMAQ_CONC_LocalTime_TotalList ,\
                            'cmaq_sim_values':CMAQ_VAR_TimeSeries_4WholePeriod_List})

    #cmaq_df.drop_duplicates('cmaq_datetimes' , inplace=True) # removes all rows with dup. data; each row all columns together

    # set date-time col. as index
    cmaq_df.index = cmaq_df['cmaq_datetimes'] # for plotting against date-time

    # define a filter for a favorite period; for error statistics
    filter_cmaq_4period = (cmaq_df.index >= dt.datetime(2014,8,10,0,0,0)) & (cmaq_df.index <= dt.datetime(2014,8,15,17,0,0))

    # only filter col.: cmaq_sim_values for our modeling period; now we have a series
    cmaq_filtered_4ModelingPeriod = cmaq_df.cmaq_sim_values[filter_cmaq_4period]

    print('-> NOTE: Filtered period for CMAQ time-series is from: <%s> to: <%s>' %(cmaq_filtered_4ModelingPeriod.index[0] , cmaq_filtered_4ModelingPeriod.index[-1]))

    # group by each hr to calculate diurnal cycle
    cmaq_mean_diurnal_timeseries = cmaq_filtered_4ModelingPeriod.groupby(cmaq_filtered_4ModelingPeriod.index.hour).mean() # based on its index=hour


# read CMAQ files
##########################################################################################


##########################################################################################
# read Observation file


# set pathes based on work directory:
obs_input_dir = '/Users/ehsan/Documents/my_Python_codes/AQ_stations/inputs'
obs_file_name = stn_name[stn_index]+'_'+StateCode[stn_index]+'_'+CountyCode[stn_index]+'_'+site_no[stn_index]+'.txt'
obs_file_full_path = os.path.join(obs_input_dir,obs_file_name)

print('-> check to see if OBS file exists...')
if (os.path.isfile(obs_file_full_path)!=True):
    print('-> ERROR: OBS input file does NOT exist at path: "%s"' %obs_file_full_path)
else:
    print('-> NOTE: OBS input file exists!')

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('-> OBS file full path is: %s' %obs_file_full_path)
    print('-> opening observation input file ...')
    obs_df = pd.read_csv(obs_file_full_path , sep=',' , dtype={'State Code':str ,
                                                         'County Code':str ,
                                                         'Site Num':str ,
                                                         'Parameter Code':str} ) # change dtype to str=obj to compare with str.

    # change column labels
    obs_df.rename(columns={'Date Local':'date_local' , '24 Hour Local':'time_local',
    'Parameter Code':'param_code', 'Sample Measurement':'obs_pol'} , inplace=True) # change obs_pol to selected pollutant name

#    print('-> head of obs_df before sort is = ' )
#    print(obs_df.head())

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # join date and time columns
    obs_df['dt_str_local'] = obs_df.date_local+' '+obs_df.time_local

#    print('-> head of OBS col: dt_str_local after joining date and time col. is:')
#    print(obs_df.dt_str_local.head(10))
#    print('-> type of elements inside col.:%s is: %s' %('dt_str_local' , obs_df.dt_str_local.dtype))

    # convert dt_str to dt_obj and add to new col.
    #obs_df['dt_obj_local'] = dt.datetime.strptime(obs_df.dt_local_str , '%Y-%m%d %H%m') # for single str, format should match dt_local_str format.
    obs_df['dt_obj_local'] = pd.to_datetime(obs_df.dt_str_local) # converts a col. to dt_obj


    # define a filter based on a pollutant=ParamCode
    filter_ParamCode = (obs_df.param_code == obs_PolCode[pol_index]) # define paramcode somewhere up here
    # apply the filter and get the filtered slice for our pollutant= obs_PolCode
    obs_df_filtered = obs_df[filter_ParamCode]

    print('-> check to see if "obs_df_filtered" has "%s" data or it is empty...'%pol_name[pol_index])
    if obs_df_filtered.obs_pol.empty == True:
        print('-> ERROR: OBS file is empty, VAR: "%s" NOT exist in the file!' %pol_name[pol_index])
    else:
        print('-> NOTE: OBS file includes pollutant: "%s"' %pol_name[pol_index])

#        print('-> head of filtered DF based on favorite VAR at "obs_df_filtered" is:')
#        print(obs_df_filtered.head())

        # sort the filtered slice based on dt_obj col.; ascending; creates new DF
        # sort the dt column and make a new DF; does sort_values make a copy of a filetered DF? check it.
        obs_sorted_df = obs_df_filtered.sort_values(by=['dt_obj_local']) # sort all DF based on "dt_obj_local" col.

#        print('-> tail of "obs_sorted_df" is:')
#        print(obs_sorted_df.tail())

        # set date-time col. as index of DF
        obs_sorted_df.index = obs_sorted_df['dt_obj_local']

        # define a filter for a favorite period; for error statistics
        filter_OBS_4period = (obs_sorted_df.index >= dt.datetime(2014,8,10,0,0,0)) & (obs_sorted_df.index <= dt.datetime(2014,8,15,17,0,0))

        # we only filter and get col.: obs_pol for our modeling period; now we have a series
        OBS_filtered_4ModelingPeriod = obs_sorted_df.obs_pol[filter_OBS_4period]

        print('-> NOTE: Filtered period for OBS time-series is from: <%s> to: <%s>' %(OBS_filtered_4ModelingPeriod.index[0] , OBS_filtered_4ModelingPeriod.index[-1]))

        # group by each hr to calculate diurnal cycle
        obs_mean_diurnal_timeseries = OBS_filtered_4ModelingPeriod.groupby(OBS_filtered_4ModelingPeriod.index.hour).mean() # based on its index=hour

# read Observation file
##########################################################################################


##########################################################################################
# statistical analysis

if (os.path.isfile(mcip_file_full_path)!=True):
    print('-> ERROR: MCIP file missing!')
elif (os.path.isfile(cmaq_CONC_file_full_path)!=True):
    print('-> ERROR: CMAQ file missing!')

elif (os.path.isfile(obs_file_full_path)!=True):
    print('-> ERROR: OBS file missing!')

elif obs_df_filtered.obs_pol.empty == True:
    print('-> ERROR: missing VAR in OBS file!')
else:
    print('-> start statistics and plotting data ...')


    print('-> Check if CMAQ and OBS have same data size:')
#    filter_obs_avail = cmaq_filtered_4ModelingPeriod.index.isin(OBS_filtered_4ModelingPeriod.index)
#    missing_obs = OBS_filtered_4ModelingPeriod[filter_obs_avail == False]
#    print(missing_obs)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    len_list = [len(cmaq_filtered_4ModelingPeriod) , len(OBS_filtered_4ModelingPeriod)]
    #print('-> len_list elements are: %s' %len_list)
    print('-> CMAQ has %s elements and OBS has %s elements!' %(len_list[0] , len_list[1]))

    if len_list[0] == len_list[1]:
        print('-> Note: both CMAQ and Obs time series have same size!')
    else:
        print('-> Note: CMAQ and Obs time series do NOT have same size, so we have to remove lines that are not available in other list!')

    loop_len = min(len_list)

    for itry in range(loop_len):
        if cmaq_filtered_4ModelingPeriod.index[itry] == OBS_filtered_4ModelingPeriod.index[itry]:
            print('-> index matches for %s' %cmaq_filtered_4ModelingPeriod.index[itry])
        else:
            while cmaq_filtered_4ModelingPeriod.index[itry] != OBS_filtered_4ModelingPeriod.index[itry]:
                cmaq_index = cmaq_filtered_4ModelingPeriod.index[itry]
                obs_index = OBS_filtered_4ModelingPeriod.index[itry]
                #print('-> cmaq index is %s, but OBS index is %s' %(cmaq_index , obs_index))

                if cmaq_index > obs_index:
                    print('-> %s is missing in CMAQ. file' %obs_index)#  strftime(obs_index))# , %Y-%m-%d %H:%M:%S))
                    del(OBS_filtered_4ModelingPeriod[OBS_filtered_4ModelingPeriod.index == obs_index])
                elif cmaq_index < obs_index:
                    print('-> %s is missing in OBS file' %cmaq_index) # strftime(cmaq_index , %Y-%m-%d %H:%M:%S)))
                    del(cmaq_filtered_4ModelingPeriod[cmaq_index])
                else:
                    print('-> something wierd is happening')
#                cmaq_filtered_4ModelingPeriod.index[itry] == OBS_filtered_4ModelingPeriod.index[itry]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print('*******************************************************************************')
    print('-> Statistics for VAR: "%s" at station: "%s":' %(VAR,stn_name[stn_index]))

    print('==> OBS max for VAR: "%s" during the modeling period is: %s %s' %(VAR,OBS_filtered_4ModelingPeriod.max(),VAR_unit))
    print('==> OBS min for VAR: "%s" during the modeling period is: %s %s' %(VAR,OBS_filtered_4ModelingPeriod.min(),VAR_unit))
    print('==> CMAQ max for VAR: "%s" during the modeling period is: %s %s' %(VAR,cmaq_filtered_4ModelingPeriod.max(),VAR_unit))
    print('==> CMAQ min for VAR: "%s" during the modeling period is: %s %s' %(VAR,cmaq_filtered_4ModelingPeriod.min(),VAR_unit))



    cmaq_bias = cmaq_filtered_4ModelingPeriod - OBS_filtered_4ModelingPeriod
    print('==> CMAQ mean BIAS (MBE) is: %s %s' %(cmaq_bias.mean(),VAR_unit))

    cmaq_rmse = ((cmaq_bias**2).mean())**0.5
    print('==> CMAQ RMSE is: %s %s' %(cmaq_rmse,VAR_unit))

    cmaq_R2 = np.corrcoef(cmaq_filtered_4ModelingPeriod , OBS_filtered_4ModelingPeriod , rowvar=False)[0,1]**2 # rowvar=False => col.=VAR; rows=obs
    print('==> CMAQ R2 is: %s' %cmaq_R2)

    print('*******************************************************************************')
    print('-> DONE - run duration of processing section is: %s seconds' %(time.time() - start_time))

# statistical analysis
##########################################################################################


##########################################################################################
# plotting- mean diurnal


#    # set min and max values based on data range
#    y_min = 0
#    y_max_obs = obs_mean_diurnal_timeseries.max()
#    y_max = y_max_obs*2.2 + y_max_obs
#    y_int = (y_max-y_min)/10
#
##    [y_min , y_max , y_int] = [0, 100, 5]           # change this range based on pollutant; O3:0-0.1; CO:0-1; NO2:0-50; SO2:0-5; PM:0-100
#    y_ticks = np.arange(y_min, y_max+y_int, y_int)  # returns a range with float range
#    y_units_label = VAR_unit
#
#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#    fig_num = 101
#    fig = plt.figure(num=fig_num,figsize=(10,5))
#    plt.clf()
#
#    #   Observation data
#    obs_dt =      obs_mean_diurnal_timeseries.index           #obs_sorted_df.dt_obj_local
#    obs_values =  obs_mean_diurnal_timeseries      #obs_sorted_df.obs_pol # param_code is the filtered col. based on my defined filter
#
#    #   CMAQ data
#    cmaq_dt =     cmaq_mean_diurnal_timeseries.index #dt_CMAQ_CONC_LocalTime_TotalList
#    cmaq_values =        cmaq_mean_diurnal_timeseries   #CMAQ_VAR_TimeSeries_4WholePeriod_List
#
#    #   plot for Observation
#    plt.plot(obs_dt, obs_values, 'k', linestyle='-', label='Obs', linewidth=2.0, markersize=6, markeredgecolor='k')
#
#    #   plot for CMAQ
#    plt.plot(cmaq_dt, cmaq_values, 'r', linestyle='-', label='CMAQ', linewidth=2.0, markersize=6, markeredgecolor='k')
#
#    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#    plt.legend(loc=2,fontsize=12,ncol=1)
##    plt.title('CMAQ vs OBS at "%s" station - Aug. 2014' %( stn_name[stn_index]), \
##         fontsize=12, x=0.5, y=1.01)
#    plt.title('24-hr diurnal cycle of CMAQ vs OBS at "%s" station - Aug. 2014' %( stn_name[stn_index]), \
#         fontsize=12, x=0.5, y=1.01)
#
#    plt.ylabel('%s [%s]' %(VAR,y_units_label),\
#               fontsize=12,labelpad=20)
#
##    plt.yticks(y_ticks)
#    plt.ylim([y_min, y_max])
#    plt.xlabel('hour [PST]',fontsize=12,labelpad=00)
##    plt.xlim([datetime_ticks[0], datetime_ticks[-1]])
##    plt.gca().xaxis.set_major_formatter(DateFormatter(datetick_format))
##    plt.xticks(datetime_ticks,visible=True)
#    plt.show()
#
#
#    # save the plot
#    #filename = 'CMAQ_timeseries_comparison_'+VAR+'_'+stn_name[stn_index]+'.png'
#    filename = 'mean_durnal_cycle_'+VAR+'_'+stn_name[stn_index]+'.png'
#    plot_name = os.path.join(figure_dir,filename)
#    plt.savefig(plot_name)
#

# plotting- mean diurnal
##########################################################################################


##########################################################################################
# plotting- time series of modeling period

    # set min and max values based on data range
    y_min = 0
    y_max_obs  = OBS_filtered_4ModelingPeriod.max()
    y_max_cmaq = cmaq_filtered_4ModelingPeriod.max()
    y_limit = max([y_max_obs,y_max_cmaq])

    y_max = y_limit + 0.2*y_limit
    y_int = (y_max-y_min)/10

#    [y_min , y_max , y_int] = [0, 100, 5]           # change this range based on pollutant; O3:0-0.1; CO:0-1; NO2:0-50; SO2:0-5; PM:0-100
    y_ticks = np.arange(y_min, y_max+y_int, y_int)  # returns a range with float range
    y_units_label = VAR_unit

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    fig_num = 101
    fig = plt.figure(num=fig_num,figsize=(10,5))
    plt.clf()

    #   Observation data
    obs_dt     =  OBS_filtered_4ModelingPeriod.index           #obs_sorted_df.dt_obj_local
    obs_values =  OBS_filtered_4ModelingPeriod                 #obs_sorted_df.obs_pol # param_code is the filtered col. based on my defined filter

    #   CMAQ data
    cmaq_dt     = cmaq_filtered_4ModelingPeriod.index   #dt_CMAQ_CONC_LocalTime_TotalList
    cmaq_values = cmaq_filtered_4ModelingPeriod         #CMAQ_VAR_TimeSeries_4WholePeriod_List

    #   plot for Observation
    plt.plot(obs_dt, obs_values, 'k', linestyle='-', label='Obs', linewidth=2.0, markersize=6, markeredgecolor='k')

    #   plot for CMAQ
    plt.plot(cmaq_dt, cmaq_values, 'r', linestyle='-', label='CMAQ', linewidth=2.0, markersize=6, markeredgecolor='k')

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    plt.legend(loc=2,fontsize=12,ncol=1)
#    plt.title('CMAQ vs OBS at "%s" station - Aug. 2014' %( stn_name[stn_index]), \
#         fontsize=12, x=0.5, y=1.01)
    plt.title('time series of CMAQ vs OBS at "%s" station - Aug. 2014' %( stn_name[stn_index]), \
         fontsize=12, x=0.5, y=1.01)

    plt.ylabel('%s [%s]' %(VAR,y_units_label),\
               fontsize=12,labelpad=20)

#    plt.yticks(y_ticks)
    plt.ylim([y_min, y_max])
    plt.xlabel('hour [PST]',fontsize=12,labelpad=00)
#    plt.xlim([datetime_ticks[0], datetime_ticks[-1]])
#    plt.gca().xaxis.set_major_formatter(DateFormatter(datetick_format))
#    plt.xticks(datetime_ticks,visible=True)
    plt.show()


    # save the plot
    #filename = 'CMAQ_timeseries_comparison_'+VAR+'_'+stn_name[stn_index]+'.png'
    filename = 'timeseries_'+VAR+'_'+stn_name[stn_index]+'.png'
    plot_name = os.path.join(figure_dir,filename)
    plt.savefig(plot_name)

print('-> run duration is: %s seconds' %(time.time() - start_time))

# plotting- time series of modeling period
##########################################################################################