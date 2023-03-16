import  subprocess
import os
import shutil


os.chdir(r'C:\Users\marti\Dropbox\research_projects\air_pollution')




# china 
#echo  27.0 112.0 10.0         >>CONTROL
#echo  32.0 117.0 10.0         >>CONTROL
#echo  27.5 112.5 10.0         >>CONTROL

# pakistan
#echo  28.0 70.0 10.0         >>CONTROL
#echo  33.0 75.0 10.0         >>CONTROL
#echo  28.5 70.5 10.0         >>CONTROL




# see here for the description of the parameters
# https://www.ready.noaa.gov/hysplitusersguide/S310.htm
# https://www.ready.noaa.gov/hysplitusersguide/S312.htm


hysplit_batch_script =  """@echo off
setLocal EnableDelayedExpansion

set WRK=%CD%

set DIR=C:/Users/marti/Dropbox/research_projects/air_pollution
set PGM=%DIR%\hysplit

cd %PGM%\working

SET SYR={source_year}
SET SMO={source_month}
SET SDA={source_day}
SET SHR={source_hour}

SET INP={input_file_name}


echo %SYR% %SMO% %SDA% %SHR% >CONTROL
echo  3                      >>CONTROL
echo  {lower_left_source_grid} {source_height}       >>CONTROL
echo  {upper_right_source_grid} {source_height}        >>CONTROL
echo  {resolution_source_grid} {source_height}        >>CONTROL
echo {length_in_hours}       >>CONTROL
echo 0                       >>CONTROL
echo 10000.0                 >>CONTROL
echo 1                       >>CONTROL
echo %DIR%/met_files/  >>CONTROL
echo RP20{source_year}{source_month}.gbl   >>CONTROL
echo 1                      >>CONTROL
echo TEST                   >>CONTROL
echo 1000000.0               >>CONTROL
echo 1.0                    >>CONTROL
echo 00 00 00 00 00         >>CONTROL
echo 1                      >>CONTROL
echo 0.0 0.0                >>CONTROL
echo 0.05 0.05              >>CONTROL
echo 30.0 30.0              >>CONTROL
echo ./                     >>CONTROL
echo %INP%                  >>CONTROL
echo 1                      >>CONTROL
echo 1200                   >>CONTROL
echo 00 00 00 00 00         >>CONTROL
echo 00 00 00 00 00         >>CONTROL
echo 00 {conc_measurement_freq_hours} 00     >>CONTROL
echo 1                      >>CONTROL
echo {part_diam} {part_dens} {part_shape}          >>CONTROL
echo {dep_vel} 0.0 0.0 0.0 0.0    >>CONTROL
echo 0.0 0.0 0.0            >>CONTROL
echo 0.0                    >>CONTROL
echo 0.0                    >>CONTROL

REM -------------------------------------------


%PGM%\exec\latlon
IF EXIST %INP% DEL %INP%
%PGM%\exec\hycs_std
"""

# explanation 
# echo 1                      >>CONTROL
# echo 0.0 0.0 0.0            (part diameter, density, shape)
# echo 0.0 0.0 0.0 0.0 0.0    (vel, mol, a-ratio, d-ratio, henry)
# echo 0.0 0.0 0.0            ()
# echo 0.0                    >>CONTROL
# echo 0.0                    >>CONTROL




#print(hysplit_batch_script.format(**params_dict))



#import subprocess
#subprocess.call([r'C:\Users\marti\Dropbox\research_projects\air_pollution\my_hysplit_test.bat'])


# dep_vel should be in meters per second

def run_hysplit(source_year='14', source_month = '10', source_day = '01', source_hour = '00',
                conc_file_name_prefix = 'cdump_', 
                conc_folder_path = 'cdumps/ten_days_length/',
                length_in_hours = '96',
                conc_measurement_freq_hours = '12',
                lower_left_source_grid = '27.0 73.5', upper_right_source_grid = '32.5 78.5',
                resolution_source_grid = '27.5 74.0', source_height = '10.0', 
                part_diam = '1.0', part_dens = '2.0', part_shape='1.0', 
                dep_vel='0.002'):
    """Runs HYSPLIT simulation with specified parameters and creates a concetration file in the specified path. 

    Args:
        file_loc (str): The file location of the spreadsheet
        print_cols (bool): A flag used to print the columns to the console
            (default is False)

    Returns:
        list: a list of strings representing the header columns
"""


    final_file_name = conc_file_name_prefix + source_year + source_month + source_day + source_hour
    final_path_name = conc_folder_path + final_file_name

    params_dict = {'source_year': source_year, 'source_month': source_month, 'source_day': source_day, 'source_hour': source_hour,
                    'input_file_name': final_path_name, 'length_in_hours': length_in_hours,
                    'conc_measurement_freq_hours': conc_measurement_freq_hours,
                    'lower_left_source_grid': lower_left_source_grid,
                    'upper_right_source_grid': upper_right_source_grid,
                    'resolution_source_grid': resolution_source_grid,
                    'source_height': source_height, 'part_diam': part_diam,
                    'part_dens':part_dens, 'part_shape':part_shape,
                    'dep_vel': dep_vel}

    hysplit_cmd_final = hysplit_batch_script.format(**params_dict)

    my_batch_file = open(r'my_hysplit_test.bat','w+')
    my_batch_file.write(hysplit_cmd_final)
    my_batch_file.close()

    subprocess.call([r'my_hysplit_test.bat'])

    shutil.copy('hysplit/working/CONC.CFG', 
                os.path.join('hysplit/working/', conc_folder_path, 'CONC.CFG'))  
    shutil.copy('hysplit/working/CONTROL', 
                os.path.join('hysplit/working/', conc_folder_path, 'CONTROL'))  
    shutil.copy('hysplit/working/SETUP.CFG', 
                os.path.join('hysplit/working/', conc_folder_path, 'SETUP.CFG'))  

    






run_hysplit(source_year = '06', source_day = '01', source_hour = '09',
                length_in_hours = '240', conc_file_name_prefix = 'cdumps/test_cdump_final_pakistan_',
                conc_measurement_freq_hours = '80')


year_list = ['06', '07', '08', '09', '10', '11', '12', '15', '16', '17']

for year in year_list:
    print(year)
#    run_hysplit(source_year = year)
# second batch
    run_hysplit(source_year = year, source_day = '12', source_hour = '09',
                length_in_hours = '96', conc_file_name_prefix = 'cdumps/india_low_res/cdump_',
                conc_measurement_freq_hours = '96', source_height = '7.5',
                part_diam = '0.0', part_dens = '0.0', part_shape='0.0',)


#07, 12


# Vietnam
year_list = ['06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17']
# region ssof interest
# [104.5,9,107,11]

for year in year_list:
    print(year)
#    run_hysplit(source_year = year)
# second batch
    run_hysplit(source_year = year, source_day = '12', source_hour = '09', source_month='02',
                length_in_hours = '96', conc_file_name_prefix = 'cdumps/vietnam/cdump_',
                conc_measurement_freq_hours = '96', source_height = '7.5',
                lower_left_source_grid = '9 104.5', upper_right_source_grid = '11.0 107.0',
                resolution_source_grid = '9.5 105.0')



# China
year_list = ['06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17']
# region of interest
# [104.5,9,107,11]

for year in year_list:
    print(year)
#    run_hysplit(source_year = year)
# second batch
    run_hysplit(source_year = year, source_day = '12', source_hour = '09', source_month='02',
                length_in_hours = '96', conc_file_name_prefix = 'cdumps/vietnam/cdump_',
                conc_measurement_freq_hours = '96', source_height = '7.5',
                lower_left_source_grid = '9 104.5', upper_right_source_grid = '11.0 107.0',
                resolution_source_grid = '9.5 105.0')











# Higher source grid resolution 
year_list = ['06', '07', '08', '09', '10', '11', '12', '13',
             '14', '15', '16', '17', '18', '19']


days_list = ['01', '07', '20', '12']

# source days: 07, 01,20 to be done:  12

for day in days_list:
    print('day: ' + day)
    for year in year_list:
        print('year: ' + year)
        run_hysplit(source_year = year, source_day = day, source_hour = '09',
                    length_in_hours = '96', conc_file_name_prefix = 'cdump_',
                    conc_folder_path = 'cdumps/india/',
                    conc_measurement_freq_hours = '96', source_height = '7.5',
                    lower_left_source_grid = '27.0 73.5', upper_right_source_grid = '32.5 78.5',
                resolution_source_grid = '27.25 73.75')



year_list = ['06', '07', '08', '09', '10', '11', '12', '15', '16', '17']

# source days: 07, 01,20 to be done:  12

days_list = ['01', '07', '20', '12']

for day in days_list:
    print('day: ' + day)
    for year in year_list:
        print('year: ' + year)
        run_hysplit(source_year = year, source_day = day, source_hour = '09',
                    length_in_hours = '96', conc_file_name_prefix = 'cdump_',
                    conc_folder_path = 'cdumps/india_low_res/',
                    conc_measurement_freq_hours = '96', source_height = '7.5')


