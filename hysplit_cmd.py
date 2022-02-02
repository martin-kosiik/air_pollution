import  subprocess

hysplit_cmd = """
set DIR=C:
set PGM=%DIR%\HYSPLIT
cd %PGM%\working

%PGM%\exec\hycs_std
"""

from subprocess import PIPE, run

command = ['echo', 'hello']
result = subprocess.run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
print(result.returncode, result.stdout, result.stderr)


subprocess.call(hysplit_cmd, shell=True)

subprocess.check_output()


myBat = open(r'C:\Users\marti\Dropbox\research_projects\air_pollution\my_test_2.bat','w+')
myBat.write('''@echo off
color a & date
pause
''')
myBat.close()



import subprocess
subprocess.call([r'+\run_exe.bat'])

import os
os.system(r"C:\Users\marti\Dropbox\research_projects\air_pollution\run_exe.bat")




params_dict = {'source_year': '14', 'source_month': '10', 'source_day': '01', 'source_hour': '00',
                'input_file_name': 'test_exec_cdump_new_2'}



hysplit_batch_script =  """@echo off
setLocal EnableDelayedExpansion

cd ..
set WRK=%CD%

set DIR=c:
set PGM=%DIR%\hysplit
set MAP=%PGM%\graphics\arlmap

cd %PGM%\working

SET SYR={source_year}
SET SMO={source_month}
SET SDA={source_day}
SET SHR={source_hour}

SET INP={input_file_name}


echo %SYR% %SMO% %SDA% %SHR% >CONTROL
echo  3                      >>CONTROL
echo  27.0 73.5 10.0         >>CONTROL
echo  32.0 78.5 10.0         >>CONTROL
echo  27.5 74.0 10.0         >>CONTROL
echo 72                      >>CONTROL
echo 0                       >>CONTROL
echo 10000.0                 >>CONTROL
echo 4                       >>CONTROL
echo C:/Users/marti/Dropbox/research_projects/air_pollution/met/  >>CONTROL
echo RP201410.gbl           >>CONTROL
echo C:/Users/marti/Dropbox/research_projects/air_pollution/met/  >>CONTROL
echo RP201310.gbl           >>CONTROL
echo C:/Users/marti/Dropbox/research_projects/air_pollution/met/  >>CONTROL
echo RP201610.gbl           >>CONTROL
echo C:/Users/marti/Dropbox/research_projects/air_pollution/met/  >>CONTROL
echo RP201710.gbl           >>CONTROL
echo 1                      >>CONTROL
echo TEST                   >>CONTROL
echo 1.0                    >>CONTROL
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
echo 00 00 00 01 00         >>CONTROL
echo 00 00 00 00 00         >>CONTROL
echo 00 12 00               >>CONTROL
echo 1                      >>CONTROL
echo 0.0 0.0 0.0            >>CONTROL
echo 0.0 0.0 0.0 0.0 0.0    >>CONTROL
echo 0.0 0.0 0.0            >>CONTROL
echo 0.0                    >>CONTROL
echo 0.0                    >>CONTROL

REM -------------------------------------------


%PGM%\exec\latlon
IF EXIST %INP% DEL %INP%
%PGM%\exec\hycs_std

"""


print(hysplit_batch_script.format(**params_dict))

hysplit_cmd_final = hysplit_batch_script.format(**params_dict)

myBat = open(r'C:\Users\marti\Dropbox\research_projects\air_pollution\my_hysplit_test.bat','w+')
myBat.write(hysplit_cmd_final)
myBat.close()



import subprocess
subprocess.call([r'C:\Users\marti\Dropbox\research_projects\air_pollution\my_hysplit_test.bat'])






def run_hysplit(source_year='14', source_month = '10', source_day = '01', source_hour = '00',
                conc_file_name_prefix = 'cdump_'):

    final_file_name = conc_file_name_prefix + source_year + source_month + source_day + source_hour
    params_dict = {'source_year': source_year, 'source_month': source_month, 'source_day': source_day, 'source_hour': source_hour,
                    'input_file_name': final_file_name}

    hysplit_cmd_final = hysplit_batch_script.format(**params_dict)

    my_batch_file = open(r'C:\Users\marti\Dropbox\research_projects\air_pollution\my_hysplit_test.bat','w+')
    my_batch_file.write(hysplit_cmd_final)
    my_batch_file.close()

    subprocess.call([r'C:\Users\marti\Dropbox\research_projects\air_pollution\my_hysplit_test.bat'])





year_list = ['13', '14', '16', '17']


for year in year_list:
    run_hysplit(source_year = year)





new_cmd = """@echo off
setLocal EnableDelayedExpansion

cd ..
set WRK=%CD%

set DIR=c:
set PGM=%DIR%\hysplit
set MAP=%PGM%\graphics\arlmap

cd %PGM%\working

REM rerun model -------------------------------------------

SET SYR=83
SET SMO=09
SET SDA=01
SET SHR=00

SET LAT=43.0
SET LON=-75.0
SET LVL=10.0

SET RUN=60
SET TOP=10000.0

SET MET=%WRK%\captex
SET DAT=RP198309.gbl
SET INP=srm.bin

echo %SYR% %SMO% %SDA% %SHR% >CONTROL
echo 3                      >>CONTROL
echo 40.0 -78.0 %LVL%       >>CONTROL
echo 45.0 -70.0 %LVL%       >>CONTROL
echo 40.5 -77.5 %LVL%       >>CONTROL
echo %RUN%                  >>CONTROL
echo 0                      >>CONTROL
echo %TOP%                  >>CONTROL
echo 1                      >>CONTROL
echo %MET%\                 >>CONTROL
echo %DAT%                  >>CONTROL
echo 1                      >>CONTROL
echo PMCH                   >>CONTROL
echo 1.0                    >>CONTROL
echo 60.0                   >>CONTROL
echo 83 09 01 00 00         >>CONTROL
echo 1                      >>CONTROL
echo 41.0 -73.0             >>CONTROL
echo 0.25 0.25              >>CONTROL
echo 15.0 25.0              >>CONTROL
echo .\                     >>CONTROL
echo %INP%                  >>CONTROL
echo 1                      >>CONTROL
echo 100                    >>CONTROL
echo 83 09 01 00 00         >>CONTROL
echo 83 09 03 12 00         >>CONTROL
echo 00 03 00               >>CONTROL
echo 1                      >>CONTROL
echo 0.0 0.0 0.0            >>CONTROL
echo 0.0 0.0 0.0 0.0 0.0    >>CONTROL
echo 0.0 0.0 0.0            >>CONTROL
echo 0.0                    >>CONTROL
echo 0.0                    >>CONTROL

REM -------------------------------------------

%PGM%\exec\latlon
IF EXIST %INP% DEL %INP%
%PGM%\exec\hycs_std
"""
