# Drakkar_correct_dates.py
#
#Corrects the dates on ORCA025-N06 runs and ORCA0083-N06 runs

import os
import glob
from datetime import datetime
from datetime import timedelta
from calendar import isleap

# os.symlink("submission_scripts/submit.sh", "submit_link.sh")
data_dir = os.path.abspath("/gws/nopw/j04/nemo_vol1/ORCA0083-N006/means/")  #ORCA025
# data_dir = os.path.abspath("/gws/nopw/j04/nemo_vol1/ORCA0083-N006/means/") #ORCA12

symlink_dir = os.path.abspath("/gws/nopw/j04/aopp/astyles/TRACMASS_DATA/DRAKKAR_SET/ORCA12/5day/")

year_dir_list = sorted(glob.glob(data_dir + "/????/"))
day0_file = sorted(glob.glob(year_dir_list[0] + "/*d05T.nc"))[0]
day0_str = day0_file[-15:-7]
datetime0 = datetime.strptime(day0_str, '%Y%m%d')

#Test if the first date occurs in a leap year
# If it is a leap year advance 1 leap year to a non-leap year
while isleap(datetime0.year):
    datetime0 = datetime0 + timedelta(days=366)

date_list = []

for n in range(73):
    date_list = date_list + [datetime0.strftime("%m%d")]
    datetime0 = datetime0 + timedelta(days=5)

for year_dir in year_dir_list:
    print(year_dir)
    file_listT = sorted(glob.glob(year_dir + "/*d05T.nc"))
    file_listU = sorted(glob.glob(year_dir + "/*d05U.nc"))
    file_listV = sorted(glob.glob(year_dir + "/*d05V.nc"))

    year_str = file_listT[0][-15:-11]
    # print(year_str)
    print(file_listU[0])

    #Year quality check
    if len(file_listT) == len(file_listU) == len(file_listV) == 73:
        print(f"Year: {year_str}")
        n = 0
        for f in file_listT:
            fsym = symlink_dir + "/" + f[-28:-11] + date_list[n] + f[-7:]
            os.symlink(f, fsym)
            n = n + 1

        n = 0
        for f in file_listU:
            fsym = symlink_dir + "/" + f[-28:-11] + date_list[n] + f[-7:]
            os.symlink(f, fsym)
            n = n + 1

        n = 0
        for f in file_listV:
            fsym = symlink_dir + "/" + f[-28:-11] + date_list[n] + f[-7:]
            os.symlink(f, fsym)
            n = n + 1
    
    
    else:
        print(f"!!! Problem with year {year_str} !!!")
        print(len(file_listT))
        print(len(file_listU))
        print(len(file_listV))

        
