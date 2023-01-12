""" 
combine_csv.py >>>>>>>>>>>>>>>>

Combines multiple csv files from simultaneous runs of TRACMASS
"""

import sys
import os
import glob
import re

import numpy as np
import dask.array as da
import dask.dataframe as dd


datadir = os.path.abspath(sys.argv[1])

ini_file_list = sorted(glob.glob(datadir + '/*_ini.csv'))
out_file_list = sorted(glob.glob(datadir + '/*_out.csv'))
run_file_list = sorted(glob.glob(datadir + '/*_run.csv'))

n_ini = len(ini_file_list)
n_out = len(out_file_list)
n_run = len(run_file_list)

if (n_ini != n_out) or (n_ini != n_run):
    print("Print number of files do not match")
    print("n_ini = ", n_ini)
    print("n_out = ", n_out)
    print("n_run = ", n_run)
    print("Not safe to proceed...")
    print("Exiting...")
    exit()

else:
    print("Print number of files match")
    print("n_ini = ", n_ini)
    print("n_out = ", n_out)
    print("n_run = ", n_run)
    print("Safe to proceed")

ntraj0 = int(0)
column_names = ['ntraj', 'x', 'y', 'z', 'subvol',
                    'time', 'boxface', 'temp', 'sal',
                    'density', 'density10' ]
problem_files = []

dtype={0: 'int64',   #Ntraj
       1: 'float64', #x
       2: 'float64', #y
       3:'float64',  #z
       4:'float64',  #subvol
       5:'float64',  #time
       6:'int64',    #boxface
       7: 'float64', #temperature
       8:'float64',  #salinity
       9:'float64',  #density
       10:'float64'} #density10

for n in range(n_run):

    print(f"{n} of {n_run} >> ntraj0 = {ntraj0}")
    ini_file = ini_file_list[n]

    try:
        df_ini = dd.read_csv(ini_file, header=None, dtype=dtype)
    except:
        print("ERROR READING FILE")
        problem_files = problem_files + [ini_file]
        continue

    df_ini.columns = column_names

    ntraj_max = df_ini.ntraj.max().compute()
    ntraj_min = df_ini.ntraj.min().compute()

    print("MAX: ", ntraj_max, " MIN: ", ntraj_min)
    

    if ntraj_min <= ntraj0:
        c = ntraj0 - ntraj_min + 1  #Calculate value that must be added to ntraj for all labels to be distinct
        ntraj0 = ntraj_max + c

    else:
        ntraj0 = ntraj_max
        c = int(0)

    print("C: ", c, end='\r')
    

    c_csv_file = ini_file[:-7] + "const.csv"
    with open(c_csv_file, 'w') as out:
        out.write(str(c))
        out.close()

if len(problem_files) > 0:
    print("PROBLEMS WITH THE FOLLOWING FILES:")
    for f in problem_files: print(f)