"""
executable.py

This is the python script used to carry out trajectory analysis.

This script also contains all settings for the analysis.

"""

import sys
import os
sys.path.append(os.path.abspath("../lib"))

import input_output as io
import main



#First argument is the directory to the trajectories
data_dir = sys.argv[1]
namelist = {
"data_dir": data_dir,
"out_label": "TJtest",
"column_names": ['ntraj', 'x', 'y', 'z', 'subvol', 'time', 'boxface', 'temp', 'sal', 'density', 'density10' ],
"year0": 2000,
"month0": 1,
"day0": 1,
"load_vent":True
}




# First load the trajectories in Dask from multipe files
df_ini, df_out, df_run = io.load(data_dir, namelist["column_names"])



main.run_analysis(df_ini, df_out, df_run, namelist)

