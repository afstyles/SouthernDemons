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

#Calendar settings
"year0": 2000,
"month0": 1,
"day0": 1,

# Loading of ventilated data
"load_vent":True,

#Seeding information
"istmin":1,  #Minimum i-index seeded
"istmax":360,  #Maximum i-index seeded
"jstmin": 1,   #Minimum j-index seeded
"jstmax": 75,  #Maximum j-index seeded
"kstmin": 1,   #Minimum k-index seeded
"kstmax": 75,  #Maximum k-index seeded

#Subdomain information
"imindom": 1,  #Input values as seen in TRACMASS
"imaxdom": 360, #namelist
"jmindom": 1,
"jmaxdom": 100,
"kmindom": 1,
"kmaxdom": 75,

#Subdomain statistics
"subdomain_statistics": True,

#WG statistics
"gyre_stats": True,
"imingyre": 225,
"imaxgyre": None,
"jmingyre": None,
"jmaxgyre": None,

#Binned 2d statistics
"b2d_stats": True,

#Binned 3d statistics
"b3d_stats": True,

}




# First load the trajectories in Dask from multipe files
df_ini, df_out, df_run = io.load(data_dir, namelist["column_names"])



main.run_analysis(df_ini, df_out, df_run, namelist)

