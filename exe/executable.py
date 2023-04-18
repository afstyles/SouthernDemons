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
"grid_dir": "/gws/nopw/j04/aopp/astyles/TRACMASS_DATA/DRAKKAR_SET/ORCA025/topo/",
"out_label": "ORCA025_ndensetest2",
"column_names": ['ntraj', 'x', 'y', 'z', 'subvol', 'time', 'boxface', 'temp', 'sal', 'density', 'density10' ],

#Combine file settings
"combine_log": False,

# Chunk settings
"blocksize": 250e6,

#Calendar settings
"year0": 2012,
"month0": 12,
"day0": 31,

# Loading of ventilated data
"load_vent":True,
"ndense_bin": False,
"ndsurf_dir": "/home/users/afstyles/SouthernDemons/neutraldensity/output/ORCA025_Dec_70/",

#Seeding information
"istmin": 1,  #Minimum i-index seeded
"istmax": 1440,  #Maximum i-index seeded
"jstmin": 1,   #Minimum j-index seeded
"jstmax": 300,  #Maximum j-index seeded
"kstmin": 1,   #Minimum k-index seeded
"kstmax": 75,  #Maximum k-index seeded
"seed_month": 12,  # = None consider all seedings, e.g. = 11 consider only seedings in November when creating df_vent

#Subdomain information
"imindom": 1,  #Input values as seen in TRACMASS
"imaxdom": 1440, #namelist
"jmindom": 1,
"jmaxdom": 400,
"kmindom": 1,
"kmaxdom": 75,

#Time window information
"timewindows"  : ["twentyeight"    , "preAug"       , "postAug"   ],
"tmin_list"    : [ None            , 122*24*60*60   , None        ],
"tmax_list"    : [ 28*365*24*60*60 , 28*365*24*60*60, 122*24*60*60],

#Region information
"regions"      : [ "all" ,  "msect1", "msect2", "msect3",], # "wg_seed", "wg_out" , "testbox_out", "rg_seed", "rg_out",],
"imin_list"    : [ None  ,   1000    , 1200    , 600    ,], #   900     , 900     , 700          ,  300     , 300     ,],
"imax_list"    : [ None  ,   1100    , 1300    , 800    ,], #  None     , None    , 900          ,  700     , 700     ,],
"jmin_list"    : [ None  ,   None    , None    , None   ,], #  None     , None    , 150          ,  None    , None    ,],
"jmax_list"    : [ None  ,   None    , None    , None   ,], #  None     , None    , 300          ,  None    , None    ,],
"filterby_list": [ None  ,   "seed"  , "seed"  , "seed" ,], # "seed"    , "out"   , "out"        , "seed"   , "out"   ,],

#Subdomain statistics
"subdomain_statistics": True,

#Binned 2d statistics
# "b2d_stats": False,
"xy_stats": True,
"yz_stats": True,
"xz_stats": False,
"yrho_stats" : True,
"xrho_stats" : True,

#Binned 3d statistics
# "b3d_stats": False,

}

main.run_analysis(namelist)

