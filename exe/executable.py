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
"out_label": "ORCA025_fwd",
"column_names": ['ntraj', 'x', 'y', 'z', 'subvol', 'time', 'boxface', 'temp', 'sal', 'density', 'density10' ],

#Combine file settings
"combine_log": True,

# Chunk settings
"blocksize": 250e6,

#Calendar settings
"tdir": 'fwd',     # Direction of time for trajectories (either "fwd" or "bwd")
"year0": 1982,
"month0": 12,
"day0": 16,

# Loading of ventilated data
"load_vent":False,
"update_vent": False,

"ndense_bin": False,
"sigma_out": True,
"ndsurf_dir": "/home/users/afstyles/SouthernDemons/neutraldensity/output/ORCA025_Dec_70_398_extend/",

"sf_zint_calc": False,
"sf_zint_dir": "/home/users/afstyles/SouthernDemons/sf_zint/output/ORCA025_398/",

"bathydepth_calc": True,

#Seeding information
"istmin": 1,  #Minimum i-index seeded
"istmax": 1440, #Maximum i-index seeded
"jstmin": 1,    #Minimum j-index seeded
"jstmax": 398,  #Maximum j-index seeded
"kstmin": 1,    #Minimum k-index seeded
"kstmax": 75,   #Maximum k-index seeded
"seed_month": None,  # = None consider all seedings, e.g. = 11 consider only seedings in November when creating df_vent

#Subdomain information
"imindom": 1,  #Input values as seen in TRACMASS
"imaxdom": 1440, #namelist
"jmindom": 1,
"jmaxdom": 400,
"kmindom": 1,
"kmaxdom": 75,

#Time window information
"timewindows"  : ["all"    , "preAug"       , "postAug"   ],
"tmin_list"    : [ None    , None           , 228*24*60*60],
"tmax_list"    : [ None    , 228*24*60*60   , None        ],

#Region information
"regions"        : [   "all"  ], # , "WG_i"   , "RG_i" , "ACC_i", "ACC_o", "WG_o" , "RG_o" , "WDP" , "SEAus" , "Shelf1000" , "Shelf2000", "msect000E", "msect030E",  "msect060E",  "msect090E",  "msect120E",  "msect150E",  "msect167E",  "msect185E",  "msect210E",  "msect240E",  "msect270E",  "msect300E",  "msect330E",  "zsect70S",  "zsect65S",  "zsect60S",  "zsect55S",  "zsect50S",  "zsect45S",  "zsect40S",  "zsect35S",  "zsect30S",], 
"imin_list"      : [   None   ], # , 900      , None   , None   , None   , 900    , None   , 650   , 200     , None        , None       , 1125       , 1248       ,  1382       ,  76         ,  210        ,  329        ,  380        ,  432        ,  529        ,  661        ,  797        ,  927        ,  1038       ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,], 
"imax_list"      : [   None   ], # , None     , 900    , None   , None   , None   , 900    , 900   , 450     , None        , None       , 1125       , 1248       ,  1382       ,  76         ,  210        ,  329        ,  380        ,  432        ,  529        ,  661        ,  797        ,  927        ,  1038       ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,], 
"jmin_list"      : [   None   ], # , None     , None   , None   , None   , None   , None   , 150   , 200     , None        , None       , None       , None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  99        ,  152       ,  195       ,  232       ,  265       ,  295       ,  322       ,  347       ,  371       ,], 
"jmax_list"      : [   None   ], # , 250      , 225    , None   , None   , 250    , 225    , 295   , 350     , 200         , 200        , None       , None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  99        ,  152       ,  195       ,  232       ,  265       ,  295       ,  322       ,  347       ,  371       ,], 
"filterby_list"  : [   None   ], # , "seed"   , "seed" , "seed" , "out"  , "out"  , "out"  , "out" , "out"   , "out"       , "out"      , "seed"     , "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"     ,  "seed"    ,  "seed"    ,  "seed"    ,  "seed"    ,  "seed"    ,  "seed"    ,  "seed"    ,  "seed"    ,  "seed"    ,], 
"bdepthmax_list" : [   None   ], # , None     , None   , None   , None   , None   , None   , None  , None    , 1000.       , 2000.      , None       , None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,], 
"sfmin_list"     : [   None   ], # , 10       , 10     , None   , None   , 10.    , 10.    , None  , None    , None        , None       , None       , None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,], 
"sfmax_list"     : [   None   ], # , None     , None   , -10    , -10    , None   , None   , None  , None    , None        , None       , None       , None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,], 
"bdepthmin_list" : [   None   ], # , None     , None   , None   , None   , None   , None   , None  , None    , None        , None       , None       , None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None       ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,  None      ,], 


#Density window information
"rhowindows"  :  ["all"  ],
"rhomin_list" :  [ None  ],   #Specify minimum density BIN you would like (counting from 0, inclusive)
"rhomax_list" :  [ None  ],   #Specify maximum density BIN you would like (counting from 0, inclusive)
                                          #Remember that -1 is the density bin for unclassified density
#Subdomain statistics
"subdomain_statistics": False,

#Subset approximate neutral density surfaces
"subset_ndsurfs": False,

#Binned 2d statistics
# "b2d_stats": False,
"xy_stats": False,
"yz_stats": False,
"xz_stats": False,
"yrho_stats" : False,
"xrho_stats" : False,
"sig_vs_nd_stats" : False,

#Binned 3d statistics
# "b3d_stats": False,

}

main.run_analysis(namelist)

