"""
main.py

This is the python script where analysis libraries are called based on namelist settings
"""

import ventilation
import os
import input_output as io

from dask.diagnostics import ProgressBar

def run_analysis(df_ini, df_out, df_run, namelist):

    out_dir = namelist["data_dir"] + '/OUTPUT.' + namelist["out_label"] + '/'
    if not os.path.exists(out_dir): os.mkdir(out_dir)

    if namelist["load_vent"] == False:
            
        #Identify the trajectories that that ventilated
        df_vent = ventilation.subset_traj(df_out, mlthresh=0.01)

        #Add information about initial position, initial tracers, seeding time, and ventilation time, 
        df_vent = ventilation.add_initial_info(df_vent, df_ini)

        #Add calendar info to df_vent
        df_vent = ventilation.add_calendar_info(df_vent, namelist)

        print("Saving ventilation file")
        with ProgressBar():
            io.save_vent_file(df_vent, out_dir)


    else:
        df_vent = io.load_vent_file(out_dir)


    


    return