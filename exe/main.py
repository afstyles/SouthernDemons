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

        #Bin by initial position
        df_vent = ventilation.bin_by_initial_pos(df_vent, namelist)

        #Bin by final position
        df_vent = ventilation.bin_by_final_pos(df_vent, namelist)

        #Bin by final position
        print("Saving df_vent to csv file")
        with ProgressBar():
            io.save_vent_file(df_vent, out_dir)


    else:
        df_vent = io.load_vent_file(out_dir)
    

    if namelist["subdomain_statistics"] == True:
        print("Calculating ventilation statistics across the entire subdomain")
        with ProgressBar():
            ny_vent, nmo_vent, ndoy_vent = ventilation.subdomain_stats(df_vent)
            io.save_subdomain_stats(ny_vent, nmo_vent, ndoy_vent, out_dir)

    if namelist["gyre_stats"] == True:
        print("Calculating ventilation statistics for trajectories starting in gyre region")
        gyre_dict = {'imin':namelist["imingyre"], 'imax':namelist["imaxgyre"], "jmin":namelist["jmingyre"], "jmax":namelist["jmaxgyre"]}
        ny_vent_gyre, nmo_vent_gyre, ndoy_vent_gyre = ventilation.subdomain_stats(df_vent, **gyre_dict)
        
        with ProgressBar():
            io.save_gyre_stats(ny_vent_gyre, nmo_vent_gyre, ndoy_vent_gyre, out_dir)

    if namelist["b2d_stats"] == True:
        print("Calculating ventilation statistics for horizontal grid cells (2d)")
        xy_vent_ini = ventilation.bin_2dcells(df_vent, initialpos=True , finalpos=False)  
        xy_vent_out = ventilation.bin_2dcells(df_vent, initialpos=False, finalpos=True)

        with ProgressBar():
            io.save_xy_vent(xy_vent_ini, xy_vent_out, out_dir)


    if namelist["b3d_stats"] == True:
        print("Calculating ventilation statistics for all grid cells")


    


    return