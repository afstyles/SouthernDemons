"""
main.py

This is the python script where analysis libraries are called based on namelist settings
"""

import ventilation
import os
import input_output as io
import time

from dask.diagnostics import ProgressBar



def run_analysis(namelist):

    data_dir = namelist["data_dir"]
    out_dir = namelist["data_dir"] + '/OUTPUT.' + namelist["out_label"] + '/'
    if not os.path.exists(out_dir): os.mkdir(out_dir)

    if namelist["combine_log"] == True:
        # First load the trajectories in Dask from multipe files
        df_ini, df_out, df_run = io.load(data_dir, namelist["column_names"], blocksize=namelist["blocksize"])

        print("Saving combined ini file:")
        io.save_combined_ini_file(df_ini, out_dir)

        print("Saving combined out file:")
        io.save_combined_out_file(df_out, out_dir)
        
        #io.save_combined_run_file(df_run, out_dir)
    
    if namelist["load_vent"] == False:
        print(" Loading combined ini data")
        start = time.time()
        df_ini = io.load_combined_ini_file(out_dir, blocksize=namelist["blocksize"])
        print(time.time() - start, "seconds")

        print(" Loading combined out data")
        start = time.time()
        df_out = io.load_combined_out_file(out_dir, blocksize=namelist["blocksize"])
        print(time.time() - start, "seconds")
        # df_run = io.load_combined_run_file(out_dir)

        df_run = None


    if namelist["load_vent"] == False:
        
        #Identify the trajectories that that ventilated
        df_vent = ventilation.subset_traj(df_out, mlthresh=0.01)
        df_vent = df_vent.repartition(partition_size=namelist["blocksize"])

        #Add information about initial position, initial tracers, seeding time, and ventilation time, 
        df_vent = ventilation.add_initial_info(df_vent, df_ini)

        #Add calendar info to df_vent
        df_vent = ventilation.add_calendar_info(df_vent, namelist)

        # if namelist["seed_month"] is not None:
        #     df_vent = df_vent.repartition(partition_size=namelist["blocksize"])

        #Bin by initial position
        df_vent = ventilation.bin_by_initial_pos(df_vent, namelist)

        #Calculate neutral density bins using pre-calculated surfaces
        if namelist["ndense_bin"] == True:
            df_vent = ventilation.bin_by_initial_ndense(df_vent, namelist)

        #Bin by final position
        df_vent = ventilation.bin_by_final_pos(df_vent, namelist)


        print(df_vent)
        print("Saving df_vent to file")
        start = time.time()
        # with ProgressBar():
        io.save_vent_file(df_vent, out_dir)
        print(time.time() - start, "seconds")


    else:
        df_vent = io.load_vent_file(out_dir)
    
    
    imin_list = namelist["imin_list"]
    imax_list = namelist["imax_list"]
    jmin_list = namelist["jmin_list"]
    jmax_list = namelist["jmax_list"]
    tmin_list = namelist["tmin_list"]
    tmax_list = namelist["tmax_list"]
    filterby_list = namelist["filterby_list"]

    itwind = 0
    for timewindow in namelist["timewindows"]:
        df_vent_twindow = df_vent.copy()

        if timewindow.lower() == "all":
            tmin = None
            tmax = None
        else:
            tmin = tmin_list[itwind]
            tmax = tmax_list[itwind]

            if tmin is None: tmin = -float("inf")
            if tmax is None: tmax = +float("inf")

            df_vent_twindow =  df_vent_twindow[ (df_vent_twindow["time_i"] - df_vent_twindow["time_o"] >= tmin) ]
            df_vent_twindow =  df_vent_twindow[ (df_vent_twindow["time_i"] - df_vent_twindow["time_o"] <= tmax) ]

        
        
        ireg = 0
        for region in namelist["regions"]:

            df_vent_tmp = df_vent_twindow.copy()
            
            if region.lower() == "all":
                imin = None
                imax = None
                jmin = None
                jmax = None
                filterby = None
                

            else:
                imin = imin_list[ireg]
                imax = imax_list[ireg]
                jmin = jmin_list[ireg]
                jmax = jmax_list[ireg]
                filterby = filterby_list[ireg] 

                if imin is None: imin = -float("inf")
                if imax is None: imax = +float("inf")
                if jmin is None: jmin = -float("inf")
                if jmax is None: jmax = +float("inf")

                if filterby.lower() == "seed": 
                    suffstr = "_i"
                elif filterby.lower() == "out":
                    suffstr = "_o"
                else:
                    print("Invalid filterby condition", filterby.lower())
                    print("filterby must be either seed or out")
                    print("Will assume seeding position")
                    suffstr = "_i"
                

                df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["binnedx" + suffstr] >= imin) ]
                df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["binnedx" + suffstr] <= imax) ]
                df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["binnedy" + suffstr] >= jmin) ]
                df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["binnedy" + suffstr] <= jmax) ]
                

            print(f" Region: {region} >>>>>>>>>")
            print(f" timewindow: {timewindow}")
            print(f" tmin: {tmin}")
            print(f" tmax: {tmax}")
            print(f" imin: {imin}")
            print(f" imax: {imax}")
            print(f" jmin: {jmin}")
            print(f" jmax: {jmax}")
            if filterby is None: print("No filter applied")
            elif filterby.lower()=="seed": print("Filtering by seeding position (filterby == seed)")
            elif filterby.lower() == "out": print("Filtering by ventilation position (filterby == out)")
            


            
            if namelist["subdomain_statistics"] == True:
                start = time.time()
                print(f"Calculating aggregated regional statistics")
                sub_stats = ventilation.subdomain_stats(df_vent_tmp)
                io.save_subdomain_stats(sub_stats, out_dir, label=f"{region}.{timewindow}")
                print((time.time() - start)/60, "minutes")
                print("Complete")

            if namelist["xy_stats"] == True:
                print("Calculating ventilation statistics for horizontal grid cells (2d)")
                start = time.time()
                # xy_n_vent_ini, xy_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, initialpos=True , finalpos=False)  
                # xy_n_vent_out, xy_vol_vent_out = ventilation.bin_2dcells(df_vent_tmp, initialpos=False, finalpos=True)
                xy_n_vent_ini, xy_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_i", ystr="binnedy_i")  
                xy_n_vent_out, xy_vol_vent_out = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_o", ystr="binnedy_o")  
                io.save_xy_vent(xy_n_vent_ini, xy_n_vent_out, xy_vol_vent_ini, xy_vol_vent_out, out_dir, label=f"{region}.{timewindow}")
                print((time.time() - start)/60, "minutes")
                print("Complete")

            if namelist["yz_stats"] == True:
                print("Calculating ventilation statistics for vertical grid cells (y,z)")
                start = time.time()
                yz_n_vent_ini, yz_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedy_i", ystr="binnedz_i")  
                yz_n_vent_out, yz_vol_vent_out = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedy_o", ystr="binnedz_o")  
                io.save_yz_vent(yz_n_vent_ini, yz_n_vent_out, yz_vol_vent_ini, yz_vol_vent_out, out_dir, label=f"{region}.{timewindow}")
                print((time.time() - start)/60, "minutes")
                print("Complete")

            if namelist["xz_stats"] == True:
                print("Calculating ventilation statistics for vertical grid cells (x,z)")
                start = time.time()
                xz_n_vent_ini, xz_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_i", ystr="binnedz_i")  
                xz_n_vent_out, xz_vol_vent_out = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_o", ystr="binnedz_o")  
                io.save_xz_vent(xz_n_vent_ini, xz_n_vent_out, xz_vol_vent_ini, xz_vol_vent_out, out_dir, label=f"{region}.{timewindow}")
                print((time.time() - start)/60, "minutes")
                print("Complete")

            if namelist["yrho_stats"] == True:
                print("Calculating ventilation statistics for density grid cells (y,rho)")
                start = time.time()
                yrho_n_vent_ini, yrho_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedy_i", ystr="nd_bin_ini")  
                io.save_yrho_vent(yrho_n_vent_ini, yrho_vol_vent_ini, out_dir, label=f"{region}.{timewindow}")
                print((time.time() - start)/60, "minutes")
                print("Complete")

            if namelist["xrho_stats"] == True:
                print("Calculating ventilation statistics for density grid cells (x,rho)")
                start = time.time()
                xrho_n_vent_ini, xrho_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_i", ystr="nd_bin_ini")  
                io.save_xrho_vent(xrho_n_vent_ini, xrho_vol_vent_ini, out_dir, label=f"{region}.{timewindow}")
                print((time.time() - start)/60, "minutes")
                print("Complete")



            # if namelist["b3d_stats"] == True:
            #     print("Calculating ventilation statistics for all grid cells")

            ireg = ireg + 1
        
        itwind = itwind + 1

        print(" >>>>>>>>>>>>>>>>>>>>>>>>")


    


    return