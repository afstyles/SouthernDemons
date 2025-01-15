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

        #Bin by final position
        df_vent = ventilation.bin_by_final_pos(df_vent, namelist)


        print(df_vent)
        print("Saving df_vent to file")
        start = time.time()
        # with ProgressBar():
        io.save_vent_file(df_vent, out_dir)
        print(time.time() - start, "seconds")


    else:

            # dd.to_csv(df_vent, out_dir + "/df_vent.csv", single_file=True, mode="wt", index=False)
        if namelist["update_vent"] == True:
            print("Saving df_vent backup")
            timestring = str(int(time.time()))
            print(f"/df_vent.parquet --> df_vent.backup.{timestring}.parquet")
            os.rename(out_dir + "/df_vent.parquet", out_dir + f"/df_vent.backup.{timestring}.parquet")
            df_vent = io.load_vent_file(out_dir, timestring=timestring)

        else:
            df_vent = io.load_vent_file(out_dir)

    
    #Calculate neutral density bins using pre-calculated surfaces
    if namelist["ndense_bin"] == True:
        df_vent = ventilation.bin_by_initial_ndense(df_vent, namelist)

    if namelist["sigma_out"] == True:
        df_vent = ventilation.bin_by_out_sigma(df_vent, namelist)

    #Calculate the value of the depth-integrated stream function at seeding
    if namelist["sf_zint_calc"] == True:
        df_vent = ventilation.bin_by_initial_sf_zint(df_vent, namelist)

    if namelist["bathydepth_calc"] == True:
        df_vent = ventilation.bin_by_bathy_depth(df_vent, namelist)

    #Save the updated version of df_vent
    if namelist["update_vent"] == True:
        print("Saving updated version of df_vent")
        io.save_vent_file(df_vent, out_dir)


    
    imin_list = namelist["imin_list"]
    imax_list = namelist["imax_list"]
    jmin_list = namelist["jmin_list"]
    jmax_list = namelist["jmax_list"]
    tmin_list = namelist["tmin_list"]
    tmax_list = namelist["tmax_list"]
    sfmin_list = namelist["sfmin_list"]
    sfmax_list = namelist["sfmax_list"]
    bdepthmin_list = namelist["bdepthmin_list"]
    bdepthmax_list = namelist["bdepthmax_list"]
    rhomin_list = namelist["rhomin_list"]
    rhomax_list = namelist["rhomax_list"]
    filterby_list = namelist["filterby_list"]

    nt = len(tmin_list)
    nreg = len(imin_list)
    nrho = len(rhomin_list)


    for timewindow, itwind in zip(namelist["timewindows"],range(nt)):
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


        for rhowindow, irho in zip(namelist["rhowindows"],range(nrho)):
            rhomin = rhomin_list[irho]
            rhomax = rhomax_list[irho]

            if rhomin is None: rhomin = -float("inf")
            if rhomax is None: rhomax = +float("inf")

            df_vent_twindow = df_vent_twindow[ (df_vent_twindow["nd_bin_ini"] >= rhomin) ]
            df_vent_twindow = df_vent_twindow[ (df_vent_twindow["nd_bin_ini"] <= rhomax) ]


            for region, ireg in zip(namelist["regions"],range(nreg)):

                df_vent_tmp = df_vent_twindow.copy()
                
                if region.lower() == "all":
                    imin = None
                    imax = None
                    jmin = None
                    jmax = None
                    filterby = None
                    sfmin = None
                    sfmax = None
                    bdepthmin = None
                    bdepthmax = None
                    

                else:
                    imin = imin_list[ireg]
                    imax = imax_list[ireg]
                    jmin = jmin_list[ireg]
                    jmax = jmax_list[ireg]
                    filterby = filterby_list[ireg] 
                    sfmin = sfmin_list[ireg]
                    sfmax = sfmax_list[ireg]
                    bdepthmin = bdepthmin_list[ireg]
                    bdepthmax = bdepthmax_list[ireg]

                    if imin is None: imin = -float("inf")
                    if imax is None: imax = +float("inf")
                    if jmin is None: jmin = -float("inf")
                    if jmax is None: jmax = +float("inf")
                    if sfmin is None: sfmin = -float("inf")
                    if sfmax is None: sfmax = +float("inf")
                    if bdepthmin is None: bdepthmin = -float("inf")
                    if bdepthmax is None: bdepthmax = +float("inf")

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
                    df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["sf_zint"] >= sfmin) ]
                    df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["sf_zint"] <  sfmax)]
                    df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["sf_zint"] >= sfmin) ]
                    df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["bathy_depth" + suffstr] <  bdepthmax) ]
                    df_vent_tmp =  df_vent_tmp[ (df_vent_tmp["bathy_depth" + suffstr] >= bdepthmin) ]

                #Write to namelist for specific region / timewindow / rho combination
                io.write_reg_tim_rho_namelist(region, timewindow, rhowindow, filterby,
                                              tmin, tmax, imin, imax, jmin, jmax, sfmin,
                                              sfmax, rhomin, rhomax, bdepthmin, bdepthmax, out_dir)

                
                if namelist["subdomain_statistics"] == True:
                    start = time.time()
                    print(f"Calculating aggregated regional statistics")
                    sub_stats = ventilation.subdomain_stats(df_vent_tmp)
                    io.save_subdomain_stats(sub_stats, out_dir, label=f"Reg_{region}.Time_{timewindow}.Rho_{rhowindow}")
                    print((time.time() - start)/60, "minutes")
                    print("Complete")

                if namelist["xy_stats"] == True:
                    print("Calculating ventilation statistics for horizontal grid cells (2d)")
                    start = time.time()
                    # xy_n_vent_ini, xy_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, initialpos=True , finalpos=False)  
                    # xy_n_vent_out, xy_vol_vent_out = ventilation.bin_2dcells(df_vent_tmp, initialpos=False, finalpos=True)
                    xy_n_vent_ini, xy_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_i", ystr="binnedy_i")  
                    xy_n_vent_out, xy_vol_vent_out = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_o", ystr="binnedy_o")  
                    io.save_xy_vent(xy_n_vent_ini, xy_n_vent_out, xy_vol_vent_ini, xy_vol_vent_out, out_dir, label=f"Reg_{region}.Time_{timewindow}.Rho_{rhowindow}")
                    print((time.time() - start)/60, "minutes")
                    print("Complete")

                if namelist["yz_stats"] == True:
                    print("Calculating ventilation statistics for vertical grid cells (y,z)")
                    start = time.time()
                    yz_n_vent_ini, yz_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedy_i", ystr="binnedz_i")  
                    yz_n_vent_out, yz_vol_vent_out = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedy_o", ystr="binnedz_o")  
                    io.save_yz_vent(yz_n_vent_ini, yz_n_vent_out, yz_vol_vent_ini, yz_vol_vent_out, out_dir, label=f"Reg_{region}.Time_{timewindow}.Rho_{rhowindow}")
                    print((time.time() - start)/60, "minutes")
                    print("Complete")

                if namelist["xz_stats"] == True:
                    print("Calculating ventilation statistics for vertical grid cells (x,z)")
                    start = time.time()
                    xz_n_vent_ini, xz_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_i", ystr="binnedz_i")  
                    xz_n_vent_out, xz_vol_vent_out = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_o", ystr="binnedz_o")  
                    io.save_xz_vent(xz_n_vent_ini, xz_n_vent_out, xz_vol_vent_ini, xz_vol_vent_out, out_dir, label=f"Reg_{region}.Time_{timewindow}.Rho_{rhowindow}")
                    print((time.time() - start)/60, "minutes")
                    print("Complete")

                if namelist["yrho_stats"] == True:
                    print("Calculating ventilation statistics for density grid cells (y,rho)")
                    start = time.time()
                    yrho_n_vent_ini, yrho_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedy_i", ystr="nd_bin_ini")  
                    io.save_yrho_vent(yrho_n_vent_ini, yrho_vol_vent_ini, out_dir, label=f"Reg_{region}.Time_{timewindow}.Rho_{rhowindow}")
                    print((time.time() - start)/60, "minutes")
                    print("Complete")

                if namelist["xrho_stats"] == True:
                    print("Calculating ventilation statistics for density grid cells (x,rho)")
                    start = time.time()
                    xrho_n_vent_ini, xrho_vol_vent_ini = ventilation.bin_2dcells(df_vent_tmp, xstr="binnedx_i", ystr="nd_bin_ini")  
                    io.save_xrho_vent(xrho_n_vent_ini, xrho_vol_vent_ini, out_dir, label=f"Reg_{region}.Time_{timewindow}.Rho_{rhowindow}")
                    print((time.time() - start)/60, "minutes")
                    print("Complete")

                if namelist["sig_vs_nd_stats"] == True:
                    print("Calculating ventilation statistics for density only grid cells (rho_i,rho_o)")
                    start = time.time()
                    sig_vs_nd_n_vent, sig_vs_nd_vol_vent = ventilation.bin_2dcells(df_vent_tmp, xstr="sigmabin_o", ystr="nd_bin_ini")
                    io.save_sig_vs_nd_vent(sig_vs_nd_n_vent, sig_vs_nd_vol_vent, out_dir, label=f"Reg_{region}.Time_{timewindow}.Rho_{rhowindow}")
                    print((time.time() - start)/60, "minutes")
                    print("Complete")


        


                print(" >>>>>>>>>>>>>>>>>>>>>>>>")

    if namelist["subset_ndsurfs"] == True:
        
        for region in namelist["regions"]:
            print("REGION is", region)
            zsurf_xmean, zsurf_ymean, vol_bin_xint_cube, vol_bin_yint_cube = ventilation.subset_ndsurfs(namelist, region)
            io.save_subset_nd(zsurf_xmean, zsurf_ymean, vol_bin_xint_cube, vol_bin_yint_cube, out_dir, label=f"Reg_{region}")

            



    


    return