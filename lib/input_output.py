"""
input_output.py

This file contains all functions related to the loading of trajectory data and the saving of data
"""

import dask.array as da
import dask.dataframe as dd
import glob
import re

def load(data_dir, column_names):

    print("io.load >>")

    ini_file_list = sorted(glob.glob(data_dir + '/*_ini.csv'))
    out_file_list = sorted(glob.glob(data_dir + '/*_out.csv'))
    run_file_list = sorted(glob.glob(data_dir + '/*_run.csv'))
    const_file_list = sorted(glob.glob(data_dir + '/*_const.csv'))

    nini   = len(ini_file_list)
    nout   = len(out_file_list)
    nrun   = len(run_file_list)
    nconst = len(const_file_list)

    if (nini == nout) and (nini == nrun) and (nini == nconst):
        print("nini = nout = nrun = nconst = ", nini)
        print("Safe to proceed with loading")

    else:
        print("io.load >>")
        print("nini, nout, nrun, and nconst do not match")
        print("nini = ", nini)
        print("nrun = ", nrun)
        print("nout = ", nout)
        print("nconst = ", nconst)
        print("Not safe to proceed with loading")
        print(">>")

        exit()

    dfs_ini = []
    dfs_out = [] 
    dfs_run = []


    for n in range(nini):
        print(n, " of ", nini, end="\r")
        ini_file = ini_file_list[n]
        out_file = out_file_list[n]
        run_file = run_file_list[n]
        const_file = const_file_list[n]

        #Test that files all match
        ini_file_base = ini_file[:-7]
        run_file_base = run_file[:-7]
        out_file_base = out_file[:-7]
        const_file_base = const_file[:-9]

        if not ( run_file_base == ini_file_base and 
                out_file_base == ini_file_base and 
                const_file_base == ini_file_base ):
            print("File bases don't match, not safe to proceed")
            print("INI: ", ini_file_base)
            print("RUN: ", run_file_base)
            print("OUT: ", out_file_base)
            print("CONST: ", const_file_base)
            exit()


        
        df_ini = dd.read_csv(ini_file, header=None, names=column_names)
        df_out = dd.read_csv(out_file, header=None, names=column_names)
        df_run = dd.read_csv(run_file, header=None, names=column_names)
        const = dd.read_csv(const_file, header=None, dtype={0:'int64'}).compute().values[0][0]

        df_ini["ntrajc"] = df_ini["ntraj"] + const
        df_out["ntrajc"] = df_out["ntraj"] + const
        df_run["ntrajc"] = df_run["ntraj"] + const

        dfs_ini = dfs_ini + [df_ini]
        dfs_out = dfs_out + [df_out]
        dfs_run = dfs_run + [df_run]


    df_ini_combined = dd.multi.concat(dfs_ini)
    df_out_combined = dd.multi.concat(dfs_out)
    df_run_combined = dd.multi.concat(dfs_run)

    print("io.load complete")

    return df_ini_combined, df_out_combined, df_run_combined


def save_vent_file(df_vent, out_dir):
    dd.to_csv(df_vent, out_dir + "/df_vent.csv", single_file=True, mode="wt", index=False)
    return

def load_vent_file(out_dir):
    df_vent = dd.read_csv(out_dir + "/df_vent.csv", header=0 )
    return df_vent

def save_gyre_stats(ny_vent_gyre, nmo_vent_gyre, ndoy_vent_gyre, out_dir):
    dd.to_csv(ny_vent_gyre,out_dir + "/ny_vent_gyre.csv" , single_file=True, mode="wt", index=True)
    dd.to_csv(nmo_vent_gyre, out_dir + "/nmo_vent_gyre.csv", single_file=True, mode="wt", index=True)
    dd.to_csv(ndoy_vent_gyre, out_dir + "/ndoy_vent_gyre.csv", single_file=True, mode="wt", index=True)
    return

def save_subdomain_stats(ny_vent, nmo_vent, ndoy_vent, out_dir):
    dd.to_csv(ny_vent,out_dir + "/ny_vent.csv" , single_file=True, mode="wt", index=True)
    dd.to_csv(nmo_vent, out_dir + "/nmo_vent.csv", single_file=True, mode="wt", index=True)
    dd.to_csv(ndoy_vent, out_dir + "/ndoy_vent.csv", single_file=True, mode="wt", index=True)
    return

def save_xy_vent(xy_vent_ini, xy_vent_out, out_dir):
    dd.to_csv(xy_vent_ini, out_dir + "xy_vent_ini.csv", single_file=True, mode="wt", index=True)
    dd.to_csv(xy_vent_out, out_dir + "xy_vent_out.csv", single_file=True, mode="wt", index=True)
    return



