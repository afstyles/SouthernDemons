"""
input_output.py

This file contains all functions related to the loading of trajectory data and the saving of data
"""
import dask
import dask.array as da
import dask.dataframe as dd
import glob
import re

def load(data_dir, column_names, blocksize=1e9):

    print("io.load >>")

    ini_file_list = sorted(glob.glob(data_dir + '/out_seed*/*_ini.csv'))
    out_file_list = sorted(glob.glob(data_dir + '/out_seed*/*_out.csv'))
    run_file_list = sorted(glob.glob(data_dir + '/out_seed*/*_run.csv'))
    const_file_list = sorted(glob.glob(data_dir + '/out_seed*/*_const.csv'))

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


        
        df_ini = dd.read_csv(ini_file, header=None, names=column_names, blocksize=blocksize)
        df_out = dd.read_csv(out_file, header=None, names=column_names, blocksize=blocksize)
        # df_run = dd.read_csv(run_file, header=None, names=column_names, blocksize=blocksize)
        const = dd.read_csv(const_file, header=None, dtype={0:'int64'}).compute().values[0][0]

        df_ini["ntrajc"] = df_ini["ntraj"] + const
        df_out["ntrajc"] = df_out["ntraj"] + const
        # df_run["ntrajc"] = df_run["ntraj"] + const

        dfs_ini = dfs_ini + [df_ini]
        dfs_out = dfs_out + [df_out]
        # dfs_run = dfs_run + [df_run]


    df_ini_combined = dd.multi.concat(dfs_ini)
    df_out_combined = dd.multi.concat(dfs_out)
    # df_run_combined = dd.multi.concat(dfs_run)

    df_ini_combined = df_ini_combined.repartition(partition_size=blocksize)
    df_out_combined = df_out_combined.repartition(partition_size=blocksize)
    # df_run_combined = df_run_combined.repartition(partition_size=blocksize)

    df_run_combined = None

    print("io.load complete")

    return df_ini_combined, df_out_combined, df_run_combined

def save_combined_ini_file(df_ini, out_dir):
    # dd.to_csv(df_ini, out_dir + "/df_ini.combined.csv", single_file=True, mode="wt", index=False)
    dd.to_parquet(df_ini, out_dir + "/df_ini.combined.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    return

def save_combined_out_file(df_out, out_dir):
    # dd.to_csv(df_out, out_dir + "/df_out.combined.csv", single_file=True, mode="wt", index=False)
    dd.to_parquet(df_out, out_dir + "/df_out.combined.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    return

def save_combined_run_file(df_run, out_dir):
    # dd.to_csv(df_run, out_dir + "/df_run.combined.csv", single_file=True, mode="wt", index=False)
    dd.to_parquet(df_run, out_dir + "/df_run.combined.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    return

def load_combined_run_file(out_dir, blocksize=1e9):
    # df_run = dd.read_csv(out_dir + "/df_run.combined.csv", header=0, blocksize=blocksize)
    df_run = dd.read_parquet(out_dir + "/df_run.combined.parquet", chunksize=blocksize)
    return df_run

def load_combined_ini_file(out_dir, blocksize=1e9):
    # df_ini = dd.read_csv(out_dir + "/df_ini.combined.csv", header=0, blocksize=blocksize)
    df_ini = dd.read_parquet(out_dir + "/df_ini.combined.parquet", chunksize=blocksize)
    return df_ini

def load_combined_out_file(out_dir, blocksize=1e9):
    # df_out = dd.read_csv(out_dir + "/df_out.combined.csv", header=0, blocksize=blocksize)
    df_out = dd.read_parquet(out_dir + "/df_out.combined.parquet", chunksize=blocksize)
    return df_out

def save_vent_file(df_vent, out_dir):
    # dd.to_csv(df_vent, out_dir + "/df_vent.csv", single_file=True, mode="wt", index=False)
    dd.to_parquet(df_vent, out_dir + "/df_vent.parquet", write_index=False, overwrite=True, write_metadata_file=True )
    return

def load_vent_file(out_dir, blocksize=1e9):
    # df_vent = dd.read_csv(out_dir + "/df_vent.parquet", header=0 , blocksize=blocksize)
    df_vent = dd.read_parquet(out_dir + "/df_vent.parquet", chunksize=blocksize)
    return df_vent

def save_gyre_stats(sub_stats, out_dir):
    
    sub_stats = [ s.to_frame() for s in sub_stats]
    sub_stats = [ s.reset_index() for s in sub_stats]

    #Trajectory number statistics
    dd.to_parquet(sub_stats[0], out_dir + "/ny_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[1], out_dir + "/nmo_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[2], out_dir + "/ndoy_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    
    #Trajector volume (at seeding) statistics
    dd.to_parquet(sub_stats[3], out_dir + "/voliy_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[4], out_dir + "/volimo_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[5], out_dir + "/volidoy_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )

    #Trajector volume (at ventilation) statistics
    dd.to_parquet(sub_stats[6], out_dir + "/voloy_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[7], out_dir + "/volomo_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[8], out_dir + "/volodoy_vent_gyre.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    return

def save_subdomain_stats(sub_stats, out_dir, label="nolabel"):

    sub_stats = [ s.to_frame() for s in sub_stats ]
    sub_stats = [ s.reset_index() for s in sub_stats]

    #Trajectory number statistics
    dd.to_parquet(sub_stats[0], f"{out_dir}/nage_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[1], f"{out_dir}/ny_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[2], f"{out_dir}/nmo_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[3], f"{out_dir}/ndoy_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    
    #Trajector volume (at seeding) statistics
    dd.to_parquet(sub_stats[4], f"{out_dir}/voliage_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[5], f"{out_dir}/voliy_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[6], f"{out_dir}/volimo_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[7], f"{out_dir}/volidoy_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )

    #Trajector volume (at ventilation) statistics
    dd.to_parquet(sub_stats[8], f"{out_dir}/voloage_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[9], f"{out_dir}/voloy_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[10], f"{out_dir}/volomo_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )
    dd.to_parquet(sub_stats[11], f"{out_dir}/volodoy_vent.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True  )

    return

def save_xy_vent(xy_n_vent_ini, xy_n_vent_out, xy_vol_vent_ini, xy_vol_vent_out, out_dir, label="nolabel"):

    dd.to_parquet(dd.from_pandas(xy_n_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/xy_n_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(xy_n_vent_out.reset_index().compute(), chunksize=1000), f"{out_dir}/xy_n_vent_out.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(xy_vol_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/xy_vol_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(xy_vol_vent_out.reset_index().compute(), chunksize=1000), f"{out_dir}/xy_vol_vent_out.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)

    return

def save_yz_vent(yz_n_vent_ini, yz_n_vent_out, yz_vol_vent_ini, yz_vol_vent_out, out_dir, label="nolabel"):

    dd.to_parquet(dd.from_pandas(yz_n_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/yz_n_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(yz_n_vent_out.reset_index().compute(), chunksize=1000), f"{out_dir}/yz_n_vent_out.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(yz_vol_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/yz_vol_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(yz_vol_vent_out.reset_index().compute(), chunksize=1000), f"{out_dir}/yz_vol_vent_out.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)

    return

def save_xz_vent(xz_n_vent_ini, xz_n_vent_out, xz_vol_vent_ini, xz_vol_vent_out, out_dir, label="nolabel"):

    dd.to_parquet(dd.from_pandas(xz_n_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/xz_n_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(xz_n_vent_out.reset_index().compute(), chunksize=1000), f"{out_dir}/xz_n_vent_out.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(xz_vol_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/xz_vol_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(xz_vol_vent_out.reset_index().compute(), chunksize=1000), f"{out_dir}/xz_vol_vent_out.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)

    return

def save_yrho_vent(yrho_n_vent_ini, yrho_vol_vent_ini, out_dir, label="nolabel"):

    dd.to_parquet(dd.from_pandas(yrho_n_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/yrho_n_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(yrho_vol_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/yrho_vol_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)

    return

def save_xrho_vent(xrho_n_vent_ini, xrho_vol_vent_ini, out_dir, label="nolabel"):

    dd.to_parquet(dd.from_pandas(xrho_n_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/xrho_n_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)
    dd.to_parquet(dd.from_pandas(xrho_vol_vent_ini.reset_index().compute(), chunksize=1000), f"{out_dir}/xrho_vol_vent_ini.{label}.parquet", write_index=False, overwrite=True, write_metadata_file=True)

    return


