"""
ventilation.py

This file contains all functions related to the analysis of deep water ventilation
"""

import datesandtime
import numpy as np
import pandas as pd

def subset_traj(df_out, mlthresh=0.01):
    """
    Identify trajectories that were terminated through ventilation
    """

    df_vent = df_out[ df_out['density10'] <= mlthresh]

    return df_vent

def add_initial_info(df_vent, df_ini):

    ntrajc_col = df_vent["ntrajc"]
    df_vent = df_vent.merge(df_ini, how="left", on="ntrajc", suffixes=('_o','_i'))

    return df_vent

def add_calendar_info(df_vent, namelist):

    timedict = {"year0":namelist["year0"], "month0":namelist["month0"],
                                  "day0":namelist["day0"] }

    year_i = df_vent["time_i"].apply(datesandtime.sec_to_year, **timedict, meta=("year_i", "int64"))
    month_i = df_vent["time_i"].apply(datesandtime.sec_to_month, **timedict, meta=("month_i", "int64"))
    day_i = df_vent["time_i"].apply(datesandtime.sec_to_day, **timedict, meta=("day_i", "int64"))
    dayofyear_i = df_vent["time_i"].apply(datesandtime.sec_to_dayofyear, **timedict, meta=("dayofyear_i", "int64"))

    year_o = df_vent["time_o"].apply(datesandtime.sec_to_year, **timedict, meta=("year_o", "int64"))
    month_o = df_vent["time_o"].apply(datesandtime.sec_to_month, **timedict, meta=("month_o", "int64"))
    day_o = df_vent["time_o"].apply(datesandtime.sec_to_day, **timedict, meta=("day_o", "int64"))
    dayofyear_o = df_vent["time_o"].apply(datesandtime.sec_to_dayofyear, **timedict, meta=("dayofyear_o", "int64"))

    df_vent["year_i"] = year_i
    df_vent["month_i"] = month_i
    df_vent["day_i"] = day_i
    df_vent["dayofyear_i"] = dayofyear_i

    df_vent["year_o"] = year_o
    df_vent["month_o"] = month_o
    df_vent["day_o"] = day_o
    df_vent["dayofyear_o"] = dayofyear_o

    if namelist["seed_month"] is not None:
        print("Subsetting df_vent: Only seedings in month ", namelist["seed_month"])
        df_vent = df_vent[ df_vent['month_i'] == namelist["seed_month"]]

    return df_vent

def bin_by_initial_pos(df_vent, namelist):
    """
    Bin the df_vent dataframe by initial position.

    Bins are defined by the grid cell edges for easy comparison with model data
    """

    xmin, xmax = namelist["istmin"], namelist["istmax"]
    ymin, ymax = namelist["jstmin"], namelist["jstmax"]
    zmin, zmax = namelist["kstmin"], namelist["kstmax"]

    xbins = np.linspace(xmin-0.5, xmax+0.5, num=xmax-xmin+2)
    ybins = np.linspace(ymin-0.5, ymax+0.5, num=ymax-ymin+2)
    zbins = np.linspace(zmin-0.5, zmax+0.5, num=zmax-zmin+2)

    xbins[0], xbins[-1] = -float("inf"), float("inf")
    ybins[0], ybins[-1] = -float("inf"), float("inf")
    zbins[0], zbins[-1] = -float("inf"), float("inf")

    xcent = np.linspace(xmin, xmax, num=xmax-xmin+1, dtype=int)
    ycent = np.linspace(ymin, ymax, num=ymax-ymin+1, dtype=int)
    zcent = np.linspace(zmin, zmax, num=zmax-zmin+1, dtype=int)

    df_vent["binnedx_i"] = df_vent["x_i"].map_partitions(pd.cut, xbins, labels=xcent, retbins=False).astype(int)
    df_vent["binnedy_i"] = df_vent["y_i"].map_partitions(pd.cut, ybins, labels=ycent, retbins=False).astype(int)
    df_vent["binnedz_i"] = df_vent["z_i"].map_partitions(pd.cut, zbins, labels=zcent, retbins=False).astype(int)

    return df_vent

def bin_by_final_pos(df_vent, namelist):
    """
    Bin the df_vent dataframe by final (ventilation) position.

    Bins are defined by the grid cell edges for easy comparison with model data
    """

    xmin, xmax = namelist["imindom"], namelist["imaxdom"]
    ymin, ymax = namelist["jmindom"], namelist["jmaxdom"]
    zmin, zmax = namelist["kmindom"], namelist["kmaxdom"]

    xbins = np.linspace(xmin-0.5, xmax+0.5, num=xmax-xmin+2)
    ybins = np.linspace(ymin-0.5, ymax+0.5, num=ymax-ymin+2)
    zbins = np.linspace(zmin-0.5, zmax+0.5, num=zmax-zmin+2)

    xbins[0], xbins[-1] = -float("inf"), float("inf")
    ybins[0], ybins[-1] = -float("inf"), float("inf")
    zbins[0], zbins[-1] = -float("inf"), float("inf")

    xcent = np.linspace(xmin, xmax, num=xmax-xmin+1, dtype=int)
    ycent = np.linspace(ymin, ymax, num=ymax-ymin+1, dtype=int)
    zcent = np.linspace(zmin, zmax, num=zmax-zmin+1, dtype=int)


    df_vent["binnedx_o"] = df_vent["x_o"].map_partitions(pd.cut, xbins, labels=xcent, retbins=False).astype(int)
    df_vent["binnedy_o"] = df_vent["y_o"].map_partitions(pd.cut, ybins, labels=ycent, retbins=False).astype(int)
    df_vent["binnedz_o"] = df_vent["z_o"].map_partitions(pd.cut, zbins, labels=zcent, retbins=False).astype(int)

    return df_vent

def subdomain_stats(df_vent, imin=None, imax=None, jmin=None, jmax=None):
    """
    Calculate ventilation statistics over a subdomain of origin
    """

    if (imin == None) and (imax == None) and (jmin == None) and (jmax == None):
        df_tmp = df_vent.copy()

    else:
        if imin == None: imin = -float("inf")
        if imax == None: imax = +float("inf")
        if jmin == None: jmin = -float("inf")
        if jmax == None: jmax = +float("inf")

        df_tmp = df_vent[ (df_vent["binnedx_i"] >= imin) ]
        df_tmp =  df_tmp[ (df_tmp["binnedx_i"] <= imax) ]
        df_tmp =  df_tmp[ (df_tmp["binnedy_i"] >= jmin) ]
        df_tmp =  df_tmp[ (df_tmp["binnedy_i"] <= jmax) ]

    #Calculate age of each ventilated trajectory
    df_tmp["age"] = df_tmp["time_i"] - df_tmp["time_o"]
    df_tmp["age_years"] = np.floor(df_tmp["age"]/(365*24*60*60)).astype(int)


    group_age_vent   = df_tmp[["age_years"   , "time_o", "subvol_o", "subvol_i"]].groupby(["age_years"])
    group_year_vent  = df_tmp[["year_o", "time_o", "subvol_o", "subvol_i"]].groupby(["year_o"])
    group_month_vent = df_tmp[["month_o","time_o", "subvol_o", "subvol_i"]].groupby(["month_o"])
    group_dayofyear_vent = df_tmp[["dayofyear_o","time_o", "subvol_o", "subvol_i"]].groupby(["dayofyear_o"])

    #Numbers of ventilated trajectories in each year, month, and dayofyear
    num_age_vent = group_age_vent.count()["time_o"]
    num_year_vent = group_year_vent.count()["time_o"]
    num_month_vent = group_month_vent.count()["time_o"]
    num_dayofyear_vent = group_dayofyear_vent.count()["time_o"]

    num_age_vent.rename("count")
    num_year_vent.rename("count")
    num_month_vent.rename("count")
    num_dayofyear_vent.rename("count")

    #Volume (at seeding) of ventilated trajectories in each year, month, and dayofyear
    vol_i_age_vent = group_age_vent.sum()["subvol_i"]
    vol_i_year_vent = group_year_vent.sum()["subvol_i"]
    vol_i_month_vent = group_month_vent.sum()["subvol_i"]
    vol_i_dayofyear_vent = group_dayofyear_vent.sum()["subvol_i"]

    vol_i_age_vent.rename("vol_i")
    vol_i_year_vent.rename("vol_i")
    vol_i_month_vent.rename("vol_i")
    vol_i_dayofyear_vent.rename("vol_i")

    #Volume (at ventilation) of ventilated trajectories in each year, month, and dayofyear
    vol_o_age_vent = group_age_vent.sum()["subvol_o"]
    vol_o_year_vent = group_year_vent.sum()["subvol_o"]
    vol_o_month_vent = group_month_vent.sum()["subvol_o"]
    vol_o_dayofyear_vent = group_dayofyear_vent.sum()["subvol_o"]

    vol_o_age_vent.rename("vol_o")
    vol_o_year_vent.rename("vol_o")
    vol_o_month_vent.rename("vol_o")
    vol_o_dayofyear_vent.rename("vol_o")

    return num_age_vent, num_year_vent, num_month_vent, num_dayofyear_vent, vol_i_age_vent, vol_i_year_vent, vol_i_month_vent, vol_i_dayofyear_vent, vol_o_age_vent, vol_o_year_vent, vol_o_month_vent, vol_o_dayofyear_vent

# def bin_2dcells(df_vent, initialpos=True, finalpos=False):
def bin_2dcells(df_vent, xstr= "binnedx_i", ystr= "binnedy_i"):
    """
    Count the number and volume of ventilated trajectories that:
       initialpos=True: Originated from a given cell
       finalpos=True: Ventilation in a given cell 
    """
    # if initialpos==True:
    #     xstr = "binnedx_i"
    #     ystr = "binnedy_i"
    # else:
    #     xstr = "binnedx_o"
    #     ystr = "binnedy_o"

    # meta={"freq": "int64"             , "dayofyear_ld": "float64", "dayofyear_lq":"float64",
    #       "dayofyear_median":"float64", "dayofyear_uq": "float64", "dayofyear_ud":"float64",
    #       "dayofyear_mean":"float64"  , "month_mean": "float64"  , "month_median":"float64",
    #       "year_ld":"float64"         , "year_lq":"float64"      , "year_median":"float64" ,
    #       "year_uq":"float64"         ,"year_ud":"float64"       , "year_mean":"float64"   ,}

 
    meta={"freq": "int64"            , "vol": "float64"        ,
          "dayofyear_mean": "float64","dayofyear_sd": "float64", "dayofyear_skew": "float64", "dayofyear_kurtosis": "float64",
          "age_mean": "float64"      ,"age_sd": "float64"      , "age_skew": "float64"     , "age_kurtosis": "float64"     }

    df_vent_tmp = df_vent.copy()
    df_vent_tmp["age"] = df_vent_tmp["time_i"] - df_vent_tmp["time_o"]

    xy_n_vent = df_vent_tmp.groupby([xstr, ystr]).apply(xy_custom_func, meta=meta, weightstr=None)



    xy_vol_vent = df_vent_tmp.groupby([xstr,ystr]).apply(xy_custom_func, meta=meta, weightstr="subvol_i")

    return xy_n_vent, xy_vol_vent


def xy_custom_func(x, weightstr=None):
        d = {}
        d['freq']             = x["time_o"].count().astype(int)

        if weightstr is not None:
            d['vol']            = x[weightstr].sum()

        d['dayofyear_mean']   = mean_custom( x, "dayofyear_o", weightstr=weightstr )
        d['dayofyear_sd']     = sd_custom(x, d, "dayofyear_o", meanstr="dayofyear_mean", weightstr=weightstr, mode365=True )
        d['dayofyear_skew']   = skew_custom(x, d, "dayofyear_o", "dayofyear_sd", meanstr="dayofyear_mean", weightstr=weightstr, mode365=True)
        d['dayofyear_kurtosis']   = kurtosis_custom(x, d, "dayofyear_o", "dayofyear_sd", meanstr="dayofyear_mean", weightstr=weightstr, mode365=True)

        d['age_mean'] = mean_custom( x, "age", weightstr=weightstr)
        d['age_sd']   = sd_custom(x, d, "age", meanstr="age_mean", weightstr=weightstr, mode365=False)
        d['age_skew']   = skew_custom(x, d, "age", "age_sd", meanstr="age_mean", weightstr=weightstr, mode365=False)
        d['age_kurtosis']   = kurtosis_custom(x, d, "age", "age_sd", meanstr="age_mean", weightstr=weightstr, mode365=False)

        series = pd.Series(d, index=["freq"          , "vol"         , 
                                     "dayofyear_mean", "dayofyear_sd", "dayofyear_skew", "dayofyear_kurtosis",
                                     "age_mean"      , "age_sd"      , "age_skew"     , "age_kurtosis"])

        return series

def mean_custom(x, valstr, weightstr=None):
    """
    valstr = label for value we are calculated weighted average of
    weightstr = label for value we are using to weight the average (= None for no weighting)
    """
    if weightstr is None: 
        return np.average(x[valstr])
    else:
        return np.average(x[valstr], weights=x[weightstr])

def sd_custom(x, d, valstr, meanstr=None, weightstr=None, mode365=False):

    if meanstr is None:
        x_mean = 0.
    else:
        x_mean = d[meanstr]

    if weightstr is None:
        return np.std(delta(x[valstr],x_mean, mode365=mode365))
    else:
        return np.sqrt(np.average(delta(x[valstr], x_mean, mode365=mode365)**2, weights=x[weightstr] ))

def skew_custom(x, d, valstr, sdstr, meanstr=None, weightstr=None, mode365=False):

    if meanstr is None:
        x_mean = 0.
    else:
        x_mean = d[meanstr]

    if weightstr is None:
        return np.average(delta(x[valstr],x_mean, mode365=mode365)**3)/(d[sdstr]**3)
    else:
        return np.average(delta(x[valstr],x_mean, mode365=mode365)**3, weights=x[weightstr])/(d[sdstr]**3)

def kurtosis_custom(x, d, valstr, sdstr, meanstr=None, weightstr=None, mode365=False):

    if meanstr is None:
        x_mean = 0.
    else:
        x_mean = d[meanstr]

    if weightstr is None:
        return np.average(delta(x[valstr],x_mean,mode365=mode365)**4)/(d[sdstr]**4)
    else:
        return np.average(delta(x[valstr],x_mean,mode365=mode365)**4, weights=x[weightstr])/(d[sdstr]**4)


def delta(xcol,mu,mode365=False):
    output = xcol - mu
    if mode365 == True:
        output = output + 365. * (output < -365./2.)
        output = output - 365. * (output > 365./2.)
    
    return output

def apply_3dfunction(df, f3d, meta=None, args=(), kwarg_dict={} ):
    """
    Apply a 3 dimensional function f3d(x,y,z) to all trajectories in a DataFrame object, df

    df:  DataFrame object to apply function to
    f3d: function of x,y, and z to apply

    """

    df_out = df.apply(function, axis=1, meta=meta, args=args, **kwarg_dict )


    return df_out

def bin_by_initial_ndense(df_vent, namelist):
    """
    Bin the trajectories by the neutral density at seeding using pre-calculated surfaces
    """
    import xarray as xr
    import os

    ndsurf_dir = os.path.abspath(namelist["ndsurf_dir"])
    grid_dir = os.path.abspath(namelist["grid_dir"])


    ndsurf = xr.open_mfdataset(ndsurf_dir + "/*.nc", chunks="auto")
    zsurf = ndsurf["ns_depth"].max(dim="sigma_ref", skipna=True)



    grid_list = xr.open_dataset(grid_dir + "/mesh_zgr.nc")
    gdept_1d = grid_list["gdept_1d"].squeeze()
    e3w = grid_list["e3w_0"].squeeze()

    df_vent["nd_bin_ini"] = df_vent.apply(find_ndense_bin, axis=1, args=(zsurf, gdept_1d, e3w), meta=("nd_bin_ini",'int64'))

    return df_vent


def find_ndense_bin(df,zsurf,gdept_1d,e3w):

    x = df["binnedx_i"].astype(int)
    y = df["binnedy_i"].astype(int)
    z = df["z_i"]

    zsurf_col = zsurf[:,y-1,x-1].compute()
    
    #If no surfaces intersect column, set to -1 
    if zsurf_col.count() == 0: 
        return -1
    
    zfloor = int(np.floor(z))
    bin_no = np.linspace(0,zsurf_col.shape[0], num=zsurf_col.shape[0], dtype=int)

    if zfloor > 0: 
        depth = gdept_1d[zfloor-1] + e3w[zfloor,y,x] * (z-zfloor)
    else:
        depth = e3w[zfloor,y,x] * (z-zfloor)

    b = zsurf_col.searchsorted(depth) + int(bin_no[np.isfinite(zsurf_col)][0])

    return b



