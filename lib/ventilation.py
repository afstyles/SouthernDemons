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
                                  "day0":namelist["day0"], "leapcheck":False}

    year_i = df_vent["time_i"].apply(datesandtime.sec_to_year_30day, **timedict, meta=("year_i", "int64"))
    month_i = df_vent["time_i"].apply(datesandtime.sec_to_month_30day, **timedict, meta=("month_i", "int64"))
    day_i = df_vent["time_i"].apply(datesandtime.sec_to_day_30day, **timedict, meta=("day_i", "int64"))
    dayofyear_i = df_vent["time_i"].apply(datesandtime.sec_to_dayofyear_30day, **timedict, meta=("dayofyear_i", "int64"))

    year_o = df_vent["time_o"].apply(datesandtime.sec_to_year_30day, **timedict, meta=("year_o", "int64"))
    month_o = df_vent["time_o"].apply(datesandtime.sec_to_month_30day, **timedict, meta=("month_o", "int64"))
    day_o = df_vent["time_o"].apply(datesandtime.sec_to_day_30day, **timedict, meta=("day_o", "int64"))
    dayofyear_o = df_vent["time_o"].apply(datesandtime.sec_to_dayofyear_30day, **timedict, meta=("dayofyear_o", "int64"))

    df_vent["year_i"] = year_i
    df_vent["month_i"] = month_i
    df_vent["day_i"] = day_i
    df_vent["dayofyear_i"] = dayofyear_i

    df_vent["year_o"] = year_o
    df_vent["month_o"] = month_o
    df_vent["day_o"] = day_o
    df_vent["dayofyear_o"] = dayofyear_o

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

    df_vent["binnedx_i"] = df_vent["x_i"].map_partitions(pd.cut, xbins, labels=xcent, retbins=False, as_index=True).astype(int)
    df_vent["binnedy_i"] = df_vent["y_i"].map_partitions(pd.cut, ybins, labels=ycent, retbins=False, as_index=True).astype(int)
    df_vent["binnedz_i"] = df_vent["z_i"].map_partitions(pd.cut, zbins, labels=zcent, retbins=False, as_index=True).astype(int)

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


    df_vent["binnedx_o"] = df_vent["x_o"].map_partitions(pd.cut, xbins, labels=xcent, retbins=False, as_index=True).astype(int)
    df_vent["binnedy_o"] = df_vent["y_o"].map_partitions(pd.cut, ybins, labels=ycent, retbins=False, as_index=True).astype(int)
    df_vent["binnedz_o"] = df_vent["z_o"].map_partitions(pd.cut, zbins, labels=zcent, retbins=False, as_index=True).astype(int)

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




    group_year_vent  = df_tmp[["year_o", "time_o"]].groupby(["year_o"])
    group_month_vent = df_tmp[["month_o", "time_o"]].groupby(["month_o"])
    group_dayofyear_vent = df_tmp[["dayofyear_o","time_o"]].groupby(["dayofyear_o"])

    #Numbers of ventilated trajectories in each year, month, and dayofyear
    num_year_vent = group_year_vent.count()["time_o"]
    num_month_vent = group_month_vent.count()["time_o"]
    num_dayofyear_vent = group_dayofyear_vent.count()["time_o"]

    num_year_vent.rename("count")
    num_month_vent.rename("count")
    num_dayofyear_vent.rename("count")


    return num_year_vent, num_month_vent, num_dayofyear_vent

def bin_2dcells(df_vent, initialpos=True, finalpos=False):
    """
    Count the number of ventilated trajectories that:
       initialpos=True: Originated from a given cell
       finalpos=True: Ventilation in a given cell 
    """
    if initialpos==True:
        xstr = "binnedx_i"
        ystr = "binnedy_i"
    else:
        xstr = "binnedx_o"
        ystr = "binnedy_o"

    meta={"freq": "int64"             , "dayofyear_ld": "float64", "dayofyear_lq":"float64",
          "dayofyear_median":"float64", "dayofyear_uq": "float64", "dayofyear_ud":"float64",
          "dayofyear_mean":"float64"  , "month_mean": "float64"  , "month_median":"float64",
          "year_ld":"float64"         , "year_lq":"float64"      , "year_median":"float64" ,
          "year_uq":"float64"         ,"year_ud":"float64"       , "year_mean":"float64"   ,}

    xy_vent = df_vent.groupby([xstr, ystr]).apply(xy_custom_func, meta=meta)

    return xy_vent

def xy_custom_func(x):
        d = {}
        d['freq']             = x["time_o"].count().astype(int)
        d['dayofyear_ld']     = x["dayofyear_o"].quantile(q=0.1)
        d['dayofyear_lq']     = x["dayofyear_o"].quantile(q=0.25)
        d['dayofyear_median'] = x["dayofyear_o"].median()
        d['dayofyear_uq']     = x["dayofyear_o"].quantile(q=0.75)
        d['dayofyear_ud']     = x["dayofyear_o"].quantile(q=0.9 )
        d['dayofyear_mean']   = x["dayofyear_o"].mean()
        d['month_mean']       = x["month_o"].mean()
        d['month_median']     = x["month_o"].median()
        d['year_ld']          = x["year_o"].quantile(q=0.1)
        d['year_lq']          = x["year_o"].quantile(q=0.25)
        d['year_median']      = x["year_o"].median()
        d['year_uq']          = x["year_o"].quantile(q=0.75)
        d['year_ud']          = x["year_o"].quantile(q=0.9)
        d['year_mean']        = x["year_o"].mean()

        series = pd.Series(d, index=["freq"        , "dayofyear_ld", "dayofyear_lq"  , "dayofyear_median",
                                     "dayofyear_uq", "dayofyear_ud", "dayofyear_mean", "month_mean"      ,
                                     "month_median", "year_ld"     , "year_lq"       , "year_median"     ,
                                     "year_uq"     , "year_ud"     , "year_mean"     ,])

        return series


