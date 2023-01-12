"""
ventilation.py

This file contains all functions related to the analysis of deep water ventilation
"""

import datesandtime

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

