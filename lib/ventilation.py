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
    group_ndense_vent = df_tmp[["nd_bin_ini","time_o", "subvol_o", "subvol_i"]].groupby(["nd_bin_ini"])

    #Numbers of ventilated trajectories in each year, month, and dayofyear
    num_age_vent = group_age_vent.count()["time_o"]
    num_year_vent = group_year_vent.count()["time_o"]
    num_month_vent = group_month_vent.count()["time_o"]
    num_dayofyear_vent = group_dayofyear_vent.count()["time_o"]
    num_ndense_vent = group_ndense_vent.count()["time_o"]

    num_age_vent.rename("count")
    num_year_vent.rename("count")
    num_month_vent.rename("count")
    num_dayofyear_vent.rename("count")
    num_ndense_vent.rename("count")

    #Volume (at seeding) of ventilated trajectories in each year, month, and dayofyear
    vol_i_age_vent = group_age_vent.sum()["subvol_i"]
    vol_i_year_vent = group_year_vent.sum()["subvol_i"]
    vol_i_month_vent = group_month_vent.sum()["subvol_i"]
    vol_i_dayofyear_vent = group_dayofyear_vent.sum()["subvol_i"]
    vol_i_ndense_vent = group_ndense_vent.sum()["subvol_i"]

    vol_i_age_vent.rename("vol_i")
    vol_i_year_vent.rename("vol_i")
    vol_i_month_vent.rename("vol_i")
    vol_i_dayofyear_vent.rename("vol_i")
    vol_i_ndense_vent.rename("vol_i")

    #Volume (at ventilation) of ventilated trajectories in each year, month, and dayofyear
    vol_o_age_vent = group_age_vent.sum()["subvol_o"]
    vol_o_year_vent = group_year_vent.sum()["subvol_o"]
    vol_o_month_vent = group_month_vent.sum()["subvol_o"]
    vol_o_dayofyear_vent = group_dayofyear_vent.sum()["subvol_o"]
    vol_o_ndense_vent = group_ndense_vent.sum()["subvol_o"]

    vol_o_age_vent.rename("vol_o")
    vol_o_year_vent.rename("vol_o")
    vol_o_month_vent.rename("vol_o")
    vol_o_dayofyear_vent.rename("vol_o")
    vol_o_ndense_vent.rename("vol_o")

    return num_age_vent, num_year_vent, num_month_vent, num_dayofyear_vent, num_ndense_vent, vol_i_age_vent, vol_i_year_vent, vol_i_month_vent, vol_i_dayofyear_vent, vol_i_ndense_vent, vol_o_age_vent, vol_o_year_vent, vol_o_month_vent, vol_o_dayofyear_vent, vol_o_ndense_vent

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
          "age_mean": "float64"      ,"age_sd": "float64"      , "age_skew": "float64"     , "age_kurtosis": "float64"       ,
          "ndense_median": "int64"   , "ndende_uq": "int64"    , "ndense_lq": "int64"      , "ndense_ud" : "int64"           ,
          "ndense_ld" : "int64"      }

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

        d['ndense_median'] = percentile_custom( x , "nd_bin_ini", 0.5 , weightstr=weightstr )
        d["ndense_uq"] =     percentile_custom( x , "nd_bin_ini", 0.75, weightstr=weightstr )
        d["ndense_lq"] =     percentile_custom( x , "nd_bin_ini", 0.25, weightstr=weightstr )
        d["ndense_ud"] =     percentile_custom( x , "nd_bin_ini", 0.9 , weightstr=weightstr )
        d["ndense_ld"] =     percentile_custom( x , "nd_bin_ini", 0.1 , weightstr=weightstr )

        series = pd.Series(d, index=["freq"          , "vol"         , 
                                     "dayofyear_mean", "dayofyear_sd", "dayofyear_skew", "dayofyear_kurtosis",
                                     "age_mean"      , "age_sd"      , "age_skew"     , "age_kurtosis"       ,
                                     "ndense_median" , "ndende_uq"   , "ndense_lq"    , "ndense_ud"          ,
                                     "ndense_ld"     ])

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

def percentile_custom(x, valstr, p, weightstr=None):

    x_sort = x.sort_values(valstr, ascending=True)

    if weightstr is not None:        
        cs = np.cumsum(x_sort[weightstr])
        s = np.sum(x_sort[weightstr])
    
    else:
        cs = np.arange(1,len(x_sort[valstr])+1, dtype=float)
        s = len(x_sort[valstr])

    ind = np.searchsorted(cs/s, p, side="left")
    output = x_sort[valstr].values[ind]

    return output

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

    istmin = namelist["istmin"] - 1  #Load seeding indices (python convention)
    istmax = namelist["istmax"] - 1
    jstmin = namelist["jstmin"] - 1
    jstmax = namelist["jstmax"] - 1  


    ndsurf_dir = os.path.abspath(namelist["ndsurf_dir"])
    grid_dir = os.path.abspath(namelist["grid_dir"])


    ndsurf = xr.open_mfdataset(ndsurf_dir + "/*.nc", chunks="auto")
    zsurf = ndsurf["ns_depth"].max(dim="sigma_ref", skipna=True)


    grid_list = xr.open_dataset(grid_dir + "/mesh_zgr.nc")
    gdept_1d = grid_list["gdept_1d"].squeeze()
    e3w = grid_list["e3w_0"].squeeze()[...,jstmin:jstmax+1, istmin:istmax+1]

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

    bin_no = np.arange(zsurf_col.shape[0], dtype=int)

    zsurfs = zsurf_col[np.isfinite(zsurf_col)]
    bins   = bin_no[np.isfinite(zsurf_col)]
    bins   = np.append(bins, int(bins[-1] + 1) )

    zfloor = int(np.floor(z))

    if zfloor > 0: 
        depth = gdept_1d[zfloor-1] + e3w[zfloor,y-1,x-1] * (z-zfloor)
    else:
        depth = e3w[zfloor,y-1,x-1] * (z-zfloor)

    b = bins[zsurfs.searchsorted(depth)]

    return b

def subset_ndsurfs(namelist, region):

    import dask.array as da
    import xarray as xr
    import os

    ireg = 0
    ndsurf_dir = os.path.abspath(namelist["ndsurf_dir"])
    grid_dir = os.path.abspath(namelist["grid_dir"])
    ndsurf = xr.open_mfdataset(ndsurf_dir + "/*.nc", chunks="auto")
    zsurf = ndsurf["ns_depth"].max(dim="sigma_ref", skipna=True)

    istmin = namelist["istmin"] - 1  #Load seeding indices (python convention)
    istmax = namelist["istmax"] - 1
    jstmin = namelist["jstmin"] - 1
    jstmax = namelist["jstmax"] - 1  

    imin_list = namelist["imin_list"]
    imax_list = namelist["imax_list"]
    jmin_list = namelist["jmin_list"]
    jmax_list = namelist["jmax_list"]

    hgrid_list = xr.open_dataset(grid_dir + "/mesh_hgr.nc", chunks="auto")
    e1t = hgrid_list["e1t"].squeeze()[jstmin:jstmax+1, istmin:istmax+1]
    e2t = hgrid_list["e2t"].squeeze()[jstmin:jstmax+1, istmin:istmax+1]
    area = ( e1t * e2t ).broadcast_like(zsurf)
    area = area.where(zsurf is not np.nan )


    if region.lower() == "all":

        #Calculate volume above density surfaces
        zvolup_xint = ( zsurf * area ).sum(axis=-1, skipna=True, min_count=1)
        zvolup_yint = ( zsurf * area ).sum(axis=-2, skipna=True, min_count=1)

        #Zonal and meridional averages of density surfaces
        zsurf_xmean = zvolup_xint / area.sum(axis=-1, skipna=True, min_count=1)
        zsurf_ymean = zvolup_yint / area.sum(axis=-2, skipna=True, min_count=1)

    else:
        imin = imin_list[ireg]
        imax = imax_list[ireg]
        jmin = jmin_list[ireg]
        jmax = jmax_list[ireg]

        if imin is not None: imin = imin - 1
        if jmin is not None: jmin = jmin - 1

        zvolup_xint = subset_ij(zsurf*area,imin,imax,jmin,jmax).sum(axis=-1, skipna=True, min_count=1)
        zvolup_yint = subset_ij(zsurf*area,imin,imax,jmin,jmax).sum(axis=-2, skipna=True, min_count=1)

        zsurf_xmean = zvolup_xint / subset_ij(area,imin,imax,jmin,jmax).sum(axis=-1, skipna=True, min_count=1)
        zsurf_ymean = zvolup_yint / subset_ij(area,imin,imax,jmin,jmax).sum(axis=-2, skipna=True, min_count=1)

    #Calculate the volume of each density bin
    zgrid_list = xr.open_dataset(grid_dir + "/mesh_zgr.nc", chunks="auto")
    mask_list = xr.open_dataset(grid_dir + "/mask.nc", chunks="auto")
    e3t = zgrid_list["e3t_0"].squeeze()[...,jstmin:jstmax+1, istmin:istmax+1]
    tmask = mask_list["tmask"].squeeze()[...,jstmin:jstmax+1, istmin:istmax+1].astype(bool)
    tmask2d = mask_list["tmaskutil"].squeeze()[jstmin:jstmax+1, istmin:istmax+1].astype(bool)
    depth = e3t.where(tmask).sum(axis=-3, skipna=True)
    Hvol_xint = (e1t*e2t*depth).where(tmask2d).sum(axis=-1, skipna=True, min_count=1)
    Hvol_yint = (e1t*e2t*depth).where(tmask2d).sum(axis=-2, skipna=True, min_count=1)
    

    @da.as_gufunc(signature="(p),()->(p),()", axes=[(-2),(),(-2),()], allow_rechunk=True, vectorize=True)
    def gufoo(zvolup, Hvol):
        
        #Indices of each density surface (not density bin)
        ind = np.arange(zvolup.shape[-1], dtype=int)
        output = np.full(zvolup.shape[-1] + 1, np.nan)

        #If more than one surface exists in the column calculate volumes
        if np.count_nonzero(np.isfinite(zvolup)) > 0:
            volups = np.append(0,zvolup[np.isfinite(zvolup)])
            volups = np.append(volups, Hvol)
            bins = ind[np.isfinite(zvolup)]
            bins = np.append(bins, bins[-1]+1)

            volumes = np.abs(np.diff(volups))
            
            output[bins] = volumes

        return output[:-1], output[-1]

    vol_bin_xint, vol_lastbin_xint = gufoo(zvolup_xint, Hvol_xint)
    vol_bin_yint, vol_lastbin_yint = gufoo(zvolup_yint, Hvol_yint)

    vol_bin_xint = da.ma.masked_invalid(da.append(vol_bin_xint, vol_lastbin_xint[None,:], axis=-2))
    vol_bin_yint = da.ma.masked_invalid(da.append(vol_bin_yint, vol_lastbin_yint[None,:], axis=-2))

    zsurf_xmean.attrs["region"] = region
    zsurf_ymean.attrs["region"] = region
    zsurf_xmean.name = "zsurf_xmean"
    zsurf_ymean.name = "zsurf_ymean"

    vol_bin_xint_cube = xr.DataArray(vol_bin_xint,
                                     dims=["rho_bin", "y"],
                                     coords={"rho_bin":np.arange(vol_bin_xint.shape[-2], dtype=int)},
                                     name="vol_bin_xint",
                                     attrs={"units":"m3", "region":region})

    vol_bin_yint_cube = xr.DataArray(vol_bin_yint,
                                     dims=["rho_bin", "x"],
                                     coords={"rho_bin":np.arange(vol_bin_yint.shape[-2], dtype=int)},
                                     name="vol_bin_yint",
                                     attrs={"units":"m3", "region":region}) 

    return zsurf_xmean, zsurf_ymean, vol_bin_xint_cube, vol_bin_yint_cube

def bin_by_initial_sf_zint(df_vent, namelist):
    import xarray as xr
    meta = ("sf_zint","float64")
    sf_zint_cube = xr.open_dataset(namelist["sf_zint_dir"] + "/sf_zint_t.nc", chunks="auto")["sf_zint_t"].compute()
    df_vent["sf_zint"] = df_vent.apply(sf_function, axis=1, meta=meta, sf_zint_cube=sf_zint_cube)

    return df_vent

def sf_function(df, sf_zint_cube=None):
    x = int(df["binnedx_i"] - 1) 
    y = int(df["binnedy_i"] - 1)
    output = float(sf_zint_cube[y,x].values)
    return output

def bin_by_bathy_depth(df_vent, namelist):
    import xarray as xr
    e3t_cube = xr.open_dataset(namelist["grid_dir"] + "/mesh_zgr.nc", chunks="auto")["e3t_0"].squeeze()
    tmask_cube = xr.open_dataset(namelist["grid_dir"] + "/mask.nc", chunks="auto")["tmask"].squeeze().astype(bool)
    
    depth_cube = e3t_cube.where(tmask_cube).sum("z", skipna=True).compute()
    df_vent["bathy_depth_i"] = df_vent.apply(bathy_depth_function, axis=1, meta=("bathy_depth_i", "float64"), depth_cube=depth_cube, xstr="binnedx_i", ystr="binnedy_i")
    df_vent["bathy_depth_o"] = df_vent.apply(bathy_depth_function, axis=1, meta=("bathy_depth_o", "float64"), depth_cube=depth_cube, xstr="binnedx_o", ystr="binnedy_o")

    return df_vent

def bathy_depth_function(df, depth_cube=None, xstr=None, ystr=None):
    x = int(df[xstr] - 1) 
    y = int(df[ystr] - 1) 
    output = float(depth_cube[y,x].values)
    return output

def bin_by_out_sigma(df_vent, namelist):
    """Bin densities at ventilation (sigma0) using the same bin intervals as the neutral density surfaces"""
    import xarray as xr
    import os
    ndsurf_dir = os.path.abspath(namelist["ndsurf_dir"])
    grid_dir = os.path.abspath(namelist["grid_dir"])
    ndsurf = xr.open_mfdataset(ndsurf_dir + "/*.nc", chunks="auto")
    sigma_bins = ndsurf.coords["sigma_ver"].to_numpy() - 1000.
    sigma_bins = np.append(sigma_bins, [+float("inf")])
    sigma_bins = np.append([-float("inf")], sigma_bins)
    labels = np.arange(len(sigma_bins)-1, dtype=int)
    df_vent["sigmabin_o"] = df_vent["density_o"].map_partitions(pd.cut, sigma_bins, labels=labels, retbins=False).astype(int)
    
    return df_vent

def subset_ij(M, imin, imax, jmin, jmax):
    if (imin is not None) and (imax is not None) and (imin > imax):
        iwrap = True
    else:
        iwrap = False

    if (jmin is not None) and (jmax is not None) and (jmin > jmax):
        jwrap = True
    else:
        jwrap = False

    if (iwrap == False) and (jwrap == False):
        output = M[..., jmin:jmax, imin:imax]
    
    elif (iwrap == True) and (jwrap == False):
        output = M[..., jmin:jmax, imin:]
        output = output[...,:imax]

    elif (iwrap == False) and (jwrap == True):
        output = M[..., jmin:, imin:imax]
        output = output[...,:jmax,:]

    else:
        output = M[..., jmin:, imin:]
        output = output[...,:jmax, :imax]

    return output








