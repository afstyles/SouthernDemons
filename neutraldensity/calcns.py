# CALCNS.py

#Calculate neutral density (omega) surfaces on the seeding dates of a TRACMASS experiment

import neutralocean
from neutralocean.grid.rectilinear import build_grid
from neutralocean.surface import omega_surf
from neutralocean.surface.isopycnal import potential_surf
from neutralocean.label import veronis_density
import xarray as xr
import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import cmocean

seed_dir = sys.argv[1]
month = str(sys.argv[2])
out_dir = sys.argv[3]
month = month[0].upper() + month[1:].lower()
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# User settings 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#Seeding area information (as defined in TRACMASS namelist) >>>>>>>>>>>>>>>
istmin = 1     # Minimum i-index seeded
istmax = 1440  # Maximum i-index seeded
jstmin = 1     # Minimum j-index seeded
jstmax = 300   # Maximum j-index seeded
iperio = True  # = True if i-index is periodic
jperio = False # = True if j-index is periodic

#Density surfaces to find >>>>>>>>>>>>>>>
ndepths  = 70     # Number of surfaces
depthmin = 100.   # Minimum depth in pinned cast

# #i and j indices for pinned column >>>>>>>>>>>>>>>
# ipin = 950
# jpin = 275

# i and j indices for pinned casts
# Order specifies priority for calculating surfaces
ipins = np.array([950, 1250, 100, 300, 500, 700 , 1100, 1300, 400, 900], dtype=int)
jpins = np.array([275, 275 , 275, 275, 100, 150 , 150 , 175 , 250, 250], dtype=int) 

# Density calculation parameters >>>>>>>>>>>>>>>
ref = float(2000.)  #Reference depth for potential density calculation
ver_ref = float(0.) #Reference depth for Veronis density calculation
eos = 'gsw'     #Equation of state ('gsw' for TEOS-10)

#Iterator settings >>>>>>>>>>>>>>>
ITER_MAX = 10  #Maximum number of iterations when calculating omega surface



if month not in months: 
    print(f"Error! {month} is not a valid month")
    print(f"Use one of")
    print(months)
    exit()

#Load the data and grid information
tdata = xr.open_dataset(seed_dir + "/" + month + "T.nc")
hgriddata = xr.open_dataset(seed_dir + "/mesh_hgr.nc" )
zgriddata = xr.open_dataset(seed_dir + "/mesh_zgr.nc" )
maskdata = xr.open_dataset(seed_dir + "/mask.nc")

so_cube = tdata["so"][0,...,jstmin-1:jstmax, istmin-1:istmax]
to_cube = tdata["thetao"][0,...,jstmin-1:jstmax, istmin-1:istmax]
deptht_cube = zgriddata["gdept_0"][0,...,jstmin-1:jstmax, istmin-1:istmax]

lon_cube = hgriddata["nav_lon"][jstmin-1:jstmax, istmin-1:istmax]
lat_cube = hgriddata["nav_lat"][jstmin-1:jstmax, istmin-1:istmax]
e1u_cube = hgriddata["e1u"][0,...,jstmin-1:jstmax, istmin-1:istmax]
e1v_cube = hgriddata["e1v"][0,...,jstmin-1:jstmax, istmin-1:istmax]
e2u_cube = hgriddata["e2u"][0,...,jstmin-1:jstmax, istmin-1:istmax]
e2v_cube = hgriddata["e2v"][0,...,jstmin-1:jstmax, istmin-1:istmax]
e1t_cube = hgriddata["e1t"][0,...,jstmin-1:jstmax, istmin-1:istmax]
e2t_cube = hgriddata["e2t"][0,...,jstmin-1:jstmax, istmin-1:istmax]
e3w_cube = zgriddata["e3w_0"][0,...,jstmin-1:jstmax, istmin-1:istmax]
tmask_cube = maskdata["tmask"][0,...,jstmin-1:jstmax, istmin-1:istmax]
tmask2d_cube = maskdata["tmaskutil"][0,...]#,jstmin-1:jstmax, istmin-1:istmax]
depth_cube = ( e3w_cube * tmask_cube ).sum(axis=-3)


#Set values of salinity and temperature beneath the sea floor
so_cube = so_cube.where(so_cube > 0.1, np.nan)
to_cube = to_cube.where(so_cube > 0.1, np.nan)

#Construct the neutralocean grid dictionary
grid = build_grid( (so_cube.shape[-2],so_cube.shape[-1]), (jperio, iperio), dyC=e1u_cube.to_numpy(), dxC=e2v_cube.to_numpy(), dyG=e1v_cube.to_numpy(), dxG=e2u_cube.to_numpy()  )
# grid = build_grid( (so_cube.shape[-2],so_cube.shape[-1]), (jperio, iperio), dxC=1, dyC=1, dxG=1, dyG=1  )

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# End of user settings
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

print("Calculating surfaces >>>>>>")

nsurf = np.zeros(lat_cube.shape, dtype=int)
rho_vermin = 0.
rho_vermax = 0.
rho_verthr = 0.
rho_ver_list = []

for ipin, jpin in zip(ipins, jpins):
    depths = np.linspace(depthmin, depth_cube[jpin, ipin].values-10, num=ndepths, dtype=float )
    pin_lon = lon_cube[jpin, ipin].values
    pin_lat = lat_cube[jpin, ipin].values
    print(f"(x) ipin = {ipin}, jpin = {jpin}")
    print(f"pin lon/lat = {pin_lon:.2f} degE,{pin_lat:.2f} degN")
    print(f"Max depth = {depths[-1]:.2f}m")


    for depth in depths:
        
        print("<<< Depth = ", depth, ">>>")
 
        # Calculate the Veronis density in the pinned cast
        rho_ver = 1/veronis_density(so_cube.values[:,jpin,ipin], to_cube.values[:,jpin,ipin], deptht_cube.values[:,jpin,ipin], depth, eos=eos, p_ref=ver_ref)
        print("rho_ver = ", rho_ver)


        if ( rho_ver >= rho_vermin - rho_verthr ) and ( rho_ver <= rho_vermax + rho_verthr ): continue  #If within an already explored veronis range, move to the next depth
        
        rho_ver_list = rho_ver_list + [rho_ver]

        # Calculate the initial potential density surface
        s,t,z,d = potential_surf(so_cube, to_cube, deptht_cube, grid=grid, vert_dim="deptht", ref=ref, pin_cast=(jpin,ipin), pin_p=depth, eos=eos)
        rho_ref = 1/d["isoval"] #Save the potential density value
        print("ref = ", ref)
        print("rho_ref = ", rho_ref)

       

        # Calculate the omega surface is ITER_MAX > 1
        if ITER_MAX >= 1:
            s,t,z,d = omega_surf(so_cube, to_cube, deptht_cube, grid, vert_dim="deptht", ref=ref, pin_cast=(jpin,ipin), pin_p=depth, eos=eos, ITER_MAX=ITER_MAX)
        
        # Count how many surfaces are in a given fluid column (Grounded surfaces are within 1m of the sea floor)
        nsurf = nsurf + (z < depth_cube - 1.).astype(int).values

        # Mask grounded surfaces (within 1m of the sea floor)
        z = z.where((z < depth_cube - 1))

        #Plot the surface depth
        plt.figure(dpi=200)
        plt.pcolormesh(tmask2d_cube.values, cmap=cmocean.cm.gray, vmin=0, vmax=1)
        plt.pcolormesh(z, cmap=cmocean.cm.haline_r, vmin=0, vmax=depth_cube.max())
        plt.colorbar()
        plt.ylabel(r"y index [-]")
        plt.xlabel(r"x index [-]")
        plt.title(fr"Surface depth, $\rho_v$ = {np.round(rho_ver,2)} kg m$^3$, $z_p$ = {int(depth)} m")
        plt.contour(depth_cube, colors='k', linewidths=0.1)
        plt.scatter([ipin], [jpin], marker='x', color='red')
        plt.savefig(out_dir + "/" + f"ns_{month}_depth{int(depth)}_ipin{ipin}_jpin_{jpin}_itermax{ITER_MAX}.png")
        plt.close()

        #Store surface as DataArray (xarray)
        shared_attrs = {"ipin":ipin, "jpin":jpin, "pin_lon":pin_lon, "pin_lat":pin_lat, "month":month, "ref": ref, "eos":eos, "ITER_MAX":ITER_MAX}

        zarray = xr.DataArray(data=z, dims=['y', 'x'], name="ns_depth", attrs={**shared_attrs, "units":"m", "long_name":"isoneutral_depth"})
        zarray = zarray.expand_dims(dim={"sigma_ref":np.array([rho_ref])})
        zarray = zarray.expand_dims(dim={"sigma_ver":np.array([rho_ver])})

        tarray = xr.DataArray(data=t, dims=['y', 'x'], name="ns_temp", attrs={**shared_attrs, "units":"degC", "long_name":"isoneutral_temp"})
        tarray = tarray.expand_dims(dim={"sigma_ref":np.array([rho_ref])})
        tarray = tarray.expand_dims(dim={"sigma_ver":np.array([rho_ver])})

        sarray = xr.DataArray(data=s, dims=['y', 'x'], name="ns_sal", attrs={**shared_attrs, "units":"g/kg", "long_name":"isoneutral_sal"})
        sarray = sarray.expand_dims(dim={"sigma_ref":np.array([rho_ref])})
        sarray = sarray.expand_dims(dim={"sigma_ver":np.array([rho_ver])})

        #Save DataSet to Netcdf
        xr.merge([zarray, tarray, sarray]).to_netcdf(out_dir + "/" + f"ns_{month}_depth{int(depth)}_ipin{ipin}_jpin_{jpin}_itermax{ITER_MAX}.nc")


    rho_vermin = np.min(rho_ver_list)
    rho_vermax = np.max(rho_ver_list)

    print("Min: ", rho_vermin)
    print("Max: ", rho_vermax)

print("rho_ver_list >>>")
print(rho_ver_list)

print("")

print("rho_ver_list (sorted)>>>")
print(sorted(rho_ver_list))

#Plot surface coverage as well
nsurf = np.ma.masked_less(nsurf, 1)
plt.figure(dpi=200)
plt.ylabel(r"y index [-]")
plt.xlabel(r"x index [-]")
plt.title(r"Number of surfaces in ocean column [-]")
plt.pcolormesh(tmask2d_cube.values, cmap=cmocean.cm.gray, vmin=0, vmax=1)
plt.pcolormesh(nsurf, cmap=cmocean.cm.haline_r, vmin=0)
plt.colorbar()
plt.contour(depth_cube, colors='k', linewidths=0.1)
plt.scatter(ipins, jpins, marker='x', color='red')
plt.savefig(out_dir + "/" + f"nsurf_{month}_itermax{ITER_MAX}.png")
plt.close()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# END
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>