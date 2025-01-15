# CALCSF.py
#Calculate the depth-integrated stream function on the seeding dates of a TRACMASS experiment

import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cmocean
import sys

seed_dir = sys.argv[1]
month = str(sys.argv[2])
out_dir = sys.argv[3]
month = month[0].upper() + month[1:].lower()
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

#Seeding area information (as defined in TRACMASS namelist) >>>>>>>>>>>>>>>
istmin = 1     # Minimum i-index seeded
istmax = 1440  # Maximum i-index seeded
jstmin = 1     # Minimum j-index seeded
jstmax = 398   # Maximum j-index seeded
iperio = True  # = True if i-index is periodic
jperio = False # = True if j-index is periodic

udata = xr.open_dataset(seed_dir + "/" + month + "U.nc")
hgriddata = xr.open_dataset(seed_dir + "/mesh_hgr.nc" )
zgriddata = xr.open_dataset(seed_dir + "/mesh_zgr.nc" )
maskdata = xr.open_dataset(seed_dir + "/mask.nc")

u_cube = udata["uo"][0,...,jstmin-1:jstmax, istmin-1:istmax]
e1u = hgriddata["e1u"][0,jstmin-1:jstmax, istmin-1:istmax].data
e2u = hgriddata["e2u"][0,jstmin-1:jstmax, istmin-1:istmax].data
e1v = hgriddata["e1v"][0,jstmin-1:jstmax, istmin-1:istmax].data
e2v = hgriddata["e2v"][0,jstmin-1:jstmax, istmin-1:istmax].data
e1t = hgriddata["e1t"][0,jstmin-1:jstmax, istmin-1:istmax].data
e2t = hgriddata["e2t"][0,jstmin-1:jstmax, istmin-1:istmax].data
e1f = hgriddata["e1f"][0,jstmin-1:jstmax, istmin-1:istmax].data
e2f = hgriddata["e2f"][0,jstmin-1:jstmax, istmin-1:istmax].data

e3u = zgriddata["e3u_0"][0,...,jstmin-1:jstmax, istmin-1:istmax].data
umask = maskdata["umask"][0,...,jstmin-1:jstmax, istmin-1:istmax].astype(bool).data
tmask2d = maskdata["tmaskutil"][0,jstmin-1:jstmax, istmin-1:istmax].astype(bool).data

u_zint = (u_cube.where(umask) * e3u).sum(dim="depthu", skipna=True)
sf_zint = -(u_zint * e2u).cumsum(dim="y") #Calculate stream function centred on F-Points

sf_zint_v = (sf_zint * e1f * e2f + (sf_zint * e1f * e2f).roll(x=1))/(2*e1v*e2v)
sf_zint_t = (sf_zint_v * e1v * e2v + (sf_zint_v * e1v * e2v).roll(y=1))/(2*e1t*e2t)
sf_zint_t = sf_zint_t.where(tmask2d)/1e6

plt.figure(dpi=200)
plt.pcolormesh(tmask2d, cmap=cmocean.cm.gray, vmin=0, vmax=1)
plt.pcolormesh(sf_zint_t, cmap=cmocean.cm.balance, norm=matplotlib.colors.TwoSlopeNorm(vcenter=0.))
plt.colorbar()
plt.ylabel(r"y index [-]")
plt.xlabel(r"x index [-]")
plt.title(fr"Depth-integrated stream function (Sv)")
plt.savefig(out_dir + "/" + f"sf_zint_t_{month}.png")
plt.close()

sf_zint_t.attrs["units"] = "Sv"
sf_zint_t.name = "sf_zint_t"
sf_zint_t.to_netcdf(out_dir + "/" + f"sf_zint_t.nc")