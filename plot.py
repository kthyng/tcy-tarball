'''
Show trajectories.
'''


import numpy as np
import matplotlib.pyplot as plt
import tracpy
import tracpy.plotting
from glob import glob
import netCDF4 as netCDF


loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

Files = glob('tracks/*')

# Plot tracks
fig = plt.figure()
ax = fig.add_subplot(111)
tracpy.plotting.background(grid, ax=ax)
for File in Files:
    d = netCDF.Dataset(File)
    xg = d.variables['xg'][:]
    yg = d.variables['yg'][:]
    xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
    ind = xg==-1  # in case drifters exit domain
    xp[ind] = np.nan
    yp[ind] = np.nan

    ax.plot(xp.T, yp.T, '0.2', alpha=0.5)

# 
