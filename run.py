'''
Script to run drifters from Galveston Bay on the day of the TCY oil spill
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
from datetime import datetime, timedelta
import glob
from tracpy.tracpy_class import Tracpy


loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

def init(name):
    '''
    Initialization for the simulation.
    '''


    time_units = 'seconds since 1970-01-01'

    # horizontal_diffusivity project showed that relative dispersion did not
    # change between nsteps=25 and 50, but does between nsteps=5 and 25, and
    # interim numbers have not been tested yet.
    nsteps = 25  # in-between tracks: 12 # old tracks: 25 

    # Number of steps to divide model output for outputting drifter location
    N = 5

    # Number of days
    ndays = 10

    # This is a forward-moving simulation
    ff = 1 

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 5. # old tracks: 5.
    av = 0. # m^2/s

    # surface drifters
    z0 = 's'
    zpar = 29

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 2

    # Flag for streamlines.
    dostream = 0

    # Initialize Tracpy class
    tp = Tracpy(loc, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps, dostream=dostream, savell=False, doperiodic=0, 
                N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, 
                time_units=time_units, usebasemap=True, grid=grid)

    llcrnrlon = -94.768091; urcrnrlon = -94.670588; llcrnrlat = 29.299671; urcrnrlat = 29.400220; # New
    xcrnrs, ycrnrs = grid['basemap']([llcrnrlon, urcrnrlon], [llcrnrlat, urcrnrlat])
    X, Y = np.meshgrid(np.arange(xcrnrs[0], xcrnrs[1], 700), 
                        np.arange(ycrnrs[0], ycrnrs[1], 700))
    lon0, lat0 = grid['basemap'](X, Y, inverse=True)
    # Eliminate points that are outside domain or in masked areas
    lon0, lat0 = tracpy.tools.check_points(lon0, lat0, grid)

    # # save starting locations for future use
    # np.savez('calcs/seeds.npz', lon0=lon0, lat0=lat0)

    return tp, lon0, lat0


def run():

    # CST to GMT. started at 12:30pm CST
    overallstartdate = datetime(2014, 3, 22, 16, 0)
    overallstopdate = datetime(2014, 3, 23, 12, 0)

    date = overallstartdate

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')

    # loop through state dates
    while date < overallstopdate:

        name = date.isoformat()[0:13]

        # Read in simulation initialization
        tp, lon0, lat0 = init(name)

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc') and \
            not os.path.exists('tracks/' + name + 'gc.nc'):

            # Run tracpy
            # Save directly to grid coordinates
            lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)

        # Increment by 4 hours for next loop, to move through more quickly
        date = date + timedelta(hours=4)


if __name__ == "__main__":
    run()    
