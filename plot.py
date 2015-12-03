'''
Show trajectories.
'''


import numpy as np
import matplotlib.pyplot as plt
import tracpy
import tracpy.plotting
from glob import glob
import netCDF4 as netCDF
import os
from matplotlib import colors
import matplotlib.patches as Patches
import cmocean
import matplotlib as mpl
from matplotlib.path import Path

mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True, llcrnrlon=-97.8, llcrnrlat=27, urcrnrlon=-94.4, urcrnrlat=29.8)

dirs = glob('tracks/best/ah_20_3p_3days/*/')
index = 3
Files = glob(dirs[index] + '*.nc')
name = dirs[index].split('/')[-2]
# Files = glob('tracks/best/ah_20_3p_3days/deflection0/*.nc')
# Files = glob('tracks/*.nc')


# Find drifters that reach boxes
paths = np.load('../shelf_transport/calcs/coastpaths.npz')
pathsg = paths['pathsg']
days = np.linspace(0, 10, 901)

# convert paths from grid space to projected xy space for this small grid
pathsxy = []
for path in pathsg:
    xtemp, ytemp, _ = tracpy.tools.interpolate2d(path.vertices[:, 0], path.vertices[:, 1], grid, 'm_ij2xy')
    verts = np.vstack((xtemp, ytemp)).T
    pathsxy.append(Path(verts))

ncrosses = 5

# Plot tracks
fig = plt.figure()
ax = fig.add_subplot(111)
tracpy.plotting.background(grid, ax=ax)
for File in Files:
    if str(16) in File:  # skip first simulations
        continue
    data = netCDF.Dataset(File)
    xg = data.variables['xg'][:]
    yg = data.variables['yg'][:]
    xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy')
    ind = xg==-1  # in case drifters exit domain
    xp[ind] = np.nan
    yp[ind] = np.nan
    ax.plot(xp.T, yp.T, '0.45', lw=0.05, alpha=0.5)

    inboxfname = 'calcs/inbox-' + name + '-' + File.split('/')[-1].split('.')[0] + '.npz'  # one for each tracks file
    if not os.path.exists(inboxfname):
        inbox = np.ones((len(pathsg), xg.shape[0], ncrosses))*np.nan # to store analysis. 
        outbox = np.ones((len(pathsg), xg.shape[0], ncrosses))*np.nan # to store analysis. 
        for i, path in enumerate(pathsg):
            # which points are inside the regions
            inside = path.contains_points(np.vstack((xg.flat, yg.flat)).T).reshape(xg.shape)
            whencross = np.diff(inside.astype(int), axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            whencross_sorted = np.sort(whencross, axis=1) # will pick out when the drifters change from outside to inside region or vice versa
            isort = np.argsort(whencross, axis=1)
            inbox[i,:,:] = days[isort[:,:ncrosses]] # allow for up to ncrosses re-entrances
            # nan out the entries that aren't actually entering
            iin = whencross_sorted[:,:ncrosses]!=-1
            inbox[i,iin] = np.nan
            inbox[i,:,:] = np.sort(inbox[i,:,:], axis=1)
            outbox[i,:,:] = days[isort[:,-ncrosses:]] # allow for up to ncrosses re-exits
            # nan out the exits that aren't actually exiting
            iout = whencross_sorted[:,-ncrosses:]!=1
            outbox[i,iout] = np.nan
            outbox[i,:,:] = np.sort(outbox[i,:,:], axis=1)
        np.savez(inboxfname, inbox=inbox, outbox=outbox)


# Loop through along-coast boxes to find which other boxes they are connected to
# and calculate times
pts = np.load('../shelf_transport/calcs/alongcoastconn/inds-in-coast-paths.npz')['pts']
ndays = 10
ndboxfname = 'calcs/ndbox-' + name + '.npz'
times = []
for File in Files:
    if str(16) in File:  # skip first simulations
        continue
    inboxfname = 'calcs/inbox-' + name + '-' + File.split('/')[-1].split('.')[0] + '.npz'
    inbox = np.load(inboxfname)['inbox']
    # number of drifters reaching each coast box (yes or no, no weighting) in day days
    ndbox = np.zeros((len(pts)))
    # Drifters that enter a coast box within day days [coast box x set of drifters]
    ind = (inbox[:,:,0]<=ndays)
    # How many drifters enter each box by day days?
    ndbox += ind.sum(axis=1)
    # arrival times for Mustang Island
    for i in np.arange(104, 114):
        inds = ~np.isnan(inbox[i, :, 0])
        times.extend(inbox[i, inds, 0])
    np.savez(ndboxfname, ndbox=ndbox)

# Plot coastline oiling
# Add on vulnerability of coastline
# Need to plot the coast boxes as patches and color them according to vulnerability level
# http://matplotlib.org/1.2.1/examples/pylab_examples/hist_colormapped.html
# we need to normalize the data to 0..1 for the full
# range of the colormap
fracs = ndbox.astype(float)/ndbox.max()  # max across domain
# norm = colors.LogNorm(fracs.min()+0.00001, fracs.max())
norm = colors.Normalize(fracs.min(), fracs.max())
# Save patches together
patches = []
for path in pathsxy:
    patches.append(Patches.PathPatch(path, facecolor='orange', lw=0, edgecolor=None, zorder=5))
# assign shades of colormap to the patches according to values, and plot
for thisfrac, thispatch in zip(fracs, patches):
    color = cmocean.cm.turb(norm(thisfrac))
    thispatch.set_facecolor(color)
    ax.add_patch(thispatch)

# Plot histogram of oil arrival times to Mustang Island
# use boxes at Mustang island to make time histogram
# add subaxis for histogram
ax1 = fig.add_axes([0.475, 0.17, 0.35, 0.2])
# ax1 = fig.add_axes([0.28, 0.677, 0.35, 0.2])
ax1.hist(times, bins=9, range=(6, 8.25), color='0.3')
ax1.tick_params(axis='both', which='major', labelsize=11, bottom='off', top='off', left='off', right='off')
ax1.set_xlabel('Time to reach Mustang Island [days]', fontsize=11, labelpad=0.05)
ax1.set_ylabel('Number of drifters', fontsize=11)

# make frame red like mustand island outline
for spine in ax1.spines.values():
    spine.set_edgecolor('red')
    spine.set_lw(2)


# Plot mustand island boxes special
xs = np.hstack((pathsxy[104].vertices[2:, 0], pathsxy[113].vertices[:2, 0], pathsxy[104].vertices[2, 0]))
ys = np.hstack((pathsxy[104].vertices[2:, 1], pathsxy[113].vertices[:2, 1], pathsxy[104].vertices[2, 1]))
ax.plot(xs, ys, 'r-', lw=1, zorder=6)

# Box where drifters were initialized
llcrnrlon = -94.738918; urcrnrlon = -94.627218; llcrnrlat = 29.332186; urcrnrlat = 29.368367; # New
xll, yll = grid['basemap'](llcrnrlon, llcrnrlat)
xlr, ylr = grid['basemap'](urcrnrlon, llcrnrlat)
xur, yur = grid['basemap'](urcrnrlon, urcrnrlat)
xul, yul = grid['basemap'](llcrnrlon, urcrnrlat)
ax.plot([xll, xlr, xur, xul, xll], [yll, ylr, yur, yul, yll], 'g-', lw=2, alpha=0.8, zorder=6)

# label locations from text
ax.text(0.65, 0.9, 'Galveston\nBay', transform=ax.transAxes, fontsize=11)
ax.text(0.3, 0.65, 'Matagorda Bay', transform=ax.transAxes, fontsize=11)
ax.text(0.18, 0.525, 'Matagorda\nIsland', transform=ax.transAxes, fontsize=11)
ax.arrow(99000, 160000, 10000, -15000, head_width=6000, head_length=8000, fc='0.1', ec='0.1', zorder=6)
ax.text(0.03, 0.35, 'Mustang\nIsland', transform=ax.transAxes, fontsize=11)
ax.arrow(33200, 102500, 20000, -18000, head_width=6000, head_length=8000, fc='0.1', ec='0.1', zorder=6)

fig.savefig('figures/tracks-' + name + '.png', bbox_inches='tight')#, dpi=300)


# # # Plot up
# # fig = plt.figure()
# # ax = fig.add_subplot(111)
# # tracpy.plotting.background(grid, ax=ax)
# # plt.plot(grid['xpsi'], grid['ypsi'], 'k', grid['xpsi'].T, grid['ypsi'].T, 'k')
# # ind = grid['mask'].astype(bool)
# # plt.plot(grid['xr'][ind], grid['yr'][ind], 'bs')
# # # for path in pathsxy:
# # #     patch = patches.PathPatch(path, facecolor='orange', lw=2, zorder=10)
# # #     ax.add_patch(patch)
# # # Mustang Island region
# # for path in pathsxy[104:114]:
# #     patch = Patches.PathPatch(path, facecolor='red', lw=2, zorder=10)
# #     ax.add_patch(patch)
