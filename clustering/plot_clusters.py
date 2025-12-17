# %% Import packages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from netCDF4 import Dataset
import sys
sys.path.append('../utils/') # Path to scripts
import compute_statistics as stat
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, LAND
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.ndimage import minimum_filter, maximum_filter, gaussian_filter
# %% Input
path_root = "/home/565/fl2086/Australian-synoptic-weather-types/"
Nclusters = 40
date_start    = '1959-01-01' # Start year-month for clustering data yyyy-mm
date_end      = '2024-01-01' # End year-month for clustering data yyyy-mm (excluded)
nvec = 2
# %% Load cluster file
with Dataset(f"{path_root}clustering/clusters_{Nclusters}.nc",'r') as nc:
  time_clusters = np.array([np.datetime64('1900-01-01')+np.timedelta64(int(t),'h') for t in nc.variables['time'][:]])
  mask_clusters = (time_clusters>=np.datetime64(date_start)) & (time_clusters<np.datetime64(date_end))
  clusterSeries = nc.variables['clusterSeries'][mask_clusters]
  time_clusters = time_clusters[mask_clusters]
  if 'SWT' in nc.variables.keys():
      data['SWT'] = nc.variables['SWT'][:]
# %% Load required variables
varnames = ['u_300hPa','v_300hPa','u_850hPa','v_850hPa','pv_315K','pv_330K','msl','tcwv']
data = {}
for varname in varnames:
  with Dataset(f"{path_root}example_data/era5_data/{varname}.nc",'r') as nc:
    data['time']    = nc.variables['time'][:]
    data['latitude']     = nc.variables['latitude'][:]
    data['longitude']     = nc.variables['longitude'][:]
    temp            = nc.variables[varname.split('_')[0]][:]
    field_sum,field_sum2,field_count,time = stat.compute_sums(temp,data['time'],clusterSeries,data['time'],Nclusters)
    data[varname]   = stat.cluster_mean(field_sum,field_count)
mask_tropical_westerlies = (((data['latitude']>-20) & (data['latitude']<=0))[:,None]) & (data['u_850hPa']>0)
data['u_TW'] = np.where(mask_tropical_westerlies,data['u_850hPa'],np.nan)
data['v_TW'] = np.where(mask_tropical_westerlies,data['v_850hPa'],np.nan)
Ujet = np.sqrt(data['u_300hPa']**2+data['v_300hPa']**2)
data['jet']  = np.where(Ujet>20,Ujet,np.nan)
# %% Function to find local extremas
def extreme_vals(var, mode='wrap', window=50, sigma=4):
    def extrema(mat,mode='wrap',window=10):
        """find the indices of local extrema (min and max)
        in the input array."""
        mn = minimum_filter(mat, size=window, mode=mode)
        mx = maximum_filter(mat, size=window, mode=mode)
        return np.nonzero(mat == mn), np.nonzero(mat == mx)
    var = gaussian_filter(var, sigma)
    local_min, local_max = extrema(var, mode=mode, window=window)
    xlows = data['Lon'][local_min]; xhighs = data['Lon'][local_max]
    ylows = data['Lat'][local_min]; yhighs = data['Lat'][local_max]
    lowvals = var[local_min]; highvals = var[local_max]
    return xlows,ylows,lowvals,xhighs,yhighs,highvals
# %% Plot panels
# Create figure
Nrows = int(np.ceil(Nclusters/5))
Ncols = 5
fig = plt.figure(figsize=(8.3,11.7/8*Nrows),dpi=300)
gs = gridspec.GridSpec(Nrows, Ncols, figure=fig,hspace=0.2)
for i in range(Nrows):
  for j in range(Ncols):
    ax = fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())
    icluster = i*5+j
    # Plot land
    ax.add_feature(LAND,facecolor='lightgrey')
    ax.coastlines(linewidths=0.4)
    # Plot msl contours
    plot_levels = np.arange(900,1060,4)
    c=ax.contour(data['longitude'],data['latitude'],data['msl'][icluster,:,:]/100,
                  levels=plot_levels,colors='black',linewidths=0.75,linestyles='-',
                  transform=ccrs.PlateCarree())
    ax.clabel(c, inline=True, levels=plot_levels[::5],fontsize=5,fmt=lambda val: f'{val:.0f}hPa')
    # _,_,_,xhighs,yhighs,highvals=extreme_vals(data['msl'],data['latitude'],data['longitude'],window=70,sigma=4)
    data['Lon'], data['Lat'] = np.meshgrid(data['longitude'], data['latitude'])
    _,_,_,xhighs,yhighs,highvals=extreme_vals(data['msl'][icluster,:,:],window=70/6,sigma=4/6)
    # xlows,ylows,lowvals,_,_,_=extreme_vals(data['msl'],data['latitude'],data['longitude'],window=50,sigma=1)
    xlows,ylows,lowvals,_,_,_=extreme_vals(data['msl'][icluster,:,:],window=50/6,sigma=1/6)
    xyplotted = []
    dmin=2
    for x,y,p in zip(xlows, ylows, lowvals):
        if x < data['longitude'].max() and x > data['longitude'].min() and y < data['latitude'].max() and y > data['latitude'].min():
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                ax.text(x,y,'L',fontsize=8,fontweight='bold',
                        ha='center',va='center',color='b')
                xyplotted.append((x,y))
    for x,y,p in zip(xhighs, yhighs, highvals):
        if x < data['longitude'].max() and x > data['longitude'].min() and y < data['latitude'].max() and y > data['latitude'].min():
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                ax.text(x,y,'H',fontsize=8,fontweight='bold',
                        ha='center',va='center',color='r')
                xyplotted.append((x,y))
    # Plot jet
    plot_levels = np.arange(20,50,2.5)
    cf=ax.contourf(data['longitude'],data['latitude'],data['jet'][icluster,:,:],#*1.94384,
                    levels=plot_levels,cmap='Blues',extend='max',
                    transform=ccrs.PlateCarree(),alpha=0.75)
    # Plot total column vertically-integrated water vapour
    plot_levels = np.arange(48,1000,500)
    cf=ax.contourf(data['longitude'],data['latitude'],data['tcwv'][icluster,:,:],
                    levels=plot_levels,cmap='Greens',extend='max',
                    transform=ccrs.PlateCarree(),alpha=0.4)
    c=ax.contour(data['longitude'],data['latitude'],data['tcwv'][icluster,:,:],
                      levels=[48],colors='limegreen',linewidths=0.8,linestyles='-',
                      transform=ccrs.PlateCarree())
    ax.clabel(c, inline=True, levels=[48],fontsize=5,fmt=lambda val: f'{val:.0f}'+r'$kg\,m^{-2}$')
    # Plot PV contours
    c=ax.contour(data['longitude'],data['latitude'],data['pv_330K'][icluster,:,:]*1e6,
                      levels=[-2],colors='magenta',linewidths=1.0,linestyles='--',
                      transform=ccrs.PlateCarree())
    c=ax.contour(data['longitude'],data['latitude'],data['pv_315K'][icluster,:,:]*1e6,
                      levels=[-2],colors='magenta',linewidths=1.25,linestyles='-',
                      transform=ccrs.PlateCarree())
    # Plot tropical Westerlies
    q=ax.quiver(data['longitude'][::nvec],data['latitude'][::nvec], data['u_TW'][icluster,::nvec,::nvec], data['v_TW'][icluster,::nvec,::nvec],
                  scale=3,scale_units='xy',width=0.003,minshaft=2,color='darkgreen',# minlength=1,
                  transform=ccrs.PlateCarree())
    # Add gridlines
    gl = ax.gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
    ax.set_xlim([data['longitude'].min(),data['longitude'].max()])
    ax.set_ylim([data['latitude'].min(),data['latitude'].max()])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = i == int(np.ceil(Nclusters/5))-1 
    gl.left_labels = j == 0
    if 'SWT' in data.keys():
      ax.set_title(data['SWT'][icluster])
    else:
      ax.set_title(f"Cluster {icluster+1}")
plt.savefig(f"{path_root}clustering/clusters_{Nclusters}.png")
