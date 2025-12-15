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
from datetime import datetime,timedelta
# %% Input
path_root = "/home/565/fl2086/Australian-synoptic-weather-types/"
Nclusters = 30
date_start    = datetime(1959,1,1) # Start year-month for clustering data yyyy-mm
date_end      = datetime(2024,1,1) # End year-month for clustering data yyyy-mm (excluded)
SWT_plots     = ['WH-A','AM-A']
months        = [12,1,2]
nvec = 2
month_labels = {1:'J',2:'F',3:'M',4:'A',5:'M',6:'J',7:'J',8:'A',9:'S',10:'O',11:'N',12:'D'}
# Function to find local extremas
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
# Loop over clusters to plot
fig = plt.figure(figsize=(8.3,11.7/8*5/2),dpi=300)
gs = gridspec.GridSpec(1, len(SWT_plots), figure=fig,hspace=0.2)
for iplot,SWT_plot in enumerate(SWT_plots):
  # Load cluster file
  with Dataset(f"{path_root}clustering/clusters_{Nclusters}.nc",'r') as nc:
    time_clusters = np.array([datetime(1900,1,1)+timedelta(hours=int(t)) for t in nc.variables['time'][:]])
    mask_clusters = (time_clusters>=date_start) & (time_clusters<date_end)
    icluster = np.argwhere(nc.variables['SWT'][:]==SWT_plot).flatten()[0]
    clusterSeries = nc.variables['clusterSeries'][mask_clusters]
    time_clusters = time_clusters[mask_clusters]
    mask_summer   = [t.month in months for t in time_clusters]
    mask_cluster  = clusterSeries == icluster+1
    mask_composite = (mask_cluster & mask_summer).flatten()
    clusterSeries = clusterSeries[mask_composite]
    time_clusters = time_clusters[mask_composite]
  # Load required variables
  varnames = ['u_300hPa','v_300hPa','u_850hPa','v_850hPa','pv_315K','pv_330K','msl','tcwv']
  data = {}
  for varname in varnames:
    with Dataset(f"{path_root}example_data/era5_data/{varname}.nc",'r') as nc:
      data['time']      = nc.variables['time'][mask_composite]
      data['latitude']  = nc.variables['latitude'][:]
      data['longitude'] = nc.variables['longitude'][:]
      data[varname]     = np.nanmean(nc.variables[varname.split('_')[0]][mask_composite,:,:],axis=0)
  mask_tropical_westerlies = (((data['latitude']>-20) & (data['latitude']<=0))[:,None]) & (data['u_850hPa']>0)
  data['u_TW'] = np.where(mask_tropical_westerlies,data['u_850hPa'],np.nan)
  data['v_TW'] = np.where(mask_tropical_westerlies,data['v_850hPa'],np.nan)
  Ujet = np.sqrt(data['u_300hPa']**2+data['v_300hPa']**2)
  data['jet']  = np.where(Ujet>20,Ujet,np.nan)
  # Plot panels
  # Create figure
  ax = fig.add_subplot(gs[iplot], projection=ccrs.PlateCarree())
  # Plot land
  ax.add_feature(LAND,facecolor='lightgrey')
  ax.coastlines(linewidths=0.4)
  # Plot msl contours
  plot_levels = np.arange(900,1060,4)
  c=ax.contour(data['longitude'],data['latitude'],data['msl'][:,:]/100,
                levels=plot_levels,colors='black',linewidths=0.75,linestyles='-',
                transform=ccrs.PlateCarree())
  ax.clabel(c, inline=True, levels=plot_levels[::5],fontsize=5,fmt=lambda val: f'{val:.0f}hPa')
  # _,_,_,xhighs,yhighs,highvals=extreme_vals(data['msl'],data['latitude'],data['longitude'],window=70,sigma=4)
  data['Lon'], data['Lat'] = np.meshgrid(data['longitude'], data['latitude'])
  _,_,_,xhighs,yhighs,highvals=extreme_vals(data['msl'][:,:],window=70/6,sigma=4/6)
  # xlows,ylows,lowvals,_,_,_=extreme_vals(data['msl'],data['latitude'],data['longitude'],window=50,sigma=1)
  xlows,ylows,lowvals,_,_,_=extreme_vals(data['msl'][:,:],window=50/6,sigma=1/6)
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
  cf=ax.contourf(data['longitude'],data['latitude'],data['jet'][:,:],#*1.94384,
                  levels=plot_levels,cmap='Blues',extend='max',
                  transform=ccrs.PlateCarree(),alpha=0.75)
  # Plot total column vertically-integrated water vapour
  plot_levels = np.arange(48,1000,500)
  cf=ax.contourf(data['longitude'],data['latitude'],data['tcwv'][:,:],
                  levels=plot_levels,cmap='Greens',extend='max',
                  transform=ccrs.PlateCarree(),alpha=0.4)
  c=ax.contour(data['longitude'],data['latitude'],data['tcwv'][:,:],
                    levels=[48],colors='limegreen',linewidths=0.8,linestyles='-',
                    transform=ccrs.PlateCarree())
  ax.clabel(c, inline=True, levels=[48],fontsize=5,fmt=lambda val: f'{val:.0f}'+r'$kg\,m^{-2}$')
  # Plot PV contours
  c=ax.contour(data['longitude'],data['latitude'],data['pv_330K'][:,:]*1e6,
                    levels=[-2],colors='magenta',linewidths=1.0,linestyles='--',
                    transform=ccrs.PlateCarree())
  c=ax.contour(data['longitude'],data['latitude'],data['pv_315K'][:,:]*1e6,
                    levels=[-2],colors='magenta',linewidths=1.25,linestyles='-',
                    transform=ccrs.PlateCarree())
  # Plot tropical Westerlies
  q=ax.quiver(data['longitude'][::nvec],data['latitude'][::nvec], data['u_TW'][::nvec,::nvec], data['v_TW'][::nvec,::nvec],
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
  gl.bottom_labels = True
  gl.left_labels = iplot == 0
  title = f"{SWT_plot} for "
  savename = ''
  for month in months:
    title+=month_labels[month]
    savename+=month_labels[month]
  ax.set_title(f"{title}, {np.sum(mask_composite):.0f} days\n{np.sum(mask_composite)/np.sum(mask_cluster)*100:.0f}% of cluster data")
plt.savefig(f"{path_root}clustering/{savename}_composites.png")