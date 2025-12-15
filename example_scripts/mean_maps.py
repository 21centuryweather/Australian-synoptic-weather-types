# %% Import packages
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature, LAND
from netCDF4 import Dataset
import numpy as np
from datetime import datetime
from matplotlib.gridspec import GridSpec
import sys
sys.path.append('../utils/')
from read_era5 import read_data
import compute_statistics as stat
# %% Input
varnames     = ['U10']
data_types   = {'u10':'single-levels','v10':'single-levels','t2m':'single-levels','tp': 'single-levels'}
titles       = {'U10': '10 meter wind speed', 't2m':'2 meter temperature', 'tp': 'total surface precipitation'}
units        = {'U10': 'm/s', 't2m':'K', 'tp':'mm'}
cmaps        = {'U10': 'Blues', 't2m':'inferno', 'tp':'Blues'}
plot_limits  = {'U10': [0,15], 't2m':[275,305], 'tp':[0,20]}
varnames_path= {'u10': '10u', 'v10': '10v', 't2m': '2t', 'tp':'tp'}
scales       = {'tp':1000}
root         = '/home/565/fl2086/Australian-synoptic-weather-types/'
date_start   = '1952-01' # Start year-month yyyy-mm
date_start   = '2020-01' # Start year-month yyyy-mm
date_end     = '2023-12' # End year-month yyyy-mm (included)
utc          = 12        # Reading sampling time (utc) data is sampled daily
plevel       = 850       # Height level used (hPa)
lat_lims     = [-5,-50]  # South and North lattitude limit of analysis box
lon_lims     = [100,165] # West and East longitude limit of analysis box
Ncoarsen     = 6         # Coarsening factor in lat and lon direction
path_data    = lambda data_type : f"/g/data/rt52/era5/{data_type}/reanalysis/" # Era5 data directory
path_out     = f"{root}results/"
cluster_file = f"{root}SWT_fields/SWT_data_v1.nc"
WR_names     = ['AM','COL','WCT','FH','TH','EH','CH','WH'] # Order of plotting WR
Ncols        = 5         # Number of columns in plot map
Nlevels      = 11        # Number of colour levels for contourf
# %% Load cluster information
with Dataset(cluster_file,'r') as nc:
  clusterSeries = nc.variables['clusterSeries'][:]
  clusterID     = nc.variables['clusterID'][:]
  SWT           = nc.variables['SWT'][:]
  WR            = nc.variables['WR'][:]
  clusterTime   = nc.variables['time'][:]
  Nclusters     = len(clusterID)
# Sort cluster statistics
idx_sort = []
idx = np.arange(Nclusters)
for wr in WR_names[::-1]:
    mask = np.char.startswith(SWT.astype(str), wr) # Find indices where SWT starts with WR name
    wr_indices = idx[mask] # Extract indices
    sorted_wr_indices = wr_indices[np.argsort(SWT[wr_indices])] # Sort them alphabetically by SWT suffix
    idx_sort.extend(sorted_wr_indices) # Append to sorting indices
SWT = SWT[idx_sort]
# %% Read era5 data
for varname in varnames:
  if varname == 'U10':
    print('Reading u10')
    u,_,_,_=read_data('u10',date_start,date_end,utc,lat_lims,lon_lims,path_data(data_types['u10']),varname_path=varnames_path['u10'],Ncoarsen=Ncoarsen,level=None,progress=True,save=False)
    print('Reading v10')
    v,time,lat,lon=read_data('v10',date_start,date_end,utc,lat_lims,lon_lims,path_data(data_types['v10']),varname_path=varnames_path['v10'],Ncoarsen=Ncoarsen,level=None,progress=True,save=False)
    data = np.sqrt(u**2+v**2)
  else:
    print(f"Reading {varname}")
    data,time,lat,lon=read_data(varname,date_start,date_end,utc,lat_lims,lon_lims,path_data(data_types[varname]),varname_path=varnames_path[varname],Ncoarsen=Ncoarsen,level=plevel,progress=True,save=False)
    if varname in scales.keys():
      data*=scales[varname]
  # Get cluster statistics
  field_sum,field_sum2,field_count,time_sum=stat.compute_sums(data,time,clusterSeries,clusterTime,Nclusters,date0=datetime(1900,1,1))
  cluster_mean        = stat.cluster_mean(field_sum,field_count)
  cluster_var         = stat.cluster_var(field_sum,field_sum2,field_count)
  cluster_daily_pert  = stat.cluster_daily_pert(field_sum,field_count)
  cluster_mean       = cluster_mean[idx_sort,:,:]
  cluster_var        = cluster_var[idx_sort,:,:]
  cluster_daily_pert = cluster_daily_pert[idx_sort,:,:]
# %% Plot cluster mean map
Nrows = int(np.ceil(Nclusters/Ncols)) # Rows required for plot
fig = plt.figure(figsize=(8.3,8))
fig.suptitle(f"Cluster mean {titles[varname]}",y=0.93,fontsize=14)
gs = GridSpec(Nrows, Ncols+1, width_ratios=[4]*Ncols+[1],hspace=0)
plot_levels = np.linspace(np.min(plot_limits[varname]),np.max(plot_limits[varname]),Nlevels,endpoint=True)
iplot = 0
for i in range(Nrows):
  for j in range(Ncols):
    if iplot>=Nclusters:
      continue
    ax = fig.add_subplot(gs[i,j],projection=ccrs.PlateCarree())
    cf=ax.contourf(lon,lat,cluster_mean[iplot,:,:],levels=plot_levels,cmap=cmaps[varname],transform=ccrs.PlateCarree())
    ax.set_title(SWT[iplot])
    gl = ax.gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
    ax.set_xlim([lon.min(),lon.max()])
    ax.set_ylim([lat.min(),lat.max()])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = (i == Nrows - 1)
    gl.left_labels = (j == 0)
    iplot+=1
ax_cbar = fig.add_subplot(gs[0:2,-1])
cax = inset_axes(ax_cbar, width="100%", height="80%", loc="center") 
plt.colorbar(cf,cax=cax,label=f"{varname} ({units[varname]})")
ax_cbar.axis("off")
plt.savefig(f"{path_out}{varname}_cluster_mean.png")

# %% Plot cluster mean daily perturbation map
fig = plt.figure(figsize=(8.3,8))
fig.suptitle(f"Cluster mean daily {titles[varname]} perturbation",y=0.93,fontsize=14)
gs = GridSpec(Nrows, Ncols+1, width_ratios=[4]*Ncols+[1],hspace=0)
plot_levels = np.linspace(-np.round(np.max(np.abs(cluster_daily_pert))),np.round(np.max(np.abs(cluster_daily_pert))),Nlevels,endpoint=True)
iplot = 0
for i in range(Nrows):
  for j in range(Ncols):
    if iplot>=Nclusters:
      continue
    ax = fig.add_subplot(gs[i,j],projection=ccrs.PlateCarree())
    ax.add_feature(LAND,facecolor='lightgrey')
    ax.coastlines(linewidths=0.4)
    cf=ax.contourf(lon,lat,cluster_daily_pert[iplot,:,:],levels=plot_levels,cmap='coolwarm',transform=ccrs.PlateCarree(),extend='both')
    ax.set_title(SWT[iplot])
    gl = ax.gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
    ax.set_xlim([lon.min(),lon.max()])
    ax.set_ylim([lat.min(),lat.max()])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = (i == Nrows - 1)
    gl.left_labels = (j == 0)
    iplot+=1
ax_cbar = fig.add_subplot(gs[0:2,-1])
cax = inset_axes(ax_cbar, width="100%", height="80%", loc="center") 
plt.colorbar(cf,cax=cax,label=f"{varname} ({units[varname]})")
ax_cbar.axis("off")
plt.savefig(f"{path_out}{varname}_cluster_daily_pert.png")

# %%
