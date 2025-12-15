# %% Import packages
import numpy as np
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.patches import Polygon
from matplotlib import cm, colors
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, LAND
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.ndimage import minimum_filter, maximum_filter, gaussian_filter
from matplotlib.widgets import CheckButtons
from matplotlib.legend_handler import HandlerPatch, HandlerLine2D, HandlerBase
# %% Input
root = '/home/565/fl2086/Australian-synoptic-weather-types/'
file_clusters = f"{root}SWT_fields/SWT_data_v1.nc"
path_era5     = f"{root}example_data/era5_data/"
K = 4 # Number of principal components used
N = 3 # Number of days to plot
nvec = 2
# beta = 1
WR_names = ['AM','COL','WCT','FH','TH','EH','CH','WH']
# %% Colors and plotting settings
figsize = np.array([11.7/0.85, 8.3])
background_color = "#e8e6e1"
title_color = "#333333"
def make_linear_cmap(color, base, name):
    colors = [base, color]
    cmap = LinearSegmentedColormap.from_list(name, colors)
    return cmap
colors = [
    "#b35c44",  # Terracotta
    "#d39c55",  # Warm ochre
    "#a3a847",  # Olive green
    "#669966",  # Soft forest green
    "#5b8ea3",  # Dusty blue
    "#7b6ea7",  # Muted violet
    "#a87c7c",  # Earthy rose
    "#84786f",  # Weathered taupe
]
WR_colors = {wr:color for wr,color in zip(WR_names,colors)} # Assign color WR
WR_cmaps = {wr:make_linear_cmap(WR_colors[wr],background_color,f"{wr}_cmap") for wr in WR_names}
# %% Load cluster data
with Dataset(file_clusters) as nc:
  time = nc.variables['time'][:]
  lat = nc.variables['latitude'][:]
  lon = nc.variables['longitude'][:]
  SWT = nc.variables['SWT'][:]
  WR  = nc.variables['WR'][:]
  clusterSeries = nc.variables['clusterSeries'][:]
  clusterU = nc.variables['clusterU'][:]
  clusterV = nc.variables['clusterV'][:]
Lon, Lat = np.meshgrid(lon, lat)
Nt = len(time) # Number of days
Nc = len(SWT)  # Number of clusters
Nlat = len(lat) # Number of latitude points
Nlon = len(lon) # Number of longitude points
idx_SWT_map = {swt : i+1 for i,swt in enumerate(SWT)}
# %% Load u and v data
with Dataset(f"{path_era5}u_full.nc",'r') as nc_u, Dataset(f"{path_era5}v_full.nc",'r') as nc_v:
  u = nc_u.variables['u'][:]
  v = nc_v.variables['v'][:]
# %% Load data for plots
varnames = ['u_850hPa','v_850hPa','u_300hPa','v_300hPa','msl','tcwv','pv_315K','pv_330K']
plot_data = {}
for varname in varnames:
  with Dataset(f"{path_era5}{varname}.nc",'r') as nc:
    plot_data[varname] = nc.variables[varname.split('_')[0]][:]
    plot_data['time'] = nc.variables['time'][:]
Ujet = np.sqrt(plot_data['u_300hPa']**2+plot_data['v_300hPa']**2)
plot_data['jet'] = np.where(Ujet>20,Ujet,np.nan)
mask_tropical_westerlies = ((lat>-20) & (lat<=0))[None,:,None] & (plot_data['u_850hPa']>0)
plot_data['u_TW'] = np.where(mask_tropical_westerlies,plot_data['u_850hPa'],np.nan)
plot_data['v_TW'] = np.where(mask_tropical_westerlies,plot_data['v_850hPa'],np.nan)
# %% Get distances for each day to each cluster
d = np.zeros((Nt,Nc))
for ic in range(Nc):
  d[:,ic] = np.sqrt(np.mean((clusterU[ic,None,:,:]-u)**2+(clusterV[ic,None,:,:]-v)**2,axis=(1,2)))
# # %% Do principal component analysis for each cluster
# pca = {swt : PCA() for swt in SWT}
# E = np.zeros((Nt,Nc))
# for ic in range(Nc):
#   mask = clusterSeries == ic+1
#   Ndays = sum(mask)
#   u_dev = u[mask,:,:] - clusterU[ic,:,:]
#   v_dev = v[mask,:,:] - clusterV[ic,:,:]
#   X = np.hstack([u_dev.reshape(Ndays,Nlat*Nlon), v_dev.reshape(Ndays,Nlat*Nlon)])
#   pca[SWT[ic]].fit(X)
# # %% Calculate construction errors
# cluster_mean = np.hstack([clusterU.reshape(Nc, Nlat*Nlon), clusterV.reshape(Nc, Nlat*Nlon)])
# u_flat = u.reshape(Nt, Nlat*Nlon)
# v_flat = v.reshape(Nt, Nlat*Nlon)
# X_all = np.hstack([u_flat, v_flat])
# for ic in range(Nc):
#     # Get PCA for cluster
#     pca_ic = pca[SWT[ic]]
#     # Subtract cluster mean from all days for this cluster
#     X_dev = X_all - cluster_mean[ic]  # shape (Nt, Nlat*Nlon*2)
#     # Project onto PCA components
#     scores = pca_ic.transform(X_dev)  # shape (Nt, Nt_components)
#     # Reconstruct using first K components
#     X_reconstructed = scores[:, :K] @ pca_ic.components_[:K, :]
#     X_reconstructed += cluster_mean[ic]
#     # Compute mean absolute error for all days at once
#     E[:, ic] = np.sqrt(np.mean((X_all - X_reconstructed)**2, axis=1))
# # %% Get mean distance and mean reconstruction error per cluster
# E_cluster = np.zeros((Nc))
# d_cluster = np.zeros((Nc))
# for ic in range(Nc):
#   mask = clusterSeries == ic+1
#   E_cluster[ic] = np.mean(E[mask,ic])
#   d_cluster[ic] = np.mean(d[mask,ic])
# # %% Create scatter plot for errors
# plt.figure()
# plt.scatter(E_cluster,d_cluster)
# plt.xlabel('Mean construction error (m/s)')
# plt.ylabel('Mean distance from centroid (m/s)')
# plt.show()
# %% Get scores
# F = beta*d+(1-beta)*E
F = d
# %% Get softmax percentages
alpha = 1/np.std(np.min(F,axis=1))
P = np.exp(-alpha*F)/np.sum(np.exp(-alpha*F),axis=1)[:,None]
# %% Get confusion matrix
M = np.zeros((Nc,Nc))
for i in range(Nc):
  mask = clusterSeries == i+1
  M[i,:] = np.mean(P[mask,:],axis=0)
# %% Sort softmax data using sorting indices
idx_sort = []
idx = np.arange(Nc)
for wr in WR_names[::-1]:
  mask = np.char.startswith(SWT.astype(str), wr) # Find indices where SWT starts with WR name
  wr_indices = idx[mask] # Extract indices
  sorted_wr_indices = wr_indices[np.argsort(SWT[wr_indices])] # Sort them alphabetically by SWT suffix
  idx_sort.extend(sorted_wr_indices) # Append to sorting indices
P = P[:,idx_sort]
WR = WR[idx_sort]
SWT = SWT[idx_sort]
M = M[np.ix_(idx_sort, idx_sort)]
# %% Plot softmax matrix
fig, ax = plt.subplots(figsize=figsize)
# Create RGB matrix, scaled by maximum off diagonal value
rgb_matrix = np.zeros((Nc, Nc, 4))  # RGBA
# Color cells
for i in range(Nc):
  cmap = WR_cmaps[WR[i]]
  for j in range(Nc):
    if i == j: # Color diagonals light grey
      rgb_matrix[i,j,:] = (0.3, 0.3, 0.3, 1.0)
    else:
      rgb_matrix[i,j,:]=cmap(M[i,j]/np.max(M[~np.eye(Nc,dtype=bool)]))
    # Add value in background color (results in automatically more contrast for higher values)
    ax.text(j, i,f"{M[i,j]*100:.0f}%",ha='center', va='center',fontsize=8,color=background_color,fontweight='medium')
im = ax.imshow(rgb_matrix, aspect='auto',picker=True)
# Draw boxes around WR
boundaries = []
for i in range(1, Nc):
    if SWT[i].split('-')[0] != SWT[i-1].split('-')[0]:
        boundaries.append(i - 0.5)
for b in boundaries:
    ax.axhline(b, color=title_color, linewidth=1)
    ax.axvline(b, color=title_color, linewidth=1)
# Add SWT labels
ax.set_xticks(np.arange(Nc))
ax.set_yticks(np.arange(Nc))
ax.set_xticklabels(SWT, rotation=90)
ax.set_yticklabels(SWT)
# Add gridlines between cells
ax.set_xticks(np.arange(-0.5, Nc, 1), minor=True)
ax.set_yticks(np.arange(-0.5, Nc, 1), minor=True)
ax.grid(which='minor', color=title_color, linestyle='-', linewidth=0.25)
ax.tick_params(axis='both', which='minor', length=0)
# Set titles, labels and color axis
plt.title("Softmax confusion matrix",color=title_color,fontsize=18)
ax.set_xlabel('Cluster options',color=title_color,fontsize=14)
ax.set_ylabel('Assigned cluster',color=title_color,fontsize=14)
ax.tick_params(axis='both', colors=title_color)
for spine in ax.spines.values():
    spine.set_color(title_color)
fig.patch.set_facecolor(background_color)
ax.set_facecolor(background_color)
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
    xlows = Lon[local_min]; xhighs = Lon[local_max]
    ylows = Lat[local_min]; yhighs = Lat[local_max]
    lowvals = var[local_min]; highvals = var[local_max]
    return xlows,ylows,lowvals,xhighs,yhighs,highvals
# Create figure
def show_plot(time_plot,P_plot):
  fig2, ax2 = plt.subplots(figsize=(16*len(time_plot)/5,3.2),ncols=len(time_plot), subplot_kw=dict(projection=ccrs.PlateCarree()))
  for i in range(len(time_plot)):
    iday = np.argwhere(plot_data['time']==time_plot[i])[0][0]
    # Plot land
    ax2[i].add_feature(LAND,facecolor='lightgrey')
    ax2[i].coastlines(linewidths=0.4)
    # Plot msl contours
    plot_levels = np.arange(900,1060,4)
    c=ax2[i].contour(lon,lat,plot_data['msl'][iday,:,:]/100,
                  levels=plot_levels,colors='black',linewidths=0.75,linestyles='-',
                  transform=ccrs.PlateCarree())
    ax2[i].clabel(c, inline=True, levels=plot_levels[::5],fontsize=8,fmt=lambda val: f'{val:.0f}hPa')
    # _,_,_,xhighs,yhighs,highvals=extreme_vals(data['msl'],data['latitude'],data['longitude'],window=70,sigma=4)
    _,_,_,xhighs,yhighs,highvals=extreme_vals(plot_data['msl'][iday,:,:],window=70/6,sigma=4/6)
    # xlows,ylows,lowvals,_,_,_=extreme_vals(data['msl'],data['latitude'],data['longitude'],window=50,sigma=1)
    xlows,ylows,lowvals,_,_,_=extreme_vals(plot_data['msl'][iday,:,:],window=50/6,sigma=1/6)
    xyplotted = []
    dmin=2
    for x,y,p in zip(xlows, ylows, lowvals):
        if x < lon.max() and x > lon.min() and y < lat.max() and y > lat.min():
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                ax2[i].text(x,y,'L',fontsize=12,fontweight='bold',
                        ha='center',va='center',color='b')
                xyplotted.append((x,y))
    for x,y,p in zip(xhighs, yhighs, highvals):
        if x < lon.max() and x > lon.min() and y < lat.max() and y > lat.min():
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                ax2[i].text(x,y,'H',fontsize=12,fontweight='bold',
                        ha='center',va='center',color='r')
                xyplotted.append((x,y))
    # Plot jet
    plot_levels = np.arange(20,50,2.5)
    cf=ax2[i].contourf(lon,lat,plot_data['jet'][iday,:,:],#*1.94384,
                    levels=plot_levels,cmap='Blues',extend='max',
                    transform=ccrs.PlateCarree(),alpha=0.75)
    # Plot total column vertically-integrated water vapour
    plot_levels = np.arange(48,1000,500)
    cf=ax2[i].contourf(lon,lat,plot_data['tcwv'][iday,:,:],
                    levels=plot_levels,cmap='Greens',extend='max',
                    transform=ccrs.PlateCarree(),alpha=0.4)
    c=ax2[i].contour(lon,lat,plot_data['tcwv'][iday,:,:],
                      levels=[48],colors='limegreen',linewidths=0.8,linestyles='-',
                      transform=ccrs.PlateCarree())
    ax2[i].clabel(c, inline=True, levels=[48],fontsize=8,fmt=lambda val: f'{val:.0f}'+r'$kg\,m^{-2}$')
    # Plot PV contours
    c=ax2[i].contour(lon,lat,plot_data['pv_330K'][iday,:,:]*1e6,
                      levels=[-2],colors='magenta',linewidths=1.0,linestyles='--',
                      transform=ccrs.PlateCarree())
    c=ax2[i].contour(lon,lat,plot_data['pv_315K'][iday,:,:]*1e6,
                      levels=[-2],colors='magenta',linewidths=1.25,linestyles='-',
                      transform=ccrs.PlateCarree())
    # Plot tropical Westerlies
    q=ax2[i].quiver(lon[::nvec], lat[::nvec], plot_data['u_TW'][iday,::nvec,::nvec], plot_data['v_TW'][iday,::nvec,::nvec],
                  scale=3,scale_units='xy',width=0.003,minshaft=2,color='darkgreen',# minlength=1,
                  transform=ccrs.PlateCarree())
    # Add gridlines
    gl = ax2[i].gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
    ax2[i].set_xlim([lon.min(),lon.max()])
    ax2[i].set_ylim([lat.min(),lat.max()])
    # xlocators=list(np.arange(-180,190,10))
    #gl.xlocator = mticker.FixedLocator(xlocators)
    # Format figure
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = i==0
    ax2[i].set_title(f"{np.datetime64('1900-01-01')+np.timedelta64(int(time_plot[i]),'h')} P = {P_plot[i]*100:.2f}%")
    plt.show(block=False)

# Make matrix interactive
def on_click(event):
  artist = event.artist
  # Plot transition composite map
  if isinstance(artist, matplotlib.image.AxesImage):
    mouse_event = event.mouseevent
    if mouse_event.inaxes != ax:
      return
    j, i = int(round(mouse_event.xdata)), int(round(mouse_event.ydata))
    if 0 <= i < Nc and 0 <= j < Nc:
      swt = SWT[i]
      mask = np.argwhere((clusterSeries == idx_SWT_map[swt]) & (time>=plot_data['time'][0])).flatten()
      idx = mask[np.argsort(P[mask,j]).flatten()][-N:]
      P_plot = P[idx,j]
      time_plot = time[idx]
      show_plot(time_plot,P_plot)
      return
cid = fig.canvas.mpl_connect('pick_event', on_click)
plt.show()
# %%
