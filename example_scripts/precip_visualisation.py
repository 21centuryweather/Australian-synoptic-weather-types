# %% Import packages
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, LAND
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.backends.backend_pdf import PdfPages
# %% Input
root      = '/home/565/fl2086/Australian-synoptic-weather-types/'
path_out  = f"{root}results/"
months    = [9,10,11]
Nsim      = 5000
alpha     = 0.05
seed      = 42
WR_names  = ['AM','COL','WCT','FH','TH','EH','CH','WH'] # Order of plotting WR
month_labels = {1:'J',2:'F',3:'M',4:'A',5:'M',6:'J',7:'J',8:'A',9:'S',10:'O',11:'N',12:'D'}
month_string = ''
for m in months: month_string+=month_labels[m]
data_file = f"{root}results/precip_cluster_data_{month_string}.nc"
# %% Load data
with Dataset(data_file,'r') as nc:
  lat = nc.variables['lat'][:]
  lon = nc.variables['lon'][:]
  SWT = nc.variables['SWT'][:].astype(str)
  time = nc.variables['time_extremes'][:]
  counts = nc.variables['cluster_counts'][:]
  Nclusters,Nlat,Nlon = np.shape(counts)
  k = np.shape(time)[0]
  probs = nc.variables['probs'][:]
  precip = nc.variables['precip_extremes'][:]
  clusters = nc.variables['clusters_extremes'][:]
  precip_mean = np.zeros((Nclusters,Nlat,Nlon))
  for ic in range(Nclusters):
    precip_mean[ic,:,:] = np.sum(precip*(clusters-1==ic),axis=0)/counts[ic,:,:][None,:,:]
  precip_mean = np.where(counts>0,precip_mean,np.nan)
# %% Do Monte Carlo simulation to determine if significantly different from expected summer probabilities
dominant_cluster = np.argmax(counts, axis=0)
rng = np.random.default_rng(seed)
sims = rng.multinomial(k, probs, size=Nsim)
pvals = np.ones((Nclusters, Nlat, Nlon))
for c in range(Nclusters):
    pvals[c,:,:] = np.mean(sims[:, c][:, None, None] >= counts[c][None, :, :], axis=0)
signif = pvals<=alpha
# %% Colors and plotting settings
figsize = np.array([11.7, 8.3])
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
WR_cmaps = {wr:make_linear_cmap(WR_colors[wr],background_color,f"{wr}_cmap") for wr in WR_names} # create linear colormap WR
SWT_colors = np.zeros((len(SWT),4))
idx_sort = []
for wr in WR_names:
    swt_idx = np.where(np.char.startswith(SWT, wr))[0]
    n_swt = len(swt_idx)
    if(n_swt==0): continue
    positions = np.linspace(0.3, 1.0, n_swt)
    cmap = WR_cmaps[wr]
    for i, pos in zip(swt_idx[np.argsort(SWT[swt_idx])], positions):
        SWT_colors[i,:] = cmap(pos)
        idx_sort.append(i)

# %% Create pixels
pixels = SWT_colors[dominant_cluster]
shading = np.take_along_axis(signif,dominant_cluster[None,:,:],axis=0)[0,:,:]
# %% Create land mask
# prepare temporary plot and create mask from rasterized map
proj = {'projection': ccrs.PlateCarree()}
fig, ax = plt.subplots(figsize=(len(lon)/100, len(lat)/100), dpi=100, subplot_kw=proj)
fig.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0)
ax.set_frame_on(False)

ax.add_feature(LAND, facecolor='black')
ax.set_xlim([np.min(lon),np.max(lon)])
ax.set_ylim([np.min(lat),np.max(lat)])
fig.canvas.draw()
mask = fig.canvas.tostring_argb()
ncols, nrows = fig.canvas.get_width_height()
plt.close(fig)

mask = np.frombuffer(mask, dtype=np.uint8).reshape(nrows, ncols, 4)
mask = mask[::-1,:,1:].mean(axis=2)==0
# %% Apply mask
pixels *= mask[:,:,None]
shading *= mask
# %% Plot results
pdf = PdfPages(f'{path_out}precip_clusters_{month_string}.pdf')
fig, ax = plt.subplots(figsize=figsize,subplot_kw=dict(projection=ccrs.PlateCarree()))
ax.imshow(pixels,extent=[lon.min(), lon.max(), lat.min(), lat.max()],origin='lower',transform=ccrs.PlateCarree())
# Overlay hatching using pcolormesh
lon2d, lat2d = np.meshgrid(lon, lat)
plt.rcParams['hatch.linewidth'] = 0.2
ax.contourf(lon2d, lat2d, shading,
            levels=[0.5, 1.5],  # anything >=0.5 will be hatched
            hatches=['///'], alpha=0.0,  # transparent fill
            transform=ccrs.PlateCarree())
ax.add_feature(LAND,facecolor='lightgrey')
ax.coastlines(linewidths=0.4)
ax.set_xlim([lon.min(),lon.max()])
ax.set_ylim([lat.min(),lat.max()])
gl = ax.gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = True
gl.left_labels = True


colors_ordered = []
names_ordered = []

for wr in WR_names[::-1]:  # preserve the WR order
    swt_idx = np.where(np.char.startswith(SWT, wr))[0]  # indices of SWT in this WR
    if len(swt_idx) == 0:
        continue
    # get SWT names and colors
    swt_names_wr = SWT[swt_idx]
    swt_colors_wr = SWT_colors[swt_idx]
    
    # sort alphabetically by SWT name
    sorted_idx = np.argsort(swt_names_wr)[::-1]
    swt_names_wr_sorted = swt_names_wr[sorted_idx]
    swt_colors_wr_sorted = swt_colors_wr[sorted_idx]
    
    # append to final lists
    names_ordered.extend(swt_names_wr_sorted.tolist())
    colors_ordered.extend(swt_colors_wr_sorted)
cmap_manual = ListedColormap(colors_ordered)
norm_manual = BoundaryNorm(np.arange(len(colors_ordered)+1), len(colors_ordered))

sm = ScalarMappable(cmap=cmap_manual, norm=norm_manual)
sm.set_array([])  # required

cax = inset_axes(ax,
                 width="3%",      # colorbar width
                 height="100%",   # span full axis height
                 loc='center left',
                 borderpad=0,
                 bbox_to_anchor=(1.02, 0., 1, 1),
                 bbox_transform=ax.transAxes)

cbar = fig.colorbar(sm, cax=cax, orientation='vertical', ticks=np.arange(len(colors_ordered))+0.5)
cbar.ax.set_yticklabels(names_ordered,fontsize=8)
# cbar.set_label('SWT clusters (WR order, SWT alphabetical)')
ax.set_title(f'Dominant SWT for 5% wettest days in {month_string}',fontsize=14)
pdf.savefig()
# %% Probabilities
fig = plt.figure(figsize=(9,9))
for i in range(Nclusters):
  iplot = idx_sort[i]
  plt.bar(i,probs[iplot]*100,0.8,color=SWT_colors[iplot])
plt.xlim(np.array([0,Nclusters])-0.5)
plt.ylim([0,10])
plt.xticks(np.arange(Nclusters),labels=SWT[idx_sort],rotation=90)
plt.yticks(np.arange(0,15,5))
plt.title(f'SWT occurance in {month_string}')
plt.ylabel('probability (%)')
plt.grid(linewidth=0.5,color='lightgrey',alpha=0.5,axis='y')
pdf.savefig()
# %% Occurance
fig = plt.figure(figsize=(9,9))
fig.suptitle(f"Significant SWT occurance in 5% wettest days in {month_string}",y=0.93,fontsize=14)
gs = GridSpec(6, 5+1, width_ratios=[4]*5+[1],hspace=0)
i = 0
for irow in range(6):
    for icol in range(5):
        if i>=Nclusters:
            continue
        iplot=idx_sort[i]
        ax = fig.add_subplot(gs[irow,icol],projection=ccrs.PlateCarree())
        cf=ax.contourf(lon2d,lat2d,np.where((signif[iplot,:,:] & mask),counts[iplot,:,:]/k*100,np.nan),levels=np.linspace(0,20,5),cmap="viridis",extend='max')
        ax.set_title(SWT[iplot])
        ax.add_feature(LAND,facecolor='lightgrey')
        ax.coastlines(linewidths=0.4)
        ax.set_xlim([lon.min(),lon.max()])
        ax.set_ylim([lat.min(),lat.max()])
        gl = ax.gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = (irow == 6 - 1)
        gl.left_labels = (icol == 0)
        i+=1
ax_cbar = fig.add_subplot(gs[0:2,-1])
cax = inset_axes(ax_cbar, width="100%", height="80%", loc="center") 
plt.colorbar(cf,cax=cax,label=f"% of days")
ax_cbar.axis("off")
pdf.savefig()
# %% Precip
fig = plt.figure(figsize=(9,9))
fig.suptitle(f"Mean precipitation for 5% wettest days in {month_string}",y=0.93,fontsize=14)
gs = GridSpec(6, 5+1, width_ratios=[4]*5+[1],hspace=0)
i = 0
for irow in range(6):
    for icol in range(5):
        if i>=Nclusters:
            continue
        iplot=idx_sort[i]
        ax = fig.add_subplot(gs[irow,icol],projection=ccrs.PlateCarree())
        # cf=ax.contourf(lon2d,lat2d,np.where(signif[iplot,:,:],counts[iplot,:,:]/k*100,np.nan),levels=np.linspace(0,50,11),cmap="inferno")
        cf=ax.contourf(lon2d,lat2d,np.where((signif[iplot,:,:] & mask),precip_mean[iplot,:,:],np.nan),levels=np.linspace(0,60,7),cmap="viridis",extend='max')
        ax.set_title(SWT[iplot])
        ax.add_feature(LAND,facecolor='lightgrey')
        ax.coastlines(linewidths=0.4)
        ax.set_xlim([lon.min(),lon.max()])
        ax.set_ylim([lat.min(),lat.max()])
        gl = ax.gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = (irow == 6 - 1)
        gl.left_labels = (icol == 0)
        i+=1
ax_cbar = fig.add_subplot(gs[0:2,-1])
cax = inset_axes(ax_cbar, width="100%", height="80%", loc="center") 
plt.colorbar(cf,cax=cax,label=f"mm")
ax_cbar.axis("off")
pdf.savefig()
pdf.close()
# %%
