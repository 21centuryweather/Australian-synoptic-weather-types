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
file_data = "/home/565/fl2086/Australian-synoptic-weather-types/reassignment_analysis/data_interactive_transition_matrix_extratropics.nc"
WR_names = ['AM','COL','WCT','FH','TH','EH','CH','WH'] # Order of plotting WR
label_margin = 0.4
# %% Colors and plotting settings
figsize = np.array([11.7/0.85, 8.3])
background_color = "#e8e6e1"
title_color = "#333333"
nvec=2 # nvec=15
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
# %% Legend definitions
class ArrowHandler(HandlerBase):
    def __init__(self, color='k'):
          super().__init__()
          self.color = color
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        shaft_thickness = height * 0.1
        head_width = height * 0.4
        head_length = width * 0.3

        arrow = mpatches.FancyArrow(
            x=x0,
            y=y0 + height / 2 - shaft_thickness / 2,
            dx=width,
            dy=0,
            width=shaft_thickness,
            length_includes_head=True,
            head_width=head_width,
            head_length=head_length,
            color=self.color
        )
        return [arrow]
def make_quiver_legend(color='k'):
    arrow._legend_color = color
    arrow = mpatches.FancyArrow(0, 0, 1, 0, color=color)
    return arrow
legend_defs = {
    'msl': mlines.Line2D([], [], color='black', linestyle='-', linewidth=0.75),
    'tcwv': mpatches.Patch(color='green', alpha=0.4),
    'jet': mpatches.Patch(color='blue', alpha=0.75),
    'z_500hPa': mlines.Line2D([], [], color='C0', linestyle='-', linewidth=0.75),
    'pv_330K': mlines.Line2D([], [], color='magenta', linestyle='--', linewidth=1),
    'pv_315K': mlines.Line2D([], [], color='magenta', linestyle='-', linewidth=1.25),
    'uv_850hPa_TW': mpatches.FancyArrow(0, 0, 1, 0, color='darkgreen'),
    'uv_850hPa': mpatches.FancyArrow(0, 0, 1, 0, color='black'),
    'uv_500hPa': mpatches.FancyArrow(0, 0, 1, 0, color='black'),
}
# %% Load transition matrix and composites
nc = Dataset(file_data,'r')
data = {var : nc.variables[var][:] for var in nc.variables}
nc.close()
Nclusters = len(data['SWT_original'])
data['Lon'], data['Lat'] = np.meshgrid(data['longitude'], data['latitude'])
# %% Sort data using sorting indices
idx_sort = []
idx = np.arange(Nclusters)
for wr in WR_names[::-1]:
    mask = np.char.startswith(data['SWT_original'].astype(str), wr) # Find indices where SWT starts with WR name
    wr_indices = idx[mask] # Extract indices
    sorted_wr_indices = wr_indices[np.argsort(data['SWT_original'][wr_indices])] # Sort them alphabetically by SWT suffix
    idx_sort.extend(sorted_wr_indices) # Append to sorting indices
for key, value in data.items():
    array = value
    for axis, dim in enumerate(array.shape):
        if dim == Nclusters:
            array = np.take(array, idx_sort, axis=axis)
    data[key] = array
# data['clusterID'] = data['clusterID'][idx_sort]
# data['SWT_original'] = data['SWT_original'][idx_sort]
# data['WR_original']  = data['WR_original'][idx_sort]
# data['Tmatrix'] = data['Tmatrix'][np.ix_(idx_sort, idx_sort)]
data['Tmatrix_perc'] = data['Tmatrix']/np.sum(data['Tmatrix'],axis=1)[:,None]*100
# %% Plot transition matrix
fig, ax = plt.subplots(figsize=figsize)
# Create RGB matrix, scaled by maximum off diagonal value
rgb_matrix = np.zeros((Nclusters, Nclusters, 4))  # RGBA
Tmatrix_perc_max = np.max(data['Tmatrix_perc'][~np.eye(Nclusters, dtype=bool)])
Tmatrix_perc_max = np.ceil(Tmatrix_perc_max/10*2)*10/2
data['Tmatrix_perc_norm'] = data['Tmatrix_perc']/Tmatrix_perc_max
# Color cells
for i in range(Nclusters):
  cmap = WR_cmaps[data['WR_original'][i]]
  for j in range(Nclusters):
    if i == j: # Color diagonals light grey
      rgb_matrix[i,j,:] = (0.3, 0.3, 0.3, 1.0)  # light grey
    else:
      rgb_matrix[i,j,:]=cmap(data['Tmatrix_perc_norm'][i,j])
    # Add value in background color (results in automatically more contrast for higher values)
    ax.text(j, i,f"{data['Tmatrix_perc'][i,j]:.0f}%",ha='center', va='center',fontsize=8,color=background_color,fontweight='medium')
im = ax.imshow(rgb_matrix, aspect='auto',picker=True)
# Draw boxes around WR
boundaries = []
for i in range(1, Nclusters):
    if data['SWT_original'][i].split('-')[0] != data['SWT_original'][i-1].split('-')[0]:
        boundaries.append(i - 0.5)
for b in boundaries:
    ax.axhline(b, color=title_color, linewidth=1)
    ax.axvline(b, color=title_color, linewidth=1)
# Add SWT labels
ax.set_xticks(np.arange(Nclusters))
ax.set_yticks(np.arange(Nclusters))
ax.set_xticklabels(data['SWT_original'], rotation=90)
ax.set_yticklabels(data['SWT_original'])
# Add gridlines between cells
ax.set_xticks(np.arange(-0.5, Nclusters, 1), minor=True)
ax.set_yticks(np.arange(-0.5, Nclusters, 1), minor=True)
ax.grid(which='minor', color=title_color, linestyle='-', linewidth=0.25)
ax.tick_params(axis='both', which='minor', length=0)
# Set titles, labels and color axis
plt.title("Synoptic weather types original \u2192 extratropics only transition matrix",color=title_color,fontsize=18)
ax.set_xlabel('Original SWT',color=title_color,fontsize=14)
ax.set_ylabel('Reassigned SWT',color=title_color,fontsize=14)
ax.tick_params(axis='both', colors=title_color)
for spine in ax.spines.values():
    spine.set_color(title_color)
fig.patch.set_facecolor(background_color)
ax.set_facecolor(background_color)
# Set clickable ticklabels
for i, label in enumerate(data['SWT_original']):
    rect = plt.Rectangle((-2.3, i - 0.5), 1.8, 1, facecolor='none',linewidth=0.25, edgecolor=title_color, picker=True, clip_on=False, zorder=10)
    ax.add_patch(rect)
    rect._label_index = i  # Store the index
    rect._label_type = 'original'
for j, label in enumerate(data['SWT_original']):
    rect = plt.Rectangle((j - 0.5, Nclusters-0.5), 1, 2.7, facecolor='none',linewidth=0.25, edgecolor=title_color, picker=True, clip_on=False, zorder=10)
    ax.add_patch(rect)
    rect._label_index = j  # Store the index
    rect._label_type = 'reassigned'
# Add clickable tickboxes for plot options
labels = ['msl','tcwv','jet','pv_315K','pv_330K','z_500hPa','uv_850hPa_TW','uv_850hPa','uv_500hPa']
labels_text = dict(msl = r'$\overline{msl}$', tcwv = r'$\overline{tcwv}$', jet = 'jet', pv_315K = r'$\overline{pv}$ 315K', pv_330K = r'$\overline{pv}$ 330K',z_500hPa = r"$\overline{z'}$ 500hPa", uv_850hPa_TW = r'$\overline{u},\overline{v}$ 850hPa TW', uv_850hPa = r'$\overline{u},\overline{v}$ 850hPa', uv_500hPa=r'$\overline{u},\overline{v}$ 500hPa')
initially_active = dict(msl = True, tcwv = True, jet = True, pv_315K = True, pv_330K = True, 
                        uv_850hPa_TW = True, uv_850hPa = False, uv_500hPa = False,z_500hPa = False)
currently_active = initially_active.copy()
label_tick = {}
hcbox = 1
for i,label in enumerate(labels):
  x0 = Nclusters-0.1
  y0 = 1.5*hcbox*i-0.4
  rect = plt.Rectangle((x0,y0), hcbox, hcbox, facecolor='none',linewidth=0.25, edgecolor=title_color, picker=True, clip_on=False, zorder=10)
  ax.add_patch(rect)
  rect._label_type = label
  ax.text(x0+hcbox+0.2, y0+hcbox/2,labels_text[label],ha='left', va='center',fontsize=12,color=title_color)
  # Draw tick (initially shown/hidden)
  label_tick[label] =  ax.text(x0 + hcbox / 2, y0 + hcbox / 2+0.05, '✓' if initially_active[label] else '',
                ha='center', va='center', fontsize=16, color=title_color, zorder=11)
i = i+1
rect = plt.Rectangle((Nclusters-0.1, 1.5*hcbox*i-0.4), 3.15, hcbox, facecolor='none',linewidth=0.25, edgecolor=title_color, picker=True, clip_on=False, zorder=10)
ax.add_patch(rect)
rect._label_type = 'deselect'
ax.text(Nclusters, 1.5*hcbox*i+hcbox/2-0.35,'Deselect all',ha='left', va='center',fontsize=12,color=title_color)
i = i+1
rect = plt.Rectangle((Nclusters-0.1, 1.5*hcbox*i-0.4), 1.6, hcbox, facecolor='none',linewidth=0.25, edgecolor=title_color, picker=True, clip_on=False, zorder=10)
ax.add_patch(rect)
rect._label_type = 'reset'
ax.text(Nclusters, 1.5*hcbox*i+hcbox/2-0.35,'Reset',ha='left', va='center',fontsize=12,color=title_color)
fig.subplots_adjust(right=0.85)
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
# %% Get data to be plotted
def get_composite_transition(i,j):
  u_850hPa = data['u_850hPa'][i,j,::nvec,::nvec]
  v_850hPa = data['v_850hPa'][i,j,::nvec,::nvec]
  u_500hPa = data['u_500hPa'][i,j,::nvec,::nvec]
  v_500hPa = data['v_500hPa'][i,j,::nvec,::nvec]  
  mask_tropical_westerlies = (((data['latitude'][::nvec]>-20) & (data['latitude'][::nvec]<=0))[:,None]) & \
                               (u_850hPa>0)
  u_TW = np.where(mask_tropical_westerlies,u_850hPa,np.nan)
  v_TW = np.where(mask_tropical_westerlies,v_850hPa,np.nan)
  Ujet = np.sqrt(data['u_300hPa'][i,j,:,:]**2+data['v_300hPa'][i,j,:,:]**2)
  jet = np.where(Ujet>20,Ujet,np.nan)
  plot_data = dict(pv_315K = data['pv_315K'][i,j,:,:], pv_330K = data['pv_330K'][i,j,:,:], msl = data['msl'][i,j,:,:], 
                   tcwv = data['tcwv'][i,j,:,:], u_TW = u_TW, v_TW = v_TW, jet = jet,z_500hPa = data['z_500hPa'][i,j,:,:],
                   u_850hPa = u_850hPa, v_850hPa = v_850hPa, u_500hPa = u_500hPa, v_500hPa = v_500hPa)
  return plot_data

def get_composite_original(i):
  u_850hPa = np.nansum(data['u_850hPa'][i,:,::nvec,::nvec]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  v_850hPa = np.nansum(data['v_850hPa'][i,:,::nvec,::nvec]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  u_500hPa = np.nansum(data['u_500hPa'][i,:,::nvec,::nvec]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  v_500hPa = np.nansum(data['v_500hPa'][i,:,::nvec,::nvec]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  u_300hPa = np.nansum(data['u_300hPa'][i,:,:,:]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  v_300hPa = np.nansum(data['v_300hPa'][i,:,:,:]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  msl      = np.nansum(data['msl'][i,:,:,:]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  z_500hPa = np.nansum(data['z_500hPa'][i,:,:,:]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  tcwv     = np.nansum(data['tcwv'][i,:,:,:]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  pv_315K  = np.nansum(data['pv_315K'][i,:,:,:]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  pv_330K  = np.nansum(data['pv_330K'][i,:,:,:]*data['Tmatrix'][i,:,None,None],axis=0)/np.sum(data['Tmatrix'][i,:])
  mask_tropical_westerlies = (((data['latitude'][::nvec]>-20) & (data['latitude'][::nvec]<=0))[:,None]) & (u_850hPa>0)
  u_TW = np.where(mask_tropical_westerlies,u_850hPa,np.nan)
  v_TW = np.where(mask_tropical_westerlies,v_850hPa,np.nan)
  Ujet = np.sqrt(u_300hPa**2+v_300hPa**2)
  jet = np.where(Ujet>20,Ujet,np.nan)
  plot_data = dict(pv_315K = pv_315K, pv_330K = pv_330K, msl = msl, tcwv = tcwv, u_TW = u_TW, v_TW = v_TW, jet = jet,
                   u_850hPa = u_850hPa, v_850hPa = v_850hPa, u_500hPa = u_500hPa, v_500hPa = v_500hPa, z_500hPa = z_500hPa)
  return plot_data

def get_composite_reassigned(j):
  u_850hPa = np.nansum(data['u_850hPa'][:,j,::nvec,::nvec]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  v_850hPa = np.nansum(data['v_850hPa'][:,j,::nvec,::nvec]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  u_500hPa = np.nansum(data['u_500hPa'][:,j,::nvec,::nvec]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  v_500hPa = np.nansum(data['v_500hPa'][:,j,::nvec,::nvec]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  u_300hPa = np.nansum(data['u_300hPa'][:,j,:,:]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  v_300hPa = np.nansum(data['v_300hPa'][:,j,:,:]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  msl      = np.nansum(data['msl'][:,j,:,:]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  z_500hPa = np.nansum(data['z_500hPa'][:,j,:,:]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  tcwv     = np.nansum(data['tcwv'][:,j,:,:]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  pv_315K  = np.nansum(data['pv_315K'][:,j,:,:]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  pv_330K  = np.nansum(data['pv_330K'][:,j,:,:]*data['Tmatrix'][:,j,None,None],axis=0)/np.sum(data['Tmatrix'][:,j])
  mask_tropical_westerlies = (((data['latitude'][::nvec]>-20) & (data['latitude'][::nvec]<=0))[:,None]) & (u_850hPa>0)
  u_TW = np.where(mask_tropical_westerlies,u_850hPa,np.nan)
  v_TW = np.where(mask_tropical_westerlies,v_850hPa,np.nan)
  Ujet = np.sqrt(u_300hPa**2+v_300hPa**2)
  jet = np.where(Ujet>20,Ujet,np.nan)
  plot_data = dict(pv_315K = pv_315K, pv_330K = pv_330K, msl = msl, tcwv = tcwv, u_TW = u_TW, v_TW = v_TW, jet = jet,
                   u_850hPa = u_850hPa, v_850hPa = v_850hPa, u_500hPa = u_500hPa, v_500hPa = v_500hPa, z_500hPa = z_500hPa)
  return plot_data
# %% Function to get and plot composite map
def show_composite_plot(plot_data,i,j,mode):
  # Create figure
  fig2, ax2 = plt.subplots(figsize=figsize*0.75, subplot_kw=dict(projection=ccrs.PlateCarree()))
  # Plot land
  ax2.add_feature(LAND,facecolor='lightgrey')
  ax2.coastlines(linewidths=0.4)
  # Plot msl contours
  if(currently_active['msl']):
    plot_levels = np.arange(900,1060,4)
    c=ax2.contour(data['longitude'],data['latitude'],plot_data['msl']/100,
                  levels=plot_levels,colors='black',linewidths=0.75,linestyles='-',
                  transform=ccrs.PlateCarree())
    ax2.clabel(c, inline=True, levels=plot_levels[::5],fontsize=8,fmt=lambda val: f'{val:.0f}hPa')
    # _,_,_,xhighs,yhighs,highvals=extreme_vals(data['msl'],data['latitude'],data['longitude'],window=70,sigma=4)
    _,_,_,xhighs,yhighs,highvals=extreme_vals(plot_data['msl'],window=70/6,sigma=4/6)
    # xlows,ylows,lowvals,_,_,_=extreme_vals(data['msl'],data['latitude'],data['longitude'],window=50,sigma=1)
    xlows,ylows,lowvals,_,_,_=extreme_vals(plot_data['msl'],window=50/6,sigma=1/6)
    xyplotted = []
    dmin=2
    for x,y,p in zip(xlows, ylows, lowvals):
        if x < data['longitude'].max() and x > data['longitude'].min() and y < data['latitude'].max() and y > data['latitude'].min():
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                ax2.text(x,y,'L',fontsize=14,fontweight='bold',
                        ha='center',va='center',color='b')
                xyplotted.append((x,y))
    for x,y,p in zip(xhighs, yhighs, highvals):
        if x < data['longitude'].max() and x > data['longitude'].min() and y < data['latitude'].max() and y > data['latitude'].min():
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0,y0 in xyplotted]
            if not dist or min(dist) > dmin:
                ax2.text(x,y,'H',fontsize=14,fontweight='bold',
                        ha='center',va='center',color='r')
                xyplotted.append((x,y))
  # Plot jet
  if(currently_active['jet']):
    plot_levels = np.arange(20,50,2.5)
    cf=ax2.contourf(data['longitude'],data['latitude'],plot_data['jet'],#*1.94384,
                    levels=plot_levels,cmap='Blues',extend='max',
                    transform=ccrs.PlateCarree(),alpha=0.75)
  # Plot total column vertically-integrated water vapour
  if(currently_active['tcwv']):
    plot_levels = np.arange(48,1000,500)
    cf=ax2.contourf(data['longitude'],data['latitude'],plot_data['tcwv'],
                    levels=plot_levels,cmap='Greens',extend='max',
                    transform=ccrs.PlateCarree(),alpha=0.4)
    c=ax2.contour(data['longitude'],data['latitude'],plot_data['tcwv'],
                      levels=[45],colors='limegreen',linewidths=0.8,linestyles='-',
                      transform=ccrs.PlateCarree())
    ax2.clabel(c, inline=True, levels=[45],fontsize=8,fmt=lambda val: f'{val:.0f}'+r'$kg\,m^{-2}$')
  # Plot geopotential height perturbations
  if(currently_active['z_500hPa']):
    plot_levels = np.arange(-5.1,5.1,0.2)
    c=ax2.contour(data['longitude'],data['latitude'],plot_data['z_500hPa']/9.80665,
                      levels=plot_levels,colors='C0',linewidths=0.75,
                      transform=ccrs.PlateCarree())
    ax2.clabel(c, inline=True, levels=plot_levels[::3],fontsize=8,fmt=lambda val: f'{val*100:.0f}cm')
  # Plot PV contours
  if(currently_active['pv_330K']):
    c=ax2.contour(data['longitude'],data['latitude'],plot_data['pv_330K']*1e6,
                      levels=[-2],colors='magenta',linewidths=1.0,linestyles='--',
                      transform=ccrs.PlateCarree())
  if(currently_active['pv_315K']):
    c=ax2.contour(data['longitude'],data['latitude'],plot_data['pv_315K']*1e6,
                      levels=[-2],colors='magenta',linewidths=1.25,linestyles='-',
                      transform=ccrs.PlateCarree())
  # Plot tropical Westerlies
  if(currently_active['uv_850hPa_TW']):
    q=ax2.quiver(data['longitude'][::nvec],data['latitude'][::nvec], plot_data['u_TW'], plot_data['v_TW'],
                  scale=3,scale_units='xy',width=0.003,minshaft=2,color='darkgreen',# minlength=1,
                  transform=ccrs.PlateCarree())
  # Plot wind at 850 hPa
  if(currently_active['uv_850hPa']):
    q=ax2.quiver(data['longitude'][::nvec],data['latitude'][::nvec], plot_data['u_850hPa'], plot_data['v_850hPa'],
                  scale=3,scale_units='xy',width=0.003,minshaft=2,color='k',# minlength=1,
                  transform=ccrs.PlateCarree())
  # Plot wind at 500 hPa
  if(currently_active['uv_500hPa']):
    q=ax2.quiver(data['longitude'][::nvec],data['latitude'][::nvec], plot_data['u_500hPa'], plot_data['v_500hPa'],
                  scale=5,scale_units='xy',width=0.003,color='k',minshaft=2,# minlength=1,
                  transform=ccrs.PlateCarree())
  # Add gridlines
  gl = ax2.gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
  ax2.set_xlim([data['longitude'].min(),data['longitude'].max()])
  ax2.set_ylim([data['latitude'].min(),data['latitude'].max()])
  # Move axis to create space for legend
  box = ax2.get_position()
  ax2.set_position([box.x0 - 0.08, box.y0, box.width, box.height])
  # Create legend
  legend_labels = []
  legend_handles = []
  legend_handlers = {}

  for key in labels:
      if currently_active.get(key):
          handle = legend_defs[key]
          legend_handles.append(handle)
          legend_labels.append(labels_text[key])
          if 'uv_' in key:       
              legend_handlers[handle] = ArrowHandler(color=handle.get_facecolor())

  ax2.legend(handles=legend_handles,
           labels=legend_labels,
           handler_map=legend_handlers,
           loc='upper left', fontsize=12, title="Plotted variables", title_fontsize=12,
           bbox_to_anchor=(1.0, 1.016))
  # xlocators=list(np.arange(-180,190,10))
  #gl.xlocator = mticker.FixedLocator(xlocators)
  # Format figure
  gl.xformatter = LONGITUDE_FORMATTER
  gl.yformatter = LATITUDE_FORMATTER
  gl.top_labels = False
  gl.right_labels = False
  gl.bottom_labels = True
  gl.left_labels = True
  if(mode=='transition'):
    ax2.set_title(f"{data['SWT_original'][i]} → {data['SWT_original'][j]}, {int(data['Tmatrix'][i,j])} days",fontsize=18)
  elif(mode=='original'):
    ax2.set_title(f"{data['SWT_original'][i]} before reassignment, {int(np.sum(data['Tmatrix'][i,:]))} days",fontsize=18)
  elif(mode=='reassigned'):
    ax2.set_title(f"{data['SWT_original'][j]} after reassignment, {int(np.sum(data['Tmatrix'][:,j]))} days",fontsize=18)
  plt.show(block=False)
# %% Make plot interactive
def on_click(event):
  artist = event.artist
  # Plot transition composite map
  if isinstance(artist, matplotlib.image.AxesImage):
    mouse_event = event.mouseevent
    if mouse_event.inaxes != ax:
      return
    j, i = int(round(mouse_event.xdata)), int(round(mouse_event.ydata))
    if 0 <= i < Nclusters and 0 <= j < Nclusters:
      plot_data = get_composite_transition(i,j)
      show_composite_plot(plot_data,i,j,'transition')
      return

  # Plot SWT composite map
  elif isinstance(artist, plt.Rectangle):
    label_type = artist._label_type
    # Plot original composite map
    if label_type == 'original': # original composite map
      i = artist._label_index
      plot_data = get_composite_original(i)
      show_composite_plot(plot_data,i,None,'original')
    # Plot reassigned composite map
    elif label_type == 'reassigned': # reassigned composite map
      j = artist._label_index
      plot_data = get_composite_reassigned(j)
      show_composite_plot(plot_data,None,j,'reassigned')
    # Update ticks
    elif label_type in labels:
      label = label_type
      currently_active[label] = not currently_active[label]
      label_tick[label].set_text('✓' if currently_active[label] else '')
      event.canvas.draw_idle()
    elif label_type == 'reset':
      for label in labels:
        currently_active[label] = initially_active[label]
        label_tick[label].set_text('✓' if currently_active[label] else '')
      event.canvas.draw_idle()
    elif label_type == 'deselect':
      for label in labels:
        currently_active[label] = False
        label_tick[label].set_text('✓' if currently_active[label] else '')
      event.canvas.draw_idle()  

cid = fig.canvas.mpl_connect('pick_event', on_click)
plt.show()