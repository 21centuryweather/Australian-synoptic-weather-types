# %% import packages
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
from matplotlib import cm, colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
# %% Input
# Paths and data
file_clusters            = "/home/565/fl2086/Australian-synoptic-weather-types/SWT_fields/SWT_data_v1.nc"
file_clusters_reassigned = "/home/565/fl2086/Australian-synoptic-weather-types/SWT_fields/SWT_data_reassigned_500_v1.nc"
path_out = "/home/565/fl2086/Australian-synoptic-weather-types/plotting/"
# %% Load cluster data
nc = Dataset(file_clusters,'r')
surface_clusters = nc.variables['clusterSeries'][:]
WR = nc.variables['WR'][:]
SWT = nc.variables['SWT'][:]
Nclusters = len(SWT)
nc.close()
nc = Dataset(file_clusters_reassigned,'r')
upper_clusters = nc.variables['clusterSeries'][:]
nc.close()
# %% Colors and plotting settings
figsize = (11.7, 8.3)
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
WR_names = ['AM','COL','WCT','FH','TH','EH','CH','WH'] # Order of plotting WR
WR_colors = {wr:color for wr,color in zip(WR_names,colors)} # Assign color WR
WR_cmaps = {wr:make_linear_cmap(WR_colors[wr],background_color,f"{wr}_cmap") for wr in WR_names} # create linear colormap WR
pdf = PdfPages(f'{path_out}cluster_transitions_500hpa.pdf') # PDF to write figures to
# %% Calculate transitions between surface and upper clusters for both WR and SWT
transition_matrix = np.zeros((Nclusters,Nclusters))
NWR = len(WR_names)
WR_to_idx = {wr : idx for idx,wr in enumerate(WR_names)}
transition_matrix_WR = np.zeros((NWR,NWR))
for isurf,iupper in zip(surface_clusters,upper_clusters):
  transition_matrix[isurf-1,iupper-1] += 1
  transition_matrix_WR[WR_to_idx[WR[isurf-1]],WR_to_idx[WR[iupper-1]]] += 1
# %% Sort clusters by name and WR
idx_sort = []
idx = np.arange(Nclusters)
for wr in WR_names[::-1]:
    mask = np.char.startswith(SWT.astype(str), wr) # Find indices where SWT starts with WR name
    wr_indices = idx[mask] # Extract indices
    sorted_wr_indices = wr_indices[np.argsort(SWT[wr_indices])] # Sort them alphabetically by SWT suffix
    idx_sort.extend(sorted_wr_indices) # Append to sorting indices
# Sort data using sorting indices
SWT = SWT[idx_sort]
WR  = WR[idx_sort]
transition_matrix = transition_matrix[np.ix_(idx_sort, idx_sort)]
transition_matrix_perc = transition_matrix/np.sum(transition_matrix,axis=1)[:,None]*100
# %% Sankey plot for WR changes
fig, ax = plt.subplots(figsize=figsize)
# Parameters
bar_width = 0.8
min_bar_height = 0.4
height_scale = 5  # scales percentage to bar height
row_spacing = 0.01  # fixed vertical spacing per cluster row
y_positions = np.arange(NWR) * row_spacing
# Change percentage per surface cluster
percent_change_WR = 100.0 - np.diag(transition_matrix_WR/np.sum(transition_matrix_WR,axis=1)[:,None]*100)
surface_bar_heights = np.clip((percent_change_WR / 100) * height_scale, min_bar_height, None)
scale_surface_heights = (1-NWR*row_spacing)/np.sum(surface_bar_heights)
surface_bar_heights_scaled = surface_bar_heights*scale_surface_heights 
# Inflow percentage per upper cluster
percent_inflow_WR = 100.0 - np.diag(transition_matrix_WR/np.sum(transition_matrix_WR,axis=0)[None,:]*100)
upper_bar_heights = np.clip((percent_inflow_WR / 100) * height_scale, min_bar_height, None)
scale_upper_heights = (1-NWR*row_spacing)/np.sum(upper_bar_heights)
upper_bar_heights_scaled = upper_bar_heights*scale_upper_heights
# Plot surface cluster bars (left)
y_pos_surface = np.cumsum(np.concatenate([[0],surface_bar_heights_scaled+row_spacing]))
for wr in WR_names:
  iwr = WR_to_idx[wr]
  y_center = y_pos_surface[iwr]+surface_bar_heights_scaled[iwr]/2
  rect = Rectangle((0, y_pos_surface[iwr]), bar_width, surface_bar_heights_scaled[iwr],
                    color=WR_colors[wr], ec=None, lw=0.5)
  ax.add_patch(rect)
  # Add WR and percent label in cluster bar
  label_text = f"{wr}\n{percent_change_WR[iwr]:.1f}%"
  ax.text(bar_width/2, y_center, label_text, va='center', ha='center',
          fontsize=14, color='white')
# Plot upper cluster bars (right)
y_pos_upper = np.cumsum(np.concatenate([[0],upper_bar_heights_scaled+row_spacing]))
for wr in WR_names:
  iwr = WR_to_idx[wr]
  y_center = y_pos_upper[iwr]+upper_bar_heights_scaled[iwr]/2
  rect = Rectangle((8, y_pos_upper[iwr]), bar_width, upper_bar_heights_scaled[iwr],
                    color=WR_colors[wr], ec=None, lw=0.5)
  ax.add_patch(rect)
  # Add WR and percent label in cluster bar
  label_text = f"{wr}\n{percent_inflow_WR[iwr]:.1f}%"
  ax.text(8+bar_width/2, y_center, label_text, va='center', ha='center',
          fontsize=14, color='white')
# Plot transition lines (exclude diagonal)
transition_width_surface = (transition_matrix_WR/np.sum(transition_matrix_WR,axis=1)[:,None]*100)/percent_change_WR[:,None]*surface_bar_heights_scaled[:,None]
transition_width_upper   = (transition_matrix_WR/np.sum(transition_matrix_WR,axis=0)[None,:]*100)/percent_inflow_WR[None,:]*upper_bar_heights_scaled[None,:]
transition_width_surface[np.arange(NWR),np.arange(NWR)] = 0
transition_width_upper[np.arange(NWR),np.arange(NWR)] = 0
# Get lower y coordinate of start transition line
yline_surface = y_pos_surface[:NWR,None]+np.cumsum(np.concatenate([np.zeros(NWR)[:,None],transition_width_surface],axis=1),axis=1)
# Get lower y coodrinate of end transition line
yline_upper   = y_pos_upper[None,:NWR]+np.cumsum(np.concatenate([np.zeros(NWR)[None,:],transition_width_upper],axis=0),axis=0)
# Get third largest values to display in transition lines
third_largest_surface = np.sort(transition_matrix_WR/np.sum(transition_matrix_WR,axis=1)[:,None]*100, axis=1)[:, -4]
third_largest_upper = np.sort(transition_matrix_WR/np.sum(transition_matrix_WR,axis=0)[None,:]*100, axis=0)[-4,:]
# x coordinate for transition lines
x = np.linspace(bar_width,8,100) 
for wr1 in WR_names:
  iwr1 = WR_to_idx[wr1]
  for wr2 in WR_names:
    iwr2 = WR_to_idx[wr2]
    if iwr1 == iwr2:
      continue  # skip same-cluster transitions
    # Plot transition line as a fill between two scaled cosine curves
    curve_top = (np.cos((x-bar_width)/(x.max()-x.min())*np.pi)+1)/2*\
                (yline_surface[iwr1,iwr2]+transition_width_surface[iwr1,iwr2]-yline_upper[iwr1,iwr2]-transition_width_upper[iwr1,iwr2])+\
                yline_upper[iwr1,iwr2]+transition_width_upper[iwr1,iwr2]
    curve_bottom = (np.cos((x-bar_width)/(x.max()-x.min())*np.pi)+1)/2*\
                (yline_surface[iwr1,iwr2]-yline_upper[iwr1,iwr2])+\
                yline_upper[iwr1,iwr2]
    x_full = np.concatenate([x, x[::-1]])
    y_full = np.concatenate([curve_top, curve_bottom[::-1]])
    ax.fill(x_full, y_full, color=WR_colors[wr1], alpha=0.5, linewidth=0, zorder=0)
    ax.plot(x,curve_top,WR_colors[wr1],linewidth=0.1,alpha=0.5)
    ax.plot(x,curve_bottom,WR_colors[wr1],linewidth=0.1)
    # Add 3 highest contributors as text labels in transition lines
    if((transition_matrix_WR/np.sum(transition_matrix_WR,axis=1)[:,None]*100)[iwr1,iwr2]>=third_largest_surface[iwr1]):
      ax.text(bar_width,yline_surface[iwr1,iwr2]+transition_width_surface[iwr1,iwr2]/2,f" {wr2} {(transition_matrix_WR/np.sum(transition_matrix_WR,axis=1)[:,None]*100)[iwr1,iwr2]:.1f}%", va='center', ha='left',color='w',fontsize=7)
    if((transition_matrix_WR/np.sum(transition_matrix_WR,axis=0)[None,:]*100)[iwr1,iwr2]>=third_largest_upper[iwr2]):
      ax.text(8,yline_upper[iwr1,iwr2]+transition_width_upper[iwr1,iwr2]/2,f" {wr1} {(transition_matrix_WR/np.sum(transition_matrix_WR,axis=0)[None,:]*100)[iwr1,iwr2]:.1f}%", va='center', ha='right',color='w',fontsize=7)
# Final layout
ax.set_xlim(0, 8+bar_width)
ax.set_ylim(0, 1)
ax.axis('off')
ax.text(bar_width/2, 1+row_spacing, "850 hPa", va='center', ha='center',
          fontsize=14, color=title_color)
ax.text(bar_width/2+8, 1+row_spacing, "500 hPa", va='center', ha='center',
          fontsize=14, color=title_color)
ax.set_title("Weather regimes 850 hPa \u2192 500 hPa transitions\n(Bar height = % of total data in regime)",
             fontsize=18,color=title_color)
fig.patch.set_facecolor(background_color)
ax.set_facecolor(background_color)
# plt.show()
pdf.savefig(fig)
# %% Plot full transition matrix as percentage from 850hPa cluster
fig, ax = plt.subplots(figsize=figsize)
# Create RGB matrix, scaled by maximum off diagonal value
rgb_matrix = np.zeros((Nclusters, Nclusters, 4))  # RGBA
transition_matrix_perc_max = np.max(transition_matrix_perc[~np.eye(Nclusters, dtype=bool)])
transition_matrix_perc_max = np.ceil(transition_matrix_perc_max/10*2)*10/2
transition_matrix_perc_norm = transition_matrix_perc/transition_matrix_perc_max
# Color cells
for i in range(Nclusters):
  cmap = WR_cmaps[WR[i]]
  for j in range(Nclusters):
    if i == j: # Color diagonals light grey
      rgb_matrix[i,j,:] = (0.3, 0.3, 0.3, 1.0)  # light grey
    else:
      rgb_matrix[i,j,:]=cmap(transition_matrix_perc_norm[i,j])
    # Add value in background color (results in automatically more contrast for higher values)
    ax.text(j, i,f"{transition_matrix_perc[i,j]:.0f}%",ha='center', va='center',fontsize=8,color=background_color,fontweight='medium')
im = ax.imshow(rgb_matrix, aspect='auto')
# Draw boxes around WR
boundaries = []
for i in range(1, Nclusters):
    if SWT[i].split('-')[0] != SWT[i-1].split('-')[0]:
        boundaries.append(i - 0.5)
for b in boundaries:
    ax.axhline(b, color=title_color, linewidth=1)
    ax.axvline(b, color=title_color, linewidth=1)
# Add SWT labels
ax.set_xticks(np.arange(Nclusters))
ax.set_yticks(np.arange(Nclusters))
ax.set_xticklabels(SWT, rotation=90)
ax.set_yticklabels(SWT)
# Add gridlines between cells
ax.set_xticks(np.arange(-0.5, Nclusters, 1), minor=True)
ax.set_yticks(np.arange(-0.5, Nclusters, 1), minor=True)
ax.grid(which='minor', color=title_color, linestyle='-', linewidth=0.25)
ax.tick_params(axis='both', which='minor', length=0)
# Set titles, labels and color axis
plt.title("Synoptic weather types 850 hPa \u2192 500 hPa transition matrix",color=title_color,fontsize=18)
ax.set_xlabel('500 hPa SWT',color=title_color,fontsize=14)
ax.set_ylabel('850 hPa SWT',color=title_color,fontsize=14)
ax.tick_params(axis='both', colors=title_color)
for spine in ax.spines.values():
    spine.set_color(title_color)
fig.patch.set_facecolor(background_color)
ax.set_facecolor(background_color)
# plt.show()
pdf.savefig(fig)
pdf.close()
# %%
