# %% Import packages
import sys
sys.path.append('../utils/') # Path to scripts
import read_era5 as era5
import compute_statistics as stat
import kmeans_clustering as kmeans
from netCDF4 import Dataset
import numpy as np
# %% Input
date_start    = '1952-01' # Start year-month for clustering data yyyy-mm 1952
date_end      = '2023-12' # End year-month for clustering data yyyy-mm (included)
utc           = 12        # Clustering sampling time (utc) data is sampled daily
plevel        = 850       # Height level used (hPa)
lat_lims      = np.array([-50,-25])  # South and North lattitude limit of analysis box
lon_lims      = np.array([100,165]) # West and East longitude limit of analysis box
Ncoarsen      = 6         # Coarsening factor in lat and lon direction
Nclusters     = 30        # Number of clusters
exclude       = ['AM']
file_clusters = "/home/565/fl2086/Australian-synoptic-weather-types/SWT_fields/SWT_data_v1.nc"
path_out      = "/home/565/fl2086/Australian-synoptic-weather-types/reassignment_analysis/" # Directory to save cluster results
path_data     = "/g/data/rt52/era5/pressure-levels/reanalysis/"         # Era5 data directory
# %% Load clustering information (SWT)
nc = Dataset(file_clusters,'r')
cluster_series = nc.variables['clusterSeries'][:]
cluster_time = nc.variables['time'][:]
lat_cluster = nc.variables['latitude'][:]
lon_cluster = nc.variables['longitude'][:]
mask_lat = np.argwhere((lat_cluster>=np.min(lat_lims)) & (lat_cluster<=np.max(lat_lims))).flatten()
mask_lon = np.argwhere((lon_cluster>=np.min(lon_lims)) & (lon_cluster<=np.max(lon_lims))).flatten()
WR = nc.variables['WR'][:]
clusterU    = nc.variables['clusterU'][:,mask_lat,mask_lon]
clusterV    = nc.variables['clusterV'][:,mask_lat,mask_lon]
for wr in exclude:
  mask_WR = np.argwhere(WR==wr)
  clusterU[mask_WR,:,:] *= np.NaN
  clusterV[mask_WR,:,:] *= np.NaN
nc.close()
# %% Load ERA5 data at specified height level hpa and calculate cluster mean
data = {}
for varname in ['u','v']:
  print(f"Read data for {varname}")
  data[varname],time,lat,lon = era5.read_data(varname,date_start,date_end,utc,lat_lims,lon_lims,path_data,varname_path=varname,Ncoarsen=Ncoarsen,level=plevel,progress=True)
  # print(f"Calculate cluster mean at {plevel}hpa for {varname}")
  # field_sum,field_sum2,field_count,time_era5 = stat.compute_sums(data[varname],time,cluster_series,cluster_time,Nclusters)
  # cluster_mean[varname] = stat.cluster_mean(field_sum,field_count)
# %% Assign days to new cluster centres
cluster_series_new = np.zeros((len(time)))
for it,t in enumerate(time):
  cluster_series_new[it] = kmeans.assign(data['u'][it,:,:],data['v'][it,:,:],clusterU,clusterV)
# %% save new cluster time series and cluster centres to output
fileout = path_out+f"SWT_data_reassigned_extratropics_v1.nc"
kmeans.save(fileout,Nclusters,time,lat,lon,clusterU,clusterV,cluster_series_new,None,None)
# %%
