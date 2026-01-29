# %% Import packages
import sys
sys.path.append('../utils/') # Path to scripts
import read_era5 as era5
import compute_statistics as stat
import kmeans_clustering as kmeans
from netCDF4 import Dataset
import numpy as np
# %% Input
root          = '/home/565/gdata-gb02/mb0427/Australian-synoptic-weather-types/'
date_start    = '1952-01' # Start year-month for clustering data yyyy-mm 1952
date_end      = '2023-12' # End year-month for clustering data yyyy-mm (included)
utc           = 12        # Clustering sampling time (utc) data is sampled daily
plevel        = 500       # Height level used (hPa)
lat_lims      = [-5,-50]  # South and North lattitude limit of analysis box
lon_lims      = [100,165] # West and East longitude limit of analysis box
Ncoarsen      = 6         # Coarsening factor in lat and lon direction
Nclusters     = 30        # Number of clusters
file_clusters = f"{root}SWT_fields/SWT_data_v1.nc"
path_out      = f"{root}reassignment_analysis/" # Directory to save cluster results
path_data     = "/g/data/rt52/era5/pressure-levels/reanalysis/"         # Era5 data directory
# %% Load clustering information (SWT)
nc = Dataset(file_clusters,'r')
cluster_series = nc.variables['clusterSeries'][:]
cluster_time = nc.variables['time'][:]
lat_cluster = nc.variables['latitude'][:]
lon_cluster = nc.variables['longitude'][:]
nc.close()
# %% Load ERA5 data at specified height level hpa and calculate cluster mean
cluster_mean = {}; data = {}
for varname in ['u','v']:
  print(f"Read data for {varname}")
  data[varname],time,lat,lon = era5.read_data(varname,date_start,date_end,utc,lat_lims,lon_lims,path_data,varname_path=varname,Ncoarsen=Ncoarsen,plevel=plevel,progress=True)
  print(f"Calculate cluster mean at {plevel}hpa for {varname}")
  field_sum,field_sum2,field_count,time_era5 = stat.compute_sums(data[varname],time,cluster_series,cluster_time,Nclusters)
  cluster_mean[varname] = stat.cluster_mean(field_sum,field_count)
# %% Assign days to new cluster centres
cluster_series_new = np.zeros((len(time)))
for it,t in enumerate(time):
  cluster_series_new[it] = kmeans.assign(data['u'][it,:,:],data['v'][it,:,:],cluster_mean['u'],cluster_mean['v'])
# %% save new cluster time series and cluster centres to output
fileout = path_out+f"SWT_data_reassigned_{plevel}_v1.nc"
kmeans.save(fileout,Nclusters,time,lat,lon,cluster_mean['u'],cluster_mean['v'],cluster_series_new,None,None)