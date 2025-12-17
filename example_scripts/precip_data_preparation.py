# %% Import packages
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta
# %% Input
root         = '/home/565/fl2086/Australian-synoptic-weather-types/'
cluster_file = f"{root}SWT_fields/SWT_data_v1.nc"
path_precip    = "/g/data/zv2/agcd/v1/precip/total/r005/01day/"
path_out     = f"{root}results/"
start_year   = 1952
end_year     = 2019
timezone     = timedelta(hours=11)
months       = [9,10,11]
Ncoarsen     = 1
month_labels = {1:'J',2:'F',3:'M',4:'A',5:'M',6:'J',7:'J',8:'A',9:'S',10:'O',11:'N',12:'D'}
month_string = ''
for m in months: month_string+=month_labels[m]
# %% Load cluster information
with Dataset(cluster_file,'r') as nc:
  clusterSeries = nc.variables['clusterSeries'][:]
  clusterID     = nc.variables['clusterID'][:]
  SWT           = nc.variables['SWT'][:]
  WR            = nc.variables['WR'][:]
  time_clusters = nc.variables['time'][:]
  Nclusters     = len(clusterID)
time_clusters = np.array([datetime(1900,1,1)+timedelta(hours=int(t))+timezone for t in time_clusters])
mask          = [t.month in months and t.year>=start_year and t.year<=end_year for t in time_clusters]
time_clusters = time_clusters[mask]
clusterSeries = clusterSeries[mask]
# %% Load AGCD data max temperature and total daily rainfall data
precip = []
time_precip = []
for year in np.arange(start_year,end_year+1):
  filename = f"{path_precip}agcd_v1_precip_total_r005_daily_{int(year)}.nc"
  with Dataset(filename,'r') as nc:
    precip.append(nc.variables['precip'][:,::Ncoarsen,::Ncoarsen])
    lat = nc.variables['lat'][::Ncoarsen]
    lon = nc.variables['lon'][::Ncoarsen]
    time_precip.append(nc.variables['time'][:])
precip = np.concatenate(precip)
time_precip = np.concatenate(time_precip).astype(int)
time_precip = np.array([datetime(1850,1,1)+timedelta(days=int(t)) for t in time_precip])
mask      = [t.month in months and t.year>=start_year and t.year<=end_year for t in time_precip]
precip      = precip[mask]
time_precip = time_precip[mask]
print('Loaded data')
# %% Get the 5% wettest days
k = max(1, int(np.ceil(0.05 * len(time_precip))))
top_idx = np.argpartition(precip, -k, axis=0)[-k:]
clusters_of_extremes = clusterSeries[top_idx]
cluster_counts = np.zeros((Nclusters, len(lat), len(lon)), dtype=int)
for c in range(Nclusters):
    cluster_counts[c] = np.sum(clusters_of_extremes == c+1, axis=0)
    probs = np.array([np.sum(clusterSeries==c+1) for c in range(Nclusters)])/len(clusterSeries)
print("Got wettest days")
# %% Save information
outfile = f"precip_cluster_data_{month_string}.nc"
with Dataset(f"{path_out}{outfile}", "w", format="NETCDF4") as nc:
    # Dimensions
    nc.createDimension("k", k)     # number of extreme days per gridcell
    nc.createDimension("lat", len(lat))
    nc.createDimension("lon", len(lon))
    nc.createDimension("clusters", Nclusters)
    # Coordinates
    lat_var = nc.createVariable("lat", "f4", ("lat",))
    lon_var = nc.createVariable("lon", "f4", ("lon",))
    lat_var[:] = lat
    lon_var[:] = lon
    # Variables
    v_precip = nc.createVariable("precip_extremes", "f4",
                               ("k", "lat", "lon"))
    v_precip[:] = np.take_along_axis(precip,top_idx,axis=0)
    v_clusters = nc.createVariable("clusters_extremes", "i4",
                               ("k", "lat", "lon"))
    v_clusters[:] = clusters_of_extremes
    v_time = nc.createVariable("time_extremes", "i8",
                               ("k", "lat", "lon"))
    v_time[:] = ((time_precip[top_idx]-datetime(1850,1,1))/timedelta(days=1)).astype("int64")
    v_time.units = 'Days since 1850-01-01'
    v_counts = nc.createVariable("cluster_counts", "i4",
                                 ("clusters", "lat", "lon"))
    v_counts[:] = cluster_counts
    v_probs = nc.createVariable("probs", "f4", ("clusters",))
    v_probs[:] = probs
    v_swt = nc.createVariable("SWT", str, ("clusters",))
    v_swt[:] = SWT
    nc.description = f"Clusters for 5% wettest days during {month_string} of {start_year}-{end_year+1}"