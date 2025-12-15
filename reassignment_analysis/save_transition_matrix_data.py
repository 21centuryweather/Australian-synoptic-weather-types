# %% Import packages
import sys
sys.path.append('../utils/') # Path to scripts
import numpy as np
from netCDF4 import Dataset
import compute_statistics as stat
# %% Paths, data and input
root          = "/home/565/fl2086/Australian-synoptic-weather-types/" 
date_start    = '1959-01-01' # Start year-month for clustering data yyyy-mm
date_end      = '2024-01-01' # End year-month for clustering data yyyy-mm (excluded)
affix         = 'extratropics' # Affix to select reassigned data: extratropics, tropics, 500hPa
file_clusters            = f"{root}SWT_fields/SWT_data_v1.nc"
file_clusters_reassigned = f"{root}reassignment_analysis/SWT_data_reassigned_{affix}_v1.nc"
path_out  = f"{root}reassignment_analysis/"
path_data = f"{root}example_data/era5_data/"
varnames = ['msl','tcwv','pv_315K','pv_330K','u_850hPa','v_850hPa','u_300hPa','v_300hPa','u_500hPa','v_500hPa','z_500hPa']
mode = {varname : 'mean' for varname in varnames}
mode['z_500hPa'] = 'pert'
Nlat = 31
Nlon = 44
# %% Load cluster data
nc = Dataset(file_clusters,'r')
surface_clusters = nc.variables['clusterSeries'][:]
WR = nc.variables['WR'][:]
SWT = nc.variables['SWT'][:]
clusterID = nc.variables['clusterID'][:]
Nclusters = len(SWT)
time = nc.variables['time'][:]
time_clusters = np.array([np.datetime64('1900-01-01')+np.timedelta64(int(t),'h') for t in time])
nc.close()
nc = Dataset(file_clusters_reassigned,'r')
upper_clusters = nc.variables['clusterSeries'][:]
nc.close()
mask_clusters = (time_clusters>=np.datetime64(date_start)) & (time_clusters<np.datetime64(date_end))
surface_clusters = surface_clusters[mask_clusters]
upper_clusters = upper_clusters[mask_clusters]
time_clusters = time_clusters[mask_clusters]
time = time[mask_clusters]
# %% Calculate transitions between surface and upper clusters for both WR and SWT
transition_matrix = np.zeros((Nclusters,Nclusters))
for isurf,iupper in zip(surface_clusters,upper_clusters):
  transition_matrix[isurf-1,iupper-1] += 1
# %% Get composites
data_out = {var : np.zeros((Nclusters,Nclusters,Nlat,Nlon)) for var in varnames}
for i in range(Nclusters):
  for j in range(Nclusters):
    # Load necessary data
    data = {}
    mask = (surface_clusters == clusterID[i]) & (upper_clusters == clusterID[j])
    for varname in varnames:
      if mode[varname] == 'mean':
        nc = Dataset(f"{path_data}{varname}.nc",'r')
        data_out[varname][i,j,:,:] = np.mean(nc.variables[varname.split('_')[0]][mask,:,:],axis=0)
        nc.close()
      elif mode[varname] == 'pert':
        nc = Dataset(f"{path_data}{varname}.nc",'r')
        temp = nc.variables[varname.split('_')[0]][:,:,:]
        temp_transitionIDs = np.where(mask,0,1)
        field_sum,field_sum2,field_count,_ = stat.compute_sums(temp,time,temp_transitionIDs,time,2)
        data_out[varname][i,j,:,:] = stat.cluster_daily_pert(field_sum,field_count)[0,:,:]
      else:
        sys.exit(f"{mode[varname]} for {varname} invalid mode")

nc = Dataset(f"{path_data}{varname}.nc",'r')
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
nc.close()
# %% Save to netcdf file
nc = Dataset(f"{path_out}data_interactive_transition_matrix_{affix}.nc",'w')
nc.createDimension('clusterID_original',Nclusters)
nc.createDimension('clusterID_reassigned',Nclusters)
nc.createDimension('latitude',len(lat))
nc.createDimension('longitude',len(lon))
nc_cID = nc.createVariable('clusterID_original','i4',('clusterID_original',)); nc_cID[:] = clusterID
nc_cID = nc.createVariable('clusterID_reassigned','i4',('clusterID_reassigned',)); nc_cID[:] = clusterID
nc_SWT = nc.createVariable('SWT_original',str,('clusterID_original',)); nc_SWT[:] = SWT
nc_WR = nc.createVariable('WR_original',str,('clusterID_original',)); nc_WR[:] = WR
nc_SWT = nc.createVariable('SWT_reassigned',str,('clusterID_reassigned',)); nc_SWT[:] = SWT
nc_WR = nc.createVariable('WR_reassigned',str,('clusterID_reassigned',)); nc_WR[:] = WR    
nc_lat  = nc.createVariable('latitude','f4',('latitude',)); nc_lat[:] = lat
nc_lon  = nc.createVariable('longitude','f4',('longitude',)); nc_lon[:] = lon
for var in varnames:
  nc_var = nc.createVariable(var,'f4',('clusterID_original','clusterID_reassigned','latitude','longitude'))
  nc_var[:,:,:,:]  = data_out[var]
nc_var = nc.createVariable('Tmatrix','i4',('clusterID_original','clusterID_reassigned')); nc_var[:] = transition_matrix
nc.close()
# %%
