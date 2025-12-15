"""
File: kmeans_clustering.py
Author: Frans Liqui Lung
Date: 2024-07-15
Description: Cluster era5 wind data.
Contains:
clusterU,clusterV,labels,inertia,silhouette = cluster(uv_stack,latitude,longitude,Nclusters,n_jobs=1)
    Description: Create clusters from wind data.
    Input:
        uv_stack: stacked u,v wind data per point (float), output of stack(u,v)
        latitude: latitude of original unstacked data (float)
        longitude: longitude of original unstacked data (float)
        Nclusters: number of clusters to form (int)
    Optional input:
        n_jobs: number of processes to use (int)
    Output:
        clusterU: u wind component for clusters (float)
        clusterV: v wind component for clusters (float)
        labels: timeseries of cluster numbers (int)
        inertia: inertia for resulting clusters (float)
        silhouette: silhouette score for resulting clusters (float)
uv_stack = stack(u,v)
    Description: stack u and v data per point
    Input:
        u: u wind component (float)
        v: v wind component (float)
    Output:
        uv_stack: stacked data (float)
save(fileout,Nclusters,time,lat,lon,clusterU,clusterV,labels,inertia,silhouette)
    Description: write cluster information to netcdf file
    Input:
        fileout: path + filename where to save data "path/filename.nc" (str)
        Nclusters: number of clusters (int)
        time: time of data used (float)
        lat: latitude of unstacked data (float)
        lon: longitude of unstacked data (float)
        clusterU: u wind component for clusters (float)
        clusterV: v wind component for clusters (float)
        labels: timeseries of cluster numbers (int)
        inertia: inertia for resulting clusters (float)
        silhouette: silhouette score for resulting clusters (float)
label = assign(u,v,clusterU,clusterV)
    Description: assigns cluster number to given wind field.
    Input:
        u: u wind component (float)
        v: v wind component (float)
        clusterU: u wind component for clusters (float) 
        clusterV: v wind component for clusters (float)
    Output:
        label: cluster number for given wind field (int)
"""

import numpy as np
from netCDF4 import Dataset
from datetime import datetime,timedelta
import warnings
warnings.filterwarnings('ignore')
from joblib import parallel_backend
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import sys

def cluster(uv_stack,latitude,longitude,Nclusters,n_jobs=1):
    # Cluster stacked uv data
    nlats = latitude.size
    nlons = longitude.size
    clusterU=np.zeros((Nclusters,nlats,nlons))
    clusterV=np.zeros((Nclusters,nlats,nlons))
    with parallel_backend('threading', n_jobs=n_jobs):
        uv_clusters = KMeans(n_clusters=Nclusters, random_state=0).fit(uv_stack)
    uv_centroids=uv_clusters.cluster_centers_
    inertia = uv_clusters.inertia_
    labels=uv_clusters.labels_+1
    # Calculate silhouette score
    silhouette = silhouette_score(uv_stack, labels)
    # Reshape clusters
    for nc in range(0,Nclusters):
        clusterU[nc,:,:]=uv_centroids[nc,0:nlats*nlons].reshape((nlats,nlons))
        clusterV[nc,:,:]=uv_centroids[nc,nlats*nlons:2*nlats*nlons].reshape((nlats,nlons))
    return clusterU,clusterV,labels,inertia,silhouette

def stack(u,v):
    # Stack u and v data per grid point
    Ntime,Nlat,Nlon = np.shape(u)
    uv_stack = np.zeros((Ntime,Nlat*Nlon*2))
    for it in range(Ntime):
        uv_stack[it,:] = np.concatenate([u[it,:,:].flatten(),v[it,:,:].flatten()])
    return uv_stack

def save(fileout,Nclusters,time,lat,lon,clusterU,clusterV,labels,inertia,silhouette):
    # Save cluster information to output
    nc = Dataset(fileout,'w')
    nc.createDimension('time',None)
    nc.createDimension('latitude',len(lat))
    nc.createDimension('longitude',len(lon))
    nc.createDimension('clusterID',Nclusters)
    nc_time = nc.createVariable('time','f4',('time',)); nc_time[:] = time
    nc_lat  = nc.createVariable('latitude','f4',('latitude',)); nc_lat[:] = lat
    nc_lon  = nc.createVariable('longitude','f4',('longitude',)); nc_lon[:] = lon
    nc_cID  = nc.createVariable('clusterID','i4',('clusterID',)); nc_cID[:] = np.arange(1,Nclusters+1,dtype=np.int16)
    nc_cIDs = nc.createVariable('clusterSeries','i4',('time',)); nc_cIDs[:] = labels
    nc_u    = nc.createVariable('clusterU','f4',('clusterID','latitude','longitude')); nc_u[:,:,:] = clusterU
    nc_v    = nc.createVariable('clusterV','f4',('clusterID','latitude','longitude')); nc_v[:,:,:] = clusterV
    nc_in   = nc.createVariable('inertia','f4'); nc_in[:] = inertia
    nc_sil  = nc.createVariable('silhouette_score','f4'); nc_sil[:] = silhouette
    # Add meta data
    nc_time.units = "hours since 1900-01-01 00:00:00.0"
    nc_time.long_name = "time"
    nc_time.calendar = "gregorian"
    nc_lon.units = "degrees_east"
    nc_lon.long_name = "longitude"
    nc_lat.units = "degrees_north"
    nc_lat.long_name = "latitude"
    nc_cIDs.long_name = "time series of cluster IDs"
    nc_in.long_name = "inertia"
    nc_sil.long_name = "silhouette score"
    nc.close()

def assign(u,v,clusterU,clusterV):
    if(np.shape(u)!=np.shape(clusterU)[1:]): 
        sys.exit('Wind velocity field not the same shape as cluster field, interpolate data to cluster grid first.')
    return np.argmin(np.sum((u[None,:,:]-clusterU)**2+(v[None,:,:]-clusterV)**2,axis=(-1,-2)))+1