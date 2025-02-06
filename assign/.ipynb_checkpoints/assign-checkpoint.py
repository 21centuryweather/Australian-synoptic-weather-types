import xarray as xr
import numpy as np
import sys
from scipy import interpolate
from pathlib import Path

def assign_grid_to_SWT(u,v,latname="latitude",lonname="longitude",
                             cluster_filename='SWT_fields/SWT_data.nc',
                             quiet=False):
    
    ''' Assign a 850hPa wind field to a Australian Synoptic Weather Type
     
    Created by Frans Liqui Lung (Monash University) and Michael A. Barnes (ARC CoE 21st Century Weather, Monash University)
    
    Parameters
    ----------
    
    u : xr.DataArray with shape (latitude,longitude) and dtype float64
        The 850hPa zonal wind field.
    v : xr.DataArray with shape (latitude,longitude) and dtype float64
        The 850hPa meridional wind field.
    
    Options
    -------
    
    latname : str (default: "latitude")
        Name of the latitude variable name in the input u and v DataArrays
    lonname : str (default: "longitude")
        Name of the longitude variable name in the input u and v DataArrays
    cluster_filename : str (default: "SWT_fields/SWT_data.nc")
        Location of the netcdf containing the clustered SWT fields.
    quiet : bool (default: False)
        If false, prints out the cluster_ID and SWT label, else simply and quietly returns the variables as specified.
    
    Returns
    -------
    label : str
        The Synoptic Weather Type code.
    cluster_ID : int
        The clusterID of the Synoptic Weather Type as per the SWT_data file. 
    '''
    
    if not Path(cluster_filename).is_file():
        sys.exit(f'The SWT data file cannot be located at {cluster_filename}.\nPoint to a different file location using the option cluster_filename')
    clusters=xr.open_dataset(cluster_filename)
    SWTs=clusters.SWT
    lat_cluster,lon_cluster = np.meshgrid(clusters.latitude, clusters.longitude, indexing='ij')

    ## Interpolate the input u-grid to the clustering grid
    if u["latitude"].values[0]>u["latitude"].values[-1]:
        u=u.reindex(latname=list(reversed(u[latname]))) ## Reverse the latitude so that it is increasing in value
    ulat=u.latitude.values; ulon=u.longitude.values
    f = interpolate.RegularGridInterpolator((ulat,ulon),u.u.values)
    u_int = f((lat_cluster,lon_cluster))

    ## Interpolate the input v-grid to the clustering grid
    if v["latitude"].values[0]>v["latitude"].values[-1]:
        v=v.reindex(latname=list(reversed(v[latname]))) ## Reverse the latitude so that it is increasing in value
    vlat=v.latitude.values; vlon=v.longitude.values
    f = interpolate.RegularGridInterpolator((vlat,vlon),v.v.values)
    v_int = f((lat_cluster,lon_cluster))

    cluster_ID = assign(u_int,v_int,clusters.clusterU.values,clusters.clusterV.values)
    label = str(SWTs.sel(clusterID=cluster_ID).values)

    if not quiet:
        print(f'Synoptic Weather Type: {label}')
        print(f'Cluster ID: {cluster_ID}')

    return label,cluster_ID

def assign(u,v,clusterU,clusterV):
    if(np.shape(u)!=np.shape(clusterU)[1:]): 
        sys.exit('Wind velocity field not the same shape as cluster field, interpolate data to cluster grid first.')
    return np.argmin(np.sum((u[None,:,:]-clusterU)**2+(v[None,:,:]-clusterV)**2,axis=(-1,-2)))+1