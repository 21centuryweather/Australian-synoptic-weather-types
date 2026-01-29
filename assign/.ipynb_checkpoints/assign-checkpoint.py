import xarray as xr
import numpy as np
import sys
from scipy import interpolate
from pathlib import Path

def assign_grid_to_SWT(u,v,latname="latitude",lonname="longitude",
                             cluster_filename='SWT_fields/SWT_data.nc',
                             interpolation=True,quiet=False,silence_warning=False):
    
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
    interpolate : bool (default: True)
        If false, does not do any interpolation on the grid. Note: if false, grids must be exactly the same.
    quiet : bool (default: False)
        If false, prints out the SWT label, else simply and quietly returns the variables as specified.
    silence_warning : bool (default: False)
        If false, prints out any defined warnings (e.g. if NaNs are present), otherwise silences them.
    
    Returns
    -------
    SWT : str
        The Synoptic Weather Type code.
    '''
    
    if not Path(cluster_filename).is_file():
        sys.exit(f'The SWT data file cannot be located at {cluster_filename}.\nPoint to a different file location using the option cluster_filename')
    clusters=xr.open_dataset(cluster_filename)
    lat_cluster,lon_cluster = np.meshgrid(clusters.latitude, clusters.longitude, indexing='ij')

    if np.isnan(u.values).any() or np.isnan(v.values).any():
        if not silence_warning:
            print('Warning! NaN-values detected! Filling missing values before assignment!')
        u = u.interpolate_na(dim=latname).interpolate_na(dim=lonname)
        v = v.interpolate_na(dim=latname).interpolate_na(dim=lonname)

    if interpolation:
        ## Interpolate the input u-grid to the clustering grid
        if u[latname].values[0]>u[latname].values[-1]:
            u=u.reindex(latname=list(reversed(u[latname]))) ## Reverse the latitude so that it is increasing in value
        ulat=u[latname].values; ulon=u[lonname].values
        f = interpolate.RegularGridInterpolator((ulat,ulon),u.values)
        u_int = f((lat_cluster,lon_cluster))
    
        ## Interpolate the input v-grid to the clustering grid
        if v[latname].values[0]>v[latname].values[-1]:
            v=v.reindex(latname=list(reversed(v[latname]))) ## Reverse the latitude so that it is increasing in value
        vlat=v[latname].values; vlon=v[lonname].values
        f = interpolate.RegularGridInterpolator((vlat,vlon),v.values)
        v_int = f((lat_cluster,lon_cluster))
    else:
        u_int= u.u.values
        v_int= v.v.values

    cluster_ID = assign(u_int,v_int,clusters.clusterU.values,clusters.clusterV.values)
    SWT = clusters.SWT.isel(SWT=cluster_ID).item()

    if not quiet:
        print(f'Synoptic Weather Type: {SWT}')

    return SWT

def assign(u,v,clusterU,clusterV):
    if(np.shape(u)!=np.shape(clusterU)[1:]): 
        sys.exit('Wind velocity field not the same shape as cluster field, interpolate data to cluster grid first.')
    return np.argmin(np.sum((u[None,:,:]-clusterU)**2+(v[None,:,:]-clusterV)**2,axis=(-1,-2)))