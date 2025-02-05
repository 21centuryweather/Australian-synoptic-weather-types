""" 
File: read_era5.py
Author: Frans Liqui Lung
Date: 2024-07-15
Description: Daily samples era5 data on GADI.
Contains:
varout,time,lat,lon = read_data(varname,date_start,date_end,utc,lat_lims,lon_lims,path_data,varname_path=None,Ncoarsen=1,plevel=None,progress=True)
    Description: Daily samples era5 data on GADI.
    Input:
        varname: name of variable to read as defined in netcdf file (str)
        date_start: start date of data to read in format 'yyyy-mm' (str)
        date_end: end date of data to read in format 'yyyy-mm' (str), end date is included
        utc: utc hour of the day at which the data is sampled (int)
        lat_lims: latitude limits of data to read [minimum,maximum] (float), both are included
        lon_lims: longitude limits of data to read [minimum,maximum] (float), both are included
        path_data: path to reanalysis data (str)
    Optional input:
        Ncoarsen: coarsening factor for data to read, data is sampled every Ncoarsen steps in space (int)
        plevel: pressure level at which to read 3d variables (int)
        varname_path: name of variable in reanalysis path, needed if different from variable name in netcdf (str)
        progress: print progress (boolean)
    Output:
        varout: data [time,lat,lon] (float)
        time: time data (float)
        lat: latitude data (float)
        lon: longitude data (float)
"""
import glob
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import sys

def read_data(varname,date_start,date_end,utc,lat_lims,lon_lims,path_data,varname_path=None,Ncoarsen=1,plevel=None,progress=True):
    # Read era5 data
    if(varname_path == None): varname_path = varname
    # Find filenames to include
    if(progress): 
        sys.stdout.write(f"\rFinding files to include")
        sys.stdout.flush()
    filenames_included = get_filenames(varname_path,date_start,date_end,path_data)
    # Determine number of timesteps from filenames (Daily sampling)
    date_start,_ = date_from_filename(filenames_included[0])
    _,date_end   = date_from_filename(filenames_included[-1])
    Ntime        = int((date_end-date_start).total_seconds()/3600/24)+1
    # Read spatial dimensions from first file and create spatial mask
    nc    = Dataset(filenames_included[0],'r')
    lat   = nc.variables['latitude'][:]
    lon   = nc.variables['longitude'][:]
    nc.close()
    mask_lat = np.argwhere((lat>=np.min(lat_lims)) & (lat<=np.max(lat_lims))).flatten()
    mask_lon = np.argwhere((lon>=np.min(lon_lims)) & (lon<=np.max(lon_lims))).flatten()
    lat      = lat[mask_lat[::Ncoarsen]]; Nlat = len(lat)
    lon      = lon[mask_lon[::Ncoarsen]]; Nlon = len(lon)
    # Preallocate output
    if(progress): 
        sys.stdout.write(f"\rPreallocate output")
        sys.stdout.flush() 
    varout   = np.zeros([Ntime,Nlat,Nlon])
    time     = np.zeros(Ntime)
    ts = 0
    # Loop over files and read data
    for filename in filenames_included:
        file_start,_ = date_from_filename(filename)
        nc    = Dataset(filename,'r')
        if(progress): 
            sys.stdout.write(f"\rReading data for {file_start.strftime('%Y/%m')}")
            sys.stdout.flush()
        if len(np.shape(nc.variables[varname]))==3: # Single level variables
            field = nc.variables[varname][utc::24,mask_lat[::Ncoarsen],mask_lon[::Ncoarsen]]
            Ntime = np.shape(field)[0]
            varout[ts:ts+Ntime,:,:] = field
        elif len(np.shape(nc.variables[varname]))==4: # Pressure level variables
            if plevel==None: sys.exit('Need to specify pressure level for 4D data')
            ilevel = np.argwhere(nc.variables['level'][:]==plevel).flatten()
            if len(ilevel)==0: sys.exit('Pressure level not available')
            field = nc.variables[varname][utc::24,ilevel[0],mask_lat[::Ncoarsen],mask_lon[::Ncoarsen]]
            Ntime = np.shape(field)[0]
            varout[ts:ts+Ntime,:,:] = field
        else: sys.exit('Unvalid data shape')
        time[ts:ts+Ntime] = nc.variables['time'][utc::24]
        ts += Ntime
        nc.close()
    if(progress): 
        sys.stdout.write(f"\rReading data done                         \n")
        sys.stdout.flush()
    return varout,time,lat,lon

def get_ncName(filename):
    # Get filename from 'path/filename'
    return filename.split('/')[-1]

def get_filenames(varname_path,date_start,date_end,path_data):
    # Find filenames in range date_start-date_end
    filenames = glob.glob(f"{path_data}{varname_path}/*/*.nc")
    filenames.sort(key=get_ncName)
    filenames_included = []
    for filename in filenames:
        file_start,file_end = filename[:-3].split('_')[-1].split('-')
        if(datetime.strptime(file_start,'%Y%m%d')>datetime.strptime(date_end,'%Y-%m') or \
           datetime.strptime(file_end,'%Y%m%d')<datetime.strptime(date_start,'%Y-%m') ): 
            continue # Don't open file if not in time range
        filenames_included.append(filename)
    return filenames_included

def date_from_filename(filename):
    # Extract date from filename
    date_start,date_end = filename[:-3].split('_')[-1].split('-')
    return datetime.strptime(date_start,'%Y%m%d'), datetime.strptime(date_end,'%Y%m%d')