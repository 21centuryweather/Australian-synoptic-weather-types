# %% Import packages
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import sys
sys.path.append('../utils/')
from read_era5 import read_data
# %% Input
root = '/home/565/fl2086/Australian-synoptic-weather-types/'
varnames   = ['u','v']
date_start = '1952-01' # Start year-month yyyy-mm
date_end   = '2023-12' # End year-month yyyy-mm (included)
utc        = 12        # Reading sampling time (utc) data is sampled daily
plevel     = 850       # Height level used (hPa)
lat_lims   = [-5,-50]  # South and North lattitude limit of analysis box
lon_lims   = [100,165] # West and East longitude limit of analysis box
Ncoarsen   = 6         # Coarsening factor in lat and lon direction
path_out   = f"{root}example_data/era5_data/" # Directory to save cluster results
path_data  = "/g/data/rt52/era5/pressure-levels/reanalysis/" # Era5 data directory
# %% Read and save era5 data
for varname in varnames:
  _,_,_,_=read_data(varname,date_start,date_end,utc,lat_lims,lon_lims,path_data,Ncoarsen=Ncoarsen,level=plevel,progress=True,save=True,file_out=f"{path_out}{varname}_full.nc")