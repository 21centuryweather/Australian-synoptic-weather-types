import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import colormaps as cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, BoundaryNorm
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, LAND
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import metpy.calc as mpcalc
from metpy.units import units


import os
import sys
sys.path.append(os.environ['WXSYSLIBDIR'])

from utils.general.nci_utils import get_GADI_ERA5_filename

def plot_synoptic_ERA5(code,indt):
    plot_extent=[80,180,-60,-5]
    figsize=(11,8)
    fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw=dict(projection=ccrs.PlateCarree())) 

    fn = get_GADI_ERA5_filename("msl",indt,stream='hourly', level_type='single-levels')
    dsMSL=xr.open_dataset(fn).sel(time=indt,latitude=slice(-5,-60),longitude=slice(80,180)).msl
    fn = get_GADI_ERA5_filename("pv",indt,stream='hourly', level_type='potential-temperature')
    dsPV315=xr.open_dataset(fn).sel(time=indt,level=315,latitude=slice(-5,-60),longitude=slice(80,180)).pv
    fn = get_GADI_ERA5_filename("pv",indt,stream='hourly', level_type='potential-temperature')
    dsPV330=xr.open_dataset(fn).sel(time=indt,level=330,latitude=slice(-5,-60),longitude=slice(80,180)).pv
    fn = get_GADI_ERA5_filename("u",indt,stream='hourly', level_type='pressure-levels')
    dsU850=xr.open_dataset(fn).sel(time=indt,level=850,latitude=slice(-5,-60),longitude=slice(80,180)).u
    fn = get_GADI_ERA5_filename("v",indt,stream='hourly', level_type='pressure-levels')
    dsV850=xr.open_dataset(fn).sel(time=indt,level=850,latitude=slice(-5,-60),longitude=slice(80,180)).v
    fn = get_GADI_ERA5_filename("u",indt,stream='hourly', level_type='pressure-levels')
    dsU=xr.open_dataset(fn).sel(time=indt,level=300,latitude=slice(-5,-60),longitude=slice(80,180)).u
    fn = get_GADI_ERA5_filename("v",indt,stream='hourly', level_type='pressure-levels')
    dsV=xr.open_dataset(fn).sel(time=indt,level=300,latitude=slice(-5,-60),longitude=slice(80,180)).v
    fn = get_GADI_ERA5_filename("tcwv",indt,stream='hourly', level_type='single-levels')
    dsTCWV=xr.open_dataset(fn).sel(time=indt,latitude=slice(-5,-60),longitude=slice(80,180)).tcwv
    
    spd300=mpcalc.wind_speed(dsU * units('m/s'), dsV  * units('m/s')).metpy.dequantify()
    ujet300=dsU.where(spd300>40*0.514444)
    vjet300=dsV.where(spd300>40*0.514444)
    jet300=spd300.where(spd300>35*0.514444)
    
    uW850=dsU850.where(dsU850>=2).sel(latitude=slice(0,-20))
    vW850=dsV850.where(dsU850>=2).sel(latitude=slice(0,-20))
    
    ax.set_extent(plot_extent,crs=ccrs.PlateCarree())
    
    ax.add_feature(LAND,facecolor='lightgrey')
    ax.coastlines(linewidths=0.4)
    
    plot_levels = list(range(900,1060,2))
    c=ax.contour(dsMSL.longitude,dsMSL.latitude,dsMSL/100,
                     levels=plot_levels,colors='black',linewidths=0.75,linestyles='-',
                     transform=ccrs.PlateCarree())

    fmt = '%i'
    ax.clabel(c, c.levels, inline=True, fmt=fmt, fontsize=8) 
    
    plot_levels = range(35,80,5)
    cf=ax.contourf(jet300.longitude,jet300.latitude,jet300*1.94384,
                   levels=plot_levels,cmap='Blues',extend='max',
                   transform=ccrs.PlateCarree(),alpha=0.5)
    
    plot_levels = range(48,1000,500)
    cf=ax.contourf(dsTCWV.longitude,dsTCWV.latitude,dsTCWV,
                   levels=plot_levels,cmap='Greens',extend='max',
                   transform=ccrs.PlateCarree(),alpha=0.4)
    
    c=ax.contour(dsTCWV.longitude,dsTCWV.latitude,dsTCWV,
                     levels=[48],colors='limegreen',linewidths=0.8,linestyles='-',
                     transform=ccrs.PlateCarree())
    
    nvec=10
    q=ax.barbs(uW850.longitude[::nvec],uW850.latitude[::nvec], 
                 uW850[::nvec,::nvec], vW850[::nvec,::nvec],color='darkgreen',
                 length=6,
                  #scale=10,width=0.0008,color='darkgreen',#pivot='tail',
                  transform=ccrs.PlateCarree())
    nvec=10
    q=ax.quiver(jet300.longitude[::nvec],jet300.latitude[::nvec], 
                 ujet300[::nvec,::nvec], vjet300[::nvec,::nvec],
                  scale=800,width=0.004,
                  transform=ccrs.PlateCarree(), alpha=0.5)
    
    c=ax.contour(dsPV330.longitude,dsPV330.latitude,dsPV330*1e6,
                     levels=[-2],colors='magenta',linewidths=2.0,linestyles='--',
                     transform=ccrs.PlateCarree())
    c=ax.contour(dsPV315.longitude,dsPV315.latitude,dsPV315*1e6,
                     levels=[-2],colors='magenta',linewidths=2.25,linestyles='-',
                     transform=ccrs.PlateCarree())
    
    label=ax.text(82,-9.5,f"{indt.strftime('%Y-%m-%d')} | {code}",ha='left',va='center',fontsize=20,
                          bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='black'),
                          zorder=10)
    
    gl = ax.gridlines(transform=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.4, color='lightgrey', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.left_labels = False

    plt.show()