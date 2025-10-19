""" 
File: statistics.py
Author: Frans Liqui Lung
Date: 2024-07-15
Description: Calculates climatological statistics.
Contains:
field_sum,field_sum2,field_count,time = compute_sums(data,time,cluster_timeseries,cluster_time,Nclusters,date0=datetime(1900,1,1))
    Description: Calculate sampled sums and sums of squares from input data, which are used to calculate the statistics.
    Input:
        data(time,lat,lon): data for which to compute sums (float)
        time: time of data (float)
        cluster_timeseries: cluster numbers for each day (int)
        cluster_time: time for cluster_timeseries (float)
        Nclusters: number of clusters (int)
    Optional input:
        date0: reference date (datetime)
    Output:
        field_sum(Nclusters,doy,lat,lon): sampled sums, sampled over cluster number, day of the year (float)
        field_sum2(Nclusters,doy,lat,lon): sampled sums of data squared, sampled over cluster number and day of the year (float)
        field_count(Nclusters,doy): number of days present in each sample (int)
        time: time epochs included in the summation (float)
daily_mean(field_sum,field_count)
    Description: calculate daily mean.
daily_var(field_sum,field_sum2,field_count)
    Description: calculate daily variance.
cluster_mean(field_sum,field_count)
    Description: calculate cluster mean
cluster_var(field_sum,field_sum2,field_count)
    Description: calculate cluster variance
cluster_daily_pert(field_sum,field_count):
    Description: calculate the cluster mean daily perturbation
cluster_frequency(field_count):
    Description: calculate cluster frequency
cluster_frequency_monthly(field_count):
    Description: calculate monthly cluster frequency
"""

import numpy as np
from datetime import datetime,timedelta
from netCDF4 import Dataset

def compute_sums(data,time,cluster_timeseries,cluster_time,Nclusters,date0=datetime(1900,1,1)):
   # Calculate sampled sums and sums of squares from input data
   (Ntime,Nlat,Nlon) = np.shape(data)
   time_intersect, it_data, it_cluster = np.intersect1d(time, cluster_time, return_indices=True)
   if(len(time_intersect)!=Ntime or len(time_intersect)!=len(cluster_time)):
      time = time[it_data]
      data = data[it_data,:,:]
      cluster_time = cluster_time[it_cluster]
      cluster_timeseries = cluster_timeseries[it_cluster]
      date_start = (date0+timedelta(hours=int(time[0]))).strftime('%Y/%m/%d')
      date_end   = (date0+timedelta(hours=int(time[-1]))).strftime('%Y/%m/%d')
      Ntime = len(time_intersect)
      print(f"Data and cluster timeseries do not span same period, analysis for common period: {date_start} - {date_end}.")
   # Preallocate variables
   field_count= np.zeros([Nclusters,366])
   field_sum  = np.zeros([Nclusters,366,Nlat,Nlon])
   field_sum2 = np.zeros([Nclusters,366,Nlat,Nlon])
   idx_dict = doy_dictionary()
   # Loop over epochs
   for it in range(Ntime):
      if(time[it]!=cluster_time[it]): exit('Mismatch between time from cluster timeseries and time from data')
      date = date0+timedelta(hours=int(time[it]))
      doy = date.strftime("%m/%d")
      idx = idx_dict[doy]
      clusterID = int(cluster_timeseries[it])
      field = data[it,:,:]
      field_sum[clusterID-1,idx,:,:]  += field
      field_sum2[clusterID-1,idx,:,:] += field**2
      field_count[clusterID-1,idx] += 1
   return field_sum,field_sum2,field_count,time

def daily_mean(field_sum,field_count):
  return np.sum(field_sum, axis=0)/np.sum(field_count,axis=0)[:,None,None]

def daily_var(field_sum,field_sum2,field_count):
  return np.sum(field_sum2,axis=0)/np.sum(field_count,axis=0)[:,None,None] - daily_mean(field_sum,field_count)**2

def cluster_mean(field_sum,field_count):
  return np.sum(field_sum,axis=1)/np.sum(field_count,axis=1)[:,None,None]

def cluster_var(field_sum,field_sum2,field_count):
  return np.sum(field_sum2,axis=1)/np.sum(field_count,axis=1)[:,None,None] - cluster_mean(field_sum,field_count)**2

def cluster_daily_pert(field_sum,field_count):
  out = cluster_mean(field_sum,field_count)
  day_mean = daily_mean(field_sum,field_count)
  for icluster in range(np.shape(field_sum)[0]):
    Nsamples = np.sum(field_count[icluster,:])
    for iday in range(366):
      out[icluster,:,:] -= field_count[icluster,iday]/Nsamples*day_mean[iday,:,:]
  return out

def cluster_frequency(field_count):
  return np.sum(field_count,axis=1)

def cluster_frequency_monthly(field_count):
  # Calculate cluster occurrancy per month
  day_in_month = [31,29,31,30,31,30,31,31,30,31,30,31]
  out   = np.zeros([np.shape(field_count)[0],12])
  idx_array    = [0,*np.cumsum(day_in_month)]
  for i in range(12):
    out[:,i] = out[:,i] + np.sum(field_count[:,idx_array[i]:idx_array[i+1]],axis=1)
  return out

def doy_dictionary():
   # Dictionary to reformat 'yyyy-mm' to day number
   leapyear = 2020
   return dict((datetime.strptime(f"{leapyear}/{iday:02d}","%Y/%j").strftime("%m/%d"), int(iday-1)) for iday in np.arange(1,367))