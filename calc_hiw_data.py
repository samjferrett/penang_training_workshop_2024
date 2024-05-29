'''
An example of calculating the pre-requisite heavy rainfall data for case_example.ipynb

Written by S Ferrett
'''

#imports and definitions
import xarray as xr
import numpy as np
import hybrid_functions as H

region = 'PM'
seas = 'SON'

#load precipitation data (can be retrieved from 
#https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDF_07/summary?keywords=%22IMERG%20final%22
gpm_dir = "/gws/nopw/j04/ncas_climate_vol1/users/bc917929/GPM_Final/GPM_day/"
pr_ncs = xr.open_mfdataset(f'{gpm_dir}gpm_day*.nc')

#mask data and calculate thresholds
PM_gpm = H.mask_gpm(pr_ncs.precipitationCal,region)
P95,area_thres,pr_hiw = H.find_areathres(PM_gpm,season=seas,pc=10)

#get into desirable format for nc file and save
pr_hiw = pr_hiw.to_dataset(name='pr_count')
pr_hiw.attrs = {'P95':P95.values,'area_thres':area_thres,'pc':10}
pr_hiw = pr_hiw.convert_calendar('standard')
pr_hiw.to_netcdf(f'./inputs/GPM_hiw_pc10_{region}_{seas}.nc')

