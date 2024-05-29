import xarray as xr
import numpy as np
import xesmf as xe
from dask.diagnostics import ProgressBar
from scipy import interpolate
import pandas as pd

def find_areathres(gpm_sub,pc=10,season='DJF'):
    '''
    Calculates the local 95th percentile rainfall threshold and the area threshold required for pc,
    thresholds based data up to 2001-2014

    Inputs
    ---
    gpm_sub    - xarray DataSet/DataArray (time, lat, lon): rainfall data that has been masked accordingly (for our analysis this is GPM)
    pc         - int: percentile of required HIW event, default 10 i.e. event occurs 10% of the time
    season     - str: season required, can be one of 'DJF', 'JJA', 'MAM', 'JJA', 'OND'

    Outputs
    ---
    pc95       - xarray DataSet/DataArray (lat,lon): gridpoint threshold of 95th percentile rainfall 
    area_thres - float: the amount of the area that needs to exceed pc95 for HIW
    '''
    
    #get the defined seasonal data (one of DJF, JJA, MAM, JJA, OND)
    if season == 'OND':
        gpm_sub = gpm_sub.sel(time = np.isin(gpm_sub.time.dt.month,[10,11,12]))
    else:
        gpm_sub = gpm_sub.sel(time=gpm_sub.time.dt.season==season)

    #calculate 95th percentile over the region
    pc95 = gpm_sub.sel(time=slice(None,'2014')).chunk({'time':-1}).quantile(0.95,dim=['time','lat','lon'])
    #optional ProgressBar
    with ProgressBar():
        pc95 = pc95.compute()
    #without ProgressBar
    #pc95 = pc95.compute()
    
    #find the area threshold required for the relevant time percentile threshold defined by pc
    gpm_bool = gpm_sub>pc95
    gpm_bool = xr.where(np.isnan(gpm_sub),np.nan,gpm_bool)
    bool_count = gpm_bool.mean(['lat','lon']).compute()
    def hiw_calc(area_thres=.15,end='2014'):
        hiw_ind = (bool_count.sel(time=slice(None,end))>area_thres)
        return hiw_ind

    hiw_pc = list(map(lambda at: hiw_calc(at).mean('time').values,np.arange(0.01,1,0.01)))
    x=np.arange(0.01,1,0.01)
    f = interpolate.interp1d(hiw_pc, x)

    area_thres=f(pc/100)

    #return the spatial threshold and the area threshold
    return pc95,area_thres,bool_count

def mask_gpm(gpm_nc,region='PM',regrid=False,res=1.):
    '''
    Mask rainfall data except the provided region. Regridding optional

    Inputs
    ---
    gpm_nc     - xarray DataSet/DataArray (time, lat, lon): rainfall data
    region     - str: region string, see documentation for region definitions
    regrid     - bool: True/False if regridding to a res x res degree grid should occur
    res        - float: resolution for regridding in degrees, default 1 x 1 grid
    
    Outputs
    ---
    GPM_masked - xarray DataSet/DataArray (time,lat,lon): rainfall data that has been masked 
                 with np.nan everywhere except land points in region
    '''
    #boxes for masking
    box_points = {'PM':[100,105,2.5,7.5],
        'EM':[110,120,2.5,7.5],
        'WI':[100,105,-10,2.5],
        'SI':[105,117,-10,-5],
        'NI':[107.5,118,-5,2.5],
        'EI':[118,125,-5,2.5],
        'SP':[120,130,5,12.5],
        'NP':[120,130,12.5,20],
        'SV':[104,110,8,11.5],
        'CV':[107,110,11.5,16.5],
        'NV':[103,108,17,23],
        'ST':[98,101,6.5,14],
        'NT':[99,105,14,19],
        'ING':[130,140,-10,0],
        'PNG':[140,150,-10,0]}[region]

    #subset rainfall data to SE region and load mask
    GPM_nc = gpm_nc.sel(lat=slice(-25,25),lon=slice(85,160))
    ws_path = '/gws/nopw/j04/ncas_climate_vol1/users/bc917929/'
    mask_nc = xr.open_dataset(f"{ws_path}GPM_Final/IMERG_land_sea_mask.nc").landseamask.sel(lat=slice(-25,25),lon=slice(85,160))

    # masks for the diagonal regions PM and WI
    if region == 'PM':
        mask_nc = xr.where(mask_nc.lat>(-.85*(mask_nc.lon-106)+4),100,mask_nc)
        mask_nc = xr.where(mask_nc.lat<(-.85*(mask_nc.lon-106)-1),100,mask_nc)
        mask_nc = xr.where(np.logical_or(mask_nc.lon<97,mask_nc.lon>106),100,mask_nc)
    elif region =='WI':
        mask_nc = xr.where(mask_nc.lat>(-.85*(mask_nc.lon-106)-1),100,mask_nc)
        mask_nc = xr.where(mask_nc.lat<(-.85*(mask_nc.lon-106)-7.5),100,mask_nc)
        mask_nc = xr.where(np.logical_or(mask_nc.lon<95,mask_nc.lon>106),100,mask_nc)

    #make attributes for lon and lat
    mask_nc.lat.attrs = {'units':'degrees_north','standard_name':'latitude'}
    mask_nc.lon.attrs = {'units':'degrees_east','standard_name':'longitude'}
    
    #if rainfall and mask not same resolution will perform a regrid
    if ~np.array_equal(mask_nc.lon,GPM_nc.lon) or ~np.array_equal(mask_nc.lon,GPM_nc.lon):
        regridder = xe.Regridder(GPM_nc,mask_nc,"conservative")
        GPM_nc = regridder(GPM_nc)
    
    #add time to mask to make masking easier
    time_ind = [i for i in GPM_nc.dims].index('time')
    mask_nc = mask_nc.expand_dims({'time':np.shape(GPM_nc)[time_ind]}, 0)
    mask_nc = mask_nc.assign_coords({'time':GPM_nc.time})

    if regrid:
        #regrid mask and rainfall
        ds_out = xr.Dataset(
            {
                "newgrid": (["lat","lon"],np.zeros([31,61])),
                "lat": (["lat"], np.arange(-10, 21, res), {"units": "degrees_north","standard_name":"latitude"}),
                "lon": (["lon"], np.arange(90, 151, res), {"units": "degrees_east","standard_name":"longitude"}),
            }
        )
        mask_regridder = xe.Regridder(mask_nc, ds_out, "conservative")
        mask_nc = mask_regridder(mask_nc)
        gpm_regridder = xe.Regridder(GPM_nc, ds_out, "conservative")
        GPM_nc = gpm_regridder(GPM_nc)
        
    dims=GPM_nc.dims
    mask_nc = mask_nc.transpose(*dims)
    GPM_masked = xr.where(mask_nc<80,GPM_nc,np.nan)
    if region != 'PM' and region!='WI':
        GPM_masked = xr.where(np.logical_or(GPM_masked.lon<box_points[0],GPM_masked.lon>box_points[1]),np.nan,GPM_masked)
        GPM_masked = xr.where(np.logical_or(GPM_masked.lat<box_points[2],GPM_masked.lat>box_points[3]),np.nan,GPM_masked).transpose(*dims)
    
    return GPM_masked

def prob_construct_waves(region='EM',seas='DJF',waves=['Kelv','R1'],pc=10):
    if len(waves)>2:
        raise Exception("Can only be calculated for 1 or 2 waves")
        
    l = {'EM':114,
        'PM':102,
        'EI':121,
        'WI':102,
        'NI':114,
        'SI':112,
        'NP':122,
        'SP':124,
        'SV':107,
        'CV':107,
        'NV':107,
        'ING':137,
    }[region]
    lon = [int(l-2.5),int(l+2.5)]

    wave_nc = xr.open_dataset('./inputs/wave_phases_era5_2020-2023.nc')
    hiw_ind = xr.open_dataset(f'./inputs/GPM_hiw_pc{pc:02d}_{region}_{seas}.nc').sel(time=slice('2001','2014'))

    nbins = 2 if len(waves)==2 else 3

    wave_df = pd.DataFrame({'time':wave_nc.time.values})
    for w,wv in zip(['Kelv','R1','WMRG'],wave_nc):
        wave_df[f'{w}_phase'] = wave_nc.sel(wave=w,lon=l).wave_phase.values
        wave_df[f'{w}_amp'] = wave_nc.sel(wave=w,lon=l).wave_amplitude.values
        
    for w in ['Kelv','R1','WMRG']:
        wave_df[f'{w}_amp_dis'] = pd.cut(wave_df[f'{w}_amp'], list(np.arange(nbins))+[ float('inf')],labels=list(np.arange(nbins).astype('float')))
        wave_df[f'{w}_comb'] = wave_df[f'{w}_phase'].astype(int).astype(str) + wave_df[f'{w}_amp_dis'].astype(int).astype(str)
        wave_df[f'{w}_comb'] = wave_df[f'{w}_comb'].where(wave_df[f'{w}_amp_dis'].astype(int)!=0,'0')
        
    s_ind = {'DJF':[12,1,2], 'MAM':[3,4,5], 'JJA':[6,7,8], 'SON':[9,10,11],'OND':[10,11,12]}[seas]
    wave_seas = wave_df.loc[wave_df.time.apply(lambda x: x.month in s_ind and 2001<=x.year<=2014)]

    hiw_ind = hiw_ind.pr_count>hiw_ind.attrs['area_thres']
    
    if len(waves)==2:
        df = pd.DataFrame({waves[0]:wave_seas.loc[:,f'{waves[0]}_comb'],waves[1]:wave_seas.loc[:,f'{waves[1]}_comb'],'hiw':hiw_ind.values})
        df['comb_phase'] = df[waves[0]].str[0]+df[waves[1]].str[0]
        grouped_df = df.groupby(['hiw',waves[0],waves[1]])['comb_phase'].count().unstack()
        df_hiw = grouped_df.iloc[grouped_df.index.get_level_values('hiw') == True]
        df_hiw = df_hiw.replace(np.nan,0.)
        df_hiw.index = df_hiw.index.get_level_values(waves[0])
        df_hiw = df_hiw/df_hiw.sum().sum()
    
        df_nohiw = grouped_df.iloc[grouped_df.index.get_level_values('hiw') == False]
        df_nohiw = df_nohiw.replace(np.nan,0.)
        df_nohiw.index = df_nohiw.index.get_level_values(waves[0])
        df_nohiw = df_nohiw/df_nohiw.sum().sum()
    else:
        df = pd.DataFrame({waves[0]:wave_seas.loc[:,f'{waves[0]}_comb'],'hiw':hiw_ind.values})
        grouped_df = df.groupby(['hiw',waves[0]])[waves[0]].count().unstack()
        df_hiw = grouped_df.loc[True]
        df_hiw = df_hiw.replace(np.nan,0.)
        df_hiw = df_hiw/df_hiw.sum().sum()
    
        df_nohiw = grouped_df.loc[False]
        df_nohiw = df_nohiw.replace(np.nan,0.)
        df_nohiw = df_nohiw/df_nohiw.sum().sum()
        
    pi1 = df['hiw'].mean()
    pi0 = 1-pi1
    pi = df_hiw*pi1/(df_hiw*pi1 + df_nohiw*pi0)
    pi = pi.to_xarray() if len(waves)==1 else pi.unstack().to_xarray()
    
    if len(waves)==1:
        pi.to_netcdf(f'./inputs/pi_{region}_{seas}_{waves[0]}_pc{pc:02d}.nc')
    else:
        pi.to_netcdf(f'./inputs/pi_{region}_{seas}_{waves[0]}_{waves[1]}_pc{pc:02d}.nc')
    
    return pi
    
#=======================================
#HybridWave class and supplemental functions
#=======================================

def retrieve_func_dual(pi,i,j):
    if pd.isnull(i) or pd.isnull(j):
        return np.nan
    else:
        return pi.loc[i][j]

def retrieve_func_single(pi,i):
    if pd.isnull(i):
        return np.nan
    else:
        return pi.loc[i]

def create_phase(phase,amp,max_amp=2):
    amp_dis = xr.where(np.isnan(amp),np.nan,np.floor(amp))
    amp_dis = xr.where(amp_dis>max_amp,max_amp,amp_dis)

    comb = (xr.where(np.isnan(phase),np.nan,phase.astype(int).astype(str)) + 
        xr.where(np.isnan(amp_dis),np.nan,amp_dis.astype(int).astype(str)))
    comb = xr.where(amp_dis==0,'0',comb)
    #comb = xr.where(np.isnan(amp_dis),'nan',comb)
    return comb
    
class HybridWave:
    '''
    Class for calculating GS6 wave-based hybrid forecast
    use eg:
    wave_obj = HybridWave(wave_lon=[100,105])
    for w in ['Kelv','R1']:
        wave_obj.load_wave(wave=w)
        wave_obj.calc_phase(wave=w)
    wave_obj.calc_prob(waves=['Kelv','R1'],,pc=10,region='PM',seas='DJF')
    
    '''
    def __init__(self,wave_lon):
        self.wave_ncs = {}
        self.wave_lon = wave_lon
        self.wave_phases = {}
        
    def load_wave(self,wave='Kelv'):
        ncs = xr.open_mfdataset(f"/gws/nopw/j04/ncas_climate_vol1/users/bc917929/GLOSEA6_wave/daily_filt/*/*/*{wave}*.nc")
        self.wave_ncs[wave] = ncs
        
    def calc_phase(self,wave='Kelv'):
        '''
        Calculate phases at given longitude range for GS6 wave forecasts and save
        '''
        if wave == 'Kelv':
            x = self.wave_ncs[wave].eastward_wind
            y = self.wave_ncs[wave].divergence
            print("Calculating wave phases")
            with ProgressBar():
                x_ind = x.sel(latitude=0,longitude=slice(*self.wave_lon)).mean('longitude').compute()                
                y_ind = y.sel(latitude=0,longitude=slice(*self.wave_lon)).mean('longitude').compute()
        elif wave == 'R1':
            x = self.wave_ncs[wave].eastward_wind
            y = self.wave_ncs[wave].northward_wind
            print("Calculating wave phases")
            with ProgressBar():
                x_ind = -1*x.sel(latitude=0,longitude=slice(*self.wave_lon)).mean('longitude').compute()
                y_ind = y.sel(latitude=8,longitude=slice(*self.wave_lon)).mean('longitude').compute()
        else:
            x = self.wave_ncs[wave].eastward_wind
            y = self.wave_ncs[wave].northward_wind
            print("Calculating wave phases")
            with ProgressBar():
                x_ind = -1*x.sel(latitude=-10,longitude=slice(*self.wave_lon)).mean('longitude').compute()
                y_ind = y.sel(latitude=0,longitude=slice(*self.wave_lon)).mean('longitude').compute()
        print(x_ind)
        print(y_ind)
        x_ind = x_ind/x_ind.isel(lead=0).std(['time','number'])
        y_ind = y_ind/y_ind.isel(lead=0).std(['time','number'])

        xr.merge(wave_phase(x_ind,y_ind,wave)).to_netcdf(f"/gws/nopw/j04/ncas_climate_vol1/users/bc917929/GLOSEA6_wave/phase/phase_{wave}_{self.wave_lon[0]}-{self.wave_lon[1]}.nc")

    def calc_prob(self,waves=['Kelv','R1'],pc=10,region='PM',seas='DJF'):
        
        hiw_ind = gpm_hiw(region,seas,pc).sel(time = slice('2001','2014'))
        pis = xr.open_dataset(f'./inputs/pi_{region}_{seas}_{waves[0]}_{waves[1]}_pc{pc:02d}.nc')
        pi1=prob_construct_single(hiw_ind,region=region,seas=seas,wave=waves[0],pc=pc)
        pi2=prob_construct_single(hiw_ind,region=region,seas=seas,wave=waves[1],pc=pc)
        pi=prob_construct_dual(hiw_ind,region=region,seas=seas,waves=waves,pc=pc)

        wave1 = xr.open_dataset(f"/gws/nopw/j04/ncas_climate_vol1/users/bc917929/GLOSEA6_wave/phase/phase_{waves[0]}_{self.wave_lon[0]}-{self.wave_lon[1]}.nc")
        wave2 = xr.open_dataset(f"/gws/nopw/j04/ncas_climate_vol1/users/bc917929/GLOSEA6_wave/phase/phase_{waves[1]}_{self.wave_lon[0]}-{self.wave_lon[1]}.nc")

        self.wave1_phase = create_phase(wave1.wave_phase,wave1.wave_amplitude,max_amp=2)
        self.wave2_phase = create_phase(wave2.wave_phase,wave2.wave_amplitude,max_amp=2)

        self.phase1 = create_phase(wave1.wave_phase,wave1.wave_amplitude,max_amp=1)
        self.phase2 = create_phase(wave2.wave_phase,wave2.wave_amplitude,max_amp=1)        

        retrieve_func = np.vectorize(retrieve_func)
        self.phase1 = xr.where(self.phase1=='nan',np.nan,self.phase1)
        self.phase2 = xr.where(self.phase2=='nan',np.nan,self.phase2)
        condprob_dual = self.phase1.copy(data=retrieve_func_dual(pi,self.phase1,self.phase2))
        condprob_dual.name = f'{waves[0]}_{waves[1]}_cond'
        
        retrieve_pi1 = np.vectorize(lambda i: retrieve_func(i,pi1))
        condprob_wave1 = self.wave1_phase.copy(data=retrieve_pi1(self.wave1_phase))
        condprob_wave1.name = f'{waves[0]}_cond'
        retrieve_pi2 = np.vectorize(lambda i: retrieve_func(i,pi2))
        condprob_wave2 = self.wave2_phase.copy(data=retrieve_pi2(self.wave2_phase))
        condprob_wave2.name = f'{waves[1]}_cond'

        self.condprobs = xr.merge([condprob_dual,condprob_wave1,condprob_wave2],compat='override')
        self.condprobs.to_netcdf(f"/gws/nopw/j04/ncas_climate_vol1/users/bc917929/GLOSEA6_wave/condprob/condprob_{waves[0]}_{waves[1]}_{self.wave_lon[0]}-{self.wave_lon[1]}_pc{pc}_{region}_{seas}.nc")

if __name__ == '__main__':
    #calculate pis and save
    #see Ferrett et al 2023 for further details on methodology
    for w in [['Kelv'],['R1'],['WMRG'],['Kelv','R1'],['Kelv','WMRG'],['R1','WMRG']]:
        print(w)
        prob_construct_waves(region='PM',seas='DJF',waves=w,pc=10)
    
    