import xarray as xr
import numpy as np
import pandas as pd

def phase_func(x,y,kind='Kelv'):
    phase_dis = xr.ones_like(x)

    #amplitude
    amplitude=np.sqrt(x**2+y**2)

    #continuous phase
    phase = np.arctan2(y,x)
    phase = np.where(np.less(phase,0),phase+2*np.pi,phase)

    #discreet phase
    phase_dis = np.where(np.greater(x,np.abs(y)),4,phase_dis) if kind=='Kelv' else np.where(np.greater(x,np.abs(y)),2,phase_dis)
    phase_dis = np.where(np.less(y,-1*np.abs(x)),3,phase_dis)
    phase_dis = np.where(np.less(x,-1*np.abs(y)),2,phase_dis) if kind=='Kelv' else np.where(np.less(x,-1*np.abs(y)),4,phase_dis)

    #handling nan values
    phase = xr.where(np.isnan(amplitude),np.nan,phase)
    phase_dis = xr.where(np.isnan(amplitude),np.nan,phase_dis)

    phase.name='phi'
    phase_dis.name='wave_phase'
    amplitude.name = 'wave_amplitude'
    
    return phase,phase_dis,amplitude
