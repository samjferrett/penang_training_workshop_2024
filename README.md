# <p align="center">WCSSP: SE Asia annual meeting 2024 Training Workshop slides - Hybrid dynamical-statistical forecasting</p>

Files needed to replicate figures presented in WCSSP: SE Asia annual meeting 2024 Training Workshop "Hybrid dynamical-statistical forecasting" plus some additional methodology functions.

## Files

### case_example.ipynb
This is the main file used to produce figures for the training session slides. Data to run is available...
The event is from the 7th December 2022 in Peninsular Malaysia.
The script produces the following figures in order:
1. GPM-IMERG rainfall anomaly contours and ERA5 equatorial wave data for 3 days (6th December 2022 - 8th December 2022)
2. ERA5 Wave hovmollers for Kelvin, R1 and WMRG waves for December 2022
3. Heavy rainfall metric for the event period
4. Wave-based hybrid forecasts
5. GloSea6 direct rainfall forecasts
6. Kelvin+WMRG wave-based forecasts compared to GPM for 2 lead times
7. Atmospheric weather pattern/regime based forecast

Before plots 4. there is a function defined called build_hybrid() that shows how to calculate the hybrid forecast given wave types to use, a region, a season and a heavy rainfall event frequency. This iteration of the method requires pre-calculated datasets and so does not work for other inputs than region='PM', season='DJF' and pc=10 currently, though these required datasets can be calculated using supplied additional functions.

### hybrid_functions.py
A set of functions that can be used to calculate conditional probabilities required for the hybrid method. Plus a class that demonstrates applying the hybrid forecast methodology to preexisting ensemble forecast data. Functions require GPM and GPM landsea mask also which can be retrieved at https://gpm.nasa.gov/data/directory. ERA5 wave phase files are available on request (s.j.ferrett@reading.ac.uk).

Functions included:
1. find_areathres() - Calculates the local 95th percentile rainfall threshold and the area threshold required for pc,
    thresholds based on data up to 2001-2014.
2. mask_gpm() - Mask rainfall data except the provided region. Regridding optional.
3. prob_construct_waves() - Calculate conditional probability of heavy rainfall for 1 or 2 waves, given region, season and heavy rainfall event frequency.
4. A class and supplemental functions demonstrating an application of the hybrid construction to ensemble wave forecast data, including wave phase calculation given filtered wave data.
   
### wave_functions.py
Supplemental wave methodology function
1. phase_func() - Calculating phase of waves given x and y as detailed in Yang et al. 2021 (https://doi.org/10.1175/WAF-D-20-0144.1)
### calc_hiw_data.py
   -
### calc_wave_event.py
   -
### calc_wave_phases.ipynb
   -

## Inputs required

## Package requirements
