import numpy as np
import pandas as pd
import pyspedas.tplot_tools as pytplot
import pyspedas
import os

def load_mag(days_of_interest, spdf = True):
    """
    Loads MAG data from MAVEN using pyspedas.

    Parameters:
    day_of_interest (list): List containing dates of interest, in the form
        "YYYY-MM-DD"

    Returns:
    mag: A dictionary containing the following keys:
        - times: An array of size n with the time (np.datetime64) of each point
                in the time series at the precision of MAG.
        - B: An (n,3) array containing the vector components (Bx, By, Bz) of the 
                magnetic field B at each point in the time series.
        - posn: An (n,3) array containing the vector components (x, y, z) of 
                MAVEN's position (in MSO coordinates) at each point in the time
                series.
    """
    for day in days_of_interest:
        pyspedas.projects.maven.mag(trange=[day, day], spdf = True) #introduce spdf=True for backwaards compatibliity? 
    mag = pyspedas.get_data('OB_B')
    # times = []
    # for nanos in mag.times:   # if in UNIX nanosecond format
    #     ts = pd.Timestamp(nanos, tz = 'US/Eastern')
    #     times.append(pd.Timestamp.to_datetime64(ts))
    times = pd.to_datetime(mag.times*1e9, unit = 'ns') #from epoch in seconds to 'datetime64[ns]'
    Bvals = mag.y[:, 0:3]
    posn = pyspedas.get_data('POSN').y
    posn = posn/3389.5 # normalized to Mars radii
    return {'times': times, 'B': Bvals, 'posn': posn}

def load_mag_sts(directory, filenames):
    """
    Loads MAG data from MAVEN with pre-existing .sts files.

    Parameters:
    directory (str): Path to directory containing .sts files.
    filenames (list): List of .sts filenames to load.

    Returns:
    mag: A dictionary containing the following keys:
        - times: An array of size n with the time (np.datetime64) of each point
                in the time series at the precision of MAG.
        - B: An (n,3) array containing the vector components (Bx, By, Bz) of the 
                magnetic field B at each point in the time series.
        - posn: An (n,3) array containing the vector components (x, y, z) of 
                MAVEN's position (in MSO coordinates) at each point in the time
                series.
    """
    for file in filenames:
        f = os.path.join(directory, file)
        pytplot.sts_to_tplot(f, prefix='',suffix='', merge=True)
    mag = pyspedas.get_data('OB_B')
    times = []
    for nanos in mag.times:   # if in UNIX nanosecond format
        ts = pd.Timestamp(nanos, tz = 'US/Eastern')
        times.append(pd.Timestamp.to_datetime64(ts))
    times = np.array(times)
    Bvals = mag.values[:, 0:3]
    posn = pyspedas.get_data('POSN').values
    posn = posn/3389.5 # normalized to Mars radii
    return {'times': times, 'B': Bvals, 'posn': posn}

def load_swea(days_of_interest):
    """
    Loads SWEA flux data from MAVEN using pyspedas.

    Parameters:
    day_of_interest (list): List containing dates of interest, in the form
        "YYYY-MM-DD"

    Returns:
    swea: A dictionary containing the following keys:
        - times: An array of size n with the time (np.datetime64) of each point
                in the time series at the precision of SWEA.
        - flux: An (n,64) array containing the eflux (eV/(eV cm^2 sr s)) of 
                electrons measured by SWEA at each of its 64 energy bands, at
                each point in the time series.
        - v: An array of size 64 containing the energy values (eV) of each of
                SWEA's measured energy bands.
    """
    for day in days_of_interest:
        pyspedas.projects.maven.swea(trange=[day, day])
    swea = pyspedas.get_data('diff_en_fluxes_svyspec')
    swea_times = swea.time.values
    swea_flux = swea.values
    swea_v = swea.v.values
    return {'times': swea_times, 'flux': swea_flux, 'v': swea_v}

def load_swia_mom(days_of_interest):
    """
    Loads SWIA ion moments from MAVEN using pyspedas.

    Parameters:
    day_of_interest (list): List containing dates of interest, in the form
        "YYYY-MM-DD"

    Returns:
    swia_mom: A dictionary containing the following keys:
        - times: An array of size n with the time (np.datetime64) of each point
                in the time series at the precision of SWIA.
        - temp: An (n,3) array containing the vector components (Tx, Ty, Tz) of
                the ion temperature at each point in the time series.
        - vel: An (n,3) array containing the vector components (v_x, v_y, v_z) of
                the ion velocity at each point in the time series.
        - density: An array of size n containing the ion density measured at 
                each point in the time series.
    """
    for day in days_of_interest:
        pyspedas.projects.maven.swia(trange=[day, day], datatype = 'onboardsvymom')
    temperature = pyspedas.get_data('temperature_mso_onboardsvymom')
    swia_time = temperature.time.values
    temp = temperature.values
    velocity = pyspedas.get_data('velocity_mso_onboardsvymom').y
    density = pyspedas.get_data('density_onboardsvymom').y
    return {'times': swia_time, 'temp': temp, 'vel': velocity, 'density': density}

def load_swia_flux(days_of_interest):
    """
    Loads SWIA flux data from MAVEN using pyspedas.

    Parameters:
    day_of_interest (list): List containing dates of interest, in the form
        "YYYY-MM-DD"

    Returns:
    swia_flux: A dictionary containing the following keys:
        - times: An array of size n with the time (np.datetime64) of each point
                in the time series at the precision of SWIA.
        - flux: An (n,48) array containing the eflux (eV/(eV cm^2 sr s)) of 
                ions (assumed to be H+ protons) measured by SWIA at each of its
                64 energy bands, at each point in the time series.
        - v: An array of size 48 containing the energy values (eV) of each of
                SWIA's measured energy bands.
    """
    for day in days_of_interest:
        pyspedas.projects.maven.swia(trange=[day, day], datatype = 'onboardsvyspec')
    swia = pyspedas.get_data('spectra_diff_en_fluxes_onboardsvyspec')
    swia_times = swia.times
    swia_flux = swia.y
    swia_v = swia.v
    return {'times': swia_times, 'flux': swia_flux, 'v': swia_v}

def load_static(days_of_interest):
    """
    Loads STATIC flux data from MAVEN using pyspedas.

    Parameters:
    day_of_interest (list): List containing dates of interest, in the form
        "YYYY-MM-DD"

    Returns:
    static: A dictionary containing the following keys:
        - times: An array of size n with the time (np.datetime64) of each point
                in the time series at the precision of STATIC.
        - flux: An (n,48) array containing the eflux (eV/(eV cm^2 sr s)) of 
                ions (assumed to be H+ protons) measured by SWIA at each of its
                64 energy bands, at each point in the time series.
        - v: An array of size 48 containing the energy values (eV) of each of
                SWIA's measured energy bands.
    """
    for day in days_of_interest:
        pyspedas.projects.maven.sta(trange=[day, day], level='l2',datatype='c6', 
                           get_support_data = True)
    eflux = pyspedas.get_data('eflux_c6-32e64m')
    static_times = pd.to_datetime(eflux.times*1e9, unit = 'ns') #from epoch in seconds to 'datetime64[ns]'
    swp_ind = pyspedas.get_data('swp_ind_c6-32e64m')
    energy = pyspedas.get_data('energy_c6-32e64m')
    mass = pyspedas.get_data('mass_arr_c6-32e64m')

    mas = []
    for i in range(len(swp_ind.y)):
        sweep = swp_ind.y[i]
        mas.append(mass[:,:,sweep])
    mas = np.asarray(mas)

    return {'times': static_times, 'flux': eflux.y,
            'sweep_index': swp_ind.y, 'energy': energy, 'mass': mas}