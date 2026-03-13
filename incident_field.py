import numpy as np
from numba import njit


#Incident fields and wavelets
@njit
def ezfield(current_time, location, time_delay, loc_delay, width, courant):
    return np.exp((-(current_time + time_delay - (location + loc_delay) / courant)**2)/width)

@njit
def rampf(current_time, location, time_delay, loc_delay, steps_wavelength, ramp_periods, courant):
    ramp_duration = (steps_wavelength * ramp_periods)/(2*courant)
    phase_ramp = (1/ramp_duration)*np.pi*(current_time + time_delay)
    phase_cos = (1/ramp_duration)*np.pi*ramp_periods*((current_time + time_delay) - (location + loc_delay))
    if current_time < 0:
        return 0
    elif(current_time < ramp_duration):
        return 0.5 * (1 - np.cos(phase_ramp))*np.cos(phase_cos)
    else:
        return np.cos(phase_cos)

@njit
def ricker(current_time, location, time_delay, loc_delay, steps_pkwavelength, delay_multiple, courant):
    exponent = ((np.pi)*((courant*(current_time + time_delay) - (location + loc_delay))/steps_pkwavelength - delay_multiple))**2
    return (1-2*exponent)*np.exp(-exponent)


#Calculates loss (conductivity) necessary for specific penetration distance at set wavelength
def adjusted_loss(ppw, penetration_dist, er, ur, courant):
    return (np.pi / ppw)*courant*((1+(ppw**2/(2 * np.pi**2 * penetration_dist**2 * er * ur)))**2 - 1)**0.5



"""CODIGO EN FASE DE ELABORACION
def transmission_c(location, incidident_f, transmitted_f, total_freq, total_time):
    transm_array = np.zeros(total_freq, dtype=np.complex64)
    for n_freq in range(0, total_freq):
        transm_array[n_freq] = np.exp(1j * 4 * np.pi * location * n_freq / total_time) * (transmitted_f[n_freq] / incidident_f[n_freq])
    
    return transm_array


incident_f = hdf5_handler.retrieve_array("/Probes/", "Probe1")
transmitted_f = hdf5_handler.retrieve_array("/Probes/", "Probe2")
trasmission_coef = df.transmission_c(275, incident_f, transmitted_f, cf.TOTAL_TIME //2 + 1, cf.TOTAL_TIME)
hdf5_handler.save_array("Coef/", "transmission", trasmission_coef)
"""

    




    