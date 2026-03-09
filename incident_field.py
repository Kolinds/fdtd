import numpy as np
from numba import njit

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
