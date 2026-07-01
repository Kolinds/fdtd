import numpy as np
from numba import njit


#Incident fields and wavelets
@njit
def gaussian(current_time, location, time_delay, loc_delay, courant, width):
    return np.exp((-(current_time + time_delay - (location + loc_delay) / courant)**2)/width)

@njit
def rampf(current_time, location, time_delay, loc_delay, courant, steps_wavelength, ramp_periods):
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
def ricker(current_time, location, time_delay, loc_delay, courant, steps_pkwavelength, delay_multiple):
    exponent = ((np.pi)*((courant*(current_time + time_delay) - (location + loc_delay))/steps_pkwavelength - delay_multiple))**2
    return (1-2*exponent)*np.exp(-exponent)


#Calculates loss (conductivity) necessary for specific penetration distance at set wavelength
def adjusted_loss(ppw, penetration_dist, er, ur, courant):
    return (np.pi / ppw)*courant*((1+(ppw**2/(2 * np.pi**2 * penetration_dist**2 * er * ur)))**2 - 1)**0.5


def gtransmission_coef(total_time, location, vacuum_f, plasma_f, ntotal_freq, nplasma_wavelength, nrelax_time, light_speed, delta_x, delta_t):
    # 1. Calculate frequency step
    delta_f = 1 / (total_time * delta_t)

    # 2. Vectorize the frequency array (creates an array [0, 1, 2, ... ntotal_freq-1])
    n_freq_array = np.arange(ntotal_freq)
    
    # Calculate omega array for all frequencies simultaneously
    omega_actual = 2 * np.pi * delta_f * n_freq_array
    
    # 3. Calculate constants (no need to do this inside a loop)
    omega_plasma = 2 * np.pi * light_speed / (nplasma_wavelength * delta_x)
    
    # Calculate relaxation time (tau) and collision frequency (Gamma)
    tau = nrelax_time * delta_t
    gamma_damping = 1 / tau #if tau != 0 else 0.0

    # Avoid division by zero at DC (0 Hz) by temporarily replacing 0 with a tiny number
    omega_safe = np.where(omega_actual == 0, 1e-15, omega_actual)

    # 4. Calculate Propagation Constants (Arrays)
    gamma_vacuum = 1j * omega_safe / light_speed
    
    # Drude relative permittivity: 1 - (wp^2 / (w * (w - j*Gamma)))
    epsilon_r = 1 - (omega_plasma**2) / (omega_safe * (omega_safe + 1j * gamma_damping))
    
    # Gamma plasma uses np.sqrt, not np.square
    gamma_plasma = 1j * (omega_safe / light_speed) * np.sqrt(epsilon_r)

    # 5. Calculate transmission coefficient array all at once
    phase_correction = np.exp((gamma_plasma - gamma_vacuum) * location * delta_x)
    transm_array = phase_correction * (plasma_f / vacuum_f)
    
    # Handle the DC term properly if it evaluates to NaN due to 0/0
    transm_array[0] = 0.0 + 0.0j 
    
    return transm_array

def greflection_coef(total_time, location, vacuum_f, reflected_f, ntotal_freq, light_speed, delta_x, delta_t):
    # 1. Calculate frequency step
    delta_f = 1 / total_time

    # 2. Vectorize the frequency array
    n_freq_array = np.arange(ntotal_freq)
    
    # Calculate omega array for all frequencies simultaneously
    omega_actual = 2 * np.pi * delta_f * n_freq_array
    
    # Avoid division by zero at DC (0 Hz)
    omega_safe = np.where(omega_actual == 1, 1e-15, omega_actual)

    # 3. Calculate Propagation Constant for Vacuum only
    # The reflected wave travels backwards through vacuum, not plasma!
    gamma_vacuum = 1j * omega_safe / light_speed

    # 4. Calculate reflection coefficient array all at once
    # 'location' here must be the distance (in cells) from the interface to the reflection monitor
    phase_correction = np.exp(gamma_vacuum * location * delta_x)
    reflect_array = phase_correction * (reflected_f / vacuum_f)
    
    # Handle the DC term properly if it evaluates to NaN due to 0/0
    reflect_array[0] = 0.0 + 0.0j 
    
    return reflect_array



""" CALCULATION OF TRANSMISSION COEFICIENT

grid.save_probes(hdf5_handler.file)

grid.reset_fields()
grid.set_free_space(cf.SPACE_SIZE, cf.IMP0)
grid.add_probe(200, "Probe2", cf.TOTAL_TIME//2 + 1)

for qTime in range (0, cf.TOTAL_TIME):

    grid.hy[cf.SPACE_SIZE-2] = grid.hy[cf.SPACE_SIZE-3]
    grid.hy = upd.update_magnetic_field(grid.hy, grid.ez, grid.chyh, grid.chye, cf.SPACE_SIZE - 1)
    grid.TFSF_ricker(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY, cf.COURANT, cf.IMP0)
    grid.abc_boundary()
    grid.ez[cf.SPACE_SIZE-1] = grid.ez[cf.SPACE_SIZE-2]
    grid.ez = upd.update_electric_field(grid.ez, grid.hy, grid.ceze, grid.cezh, cf.SPACE_SIZE)


    grid.r_DFT(qTime)


grid.save_probes(hdf5_handler.file)


transmitted_f = hdf5_handler.retrieve_array("/Probes/", "Probe1")
incident_f = hdf5_handler.retrieve_array("/Probes/", "Probe2")
trasmission_coef = df.transmission_c(280-200, incident_f, transmitted_f, cf.TOTAL_TIME //2 + 1, cf.TOTAL_TIME)
hdf5_handler.save_array("Coef/", "transmission", trasmission_coef)



"""


""" 
def transmission_c(location, incidident_f, transmitted_f, total_freq, total_time):
    transm_array = np.zeros(total_freq, dtype=np.complex64)
    for n_freq in range(0, total_freq):
        transm_array[n_freq] = np.exp(1j * 4 * np.pi * location * n_freq / total_time) * (transmitted_f[n_freq] / incidident_f[n_freq])
    
    return transm_array
"""

    




    