import numpy as np
from numba import njit
import h5py


def hpf5_maxValue(file_name, dataset_name):
    DIM_JUMP = 10
    total_max = 0
    with h5py.File(file_name, "r") as f:
        dset = f[dataset_name]
        for actual_line in range(0, dset.shape[0], DIM_JUMP):
            jump_matrix = dset[actual_line:actual_line+DIM_JUMP, ...]
            actual_max = np.max(jump_matrix)
            if (abs(total_max) < (actual_max)): total_max = actual_max
    return total_max

def hpf5_normalization(file_name, dataset_name, maxValue):
    DIM_JUMP = 10
    with h5py.File(file_name, "r+") as f:
        dset = f[dataset_name]
        for actual_line in range(0, dset.shape[0], DIM_JUMP):
            jump_matrix = dset[actual_line:actual_line+DIM_JUMP, ...]
            norm_matrix = jump_matrix / maxValue
            dset[actual_line:actual_line+DIM_JUMP, ...] = norm_matrix


def update_magnetic_field(hy, ez, imp0, SPACE_SIZE):
    for m in range(0, SPACE_SIZE - 1):
        hy[m] = hy[m] + (1.0 / imp0) * (ez[m + 1] - ez[m])
        
    return hy


def update_electric_field(ez, hy, imp0, SPACE_SIZE):
    for m in range(1, SPACE_SIZE):
        ez[m] = ez[m] + imp0 * (hy[m] - hy[m - 1])
        
    return ez



SPACE_SIZE = 200
TOTAL_TIME = 5000
TIME_BUFFER = 100

ez_array= np.zeros(SPACE_SIZE)
hy_array= np.zeros(SPACE_SIZE)
ez_buffer= np.zeros((TIME_BUFFER, SPACE_SIZE))
hy_buffer= np.zeros((TIME_BUFFER, SPACE_SIZE))

imp0 = 337.0

with h5py.File("wave_data.hdf5", "w") as f:
    hy_dset = f.create_dataset("mag_fdata", (TOTAL_TIME, SPACE_SIZE))
    ez_dset = f.create_dataset("elec_fdata", (TOTAL_TIME, SPACE_SIZE))

    # Time loop
    for qTime in range (0, TOTAL_TIME):
        hy_array = update_magnetic_field(hy_array, ez_array, imp0, SPACE_SIZE)
        ez_array = update_electric_field(ez_array, hy_array, imp0, SPACE_SIZE)
        
        #Buffer stuff
        buffer_index = qTime % TIME_BUFFER
        ez_buffer[buffer_index, :] = ez_array
        hy_buffer[buffer_index, :] = hy_array

        if (buffer_index == 0)  and (qTime != 0):
            ez_dset[qTime - TIME_BUFFER:qTime, 0:SPACE_SIZE] = ez_buffer
            hy_dset[qTime - TIME_BUFFER:qTime, 0:SPACE_SIZE] = hy_buffer

        ez_array[0] = np.exp((-(qTime - 30)**2)/100.)

"""
#Normalization of the values
max_mag = hpf5_maxValue("wave_data.hdf5", "mag_fdata")
hpf5_normalization("wave_data.hdf5", "mag_fdata", max_mag)
max_elec = hpf5_maxValue("wave_data.hdf5", "elec_fdata")
hpf5_normalization("wave_data.hdf5", "elec_fdata", max_elec)
"""

# holi uwu

