import numpy as np
import numba 
import h5py

SPACE_SIZE = 200
TOTAL_TIME = 500

ez = np.zeros(SPACE_SIZE)
hy = np.zeros(SPACE_SIZE)
imp0 = 337.0


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




with h5py.File("wave_data.hdf5", "w") as f:
    hy_data = f.create_dataset("mag_fdata", (500, SPACE_SIZE))
    ez_data = f.create_dataset("elec_fdata", (500, SPACE_SIZE))

    # Time loop
    for qTime in range (0, TOTAL_TIME):
        # Update all magnetic field values half a time step
        for m in range (0, SPACE_SIZE - 2):
            hy[m] = hy[m] + (1/imp0) * (ez[m + 1] - ez[m])

        # Update all electric field values half a time step
        for m in range (1, SPACE_SIZE - 1):
            ez[m] = ez[m] + imp0 * (hy[m] - hy[m - 1])
        
        
        ez_data[qTime, 0:SPACE_SIZE] = ez[:]
        hy_data[qTime, 0:SPACE_SIZE] = hy[:]

        ez[0] = np.exp((-(qTime - 30)**2)/100.)


#Normalization of the values
max_mag = hpf5_maxValue("wave_data.hdf5", "mag_fdata")
hpf5_normalization("wave_data.hdf5", "mag_fdata", max_mag)
max_elec = hpf5_maxValue("wave_data.hdf5", "elec_fdata")
hpf5_normalization("wave_data.hdf5", "elec_fdata", max_elec)




