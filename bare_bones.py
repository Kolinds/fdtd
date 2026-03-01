import numpy as np
import numba 
import h5py

CONS_SIZE = 200

ez = np.zeros(CONS_SIZE)
hy = np.zeros(CONS_SIZE)
imp0 = 337.0
maxTime = 250

with h5py.File("wave_data.hdf5", "w") as f:
    hy_data = f.create_dataset("mag_fdata", (500, 250))
    ez_data = f.create_dataset("elec_fdata", (500, 250))

    # Time loop
    for qTime in range (0, maxTime):
        # Update all magnetic field values half a time step
        for m in range (0, CONS_SIZE - 2):
            hy[m] = hy[m] + (1/imp0) * (ez[m + 1] - ez[m])

        # Update all electric field values half a time step
        for m in range (1, CONS_SIZE - 1):
            ez[m] = ez[m] + imp0 * (hy[m] - hy[m - 1])
        
        
        ez_data[qTime, 0:CONS_SIZE] = ez[:]
        hy_data[qTime, 0:CONS_SIZE] = hy[:]

        ez[0] = np.exp((-(qTime - 30)**2)/100.)







