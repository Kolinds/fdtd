import numpy as np
import numba 
import h5py

CONS_SIZE = 200

ez = np.random.rand(CONS_SIZE)
hy = np.random.rand(CONS_SIZE)
imp0 = 337.0
maxTime = 250

f = h5py.File("wave_data.hdf5", "w")
dset = f.create_dataset("fields_values", (1000, 200))
dset[0, :] = ez
dset[1, :] = hy
dset[2, :] = hy

print(dset[0, 3])