# HDF5 File Configuration
FILE_NAME = "wave_data.hdf5"
EDSET_NAME = "elec_fdata"
HDSET_NAME = "mag_fdata"

# Buffer & Memory
TIME_BUFFER = 100
BUFFER_JUMP = 100

# Physical & Material Parameters
IMP0 = 337.0
DIELECTRIC_LAYER = 200
LOSS_LAYER = 300
LOSS = 0.02

# Simulation Grid & Time
SPACE_SIZE = 2500
TOTAL_TIME = 100000
COURANT = 1
TFSF_BOUNDARY = 75

#Incident fields and wavelets
STEPS_WAVELENGTH = 40
RAMP_PERIODS = 10

RICKER_DELAY = 1

