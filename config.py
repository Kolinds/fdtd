# HDF5 File Configuration
FILE_NAME = "wave_data.hdf5"
EDSET_NAME = "elec_fdata"
HDSET_NAME = "mag_fdata"

# Buffer & Memory
TIME_BUFFER = 100
BUFFER_JUMP = 100

# Physical & Material Parameters
IMP0 = 337.0
DIELECTRIC_LAYER = 100
LOSS_LAYER = 150
LOSS = 0.02

# -> For plasma slab
E_CONDUCTIVITY = 0.0
H_CONDUCTIVITY = 0.0
RELAX_TIME_STEPS = 2000
PLASMA_WAVELENGTH_STEPS = 40
PERMITIVITY_INF = 1.0

# Simulation Grid & Time
TOTAL_TIME = 2000
COURANT = 0.99
TFSF_BOUNDARY = 75

DELTA_T = 1.0

#Incident fields and wavelets
STEPS_WAVELENGTH = 30
RAMP_PERIODS = 10

RICKER_DELAY = 1

