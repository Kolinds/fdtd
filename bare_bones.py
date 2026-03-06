import numpy as np
from numba import njit
import h5py


def hpf5_maxValue(file_name, dataset_name):
    DIM_JUMP = 100
    total_max = 0
    with h5py.File(file_name, "r") as f:
        dset = f[dataset_name]
        for actual_line in range(0, dset.shape[0], DIM_JUMP):
            jump_matrix = dset[actual_line:actual_line+DIM_JUMP, ...]
            actual_max = np.max(np.abs(jump_matrix))
            if (total_max < actual_max): total_max = actual_max
    return total_max


def hpf5_normalization(file_name, dataset_name, maxValue):
    DIM_JUMP = 100
    with h5py.File(file_name, "r+") as f:
        dset = f[dataset_name]
        for actual_line in range(0, dset.shape[0], DIM_JUMP):
            jump_matrix = dset[actual_line:actual_line+DIM_JUMP, ...]
            norm_matrix = jump_matrix / maxValue
            dset[actual_line:actual_line+DIM_JUMP, ...] = norm_matrix

@njit
def update_magnetic_field(hy, ez, chyh, chye, SPACE_SIZE):
    for m in range(0, SPACE_SIZE - 1):
        hy[m] = chyh[m] * hy[m] + chye[m] * (ez[m + 1] - ez[m])  
    return hy

@njit
def update_electric_field(ez, hy, ceze, cezh, SPACE_SIZE):
    for m in range(1, SPACE_SIZE - 1):
        ez[m] = ceze[m] * ez[m] + cezh[m] * (hy[m] - hy[m - 1])
    return ez



SPACE_SIZE = 240
TOTAL_TIME = 1000
TIME_BUFFER = 100
IMP0 = 337.0
LOSS_LAYER = 180
LOSS = 0.02

#Arrays for saving electric and E and H field values, and buffers
ez_array= np.zeros(SPACE_SIZE)
hy_array= np.zeros(SPACE_SIZE - 1)
ez_buffer= np.zeros((TIME_BUFFER, SPACE_SIZE))
hy_buffer= np.zeros((TIME_BUFFER, SPACE_SIZE - 1))

#Arrays for permitivity, permeability and loss of materials 
ceze_array = np.ones(SPACE_SIZE)
cezh_array = np.ones(SPACE_SIZE)
chyh_array = np.ones(SPACE_SIZE - 1)
chye_array = np.ones(SPACE_SIZE - 1)

for m in range(0, SPACE_SIZE - 1):
    if (m < LOSS_LAYER): #Free space
        chyh_array[m] = 1.0
        chye_array[m] = 1 / IMP0
    else:
        chyh_array[m] = (1.0 - LOSS) / (1.0 + LOSS)
        chye_array[m] = (1.0 / IMP0) / (1.0 + LOSS)

for m in range(0, SPACE_SIZE):
    if (m < 100): #Free space
        ceze_array[m] = 1.0
        cezh_array[m] = IMP0
    elif (m < LOSS_LAYER):   #Dielectric
        ceze_array[m] = 1.0
        cezh_array[m] = IMP0 / 9.0
    else: #Lossy Boundary Layer
        ceze_array[m] = (1.0 - LOSS) / (1.0 + LOSS)
        cezh_array[m] = (IMP0 / 9.0) / (1.0 + LOSS)



with h5py.File("wave_data.hdf5", "w") as f:
    hy_dset = f.create_dataset("mag_fdata", (TOTAL_TIME, SPACE_SIZE - 1))
    ez_dset = f.create_dataset("elec_fdata", (TOTAL_TIME, SPACE_SIZE))

    # Time loop
    for qTime in range (0, TOTAL_TIME):

        """
        Absorbing boundary condition (ABC) mag, deprecated
        hy_array[SPACE_SIZE-1] = hy_array[SPACE_SIZE-2]
        """
        hy_array = update_magnetic_field(hy_array, ez_array, chyh_array, chye_array, SPACE_SIZE - 1)
        #Correcion TFSF para quitar el efecto del incidente
        hy_array[49] -= np.exp((-(qTime - 30)**2)/100.) / IMP0

        #Absorbing boundary condition (ABC) elec
        ez_array[0] = ez_array[1]
        """
        This ABC does not work inside dielectric
        ez_array[SPACE_SIZE-1] = ez_array[SPACE_SIZE-2]
        """
        ez_array = update_electric_field(ez_array, hy_array, ceze_array, cezh_array, SPACE_SIZE)
        #Correcion TFSF para añadir el efecto del incidente
        ez_array[50] += np.exp((-(qTime + 0.5 - (-0.5) - 30)**2)/100.)


        #Buffer stuff
        buffer_index = qTime % TIME_BUFFER

        if (buffer_index == 0)  and (qTime != 0):
            ez_dset[qTime - TIME_BUFFER:qTime, 0:SPACE_SIZE] = ez_buffer
            hy_dset[qTime - TIME_BUFFER:qTime, 0:SPACE_SIZE - 1] = hy_buffer

        ez_buffer[buffer_index, :] = ez_array
        hy_buffer[buffer_index, :] = hy_array
        
        """ 
        Only equal: = -> Hardwired source
        Summed: += -> Source 
        ez_array[50] += np.exp((-(30 - 30)**2)/100)
        """
        

        # ABC's in 1D with Courant Limit 




#Normalization of the values

max_mag = hpf5_maxValue("wave_data.hdf5", "mag_fdata")
hpf5_normalization("wave_data.hdf5", "mag_fdata", max_mag)
max_elec = hpf5_maxValue("wave_data.hdf5", "elec_fdata")
hpf5_normalization("wave_data.hdf5", "elec_fdata", max_elec)



# holi uwu

