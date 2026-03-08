import hdf5_handler as h5h
import grid_struct as gr
import maxwell_update as upd
import numpy as np

FILE_NAME="wave_data.hdf5"
EDSET_NAME = "elec_fdata"
HDSET_NAME = "mag_fdata"
SPACE_SIZE = 240
TOTAL_TIME = 1000
TIME_BUFFER = 100
BUFFER_JUMP = 100
IMP0 = 337.0
LOSS_LAYER = 180
LOSS = 0.02


hdf5_file = h5h.HDF5Writer(FILE_NAME, TIME_BUFFER, TOTAL_TIME, SPACE_SIZE)
grid = gr.Grid(SPACE_SIZE)


hdf5_file.open_file(EDSET_NAME, HDSET_NAME)
grid.place_materials(SPACE_SIZE, LOSS, LOSS_LAYER, IMP0)
for qTime in range (0, TOTAL_TIME):

    grid.hy = upd.update_magnetic_field(grid.hy, grid.ez, grid.chyh, grid.chye, SPACE_SIZE - 1)
    grid.abc_boundary()
    grid.TFSF_boundary(qTime, 30, 0.5, -0.5, 0, 0, 100, 1, IMP0)
    grid.ez = upd.update_electric_field(grid.ez, grid.hy, grid.ceze, grid.cezh, SPACE_SIZE)
    hdf5_file.update_file(qTime, grid.ez, grid.hy)


h5h.normalization(FILE_NAME, HDSET_NAME, BUFFER_JUMP, h5h.maxValue(FILE_NAME, HDSET_NAME, 100)) #Normalization of the H-field
h5h.normalization(FILE_NAME, EDSET_NAME, BUFFER_JUMP, h5h.maxValue(FILE_NAME, EDSET_NAME, 100)) #Normalization of the E-field



# holi uwu

