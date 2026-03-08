import hdf5_handler as h5h
import grid_struct as gr
import maxwell_update as upd


SPACE_SIZE = 240
TOTAL_TIME = 1000
TIME_BUFFER = 100
IMP0 = 337.0
LOSS_LAYER = 180
LOSS = 0.02


hdf5_file = h5h.HDF5Writer()
grid = gr.Grid()

for qTime in range (0, TOTAL_TIME):

    grid.hy = upd.update_magnetic_field(grid.hy, grid.ez, grid.chyh, grid.chye, SPACE_SIZE - 1)
    grid.abc_boundary()
    grid.TFSF_boundary(qTime, 30, 0.5, -0.5, 0, 0, 1001, IMP0)
    grid.ez = upd.update_electric_field(grid.ez, grid.hy, grid.ceze, grid.cezh, SPACE_SIZE)
     


h5h.normalization("wave_data.hdf5", "mag_fdata", h5h.maxValue("wave_data.hdf5", "mag_fdata")) #Normalization of the H-field
h5h.normalization("wave_data.hdf5", "elec_fdata", h5h.maxValue("wave_data.hdf5", "elec_fdata")) #Normalization of the E-field



# holi uwu

