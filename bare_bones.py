import hdf5_handler as h5h
import grid_struct as gr
import maxwell_update as upd
import config as cf


hdf5_file = h5h.HDF5Writer(cf.FILE_NAME, cf.TIME_BUFFER, cf.TOTAL_TIME, cf.SPACE_SIZE)
grid = gr.Grid(cf.SPACE_SIZE)

hdf5_file.open_file(cf.EDSET_NAME, cf.HDSET_NAME)
grid.place_materials(cf.SPACE_SIZE, cf.DIELECTRIC_LAYER, cf.LOSS, cf.LOSS_LAYER, cf.IMP0)

for qTime in range (0, cf.TOTAL_TIME):

    grid.hy = upd.update_magnetic_field(grid.hy, grid.ez, grid.chyh, grid.chye, cf.SPACE_SIZE - 1)
    grid.TFSF_rampcos(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, cf.STEPS_WAVELENGTH, cf.RAMP_PERIODS, cf.COURANT, cf.IMP0)
    grid.abc_boundary()
    grid.ez = upd.update_electric_field(grid.ez, grid.hy, grid.ceze, grid.cezh, cf.SPACE_SIZE)
    
    hdf5_file.update_file(qTime, grid.ez, grid.hy)

h5h.normalization(cf.FILE_NAME, cf.HDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.HDSET_NAME, 100)) #Normalization of the H-field
h5h.normalization(cf.FILE_NAME, cf.EDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.EDSET_NAME, 100)) #Normalization of the E-field

hdf5_file.close_file()

# holi uwu
#grid.gaussian(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, 100, cf.COURANT, cf.IMP0)
#TFSF_ricker(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY, cf.COURANT, cf.IMP0)
#grid.TFSF_rampcos(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, cf.STEPS_WAVELENGTH, cf.RAMP_PERIODS, cf.COURANT, cf.IMP0)
