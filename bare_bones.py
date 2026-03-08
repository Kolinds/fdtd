import hdf5_handler as h5h
import grid_struct as gr
import maxwell_update as upd
import config as cg


hdf5_file = h5h.HDF5Writer(cg.FILE_NAME, cg.TIME_BUFFER, cg.TOTAL_TIME, cg.SPACE_SIZE)
grid = gr.Grid(cg.SPACE_SIZE)

hdf5_file.open_file(cg.EDSET_NAME, cg.HDSET_NAME)
grid.place_materials(cg.SPACE_SIZE, cg.LOSS, cg.LOSS_LAYER, cg.IMP0)

for qTime in range (0, cg.TOTAL_TIME):

    grid.hy = upd.update_magnetic_field(grid.hy, grid.ez, grid.chyh, grid.chye, cg.SPACE_SIZE - 1)
    grid.TFSF_boundary(qTime, 30, 0.5, -0.5, 0, 0, 100, 1, cg.IMP0)
    grid.abc_boundary()
    grid.ez = upd.update_electric_field(grid.ez, grid.hy, grid.ceze, grid.cezh, cg.SPACE_SIZE)
    
    hdf5_file.update_file(qTime, grid.ez, grid.hy)

h5h.normalization(cg.FILE_NAME, cg.HDSET_NAME, cg.BUFFER_JUMP, h5h.maxValue(cg.FILE_NAME, cg.HDSET_NAME, 100)) #Normalization of the H-field
h5h.normalization(cg.FILE_NAME, cg.EDSET_NAME, cg.BUFFER_JUMP, h5h.maxValue(cg.FILE_NAME, cg.EDSET_NAME, 100)) #Normalization of the E-field

hdf5_file.close_file()

# holi uwu

