import hdf5_handler as h5h
import grid_struct as gr
import config as cf
import incident_field as incf


hdf5_handler = h5h.HDF5Writer(cf.FILE_NAME, cf.TIME_BUFFER, cf.TOTAL_TIME, cf.SPACE_SIZE)
hdf5_handler.open_file(cf.EDSET_NAME, cf.HDSET_NAME)

grid = gr.Grid(cf.SPACE_SIZE, cf.TOTAL_TIME, cf.COURANT)
grid.set_dielectric(cf.DIELECTRIC_LAYER)
grid.initiate_abc()
#grid.add_probe(200, "Probe1", cf.TOTAL_TIME//2 + 1)


for qTime in range (0, cf.TOTAL_TIME):
    grid.update_Hyfield()
    grid.apply_hyTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime, 50, 0, 0, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)

    grid.update_Ezfield()
    grid.apply_ezTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime, 50, 0.5, -0.5, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)
    grid.update_abc_2order()
    

    #grid.r_DFT(qTime)

    hdf5_handler.update_file(qTime, grid.ez, grid.hy)

#grid.save_probes(hdf5_handler.file)

#h5h.normalization(cf.FILE_NAME, cf.HDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.HDSET_NAME, 100)) #Normalization of the H-field
#h5h.normalization(cf.FILE_NAME, cf.EDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.EDSET_NAME, 100)) #Normalization of the E-field

hdf5_handler.close_file()

# holi uwu
