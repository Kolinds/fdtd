import hdf5_handler as h5h
import grid_struct as gr
import config as cf
import incident_field as incf


grid = gr.Grid(cf.TOTAL_TIME, cf.COURANT)

grid.initiate_materials()
grid.materials.add_free_space(200)
grid.materials.plasma_slab_ztransf(300, cf.DELTA_T, cf.CONDUCTIVITY, cf.RELAX_TIME_STEPS,
                               cf.PLASMA_WAVELENGTH_STEPS, cf.PERMITIVITY_INF)
grid.materials.add_free_space(200)
grid.confirm_materials()


grid.initiate_abc()

#grid.add_probe(250, "Probe1", cf.TOTAL_TIME//2 + 1)



hdf5_handler = h5h.HDF5Writer(cf.FILE_NAME, cf.TIME_BUFFER, cf.TOTAL_TIME, grid.space_size)
hdf5_handler.open_file(cf.EDSET_NAME, cf.HDSET_NAME)

for qTime in range (0, cf.TOTAL_TIME):
    grid.update_Hyfield()
    grid.apply_hyTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime, 50, 0, 0, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)

    grid.update_Ezfield()
    grid.apply_ezTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime, 50, 0.5, -0.5, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)
    grid.abc.first_order()
    

    grid.r_DFT(qTime)

    hdf5_handler.update_file(qTime, grid.ez, grid.hy)

grid.save_probes(hdf5_handler.file)

#h5h.normalization(cf.FILE_NAME, cf.HDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.HDSET_NAME, 100)) #Normalization of the H-field
#h5h.normalization(cf.FILE_NAME, cf.EDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.EDSET_NAME, 100)) #Normalization of the E-field

hdf5_handler.close_file()

# holi uwu
