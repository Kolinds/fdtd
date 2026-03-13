import hdf5_handler as h5h
import grid_struct as gr
import maxwell_update as upd
import config as cf
import incident_field as df


hdf5_handler = h5h.HDF5Writer(cf.FILE_NAME, cf.TIME_BUFFER, cf.TOTAL_TIME, cf.SPACE_SIZE)
grid = gr.Grid(cf.SPACE_SIZE, cf.TOTAL_TIME)

hdf5_handler.open_file(cf.EDSET_NAME, cf.HDSET_NAME)
grid.place_materials(cf.SPACE_SIZE, cf.DIELECTRIC_LAYER, cf.LOSS, cf.LOSS_LAYER, cf.IMP0)

grid.add_probe(20, "Probe1", cf.TOTAL_TIME//2 + 1)
grid.add_probe(275, "Probe2", cf.TOTAL_TIME//2 + 1)

for qTime in range (0, cf.TOTAL_TIME):

    grid.hy = upd.update_magnetic_field(grid.hy, grid.ez, grid.chyh, grid.chye, cf.SPACE_SIZE - 1)
    grid.TFSF_rampcos(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, cf.STEPS_WAVELENGTH, cf.RAMP_PERIODS, cf.COURANT, cf.IMP0)
    grid.abc_boundary()
    grid.ez = upd.update_electric_field(grid.ez, grid.hy, grid.ceze, grid.cezh, cf.SPACE_SIZE)
    
    grid.r_DFT(qTime)

    hdf5_handler.update_file(qTime, grid.ez, grid.hy)

grid.save_probes(hdf5_handler.file)

incident_f = hdf5_handler.retrieve_array("/Probes/", "Probe1")
transmitted_f = hdf5_handler.retrieve_array("/Probes/", "Probe2")
trasmission_coef = df.transmission_c(275, incident_f, transmitted_f, cf.TOTAL_TIME //2 + 1, cf.TOTAL_TIME)
hdf5_handler.save_array("Coef/", "transmission", trasmission_coef)

h5h.normalization(cf.FILE_NAME, cf.HDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.HDSET_NAME, 100)) #Normalization of the H-field
h5h.normalization(cf.FILE_NAME, cf.EDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.EDSET_NAME, 100)) #Normalization of the E-field

hdf5_handler.close_file()

# holi uwu
#grid.TFSF_gaussian(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, 100, cf.COURANT, cf.IMP0)
#grid.TFSF_ricker(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY, cf.COURANT, cf.IMP0)
#grid.TFSF_rampcos(qTime, cf.TFSF_BOUNDARY, 50, 0.5, -0.5, 0, 0, cf.STEPS_WAVELENGTH, cf.RAMP_PERIODS, cf.COURANT, cf.IMP0)
