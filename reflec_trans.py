import hdf5_handler as h5h
import grid_struct as gr
import config as cf
import incident_field as incf
import h5py
import numpy as np

OUTPOUT_FNAME = "reflec_trans6.hdf5"

#Creamos el grid
grid1 = gr.Grid(cf.TOTAL_TIME, cf.COURANT, cf.MAX_FREQ, cf.MIN_RPERMITIVITY, cf.MIN_RPERMEABILITY)

mapped_frequencies = int(cf.MAX_FREQ / grid1.delta_f)


loc_probe_trans = 4000
loc_probe_ref = 900

#Creamos el grid
grid1 = gr.Grid(cf.TOTAL_TIME, cf.COURANT, cf.MAX_FREQ, cf.MIN_RPERMITIVITY, cf.MIN_RPERMEABILITY)

#Establecemos los materiales
grid1.initiate_materials()

grid1.materials.add_free_space(3000)
grid1.materials.implicit_plasma_ADE(300, cf.RELAX_TIME_STEPS, 
                                   cf.PLASMA_WAVELENGTH_STEPS, cf.E_CONDUCTIVITY, cf.PERMITIVITY_INF)
#grid1.materials.eplasma_slab_ztransf(150, cf.E_CONDUCTIVITY, cf.RELAX_TIME_STEPS, cf.PLASMA_WAVELENGTH_STEPS, cf.PERMITIVITY_INF)
#grid1.materials.add_free_mag(150)
grid1.materials.add_free_space(3000)

grid1.confirm_materials()

grid1.initiate_abc()

#Abrimos el gestor del archivo hdf5 y ejecutamos el civlo principal
hdf5_handler1 = h5h.HDF5Writer("plasma_data.hdf5", cf.TIME_BUFFER, cf.TOTAL_TIME, grid1.space_size)
hdf5_handler1.open_file(cf.EDSET_NAME, cf.HDSET_NAME)

grid1.add_probe(loc_probe_trans, "Probe1", mapped_frequencies)
grid1.add_probe(loc_probe_ref, "Probe2", mapped_frequencies)


for qTime in range(0, cf.TOTAL_TIME):
    grid1.ez_old[:] = grid1.ez[:]
    grid1.hy_old[:] = grid1.hy[:]

    grid1.update_Ezfield()
    grid1.apply_ezTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime, 50, 0.5, -0.5, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)


    grid1.abc.second_order()

    grid1.update_Hyfield()
    grid1.apply_hyTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime + 1, 50, 0, 0, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)

    grid1.r_DFT(qTime)
    hdf5_handler1.update_file(qTime, grid1.ez, grid1.hy)

print("ACheck2")
grid1.save_probes(hdf5_handler1.file)

print("Check1")


vacuum_fref = np.zeros(mapped_frequencies, dtype = np.complex128)
for qTime in range(0, cf.TOTAL_TIME):
    ez_t = incf.ricker(qTime, loc_probe_ref, 0, 0,
                                              grid1.courant, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)
    
    freq_array = np.arange(mapped_frequencies)
    angles_array = (2 * np.pi * freq_array * qTime) / cf.TOTAL_TIME
    kernel = np.exp(-1j * angles_array)
    vacuum_fref += (ez_t) * kernel

vacuum_fref = (1/grid1.space_size) * vacuum_fref

plasma_ftrans = hdf5_handler1.retrieve_array("/Probes/", "Probe1")
transmission_array = incf.gtransmission_coef(cf.TOTAL_TIME, 100, vacuum_fref, plasma_ftrans, mapped_frequencies, 
                                       cf.PLASMA_WAVELENGTH_STEPS, cf.RELAX_TIME_STEPS, grid1.light_speed, grid1.delta_x, grid1.delta_t)



plasma_fref = hdf5_handler1.retrieve_array("/Probes/", "Probe2")

reflection_array = np.zeros(mapped_frequencies, dtype = np.complex128)
reflection_array[0] = 0
#transmission_array = np.zeros(mapped_frequencies, dtype = np.complex128)
#reflection_array[0] = 0
for i in range(1, plasma_fref.size):
    reflection_array[i] = plasma_fref[i] / vacuum_fref[i]
#    transmission_array[i] = plasma_ftrans[i] / vacuum_fref[i]

sum_array = abs(reflection_array)**2 + abs(transmission_array)**2


with h5py.File(f"{OUTPOUT_FNAME}", "w") as f:
    f.create_dataset("transmission", data=transmission_array)
    f.create_dataset("reflection", data=reflection_array)
    f.create_dataset("sum_coef", data=sum_array)
    f.create_dataset("vacuum_fref", data=vacuum_fref)
    f.create_dataset("plasma_fref", data=plasma_fref)
    f.create_dataset("plasma_ftrans", data=plasma_fref)



#h5h.normalization(cf.FILE_NAME, cf.HDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.HDSET_NAME, 100)) #Normalization of the H-field
#h5h.normalization(cf.FILE_NAME, cf.EDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.EDSET_NAME, 100)) #Normalization of the E-field
#h5h.normalization(cf.FILE_NAME, "Probes/Probe1", cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, "Probes/Probe1", 100))
print("Check2")
hdf5_handler1.close_file()
print("Simulación terminada de manera segura.")

# holi uwu
