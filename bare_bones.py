import hdf5_handler as h5h
import grid_struct as gr
import config as cf
import incident_field as incf


#Creamos el grid
grid = gr.Grid(cf.TOTAL_TIME, cf.COURANT, cf.MAX_FREQ, cf.MIN_RPERMITIVITY, cf.MIN_RPERMEABILITY)

#Establecemos los materiales
grid.initiate_materials()

grid.materials.add_free_space(3000)
#grid.materials.implicit_plasma_ADE(300, cf.RELAX_TIME_STEPS, 
#                                   cf.PLASMA_WAVELENGTH_STEPS, cf.E_CONDUCTIVITY, cf.PERMITIVITY_INF)
grid.materials.eplasma_slab_ADE(10000, cf.E_CONDUCTIVITY, cf.RELAX_TIME_STEPS, cf.PLASMA_WAVELENGTH_STEPS, cf.PERMITIVITY_INF)
grid.materials.add_free_mag(10000)
grid.materials.add_free_space(1)



grid.confirm_materials()

#Inicializamos las condiciones de contorno
grid.initiate_abc()

#Añadimos probes de medición en ciertos puntos
mapped_frequencies = int(cf.MAX_FREQ / grid.delta_f)
grid.add_probe(3450, "Transmitted", mapped_frequencies)
grid.add_probe(900, "Reflected", mapped_frequencies)


#Abrimos el gestor del archivo hdf5 y ejecutamos el civlo principal
hdf5_handler = h5h.HDF5Writer(cf.FILE_NAME, cf.TIME_BUFFER, cf.TOTAL_TIME, grid.space_size)
hdf5_handler.open_file(cf.EDSET_NAME, cf.HDSET_NAME)

# Bucle temporal corregido con resincronización TFSF
for qTime in range(0, cf.TOTAL_TIME):
    # Guardar el pasado para el ADE implícito
    grid.ez_old[:] = grid.ez[:]
    grid.hy_old[:] = grid.hy[:]

    # 1. Actualización Eléctrica Implícita (Calcula Ez^(n+1))
    grid.update_Ezfield()
    
    # E-TFSF se mantiene con qTime (inyecta H_inc^(n+1/2) en Ez^(n+1))
    grid.apply_ezTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime, 50, 0.5, -0.5, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)
    
    # Absorción en fronteras externas
    grid.abc.second_order()

    # 2. Actualización Magnética Explícita (Calcula Hy^(n+3/2))
    grid.update_Hyfield()
    
    # ¡SOLUCIÓN! Como Hy avanzó a n+3/2, necesitamos corregir usando E_inc^(n+1).
    # Pasamos qTime + 1 para desplazar la forma de onda analítica al instante correcto.
    grid.apply_hyTFSF(incf.ricker, cf.TFSF_BOUNDARY, qTime + 1, 50, 0, 0, cf.STEPS_WAVELENGTH, cf.RICKER_DELAY)
    
    # Diagnósticos y guardado
    grid.r_DFT(qTime)
    hdf5_handler.update_file(qTime, grid.ez, grid.hy)

print("Check1")
grid.save_probes(hdf5_handler.file)


#h5h.normalization(cf.FILE_NAME, cf.HDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.HDSET_NAME, 100)) #Normalization of the H-field
#h5h.normalization(cf.FILE_NAME, cf.EDSET_NAME, cf.BUFFER_JUMP, h5h.maxValue(cf.FILE_NAME, cf.EDSET_NAME, 100)) #Normalization of the E-field
norma_const = h5h.maxValue(cf.FILE_NAME, "Probes/Reflected", 100)
h5h.normalization(cf.FILE_NAME, "Probes/Reflected", cf.BUFFER_JUMP, norma_const)
h5h.normalization(cf.FILE_NAME, "Probes/Transmitted", cf.BUFFER_JUMP, norma_const)
print("Check2")
hdf5_handler.close_file()
print("Simulación terminada de manera segura.")

# holi uwu
