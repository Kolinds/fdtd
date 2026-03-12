import numpy as np
import incident_field as incf
import h5py


class Grid():
    def __init__(self, space_size, total_time):
        self.space_size = space_size
        self.total_time = total_time

        #Arrays for E, H-fields
        self.ez= np.zeros(space_size)
        self.hy= np.zeros(space_size - 1)

        #Arrays for permitivity, permeability and loss of materials
        self.ceze = np.ones(space_size)
        self.cezh = np.ones(space_size)
        self.chyh = np.ones(space_size - 1)
        self.chye = np.ones(space_size - 1)

        #List for storing probing arrays
        self.stored_probes = []

    def add_probe(self, location, array_name, array_size):
        new_probe = np.zeros(array_size, dtype=np.complex64)
        probe_data = {
            "name": array_name,
            "size": array_size,
            "location": location,
            "array": new_probe
        }
        self.stored_probes.append(probe_data)

    def save_probes(self, file: h5py.File):
        for actual_probe in self.stored_probes:
            group = file.require_group("/probes/")
            group.create_dataset(actual_probe["name"], data=actual_probe["array"])
            
        
    def place_materials(self, space_size, dielectric_layer, loss, loss_layer, imp):
        #Magnetic field material properties
        for m in range(0, space_size - 1):
            if (m < loss_layer): #Free space
                self.chyh[m] = 1.0
                self.chye[m] = 1 / imp
            else:
                self.chyh[m] = (1.0 - loss) / (1.0 + loss)
                self.chye[m] = (1.0 / imp) / (1.0 + loss)

        #Electric field material properties
        for m in range(0, space_size):
            if (m < dielectric_layer): #Free space
                self.ceze[m] = 1.0
                self.cezh[m] = imp
            elif (m < loss_layer): #Dielectric
                self.ceze[m] = 1.0
                self.cezh[m] = imp / 9.0
            else: #lossy Boundary layer
                self.ceze[m] = (1.0 - loss) / (1.0 + loss)
                self.cezh[m] = (imp / 9.0) / (1.0 + loss)

    def abc_boundary(self):
        # ABC's in 1D with Courant Limit 
        self.ez[0] = self.ez[1]

        """
        Absorbing boundary condition (ABC) mag, deprecated
        hy[SPACE_SIZE-1] = hy[SPACE_SIZE-2]

        This ABC does not work inside dielectric
        ez[SPACE_SIZE-1] = ez[SPACE_SIZE-2]
        """
    
    def TFSF_gaussian(self, current_time, tfsf_boundary, location, etime_delay, eloc_delay, htime_delay, hloc_delay, width, courant, imp):
        self.hy[tfsf_boundary - 1] -= incf.ezfield(current_time, location, htime_delay, hloc_delay, width, courant) / imp #Correcion TFSF para quitar el efecto del incidente
        self.ez[tfsf_boundary] += incf.ezfield(current_time, location, etime_delay, eloc_delay, width, courant) #Correcion TFSF para añadir el efecto del incidente

    def TFSF_rampcos(self, current_time, tfsf_boundary, location, etime_delay, eloc_delay, htime_delay, hloc_delay, steps_wavelength, ramp_periods, courant, imp):
        self.hy[tfsf_boundary - 1] -= incf.rampf(current_time, location, htime_delay, hloc_delay, steps_wavelength, ramp_periods, courant) / imp #Correcion TFSF para quitar el efecto del incidente
        self.ez[tfsf_boundary] += incf.rampf(current_time, location, etime_delay, eloc_delay, steps_wavelength, ramp_periods, courant) #Correcion TFSF para añadir el efecto del incidente

    def TFSF_ricker(self, current_time, tfsf_boundary, location, etime_delay, eloc_delay, htime_delay, hloc_delay, steps_pkwavelength, delay_multiple, courant, imp):
        self.hy[tfsf_boundary - 1] -= incf.ricker(current_time, location, htime_delay, hloc_delay, steps_pkwavelength, delay_multiple, courant) / imp #Correcion TFSF para quitar el efecto del incidente
        self.ez[tfsf_boundary] += incf.ricker(current_time, location, etime_delay, eloc_delay, steps_pkwavelength, delay_multiple, courant) #Correcion TFSF para añadir el efecto del incidente

    """ 
        Only equal: = -> Hardwired source
        Summed: += -> Source 
        ez[50] += np.exp((-(30 - 30)**2)/100)
    """

    #Running discrete Fourier Transform
    def r_DFT(self, current_time):
        for actual_probe in self.stored_probes:
            freq_array = np.arange(actual_probe["size"])
            angles_array = (2*np.pi * freq_array * current_time)/self.total_time
            actual_probe["array"].real += self.ez[actual_probe["location"]]*np.cos(angles_array)
            actual_probe["array"].imag -= self.ez[actual_probe["location"]]*np.sin(angles_array)





