import numpy as np
import incident_field as incf



class Grid():
    def __init__(self, space_size):
        #Arrays for E, H-fields
        self.ez= np.zeros(space_size)
        self.hy= np.zeros(space_size - 1)

        #Arrays for permitivity, permeability and loss of materials
        self.ceze = np.ones(space_size)
        self.cezh = np.ones(space_size)
        self.chyh = np.ones(space_size - 1)
        self.chye = np.ones(space_size - 1)

    def place_materials(self, space_size, loss, loss_layer, imp):
        #Establishes different materials
        for m in range(0, space_size - 1):
            if (m < loss_layer): #Free space
                self.chyh[m] = 1.0
                self.chye[m] = 1 / imp
            else:
                self.chyh[m] = (1.0 - loss) / (1.0 + loss)
                self.chye[m] = (1.0 / imp) / (1.0 + loss)

        for m in range(0, space_size):
            if (m < 100): #Free space
                self.ceze[m] = 1.0
                self.cezh[m] = imp
            elif (m < loss_layer):   #Dielectric
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
    
    def TFSF_boundary(self, current_time, location, etime_delay, eloc_delay, htime_delay, hloc_delay, width, courant, imp):
        self.hy[49] -= incf.ezfield(current_time, location, htime_delay, hloc_delay, width, courant) / imp #Correcion TFSF para quitar el efecto del incidente
        self.ez[50] += incf.ezfield(current_time, location, etime_delay, eloc_delay, width, courant) #Correcion TFSF para añadir el efecto del incidente

        """ 
        Only equal: = -> Hardwired source
        Summed: += -> Source 
        ez[50] += np.exp((-(30 - 30)**2)/100)
        """



