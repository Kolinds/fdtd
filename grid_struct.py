import numpy as np
import h5py
import maxwell_update as upd



class Grid():
    def __init__(self, space_size, total_time, courant):
        self.space_size = space_size
        self.total_time = total_time
        self.courant = courant
        self.imp0 = 337.0

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
    
        
    def set_free_space(self):
        #Magnetic field material properties
        for m in range(0, self.space_size - 1):
                self.chyh[m] = 1.0
                self.chye[m] = 1 / self.imp0

        #Electric field material properties
        for m in range(0, self.space_size): 
                self.ceze[m] = 1.0
                self.cezh[m] = self.imp0    
    
    def place_materials(self, dielectric_layer, loss, loss_layer):
        #Magnetic field material properties
        for m in range(0, self.space_size - 1):
            if (m < loss_layer): #Free space
                self.chyh[m] = 1.0
                self.chye[m] = 1 / self.imp0
            else:
                self.chyh[m] = (1.0 - loss) / (1.0 + loss)
                self.chye[m] = (1.0 / self.imp0) / (1.0 + loss)

        #Electric field material properties
        for m in range(0, self.space_size):
            if (m < dielectric_layer): #Free space
                self.ceze[m] = 1.0
                self.cezh[m] = self.imp0
            elif (m < loss_layer): #Dielectric
                self.ceze[m] = 1.0
                self.cezh[m] = self.imp0 / 9.0
            else: #lossy Boundary layer
                self.ceze[m] = (1.0 - loss) / (1.0 + loss)
                self.cezh[m] = (self.imp0 / 9.0) / (1.0 + loss)



    def reset_fields(self):
        self.ez= np.zeros(self.space_size)
        self.hy= np.zeros(self.space_size - 1)

    def update_Hyfield(self):
        self.hy = upd.update_magnetic_field(self.hy, self.ez, self.chyh, self.chye, self.space_size)

    def update_Ezfield(self):
        self.ez = upd.update_electric_field(self.ez, self.hy, self.ceze, self.cezh, self.space_size)



    def initiate_abc(self):
        self.abc = Abc_conditions(self)
    
    def update_abc_1order(self):
        self.abc.first_order()

    def update_abc_2order(self):
        pass

    """
        ABC's in 1D with Courant Limit 
        self.ez[0] = self.ez[1]

        Absorbing boundary condition (ABC), deprecated. 
        hy[SPACE_SIZE-1] = hy[SPACE_SIZE-2]
        ez[SPACE_SIZE-1] = ez[SPACE_SIZE-2] This ABC does not work inside dielectric
    """
    

    def apply_hyTFSF(self, inc_func, tfsf_boundary, current_time, location, time_delay, loc_delay, *func_args):        
        self.hy[tfsf_boundary - 1] -= inc_func(current_time, location, time_delay, loc_delay, self.courant, *func_args) / self.imp0

    def apply_ezTFSF(self, inc_func, tfsf_boundary, current_time, location, time_delay, loc_delay, *func_args):
        self.ez[tfsf_boundary] += inc_func(current_time, location, time_delay, loc_delay, self.courant, *func_args)

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
            group = file.require_group("/Probes/")
            group.create_dataset(actual_probe["name"], data=(1/self.space_size)*actual_probe["array"])
        self.stored_probes.clear()





class Abc_conditions():
    def __init__(self, grid: Grid):
        self.grid = grid

        self.ezLeftOld = self.grid.ez[1]
        self.ezRightOld = self.grid.ez[self.grid.space_size - 2]

        temp1 = np.sqrt(self.grid.cezh[0] * self.grid.chye[0])
        self.abcCoefLeft = (temp1 - 1.0) / (temp1 + 1.0)
        temp2 = np.sqrt(self.grid.cezh[self.grid.space_size - 1] * self.grid.chye[self.grid.space_size - 2])
        self.abcCoefRight = (temp2 - 1.0) / (temp2 + 1.0)

    def first_order(self):
        self.grid.ez[0] = self.ezLeftOld + self.abcCoefLeft * (self.grid.ez[1] - self.grid.ez[0]) #El ez[0] aun no se ha actualizado
        self.ezLeftOld = self.grid.ez[1]

        self.grid.ez[self.grid.space_size - 1] = self.ezRightOld + self.abcCoefRight * (self.grid.ez[self.grid.space_size - 2] - self.grid.ez[self.grid.space_size - 1]) #El ez[0] aun no se ha actualizado
        self.ezRightOld = self.grid.ez[self.grid.space_size - 2]


