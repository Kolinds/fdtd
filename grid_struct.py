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
        self.abc.second_order()

    

    def apply_hyTFSF(self, inc_func, tfsf_boundary, current_time, location, time_delay, loc_delay, *func_args):        
        self.hy[tfsf_boundary - 1] -= inc_func(current_time, location, time_delay, loc_delay, self.courant, *func_args) / self.imp0

    def apply_ezTFSF(self, inc_func, tfsf_boundary, current_time, location, time_delay, loc_delay, *func_args):
        self.ez[tfsf_boundary] += inc_func(current_time, location, time_delay, loc_delay, self.courant, *func_args)


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

        #First Order Initial Conditions
        self.abc1ezLeftOld = 0
        self.abc1ezRightOld = 0

        temp1 = np.sqrt(self.grid.cezh[0] * self.grid.chye[0])
        self.abc1CoefLeft = (temp1 - 1.0) / (temp1 + 1.0)
        temp2 = np.sqrt(self.grid.cezh[self.grid.space_size - 1] * self.grid.chye[self.grid.space_size - 2])
        self.abc1CoefRight = (temp2 - 1.0) / (temp2 + 1.0)


        #Second Order Initial Conditions
        self.abc2CoefLeft = np.zeros(3)
        self.abc2CoefRight = np.zeros(3)
        temp3 = 1.0 / temp1 + 2.0 + temp1

        self.abc2CoefLeft[0] = -(1.0 / temp1 - 2.0 + temp1) / temp3
        self.abc2CoefLeft[1] = -2.0 * (temp1 - 1.0 / temp1) / temp3
        self.abc2CoefLeft[2] = 4.0 * (temp1 + 1.0 / temp1) / temp3

        self.abc2CoefRight[0] = -(1.0 / temp2 - 2.0 + temp2) / temp3
        self.abc2CoefRight[1] = -2.0 * (temp2 - 1.0 / temp2) / temp3
        self.abc2CoefRight[2] = 4.0 * (temp2 + 1.0 / temp2) / temp3

        self.abc2ezOld2Left = np.zeros(3)
        self.abc2ezOld1Left = np.zeros(3)
        self.abc2ezOld2Right = np.zeros(3)
        self.abc2ezOld1Right = np.zeros(3)


    def first_order(self):
        self.grid.ez[0] = self.abc1ezLeftOld + self.abc1CoefLeft * (self.grid.ez[1] - self.grid.ez[0]) #El ez[0] aun no se ha actualizado
        self.grid.ez[self.grid.space_size - 1] = self.abc1ezRightOld + self.abc1CoefRight * (self.grid.ez[self.grid.space_size - 2] - self.grid.ez[self.grid.space_size - 1]) #El ez[0] aun no se ha actualizado
        
        self.abc1ezLeftOld = self.grid.ez[1]
        self.abc1ezRightOld = self.grid.ez[self.grid.space_size - 2]


    def second_order(self):
        self.grid.ez[0] = (self.abc2CoefLeft[0] * (self.grid.ez[2] + self.abc2ezOld2Left[0]) 
                          +self.abc2CoefLeft[1] * (self.abc2ezOld1Left[0] + self.abc2ezOld1Left[2] - self.grid.ez[1] -self.abc2ezOld2Left[1])
                          +self.abc2CoefLeft[2] * (self.abc2ezOld1Left[1]) - self.abc2ezOld2Left[2])

        self.grid.ez[self.grid.space_size - 1] = (self.abc2CoefRight[0] * (self.grid.ez[self.grid.space_size - 3] + self.abc2ezOld2Right[0]) 
                                                +self.abc2CoefRight[1] * (self.abc2ezOld1Right[0] + self.abc2ezOld1Right[2] - self.grid.ez[self.grid.space_size - 2] -self.abc2ezOld2Right[1])
                                                +self.abc2CoefRight[2] * (self.abc2ezOld1Right[1]) - self.abc2ezOld2Right[2])

        for m in range(0, 3):
            self.abc2ezOld2Left[m] = self.abc2ezOld1Left[m]
            self.abc2ezOld1Left[m] = self.grid.ez[m]

            self.abc2ezOld2Right[m] = self.abc2ezOld1Right[m]
            self.abc2ezOld1Right[m] = self.grid.ez[self.grid.space_size - 1 - m]



