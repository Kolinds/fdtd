import numpy as np
import h5py
import maxwell_update as upd
import incident_field as incf



class Grid():
    def __init__(self, total_time, courant):
        self.total_time = total_time
        self.courant = courant
        self.imp0 = 377.0
        self.permitivity0 = 8.854e-12

        #List for storing probing arrays
        self.stored_probes = []
    

    def initiate_materials(self):
        self.materials = Material_placement(self)

    def confirm_materials(self):
        self.space_size = 0
        for action in self.materials.ez_action_sequences:
            self.space_size += action[1]

        self.ez= np.zeros(self.space_size)
        self.hy= np.zeros(self.space_size - 1)


    def initiate_abc(self):
        self.abc = Abc_conditions(self)


    def apply_hyTFSF(self, inc_func, tfsf_boundary, current_time, location, time_delay, loc_delay, *func_args):        
        self.hy[tfsf_boundary - 1] -= inc_func(current_time, location, time_delay, loc_delay, self.courant, *func_args) / self.imp0

    def apply_ezTFSF(self, inc_func, tfsf_boundary, current_time, location, time_delay, loc_delay, *func_args):
        self.ez[tfsf_boundary] += inc_func(current_time, location, time_delay, loc_delay, self.courant, *func_args)




    def reset_fields(self):
        self.ez= np.zeros(self.space_size)
        self.hy= np.zeros(self.space_size - 1)


    def update_Hyfield(self):
        initial_step = 0
        final_step = 0
        for action, width, arguments in self.materials.hy_action_sequences:
            final_step += width

            field_updf = self.materials.function_map[action]
            field_updf(initial_step, final_step, self.hy, self.ez, **arguments)
                #miauuuu :3
            initial_step += width
            
    def update_Ezfield(self):
        initial_step = 1
        final_step = 1
        for action, width, arguments in self.materials.ez_action_sequences:
            final_step += width

            field_updf = self.materials.function_map[action]
            field_updf(initial_step, final_step, self.hy, self.ez, **arguments)
                #miauuuu :3
            initial_step += width
            
  


    #Running discrete Fourier Transform
    def r_DFT(self, current_time):
        for actual_probe in self.stored_probes:
            """
            incf.running_DFT(actual_probe["array"], actual_probe["location"], self.ez, self.total_time, current_time)
            """
            freq_array = np.arange(actual_probe["size"])
            angles_array = (2*np.pi * freq_array * current_time)/self.total_time
            actual_probe["array"].real += self.ez[actual_probe["location"]]*np.cos(angles_array)
            actual_probe["array"].imag -= self.ez[actual_probe["location"]]*np.sin(angles_array)


    #Set probes 
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



class Material_placement():
    def __init__(self, grid: Grid):
        self.grid = grid

        self.function_map = {"ez_basic": upd.update_elec_field, 
                        "hy_basic": upd.update_mag_field,
                        "ez_dispersive_ADE": upd.update_disp_elec_ADE,
                        "ez_dispersive_PLRC": upd.update_disp_elec_PLRC,
                        "ez_dispersive_ztransf": upd.update_disp_elec_ztransf}
        
        self.hy_action_sequences=[]
        self.ez_action_sequences=[]

    def clear_space(self):
        self.hy_action_sequences.clear()
        self.ez_action_sequences.clear()


    def add_free_space(self, width):
        #Magnetic material
        chyh = 1.0
        chye = 1 / self.grid.imp0
        dictionary = {"chyh": chyh, "chye": chye}
        self.hy_action_sequences.append(("hy_basic", width, dictionary))

        #Electric material
        ceze = 1.0
        cezh = self.grid.imp0 
        dictionary = {"ceze": ceze, "cezh": cezh} 
        self.ez_action_sequences.append(("ez_basic", width, dictionary))

    def add_dielectric(self, width):
        #Magnetic material
        chyh = 1.0
        chye = 1 / self.grid.imp0
        dictionary = {"chyh": chyh, "chye": chye}
        self.hy_action_sequences.append(("hy_basic", width, dictionary))

        ceze = 1.0
        cezh = self.grid.imp0 / 9.0  
        dictionary = {"ceze": ceze, "cezh": cezh}
        self.ez_action_sequences.append(("ez_basic", width, dictionary))
    
    def add_lossy_material(self, loss, width):
        #Magnetic material
        chyh = (1.0 - loss) / (1.0 + loss)
        chye = (1.0 / self.grid.imp0) / (1.0 + loss)
        dictionary = {"chyh": chyh, "chye": chye}
        self.hy_action_sequences.append(("hy_basic", width, dictionary))

        #Electric field material properties
        ceze = (1.0 - loss) / (1.0 + loss)
        cezh = (self.grid.imp0 / 9.0) / (1.0 + loss)
        dictionary = {"ceze": ceze, "cezh": cezh}
        self.ez_action_sequences.append(("ez_basic", width, dictionary))
    
    def plasma_slab_ADE(self, width, delta_t, conductivity, relax_time, plasma_wavelength, permitivity_inf):
        #Magnetic material
        chyh = 1.0
        chye = 1 / self.grid.imp0
        dictionary = {"chyh": chyh, "chye": chye}
        self.hy_action_sequences.append(("hy_basic", width, dictionary))

        #Electric field material properties
        #-> Plasma slab
        ez_temp = np.zeros(width)
        pol_current = np.zeros(width)

        coef_jj = (1 - 1/(2*relax_time)) / (1 + 1/(2*relax_time))
        coef_je = (1 / (1 + 1/(2*relax_time))) * ((2 * (np.pi)**2 *self.grid.courant) / (self.grid.imp0 * plasma_wavelength**2))

        c_1 = (conductivity * delta_t) / (2* permitivity_inf * self.grid.permitivity0)
        c_2 = (coef_je * self.grid.imp0 * self.grid.courant) / (2 * permitivity_inf)

        cpez_e = (1 - c_1 - c_2) / (1 + c_1 + c_2)
        cpez_fp = ((self.grid.imp0 * self.grid.courant) / permitivity_inf) / (1 + c_1 + c_2)
        cpez_dp = 0.5 * (1 + coef_jj)

        dictionary = {"ez_temp": ez_temp, "pol_current": pol_current, "coef_jj": coef_jj, "coef_je":coef_je,
                       "cpez_e": cpez_e, "cpez_fp": cpez_fp, "cpez_dp":cpez_dp}

        self.ez_action_sequences.append(("ez_dispersive_ADE", width, dictionary))


    def plasma_slab_PLRC(self, width, delta_t, conductivity, relax_time, plasma_wavelength, permitivity_inf):
        #Magnetic material
        chyh = 1.0
        chye = 1 / self.grid.imp0
        dictionary = {"chyh": chyh, "chye": chye}
        self.hy_action_sequences.append(("hy_basic", width, dictionary))

        #Electric field material properties
        #-> Plasma slab
        ez_temp = np.zeros(width)
        rec_accumulator = np.zeros(width)

        c_rec = np.exp(-1/relax_time)
        common_factor = ((2 * np.pi * self.grid.courant)/(plasma_wavelength))**2

        chi0 = relax_time * common_factor * (1 - relax_time * (1 -c_rec))
        xi0 = relax_time * common_factor * (0.5 - relax_time**2 * (1 - (1 + 1/relax_time) * c_rec))
        delta_chi0 = - relax_time**2 * common_factor * (1 - c_rec)**2
        delta_xi0 = - relax_time**3 * common_factor * (1 - (1 + 1/relax_time)*c_rec)*(1 - c_rec)

        denominator = permitivity_inf + chi0 - xi0

        cez_ez = (permitivity_inf - xi0) / denominator
        cez_hy = (self.grid.courant * self.grid.imp0) / denominator
        cez_accum = 1 / denominator

        caccum_ezf = delta_chi0 - delta_xi0
        caccum_ezp = delta_xi0
        caccum_accum = c_rec

        dictionary = {"ez_temp": ez_temp, "rec_accumulator": rec_accumulator, 
                      "cez_ez": cez_ez, "cez_hy":cez_hy, "cez_accum":cez_accum,
                       "caccum_ezf": caccum_ezf, "caccum_ezp": caccum_ezp, "caccum_accum":caccum_accum}

        self.ez_action_sequences.append(("ez_dispersive_PLRC", width, dictionary))


    def plasma_slab_ztransf(self, width, delta_t, conductivity, relax_time, plasma_wavelength, permitivity_inf):
        #Magnetic material
        chyh = 1.0
        chye = 1 / self.grid.imp0
        dictionary = {"chyh": chyh, "chye": chye}
        self.hy_action_sequences.append(("hy_basic", width, dictionary))

        #Electric field material properties
        #-> Plasma slab
        integrator = np.zeros(width)
        low_pass = np.zeros(width)
        d_field = np.zeros(width)

        cezd = (1 / permitivity_inf)
        cezi = ((2 * np.pi * self.grid.courant)/(plasma_wavelength))**2 * (1/permitivity_inf) * relax_time
        cezl = - cezi * np.exp(-1/relax_time)

        clows = np.exp(-1/relax_time)
     
        dictionary = {"d_field": d_field, "integrator": integrator,  "low_pass": low_pass,
                     "cezd":cezd, "cezi":cezi, "cezl": cezl, "clows": clows, 
                     "imp0": self.grid.imp0, "courant": self.grid.courant}

        self.ez_action_sequences.append(("ez_dispersive_ztransf", width, dictionary))




class Abc_conditions():
    def __init__(self, grid: Grid):
        self.grid = grid

        #First Order Initial Conditions
        self.abc1ezLeftOld = 0
        self.abc1ezRightOld = 0

        cezh_first = self.grid.materials.ez_action_sequences[0][2]["cezh"]
        chye_first = self.grid.materials.hy_action_sequences[0][2]["chye"]
        cezh_last = self.grid.materials.ez_action_sequences[-1][2]["cezh"]
        chye_last = self.grid.materials.hy_action_sequences[-1][2]["chye"]

        temp1 = np.sqrt(cezh_first * chye_first)
        self.abc1CoefLeft = (temp1 - 1.0) / (temp1 + 1.0)
        temp2 = np.sqrt(cezh_last * chye_last)
        self.abc1CoefRight = (temp2 - 1.0) / (temp2 + 1.0)


        #Second Order Initial Conditions
        self.abc2CoefLeft = np.zeros(3)
        self.abc2CoefRight = np.zeros(3)
        temp3 = 1.0 / temp1 + 2.0 + temp1
        temp4 = 1.0 / temp2 + 2.0 + temp2

        self.abc2CoefLeft[0] = -(1.0 / temp1 - 2.0 + temp1) / temp3
        self.abc2CoefLeft[1] = -2.0 * (temp1 - 1.0 / temp1) / temp3
        self.abc2CoefLeft[2] = 4.0 * (temp1 + 1.0 / temp1) / temp3

        self.abc2CoefRight[0] = -(1.0 / temp2 - 2.0 + temp2) / temp4
        self.abc2CoefRight[1] = -2.0 * (temp2 - 1.0 / temp2) / temp4
        self.abc2CoefRight[2] = 4.0 * (temp2 + 1.0 / temp2) / temp4

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



