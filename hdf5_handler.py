import numpy as np
import h5py

class HDF5Writer:
    def __init__(self, file_name, time_buffer, total_time, space_size):
        self.file_name = file_name
        self.time_buffer = time_buffer
        self.total_time = total_time
        self.space_size = space_size
        self.file = None
        self.ez_dset = None
        self.hy_dset = None
        self.ez_buffer= None
        self.hy_buffer= None

    def open_file(self, edset_name, hdset_name):
        #Open file and initialize datasets and buffers, for E, H-fields
        self.file = h5py.File(self.file_name, "w")

        self.ez_dset = self.file.create_dataset(edset_name, (self.total_time, self.space_size))
        self.hy_dset = self.file.create_dataset(hdset_name, (self.total_time, self.space_size - 1))
        self.ez_buffer= np.zeros((self.time_buffer, self.space_size))
        self.hy_buffer= np.zeros((self.time_buffer, self.space_size - 1))


    def update_file(self, actual_time, ez_array, hy_array):
        buffer_index = actual_time % self.time_buffer
        if (buffer_index == 0)  and (actual_time != 0):
            self.ez_dset[actual_time - self.time_buffer:actual_time, 0:self.space_size] = self.ez_buffer
            self.hy_dset[actual_time - self.time_buffer:actual_time, 0:self.space_size - 1] = self.hy_buffer

        self.ez_buffer[buffer_index, :] = ez_array
        self.hy_buffer[buffer_index, :] = hy_array

    def save_array(self, location_group, dataset_name, array_data):
        if self.file is None:
            print("Error: El archivo no está abierto.")
            return
        group = self.file.require_group(location_group)
        group.create_dataset(dataset_name, data=array_data)

    def retrieve_array(self, location_group: str, dataset_name: str):
        if self.file is None:
            print("Error: El archivo no está abierto.")
            return
        
        with h5py.File(self.file_name, "r") as f:
            full_path = f"{location_group.rstrip('/')}/{dataset_name}"
            
            if full_path in f:
                return f[full_path][:]
            else:
                print(f"Error: No se encontró '{full_path}' en {self.file_name}")
                return None

    def close_file(self):
        if self.file is not None:
            self.file.close()



def maxValue(file_name, dataset_name, jump):
    total_max = 0
    with h5py.File(file_name, "r") as f:
        dset = f[dataset_name]
        for actual_line in range(0, dset.shape[0], jump):
            jump_matrix = dset[actual_line:actual_line+jump, ...]
            actual_max = np.max(np.abs(jump_matrix))
            if (total_max < actual_max): total_max = actual_max
    return total_max


def normalization(file_name, dataset_name, jump, max_val):
    with h5py.File(file_name, "r+") as f:
        dset = f[dataset_name]
        for actual_line in range(0, dset.shape[0], jump):
            jump_matrix = dset[actual_line:actual_line+jump, ...]
            norm_matrix = jump_matrix / max_val
            dset[actual_line:actual_line+jump, ...] = norm_matrix