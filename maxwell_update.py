from numba import njit
import numpy as np

@njit
def update_magnetic_field(hy, ez, SPACE_SIZE, chyh, chye):
    for m in range(0, SPACE_SIZE - 1):
        hy[m] = chyh[m] * hy[m] + chye[m] * (ez[m + 1] - ez[m])  
    return hy

@njit
def update_electric_field(hy, ez, SPACE_SIZE, ceze, cezh):
    for m in range(1, SPACE_SIZE - 1):
        ez[m] = ceze[m] * ez[m] + cezh[m] * (hy[m] - hy[m - 1])
    return ez



@njit
def update_mag_field(m, hy, ez, chyh, chye):
        hy[m] = chyh * hy[m] + chye * (ez[m + 1] - ez[m])  

@njit
def update_elec_field(m, hy, ez, ceze, cezh):
        ez[m] = ceze * ez[m] + cezh * (hy[m] - hy[m - 1])



