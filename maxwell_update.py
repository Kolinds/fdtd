from numba import njit
import numpy as np

@njit
def update_mag_field(start, end, hy, ez, chyh, chye):
        for m in range(start, end):
            hy[m] = chyh * hy[m] + chye * (ez[m + 1] - ez[m])

@njit
def update_elec_field(start, end, hy, ez, ceze, cezh):
        for m in range(start, end):
            ez[m] = ceze * ez[m] + cezh * (hy[m] - hy[m - 1])


def update_disp_elec(start, end, hy, ez, ez_temp, pol_current, coef_jj, coef_je, cpez_e, cpez_fp, cpez_dp):
        for m in range(start, end):
            ez_temp[m] = ez[m]
            ez[m] = cpez_e * ez[m] + cpez_fp * ((hy[m] - hy[m - 1]) - cpez_dp * pol_current[m])
            pol_current[m] = coef_jj * pol_current[m] + coef_je * (ez_temp[m] + ez[m])



"""
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
"""
