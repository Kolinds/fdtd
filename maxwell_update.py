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

@njit
def update_disp_elec_ADE(start, end, hy, ez, ez_temp, pol_current, coef_jj, coef_je, cpez_e, cpez_fp, cpez_dp):
        for m in range(start, end):
            ez_temp[m] = ez[m]
            ez[m] = cpez_e * ez[m] + cpez_fp * ((hy[m] - hy[m - 1]) - cpez_dp * pol_current[m])
            pol_current[m] = coef_jj * pol_current[m] + coef_je * (ez_temp[m] + ez[m])

@njit
def update_disp_elec_PLRC(start, end, hy, ez, ez_temp, rec_accumulator, cez_ez, cez_hy, cez_accum, caccum_ezf, caccum_ezp, caccum_accum):
        for m in range(start, end):
            ez_temp[m] = ez[m]
            ez[m] = cez_ez * ez[m] + cez_hy * (hy[m] - hy[m - 1]) + cez_accum * rec_accumulator[m]
            rec_accumulator[m] = caccum_ezf * ez[m] + caccum_ezp * ez_temp[m] + caccum_accum * rec_accumulator[m]

@njit
def update_disp_elec_ztransf(start, end, hy, ez, d_field, integrator, low_pass, cezd, cezi, cezl, clows, imp0, courant):
        for m in range(start, end):
            ez[m] = ez[m] * (1/imp0)
            d_field[m] = d_field[m] + courant * (hy[m] - hy[m - 1])
            ez[m] = cezd * d_field[m] - cezi * integrator[m] - cezl * low_pass[m]
            integrator[m] = integrator[m] + ez[m]
            low_pass[m] = clows * low_pass[m] + ez[m]
            ez[m] = ez[m] * imp0

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
