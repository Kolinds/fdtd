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
            m_retarded = m - start
            ez_temp[m_retarded] = ez[m]
            ez[m] = cpez_e * ez[m] + cpez_fp * ((hy[m] - hy[m - 1]) - cpez_dp * pol_current[m_retarded])
            pol_current[m_retarded] = coef_jj * pol_current[m_retarded] + coef_je * (ez_temp[m_retarded] + ez[m])

@njit
def update_disp_elec_PLRC(start, end, hy, ez, ez_temp, rec_accumulator, cez_ez, cez_hy, cez_accum, caccum_ezf, caccum_ezp, caccum_accum):
        for m in range(start, end):
            m_retarded = m - start
            ez_temp[m_retarded] = ez[m]
            ez[m] = cez_ez * ez[m] + cez_hy * (hy[m] - hy[m - 1]) + cez_accum * rec_accumulator[m_retarded]
            rec_accumulator[m_retarded] = caccum_ezf * ez[m] + caccum_ezp * ez_temp[m_retarded] + caccum_accum * rec_accumulator[m_retarded]

@njit
def update_disp_elec_ztransf(start, end, hy, ez, d_field, integrator, low_pass, cezd, cezi, cezl, clows, imp0, courant):
        for m in range(start, end):
            m_retarded = m - start
            ez[m] = ez[m] * (1/imp0)
            d_field[m_retarded] = d_field[m_retarded] + courant * (hy[m] - hy[m - 1])
            ez[m] = cezd * d_field[m_retarded] - cezi * integrator[m_retarded] - cezl * low_pass[m_retarded]
            integrator[m_retarded] = integrator[m_retarded] + ez[m]
            low_pass[m_retarded] = clows * low_pass[m_retarded] + ez[m]
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
