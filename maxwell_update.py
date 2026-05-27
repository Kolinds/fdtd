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


@njit
def update_disp_mag_ADE(start, end, hy, ez, hy_temp, pol_current, coef_jj, coef_jh, cphy_h, cphy_fp, cphy_dp):
        for m in range(start, end):
            m_retarded = m - start
            hy_temp[m_retarded] = hy[m]
            hy[m] = cphy_h * hy[m] + cphy_fp * ((ez[m + 1] - ez[m]) - cphy_dp * pol_current[m_retarded])
            pol_current[m_retarded] = coef_jj * pol_current[m_retarded] + coef_jh * (hy_temp[m_retarded] + hy[m])

@njit
def update_imp_ez_ADE(start, end, hy, ez, ez_old, coef_d, c_a, c_b, d_a, coef_jj, coef_je, calc_denom, coef_a, coef_c, coef_b, pol_current, delta_x):
        width = end - start
        
        # 1. Nodo de frontera izquierda (m = start)
        # Corregido: hy[start] - hy[start - 1]
        coef_d[0] = ((c_b / 2) * (1 + d_a) * hy[start] 
                     - (c_b / 2) * (1 + d_a) * hy[start - 1] 
                     + (c_a + coef_c[0]) * ez_old[start] 
                     - coef_c[0] * ez_old[start + 1] 
                     - c_b * ((1 + coef_jj) / 2) * pol_current[0])
        coef_d[0] = coef_d[0] / calc_denom[0]

        # 2. Nodos internos de la matriz 
        for m in range(start + 1, end):
            m_retarded = m - start
            # Corregido: hy[m] - hy[m - 1]
            rhs = ((c_b / 2) * (1 + d_a) * hy[m] 
                   - (c_b / 2) * (1 + d_a) * hy[m - 1] 
                   - coef_a[m_retarded] * ez_old[m - 1] 
                   + (c_a + coef_a[m_retarded] + coef_c[m_retarded]) * ez_old[m] 
                   - coef_c[m_retarded] * ez_old[m + 1] 
                   - c_b * ((1 + coef_jj) / 2) * pol_current[m_retarded])
            coef_d[m_retarded] = (rhs - coef_a[m_retarded] * coef_d[m_retarded - 1]) / calc_denom[m_retarded]

        # 3. Nodo de frontera derecha (m = end)
        # Corregido: hy[end] - hy[end - 1]
        coef_d[width] = ((c_b / 2) * (1 + d_a) * hy[end] 
                         - (c_b / 2) * (1 + d_a) * hy[end - 1] 
                         - coef_a[width] * ez_old[end - 1] 
                         + (c_a + coef_a[width]) * ez_old[end] 
                         - c_b * ((1 + coef_jj) / 2) * pol_current[width])
        coef_d[width] = (coef_d[width] - coef_a[width] * coef_d[width - 1]) / calc_denom[width]

        # 4. Sustitución hacia atrás (Thomas)
        ez[end] = coef_d[width]
        for m in range(end - 1, start - 1, -1):
            m_retarded = m - start
            ez[m] = coef_d[m_retarded] - coef_c[m_retarded] * ez[m + 1]
        
        # 5. Actualizar la corriente de polarización dispersiva
        for m in range(start, end + 1):
            m_retarded = m - start
            pol_current[m_retarded] = coef_jj * pol_current[m_retarded] + coef_je * (ez[m] + ez_old[m])

@njit
def update_imp_hy(start, end, hy, ez, ez_old, chyh, chye):
        # Corregido: Signo de la derivada espacial (ez[m+1] - ez[m]) e índices mapeados al espacio real
        for m in range(start, end):
            hy[m] = chyh * hy[m] + chye * (ez[m + 1] - ez[m] + ez_old[m + 1] - ez_old[m])   








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
