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
        # El último índice válido del array implícito no es 'width', es 'last_idx'
        last_idx = end - start - 1 
        
        # Factor de acoplamiento de frontera
        alpha_b = coef_c[1] 

        # 1. Nodo de frontera izquierda (m = start)
        # Acoplado con ez[start - 1] (Explícito)
        coef_d[0] = ((c_b / 2) * (1 + d_a) * hy[start] 
                     - (c_b / 2) * (1 + d_a) * hy[start - 1] 
                     - alpha_b * ez[start - 1] 
                     - alpha_b * ez_old[start - 1]
                     + (c_a + 2 * alpha_b) * ez_old[start] 
                     - alpha_b * ez_old[start + 1] 
                     - c_b * ((1 + coef_jj) / 2) * pol_current[0])
        coef_d[0] = coef_d[0] / calc_denom[0]

        # 2. Nodos internos (Desde start + 1 HASTA end - 2)
        for m in range(start + 1, end - 1):
            m_retarded = m - start
            rhs = ((c_b / 2) * (1 + d_a) * hy[m] 
                   - (c_b / 2) * (1 + d_a) * hy[m - 1] 
                   - coef_a[m_retarded] * ez_old[m - 1] 
                   + (c_a + coef_a[m_retarded] + coef_c[m_retarded]) * ez_old[m] 
                   - coef_c[m_retarded] * ez_old[m + 1] 
                   - c_b * ((1 + coef_jj) / 2) * pol_current[m_retarded])
            coef_d[m_retarded] = (rhs - coef_a[m_retarded] * coef_d[m_retarded - 1]) / calc_denom[m_retarded]

        # 3. Nodo de frontera derecha (m = end - 1)
        # CORREGIDO: Acoplado dinámicamente con ez[end] (Explícito)
        # Sus vecinos magnéticos son hy[end - 1] y hy[end - 2] (Ambos implícitos)
        coef_d[last_idx] = ((c_b / 2) * (1 + d_a) * hy[end - 1] 
                         - (c_b / 2) * (1 + d_a) * hy[end - 2] 
                         - alpha_b * ez[end]
                         - alpha_b * ez_old[end]
                         - alpha_b * ez_old[end - 2] 
                         + (c_a + 2 * alpha_b) * ez_old[end - 1] 
                         - c_b * ((1 + coef_jj) / 2) * pol_current[last_idx])
        coef_d[last_idx] = (coef_d[last_idx] - coef_a[last_idx] * coef_d[last_idx - 1]) / calc_denom[last_idx]

        # 4. Sustitución hacia atrás (Thomas) 
        # Empezamos a despejar desde end - 1 hacia atrás
        ez[end - 1] = coef_d[last_idx]
        for m in range(end - 2, start - 1, -1):
            m_retarded = m - start
            ez[m] = coef_d[m_retarded] - coef_c[m_retarded] * ez[m + 1]
        
        # 5. Actualizar la corriente de polarización dispersiva
        # El bloque de plasma tiene exactamente un tamaño de (end - start) nodos
        for m in range(start, end):
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
