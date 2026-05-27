        coef_d[last_idx] = ((c_b / 2) * (1 + d_a) * hy[end - 1] 
                         - (c_b / 2) * (1 + d_a) * hy[end - 2] 
                         - alpha_b * ez[end]
                         - alpha_b * ez_old[end]
                         - alpha_b * ez_old[end - 2] 
                         + (c_a + 2 * alpha_b) * ez_old[end - 1] 
                         - c_b * ((1 + coef_jj) / 2) * pol_current[last_idx])
        coef_d[last_idx] = (coef_d[last_idx] - coef_a[last_idx] * coef_d[last_idx - 1]) / calc_denom[last_idx]
