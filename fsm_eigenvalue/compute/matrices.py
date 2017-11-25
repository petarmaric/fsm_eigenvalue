from beam_integrals.characteristic_equation_solvers import find_best_root
import numpy as np

from .integral_db import get_scaled_integral
from .utils import assemble_local_matrix


def get_stiffness_matrix(I1, I2, I5, I6, I7, I8, I21, I22, I23, I24, I25, a_mu, b, t, K_x, K_y, K_1, K_xy):
    # As per eq. 4.28,4.29 from [Milasinovic1997]
    K_uu11 = ( K_x *I1/b  + K_xy*I2*b/3.) * t
    K_uu12 = (-K_x *I1/b  + K_xy*I2*b/6.) * t
    K_uu13 = (-K_1 *I5/2. - K_xy*I6/2.  ) * t * a_mu
    K_uu14 = (-K_1 *I5/2. + K_xy*I6/2.  ) * t * a_mu
    K_uu33 = ( K_xy*I8/b  + K_y *I7*b/3.) * t * a_mu**2
    K_uu34 = (-K_xy*I8/b  + K_y *I7*b/6.) * t * a_mu**2
    K_uu22 =   K_uu11
    K_uu23 =  -K_uu14
    K_uu24 =  -K_uu13
    K_uu44 =   K_uu33

    # As per eq. 4.18,4.19 from [Milasinovic1997]
    t_3  = t**3 / 12.
    D_11 = K_x  * t_3
    D_22 = K_y  * t_3
    D_12 = K_1  * t_3
    D_66 = K_xy * t_3
    K_ww11 =  12.*D_11*I21/b**3 -  6./5. *D_12*I22/b - 6./5. *D_12*I23/b + 13./35. *D_22*I24*b    + 24./5. *D_66*I25/b
    K_ww12 =   6.*D_11*I21/b**2 - 11./10.*D_12*I22   - 1./10.*D_12*I23   + 11./210.*D_22*I24*b**2 +  2./5. *D_66*I25
    K_ww13 = -12.*D_11*I21/b**3 +  6./5. *D_12*I22/b + 6./5. *D_12*I23/b + 18./140.*D_22*I24*b    - 24./5. *D_66*I25/b
    K_ww14 =   6.*D_11*I21/b**2 -  1./10.*D_12*I23   - 1./10.*D_12*I22   - 26./840.*D_22*I24*b**2 +  2./5. *D_66*I25
    K_ww22 =   4.*D_11*I21/b    -  2./15.*D_12*I22*b - 2./15.*D_12*I23*b +  2./210.*D_22*I24*b**3 +  8./15.*D_66*I25*b
    K_ww24 =   2.*D_11*I21/b    +  2./60.*D_12*I22*b + 2./60.*D_12*I23*b -  6./840.*D_22*I24*b**3 -  2./15.*D_66*I25*b
    K_ww34 =  -6.*D_11*I21/b**2 + 22./20.*D_12*I22   + 6./60.*D_12*I23   - 22./420.*D_22*I24*b**2 -  2./5. *D_66*I25
    K_ww23 = -K_ww14
    K_ww33 =  K_ww11
    K_ww44 =  K_ww22

    return assemble_local_matrix(
        X_uu=[K_uu11, K_uu12, K_uu13, K_uu14, K_uu22, K_uu23, K_uu24, K_uu33, K_uu34, K_uu44],
        X_ww=[K_ww11, K_ww12, K_ww13, K_ww14, K_ww22, K_ww23, K_ww24, K_ww33, K_ww34, K_ww44],
    )

def get_stress_matrix(I2, I7, I25, b, c):
    # As per eq. 6.78-6.80 from [Milasinovic1997]
    K_uu11 = (3. +    c)/24.*I2*b
    K_uu12 = (1. +    c)/24.*I2*b
    K_uu13 =  0.
    K_uu14 =  0.
    K_uu22 = (1. + 3.*c)/24.*I2*b
    K_uu23 =  0.
    K_uu24 =  0.
    K_uu33 = (3. +    c)/24.*I7*b
    K_uu34 = (1. +    c)/24.*I7*b
    K_uu44 = (1. + 3.*c)/24.*I7*b

    # As per eq. 6.74-6.77 from [Milasinovic1997]
    K_ww11 = (10. +  3.*c)/70.  *I25*b
    K_ww12 = (15. +  7.*c)/840. *I25*b**2
    K_ww13 = ( 9. +  9.*c)/280. *I25*b
    K_ww14 = (-7. -  6.*c)/840. *I25*b**2
    K_ww22 = ( 5. +  3.*c)/1680.*I25*b**3
    K_ww23 = ( 6. +  7.*c)/840. *I25*b**2
    K_ww24 = (-1. -     c)/560. *I25*b**3
    K_ww33 = ( 3. + 10.*c)/70.  *I25*b
    K_ww34 = (-7. - 15.*c)/840. *I25*b**2
    K_ww44 = ( 3. +  5.*c)/1680.*I25*b**3

    return assemble_local_matrix(
        X_uu=[K_uu11, K_uu12, K_uu13, K_uu14, K_uu22, K_uu23, K_uu24, K_uu33, K_uu34, K_uu44],
        X_ww=[K_ww11, K_ww12, K_ww13, K_ww14, K_ww22, K_ww23, K_ww24, K_ww33, K_ww34, K_ww44],
    )

def get_mass_matrix(I1, I8, I21, b, t, ro):
    # As per eq. 6.36 from [Milasinovic1997], will multiply by ``t * ro`` below
    M_uu11 = I1*b/3.
    M_uu12 = M_uu11/2.
    M_uu13 = 0.
    M_uu14 = 0.
    M_uu22 = M_uu11
    M_uu23 = 0.
    M_uu24 = 0.
    M_uu33 = I8*b/3.
    M_uu34 = M_uu33/2.
    M_uu44 = M_uu33

    # As per eq. 6.31 from [Milasinovic1997], will multiply by ``t * ro * I21`` below
    M_ww11 =  13./35. *b
    M_ww12 =  11./210.*b**2
    M_ww13 =   9./70. *b
    M_ww14 = -13./420.*b**2
    M_ww22 =   1./105.*b**3
    M_ww23 = -M_ww14
    M_ww24 = - 3./420.*b**3
    M_ww33 =  M_ww11
    M_ww34 = -M_ww12
    M_ww44 =  M_ww22

    return assemble_local_matrix(
        X_uu=np.array([M_uu11, M_uu12, M_uu13, M_uu14, M_uu22, M_uu23, M_uu24, M_uu33, M_uu34, M_uu44]) * t * ro,
        X_ww=np.array([M_ww11, M_ww12, M_ww13, M_ww14, M_ww22, M_ww23, M_ww24, M_ww33, M_ww34, M_ww44]) * t * ro * I21,
    )

def compute_global_matrices(integral_db, beam_type, strip_data, materials, astiff_shape, a, t_b, m):
    def get_integral(integral_id):
        return get_scaled_integral(integral_db, integral_id, a, m=m, n=m)

    I1 = I21 = get_integral(1)
    I2 = I6  = I8 = I25 = get_integral(2)
    I3 = I22 = get_integral(3)
    I5 = I23 = get_integral(5)
    I7 = I24 = get_integral(7)

    mu_m = float(find_best_root(beam_type, mode=m))
    a_mu = a / mu_m

    K_hat   = np.asmatrix(np.zeros(astiff_shape)) # global stiffness matrix
    K_sigma = np.asmatrix(np.zeros(astiff_shape)) # global stress matrix
    M       = np.asmatrix(np.zeros(astiff_shape)) # global mass matrix

    edge_data_keys = 'material_id, b, R, astiff_fill_indices'.split(', ')
    material_data_keys = 't_s, ro, c, K_x, K_y, K_1, K_xy'.split(', ')
    for _, _, edge_data in strip_data:
        material_id, b, R, astiff_fill_indices = (edge_data[k] for k in edge_data_keys)

        material = materials[material_id]
        t_s, ro, c, K_x, K_y, K_1, K_xy = (material[k] for k in material_data_keys)

        t = t_b * t_s # [mm] strip thickness

        K_hat_strip = get_stiffness_matrix(
            I1, I2, I5, I6, I7, I8, I21, I22, I23, I24, I25,
            a_mu, b, t, K_x, K_y, K_1, K_xy
        )
        K_sigma_strip = get_stress_matrix(I2, I7, I25, b, c)
        M_strip = get_mass_matrix(I1, I8, I21, b, t, ro)

        # As per eq. 3.62 from [Milasinovic1997]
        K_hat_strip   = R.T * K_hat_strip   * R
        K_sigma_strip = R.T * K_sigma_strip * R
        M_strip       = R.T * M_strip       * R

        # Deduced from Fortran block 83:94
        for astiff_indices, segment_indices in astiff_fill_indices:
            K_hat  [astiff_indices] += K_hat_strip  [segment_indices]
            K_sigma[astiff_indices] += K_sigma_strip[segment_indices]
            M      [astiff_indices] += M_strip      [segment_indices]

    return K_hat, K_sigma, M
