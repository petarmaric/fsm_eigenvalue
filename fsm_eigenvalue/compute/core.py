import numpy as np

from .matrices import compute_global_matrices
from .utils import clip_small_eigenvalues, get_relative_error


def approximate_natural_frequency_from_stress(m, a, sigma, ro):
    return m * np.pi / a * np.sqrt(sigma / ro)

def approximate_stress_from_natural_frequency(m, a, omega, ro):
    return (omega * a / (m * np.pi))**2 * ro

def solve_eigenvalue_problem(inv_G, A, normalize_eigenvalues=None):
    # As per eq. 6.48 from [Milasinovic1997]
    H = inv_G * A * inv_G.T
    eigenvalues, eigenvectors = np.linalg.eigh(H)

    # Clip the extremely small eigenvalues
    clip_small_eigenvalues(eigenvalues)

    if normalize_eigenvalues:
        eigenvalues = normalize_eigenvalues(eigenvalues)

    # As per eq. 6.45 from [Milasinovic1997]
    mode_shapes = inv_G.T * eigenvectors

    # According to Milasinovic the minimal eigenvalue is the one closest to 0
    min_idx = np.argmin(eigenvalues)
    eigenvalue_min = eigenvalues[min_idx]
    mode_shape_min = mode_shapes[:,min_idx].A1

    return eigenvalue_min, mode_shape_min

def perform_iteration(integral_db, beam_type, strip_data, materials, astiff_shape, a, t_b, m):
    K_hat, K_sigma, M = compute_global_matrices(
        integral_db, beam_type, strip_data, materials, astiff_shape, a, t_b, m
    )

    # As per eq. 6.40,6.41 from [Milasinovic1997]
    # ``G`` is the lower triangle matrix factorized from ``K_hat = G * G.T``
    inv_G = np.linalg.cholesky(K_hat).I

    # As per eq. 6.22,6.39,6.48 from [Milasinovic1997]
    # ``omega`` [rad/s] is the natural frequency, and ``Phi_omega`` is its mode shape
    omega, Phi_omega = solve_eigenvalue_problem(
        inv_G, M, normalize_eigenvalues=lambda x: np.sqrt(1./x)
    )

    # As per eq. 6.48,6.63,6.82 from [Milasinovic1997]
    # ``sigma_cr`` [MPa] is the critical buckling stress, and ``Phi_sigma_cr`` is its mode shape
    N_cr, Phi_sigma_cr = solve_eigenvalue_problem(
        inv_G, K_sigma, normalize_eigenvalues=lambda x: 1./x
    )
    sigma_cr = N_cr / (2*t_b)

    ro = float(np.mean([mat['ro'] for mat in materials.values()]))

    omega_approx = approximate_natural_frequency_from_stress(m, a, sigma_cr, ro)
    omega_rel_err = get_relative_error(omega, omega_approx)

    sigma_cr_approx = approximate_stress_from_natural_frequency(m, a, omega, ro)
    sigma_cr_rel_err = get_relative_error(sigma_cr, sigma_cr_approx)

    Phi_rel_err = get_relative_error(Phi_omega, Phi_sigma_cr)

    return (
        omega, omega_approx, omega_rel_err,
        sigma_cr, sigma_cr_approx, sigma_cr_rel_err,
        Phi_omega, Phi_sigma_cr, Phi_rel_err,
    )
