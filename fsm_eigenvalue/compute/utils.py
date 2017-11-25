import numpy as np


MIN_EIGENVALUE = 10**-12


def symmetrize_matrix(mat):
    return mat + mat.T - np.diag(np.diag(mat))

def assemble_local_matrix(X_uu, X_ww):
    X_uu11, X_uu12, X_uu13, X_uu14, X_uu22, X_uu23, X_uu24, X_uu33, X_uu34, X_uu44 = X_uu
    X_ww11, X_ww12, X_ww13, X_ww14, X_ww22, X_ww23, X_ww24, X_ww33, X_ww34, X_ww44 = X_ww

    X = np.asmatrix([
        #     0       1       2       3       4       5       6       7
        [X_uu11, X_uu13,      0,      0, X_uu12, X_uu14,      0,      0], # 0
        [     0, X_uu33,      0,      0, X_uu23, X_uu34,      0,      0], # 1
        [     0,      0, X_ww11, X_ww12,      0,      0, X_ww13, X_ww14], # 2
        [     0,      0,      0, X_ww22,      0,      0, X_ww23, X_ww24], # 3
        [     0,      0,      0,      0, X_uu22, X_uu24,      0,      0], # 4
        [     0,      0,      0,      0,      0, X_uu44,      0,      0], # 5
        [     0,      0,      0,      0,      0,      0, X_ww33, X_ww34], # 6
        [     0,      0,      0,      0,      0,      0,      0, X_ww44], # 7
    ])

    return symmetrize_matrix(X)

def clip_small_eigenvalues(eigenvalues):
    eigenvalues[np.abs(eigenvalues)<=MIN_EIGENVALUE] = MIN_EIGENVALUE

def get_relative_error(v, v_approx):
    return np.abs(1 - v_approx / v)
