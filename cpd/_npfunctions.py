import numpy as np


def sum_x_n_calc(x, z_n, k_n_1):
    """Eq 54"""
    return np.sum(z_n / (1 + k_n_1 * x) * k_n_1)
