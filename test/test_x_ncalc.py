import numpy as np

import cpd._nb_functions as nbf
import cpd._npfunctions as npf
from scipy.optimize import brentq


z_n = np.array([8.56188218e-01,   1.11672986e-01,   2.33239597e-02,
                6.09094693e-03,   1.81353598e-03,   5.88702032e-04,
                2.03177693e-04,   7.34035256e-05,   2.74758011e-05,
                1.05795032e-05,   4.16877067e-06,   1.67453239e-06,
                6.83642360e-07,   2.83009103e-07,   1.18576375e-07,
                5.02071201e-08,   2.14566276e-08,   9.24563629e-09,
                4.01340399e-09,   1.75375231e-09,   0.00000000e+00])

k_n_1 = np.array([6.68421120e+00,  -9.29405211e-01,  -9.98389058e-01,
                  -9.99939876e-01,  -9.99996899e-01,  -9.99999798e-01,
                  -9.99999984e-01,  -9.99999999e-01,  -1.00000000e+00,
                  -1.00000000e+00,  -1.00000000e+00,  -1.00000000e+00,
                  -1.00000000e+00,  -1.00000000e+00,  -1.00000000e+00,
                  -1.00000000e+00,  -1.00000000e+00,  -1.00000000e+00,
                  -1.00000000e+00,  -1.00000000e+00,   9.48172221e+03])


def test_x_n_calc():
    nbf.x_n_calc(0.5, z_n, k_n_1)


def test_results():
    """Test if the results are the same"""
    x = 0.5
    np.testing.assert_almost_equal(nbf.sum_x_n_calc(
        x, z_n, k_n_1), npf.sum_x_n_calc(x, z_n, k_n_1))


def test_brentq():
    np.testing.assert_almost_equal(
        brentq(npf.sum_x_n_calc, 0, 0.999, args=(z_n, k_n_1)),
        brentq(nbf.sum_x_n_calc, 0, 0.999, args=(z_n, k_n_1)))
