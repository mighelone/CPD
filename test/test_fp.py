import cpd._np_functions as npf
import cpd._nb_functions as nbf
import numpy as np
from scipy.optimize import newton

sigma = 4.0164
p_threasold = 0.2489778051353545
p = 0.5395


def test_fpp():
    np.testing.assert_almost_equal(npf.fp(p, sigma), nbf.fp(p, sigma))


def test_pstar():
    fpp = npf.fp(p, sigma)
    np.testing.assert_almost_equal(
        nbf.pstar_f(p_threasold * 0.5, sigma, fpp),
        npf.pstar_f(p_threasold * 0.5, sigma, fpp))


def test_pstar():
    fpp = npf.fp(p, sigma)
    np.testing.assert_almost_equal(
        nbf.pstar_f(p_threasold * 0.5, sigma, fpp),
        npf.pstar_f(p_threasold * 0.5, sigma, fpp))


def test_newton():
    fpp = npf.fp(p, sigma)
    pstar_np = newton(npf.pstar_f, p_threasold * 0.5,
                      args=(sigma, fpp))
    pstar_nb = newton(nbf.pstar_f, p_threasold * 0.5,
                      args=(sigma, fpp))

    np.testing.assert_almost_equal(pstar_nb, pstar_np)
