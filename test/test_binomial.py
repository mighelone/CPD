import numpy as np
from scipy import special

from cpd.binomial import bpmfln, bpmf, combinln


def combnp(n, k):
    return (special.gammaln(n + 1) - (special.gammaln(k + 1) +
                                      special.gammaln(n - k + 1)))


p = 0.63
nb = np.array([5,  9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65,
               69, 73, 77, 81])

ns = np.arange(20)
sigma = 5


def test_bnm():
    np.testing.assert_allclose(bpmfln(ns, nb, p),
                               bpmf(ns, nb, p))


def test_combin():
    c = np.array([combinln(nbi, p) for nbi in nb])
    np.testing.assert_allclose(c, combnp(nb, p))
