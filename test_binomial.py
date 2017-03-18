import numpy as np
from cpd.binomial import bpmfln
import numba
import math


@numba.jit(nopython=True)
def combined(n, k):
    return math.lgamma(n + 1) - (math.lgamma(k + 1) +
                                 math.lgamma(n - k + 1))


@numba.jit(nopython=True)
def binomial(k, n, p):
    bnm = np.empty_like(n)
    logp = math.log(p)
    one_logp = math.log(1 - p)
    for i in range(len(k)):
        bnm[i] = math.exp(combined(n[i], k[i]) + k[i] *
                          logp + (n[i] - k[i]) * one_logp)
    return bnm


n_bridges = np.array([5.01642226,   9.03284453,  13.04926679,  17.06568906,
                      21.08211132,  25.09853359,  29.11495585,  33.13137812,
                      37.14780038,  41.16422265,  45.18064491,  49.19706718,
                      53.21348944,  57.22991171,  61.24633397,  65.26275624,
                      69.2791785,  73.29560077,  77.31202303,  81.3284453])

p = 0.53956619110916559

k = np.arange(20)


# print(bpmfln(k, n_bridges, p))
# print(binomial(k=k, n=n_bridges, p=p))

np.testing.assert_allclose(bpmfln(k, n_bridges, p),
                           binomial(k=k, n=n_bridges, p=p))
