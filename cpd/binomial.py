'''
implementing binomial pdf with gammaln

see: Robert Kerns discussion in http://groups.google.ca/group/comp.lang.python/browse_thread/thread/839009574397dc37

'''
import numpy as np
from scipy.misc import comb
from scipy import special


def bpmf(k, n, p):
    '''
    Binomial distribution using comb function in scipy

    Parameters
    ----------

    '''
    # this does not work for large n
    return comb(n, k) * (p**k) * ((1 - p)**(n - k))


# proposed version using gammaln
def combinln(n, k):
    return (special.gammaln(n + 1) - (special.gammaln(k + 1) +
                                      special.gammaln(n - k + 1)))


def bpmfln(k, n, p):
    return np.exp(combinln(n, k) + k * np.log(p) + (n - k) * np.log(1 - p))


# n = 10
# p = 1.0e-5

# print("using gammaln")
# print(bpmfln(np.arange(11), 10, p))
# print("using comb")
# print(bpmf(np.arange(11), 10, p))
# print("difference")
# print(bpmfln(np.arange(11), 10, p) - bpmf(np.arange(11), 10, p))


# pmfnln = bpmfln(np.arange(5001), 5000, 0.99)
# print('n = 5000')
# print('nans', np.sum(np.isnan(pmfnln)))
# print('sum', np.sum(pmfnln))
# print('sum (repr)', repr(np.sum(pmfnln)))
