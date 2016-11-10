'''
CPD python implementation

Michele Vascellari
'''
from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals


import pkp
import pkp.cpd
import numpy as np
from autologging import logged
from scipy.integrate import ode
from scipy.stats import norm

Rgas = 1.987  # cal/mol-K


def invernorm(y):
    y_min = 0.0003
    y_max = 0.9997
    if y < y_min:
        y = y_min
    elif y > y_max:
        y = y_max
    return norm.ppf(y)


@logged
class CPD(pkp.cpd.CPD):

    def dydt(self, t, y):
        l, delta, c = y
        p = l + c
        T = self.T(t)
        kb, rho, kg = self.rates(T, y)
        dldt = -kb * l
        f = 1 / (1 + rho)
        ddeldt = 2 * rho * kb * f * l - kg * delta
        dcdt = kb * f * l
        return [dldt, ddeldt, dcdt]

    def _bridge_evolution(self):
        '''
        '''

        # variables are [l, d, c]
        backend = 'dopri5'
        t0 = self.operating_conditions[0, 0]
        solver = ode(self.dydt).set_integrator(backend, nsteps=1,
                                               first_step=1e-6,
                                               max_step=1e-4,
                                               verbosity=1)
        solver._integrator.iwork[2] = -1
        y0 = [self.p0 - self.c0,
              2 * (1 - self.p0),
              self.c0]
        solver.set_initial_value(y0, t0)
        solver._integrator.iwork[2] = -1
        time_end = self.operating_conditions[-1, 0]

        t = [t0]
        y = [y0]
        while solver.t < time_end:
            # print(solver.t)
            # self.__log.debug('t=%s', solver.t)
            # print(solver.t, solver.y, self.dydt(
            #    solver.t, solver.y), self.T(solver.t))
            solver.integrate(time_end, step=True)
            t.append(solver.t)
            y.append(solver.y)

        t = np.array(t)
        y = np.array(y)

        return t, y

    def rates(self, T, y):
        '''
        Calculate rates for the given temperature
        '''
        # TODO update to account for sigma
        l, delta, c = y
        p = l + c
        f = 1 - p
        g1 = 2 * f - delta
        g2 = 2 * (c - self.c0)
        g = g1 + g2
        RT = T * Rgas
        eb = self.eb + invernorm(1 - l / (self.p0 - self.c0))
        self.__log.debug('Eb %s Eb0 %s', eb, self.eb)
        kb = self.ab * np.exp(-eb / RT)
        kc = self.ac * np.exp(-self.ec / RT)
        eg = self.eg + invernorm(0.5 * g / (1 - self.c0)) * self.egsig
        kg = self.ag * np.exp(-eg / RT)
        return kb, kc, kg
