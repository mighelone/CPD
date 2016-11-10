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

Rgas = 1.987  # cal/mol-K


@logged
class CPD(pkp.cpd.CPD):

    def dydt(self, t, y):
        l, delta, c = y
        p = l + c
        T = self.T(t)
        kb, rho, kg = self.rates(T)
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

    def rates(self, T):
        '''
        Calculate rates for the given temperature
        '''
        # TODO update to account for sigma
        RT = T * Rgas
        kb = self.ab * np.exp(-self.eb / RT)
        kc = self.ac * np.exp(-self.ec / RT)
        kg = self.ag * np.exp(-self.eg / RT)
        return kb, kc, kg
