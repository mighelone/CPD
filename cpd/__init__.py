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
from scipy.stats import norm, binom
from scipy.optimize import brentq

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

    def _dydt(self, t, y):
        l, delta, c = y
        T = self.T(t)
        kb, rho, kg = self._rates(T, y)
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
        solver = ode(self._dydt).set_integrator(backend, nsteps=1,
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
        f_solid, f_gas, f_tar, f_meta, f_cross = 1, 0, 0, 0, 0
        f = [[f_solid,
              f_gas,
              f_tar,
              f_meta,
              f_cross]]
        while solver.t < time_end:
            # print(solver.t)
            # self.__log.debug('t=%s', solver.t)
            # print(solver.t, solver.y, self._dydt(
            #    solver.t, solver.y), self.T(solver.t))
            solver.integrate(time_end, step=True)
            dt = solver.t - t[1]
            T = self.T(solver.t)
            t.append(solver.t)
            y.append(solver.y)

            if f_meta > 1e-4:
                rate_cross = self._crosslinking(f_meta, T, dt)
                fract = 1 - rate_cross / f_meta
                f_meta -= rate_cross
                f_cross += rate_cross
            else:
                rate_cross = 0
                fract = 1
            percolation = self._percolation(y, f_tar, in_tar=True)

        t = np.array(t)
        y = np.array(y)

        return t, y

    def _rates(self, T, y):
        '''
        Calculate _rates for the given temperature
        '''
        l, delta, c = y
        g1, g2 = self.gas(y)
        g = g1 + g2
        RT = T * Rgas
        eb = self.eb + invernorm(1 - l / (self.p0 - self.c0))
        self.__log.debug('Eb %s Eb0 %s', eb, self.eb)
        kb = self.ab * np.exp(-eb / RT)
        kc = self.ac * np.exp(-self.ec / RT)
        eg = self.eg + invernorm(0.5 * g / (1 - self.c0)) * self.egsig
        kg = self.ag * np.exp(-eg / RT)
        return kb, kc, kg

    def _percolation(self, y, f_tar=0, in_tar=True):
        def pstar_eq(pstar):
            '''Eq. 6 in CPD summary'''
            def f(p0):
                return p0 * (1 - p0) ** (self.sig - 1)
            return f(pstar) - f(p)
        l, delta, c = y
        p = l + c
        f = 1 - p
        g1, g2 = self.gas(y)
        g = g1 + g2
        # TODO this values are global
        sigma = self.sig - 1
        # average mass of the fused ring site
        ma = self.mw - self.sig * self.mdel
        # mass of bridges
        mb = 2 * self.mdel
        rba = mb / ma

        if in_tar:
            # Eq. 36, 37
            delta_fac = delta / (1 - p) if p < 0.9999 else 1
            a = 1 + rba * (l / p + (sigma - 1) * 0.25 * delta_fac)
            b = 0.5 * delta_fac - l / p
            # p_threasold is the minimum value of p for which there are
            # nomore infinite fragments (only finite fragments)
            p_threasold = 1. / sigma + 1e-4
            if p > 0.999:
                pstar = 1
            elif p > p_threasold:
                pstar = brentq(pstar_eq, p_threasold, 1)
            else:
                pstar = p
            sfac = self.sig / (sigma - 1)
            # Eq. (5) fraction of finite fragments
            Fp = (pstar / p) ** sfac
            Kp = Fp * (1 - self.sig * 0.5 * pstar)
            # Eq. (39) mass fraction of finite fragments
            f_frag = 2 * (a * Fp + b * Kp) / \
                (2 + rba * (1 - self.c0) * self.sig)

        # mass of gas released at time t Eq (31)
        # gas is produced only considering the remaining fragments in
        # the metaplast
        mgas = rba * ma * g * self.sig * 0.25 * (1 - f_tar)
        mtot = ma + rba * ma * self.sig * 0.5 * (1 - self.c0)
        f_gas = mgas / mtot
        f_char = 1 - f_tar - f_gas

        #
        n = np.arange(1, 20)  # number of clusters in a fragment
        # broken bridges per cluster of size n
        tau = n * (sigma - 1) + 2
        s = n - 1  # intact bridges per cluster of size n
        n_bridges = tau + s
        # Eq. 32 mass of a finite fragment of size n
        m_frag_n = (n * ma + (n - 1) * mb * l / p +
                    tau * mb * delta_fac * 0.25)
        # Eqs (1-4)
        Qn = self.sig / n_bridges / n * binom.pmf(s, n_bridges, p)
        # Eq. (33) total mass of fragments of size
        f_frag_n = m_frag_n * Qn
        f_frag_n /= f_frag_n.sum() * f_frag
        return {'f_gas': f_gas,
                'f_tar': f_tar,
                'f_char': f_char,
                'f_frag': f_frag,
                'f_frag_n': f_frag_n,
                'm_frag_n': m_frag_n,
                'pstar': pstar}

    def _tar_distribution(self):
        pass

    def intact_bridges(self, y):
        l, _, c = y
        return l + c

    def broken_bridges(self, y):
        return 1 - self.intact_bridges(y)

    def gas(self, y):
        f = self.broken_bridges(y)
        _, delta, c = y
        g1 = 2 * f - delta
        g2 = 2 * (c - self.c0)
        return g1, g2

    def _crosslinking(self, f_meta, T, dt):
        return self.Acr * np.exp(-self.Ecr / Rgas / T) * f_meta * dt

    def flash_distillation(self, df_gas, mw_gas, df_n, meta_n, mw_n, fracr, T):

        def funct(x):
            '''Eq 54'''
            return np.sum(x_n_calc(x) * (k_n - 1))

        def x_n_calc(x):
            '''
            Eq. 52
            '''
            return z_n / (1 + (k_n - 1) * x)

        a = 87058.0
        b = 299.0
        g = 0.5903
        # mole fraction of n-mers contained in the metaplast
        self.__log.debug('Tot mass of fragments %s',
                         (df_n + meta_n * fracr).sum())
        F_n = np.append((df_n + meta_n * fracr) / mw_n, df_gas / mw_gas)
        F = F_n.sum()
        mw = np.append(mw_n, mw_gas)
        self.__log.debug('MW %s', mw)
        p_vap = a * np.exp(-b * mw ** g / T)
        self.__log.debug('p_vap %s', p_vap)
        k_n = p_vap * 101325 / self.pressure
        self.__log.debug('kn %s', k_n)
        z_n = F_n / F
        self.__log.debug('zn %s', z_n)
        fract_v = brentq(funct, 0, 0.999)
        self.__log.debug('F/V = %s', fract_v)
        V = fract_v * F  # moles of tar
        L = F - V
        # mole fraction of n-mers in the metaplast
        x_n = x_n_calc(fract_v)
        # mole fraction of n-mers released as tar
        y_n = k_n * x_n
        self.__log.debug('metaplast %s', x_n)
        self.__log.debug('tar %s', y_n)

        meta_n_new = x_n * L * mw
        tar_n_new = y_n * V * mw

        self.__log.debug('Mass metaplast %s', meta_n_new.sum())
        self.__log.debug('Mass tar %s', tar_n_new.sum())

        return tar_n_new[:-1], meta_n_new[:-1]
