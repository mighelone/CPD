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
from scipy.optimize import brentq, newton

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

    def _set_NMR_parameters(self, nmr_parameters=None):
        super(CPD, self)._set_NMR_parameters(nmr_parameters)
        self.sigma = self.sig - 1
        # average mass of the fused ring site
        self.ma = self.mw - self.sig * self.mdel
        # mass of bridges
        self.mb = 2 * self.mdel
        self.rba = self.mb / self.ma

    def _dydt(self, t, y):
        l, delta, c = y
        T = self.T(t)
        kb, rho, kg = self._rates(T, y)
        dldt = -kb * l
        f = 1 / (1 + rho)
        ddeldt = 2 * rho * kb * f * l - kg * delta
        dcdt = kb * f * l
        return [dldt, ddeldt, dcdt]

    def _bridge_evolution(self, n_frag=20, time_end=None):
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
        if not time_end:
            time_end = self.operating_conditions[-1, 0]
        else:
            assert time_end < self.operating_conditions[-1, 0], (
                'Set a time_end < operating_conditions[-1, 0]')

        t = [t0]
        y = [y0]
        f_solid, f_gas, f_tar, f_meta, f_cross = 1, 0, 0, 0, 0
        f = [[f_solid,
              f_gas,
              f_tar,
              f_meta,
              f_cross]]
        meta_n = np.zeros(n_frag)    # init metaplast to zeros
        f_frag_n = np.zeros(n_frag)  # init fragments to zeros
        while solver.t < time_end:
            # print(solver.t)
            # self.__log.debug('t=%s', solver.t)
            # print(solver.t, solver.y, self._dydt(
            #    solver.t, solver.y), self.T(solver.t))
            solver.integrate(time_end, step=True)
            self.__log.info(
                '\n\nStart new time step\ntime=%s y=%s\n', solver.t, solver.y)
            self.__log.debug('Gas bridges=%s', self.gas(solver.y))
            dt = solver.t - t[-1]
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
            self.__log.debug(
                'Crosslinking rate: %s / %s', rate_cross, fract)
            percolation = self._percolation(solver.y, f_tar,
                                            n_frag=n_frag, in_tar=True)
            # gas formed in the last step
            df_gas = max(percolation['f_gas'] - f_gas, 0)
            self.__log.debug('gas=%s, df_gas=%s', f_gas, df_gas)
            mw_n = percolation['mw_frag_n']
            # fragments formed in the last step
            df_n = percolation['f_frag_n'] - f_frag_n
            self.__log.debug('df_n=%s', df_n)
            mw_gas = 28.0  # TODO define better this value

            tar_n, meta_n = self._flash_distillation(
                df_gas=df_gas, mw_gas=mw_gas, df_n=df_n, meta_n=meta_n,
                mw_n=mw_n, fracr=fract, T=T)

            f_frag_n = percolation['f_frag_n']
            f_gas = percolation['f_gas']
            f_tar += tar_n.sum()
            f_meta = meta_n.sum()
            f_solid = 1 - f_tar - f_gas
            f.append([f_solid, f_gas, f_tar, f_meta, f_cross])
            self.__log.debug('F=%s', f[-1])

        t = np.array(t)
        y = np.array(y)
        f = np.array(f)

        return t, y, f

    def _rates(self, T, y):
        '''
        Calculate _rates for the given temperature
        '''
        l, delta, c = y
        g1, g2 = self.gas(y)
        g = g1 + g2
        RT = T * Rgas
        eb = self.eb + invernorm(1 - l / (self.p0 - self.c0))
        # self.__log.debug('Eb %s Eb0 %s', eb, self.eb)
        kb = self.ab * np.exp(-eb / RT)
        kc = self.ac * np.exp(-self.ec / RT)
        eg = self.eg + invernorm(0.5 * g / (1 - self.c0)) * self.egsig
        kg = self.ag * np.exp(-eg / RT)
        return kb, kc, kg

    def _percolation(self, y, f_tar=0, n_frag=20, in_tar=True):
        '''
        Percolation statistic calculation

        Parameters
        ----------
        y: list
            List of bridges parameters (l, delta, c)
        f_tar: float
            Fraction of tar already released
        n_frag: int, default=20
            Number of fragments calculated
        in_tar: bool

        Returns
        -------
        percolation, dict:
            Percolation dictionary {'f_gas', 'f_solid', 'f_frag', 'f_frag_n', 'm_frag_n', 'pstar'}
        '''
        self.__log.debug('\n\nStart Percolation\n')

        l, delta, c = y
        p = l + c
        g1, g2 = self.gas(y)
        g = g1 + g2

        if in_tar:
            # Eq. 36, 37
            # Phi->a
            # Omega->b
            delta_fac = delta / (1 - p) if p < 0.9999 else 1
            self.__log.debug('delta/(1-p)=%s', delta_fac)
            a = 1 + self.rba * (l / p + (self.sigma - 1)
                                * 0.25 * delta_fac)
            b = 0.5 * delta_fac - l / p
            # p_threasold is the maxiumum value of pstar_eq
            # pstar is search from 0 to p_threasold
            p_threasold = 1. / self.sigma
            self.__log.debug('p thresold %s', p_threasold)
            if p > 0.999:
                pstar = 1
            elif p > p_threasold:
                fp = lambda x: x * (1 - x)**(self.sigma - 1)
                fpp = fp(p)
                pstar_f = lambda x: fp(x) - fpp
                # pstar = brentq(pstar_f, 0, p_threasold)
                pstar = newton(pstar_f, p_threasold * 0.5)
                self.__log.debug('Calc pstar with brentq %s', pstar)
            else:
                pstar = p
            self.__log.debug('p %s, pstar %s', p, pstar)
            sfac = self.sig / (self.sigma - 1)
            # Eq. (5) fraction of bridges in finite fragments
            Fp = (pstar / p) ** sfac
            self.__log.debug(
                'Fraction of bridges in finite fragments=%s', Fp)
            Kp = Fp * (1 - self.sig * 0.5 * pstar)
            # Eq. (39) mass fraction of finite fragments
            f_frag = 2 * (a * Fp + b * Kp) / \
                (2 + self.rba * (1 - self.c0) * self.sig)
            self.__log.debug(
                'Mass fraction of finite fragments=%s', f_frag)

        # mass of gas released at time t Eq (31)
        # gas is produced only considering the remaining fragments in
        # the metaplast
        # gas is corrected with fraction of tar already released in the
        # gas phase. f_tar is from the previous time step
        mgas = self.mb * g * self.sig * 0.25 * (1 - f_tar)
        mtot = self.ma + self.mb * self.sig * 0.5 * (1 - self.c0)
        f_gas = mgas / mtot
        self.__log.debug(
            'fraction of gas (corrected with tar released) %s (%s)', f_gas, f_tar)
        f_solid = 1 - f_tar - f_gas
        self.__log.debug(
            'fraction of remaining solid (includes finite and inf. fragments) %s', f_solid)

        #
        n = np.arange(1, n_frag + 1)  # number of clusters in a fragment
        # broken bridges per cluster of size n
        tau = n * (self.sigma - 1) + 2
        s = n - 1  # intact bridges per cluster of size n
        n_bridges = tau + s
        # Eq. 32 mass of a finite fragment of size n
        mw_frag_n = (n * self.ma + (n - 1) * self.mb * l / p +
                     tau * self.mb * delta_fac * 0.25)

        # Eqs (1-4)
        Qn = self.sig / n_bridges / n * binom.pmf(s, n_bridges, p)
        # Eq. (33) total mass of fragments of size
        m_frag_n = mw_frag_n * Qn
        self.__log.debug(
            'mass weight of finite fragments %s', mw_frag_n)

        # Eq 35 total mass of finite fragments
        m_frag = a * self.ma * Fp + b * self.mb * Kp

        # Eq 38 fraction of finite fragments
        # TODO this has to be corrected similarly to the f_gas equation
        f_frag = m_frag / mtot
        f_frag_n = m_frag_n / mtot
        self.__log.debug(
            'mass fraction of finite fragments %s', f_frag_n)

        self.__log.debug(
            'Total fraction of fragments sum %s / Eq.35 %s', f_frag_n.sum(), f_frag)

        return {'f_gas': f_gas,
                'f_solid': f_solid,
                'f_frag': f_frag,
                'f_frag_n': f_frag_n,
                'mw_frag_n': mw_frag_n,
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

    def _flash_distillation(self, df_gas, mw_gas, df_n, meta_n, mw_n,
                            fracr, T):

        def funct(x):
            '''Eq 54'''
            return np.sum(x_n_calc(x) * (k_n - 1))

        def x_n_calc(x):
            '''
            Eq. 52
            '''
            return z_n / (1 + (k_n - 1) * x)

        self.__log.debug('\n\nStart flash_distillation\n')

        a = 87058.0
        b = 299.0
        g = 0.5903
        # mole fraction of n-mers contained in the metaplast
        self.__log.debug('Increment of fragments %s', df_n)
        self.__log.debug('Previous metaplast %s', meta_n)
        self.__log.debug('Cross-linking correction %s', fracr)
        F_n = np.append((df_n + meta_n * fracr) / mw_n, df_gas / mw_gas)
        self.__log.debug('F_n (mole) %s', F_n)
        F = F_n.sum()
        mw = np.append(mw_n, mw_gas)
        # self.__log.debug('MW %s', mw)
        p_vap = a * np.exp(-b * mw ** g / T)
        self.__log.debug('p_vap %s', p_vap)
        k_n = p_vap * 101325 / self.pressure
        # self.__log.debug('kn %s', k_n)
        z_n = F_n / F
        # self.__log.debug('zn %s', z_n)
        if funct(0) * funct(0.999) > 0:
            self.__log.debug('No vapor')
            fract_v = 0
            V = 0
            L = F
            x_n = F_n / F
            y_n = np.zeros_like(x_n)
        else:
            fract_v = brentq(funct, 0, 0.999)
            self.__log.debug('V/F = %s', fract_v)
            V = fract_v * F  # moles of tar
            L = F - V
            # mole fraction of n-mers in the metaplast
            x_n = x_n_calc(fract_v)
            # mole fraction of n-mers released as tar
            y_n = k_n * x_n if V > 0 else np.zeros_like(x_n)

        meta_n_new = x_n[:-1] * L * mw_n
        tar_n_new = y_n[:-1] * V * mw_n

        assert np.allclose(
            meta_n_new + tar_n_new, F_n[:-1] * mw_n), 'Sum of xn+yn should be equal to Fn'
        self.__log.debug('metaplast fraction %s', x_n)
        self.__log.debug('tar fraction %s', y_n)

        self.__log.debug('Mass metaplast %s', meta_n_new.sum())
        self.__log.debug('Mass tar %s', tar_n_new.sum())

        return tar_n_new, meta_n_new
