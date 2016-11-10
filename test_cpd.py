import cpd
import matplotlib.pyplot as plt
import logging

plt.style.use(['mystyle', 'mystyle-vega'])

ua = {'C': 80, 'H': 8, 'O': 12, 'N': 0, 'S': 0}
pa = {'FC': 45.1, 'VM': 50.6, 'Ash': 4.3, 'Moist': 19.0}

operating_conditions = [[0, 300], [0.01, 1200], [0.02, 1200]]


def create_cpd():
    c = cpd.CPD(ultimate_analysis=ua, proximate_analysis=pa,
                pressure=101325, name='Test')
    c.operating_conditions = operating_conditions
    return c


def test_ode():
    coal = create_cpd()
    t, y = coal._bridge_evolution()
    ax = plt.subplot(111)
    l = y[:, 0]
    delta = y[:, 1]
    c = y[:, 2]
    p = l + c
    f = 1 - p
    g1 = 2 * f - delta
    g2 = 2 * (c - c[0])
    ax.plot(t, l, label='$l$')
    ax.plot(t, delta / 2, label='$\delta/2$')
    ax.plot(t, c, label='$c$')
    ax.plot(t, p, label='$p$')
    ax.plot(t, g1 / 2, label='$g_1/2$')
    ax.plot(t, g2 / 2, label='$g_2/2$')
    ax.legend()

    ax.figure.save('cpd_0.png')
    plt.close(ax.figure)


def test_percolation():
    coal = create_cpd()
    y = [0.01, 0.3, 0.05]
    f_tar = 0
    percolation = coal._percolation(y, f_tar, True)


def test_flash():
    coal = create_cpd()
    y = [0.01, 0.3, 0.05]
    percolation = coal._percolation(y, 0, True)
    df_gas = percolation['f_gas'] * 0.0001
    f = 0.01
    f_frag_n = percolation['f_frag_n']
    df_n = f_frag_n * f
    mw_gas = 28.0
    mw_n = percolation['m_frag_n']
    meta_n = f_frag_n * (1 - f)
    T = 800.0
    fracr = 0.0
    coal.flash_distillation(df_gas=df_gas, mw_gas=mw_gas, df_n=df_n, meta_n=meta_n,
                            mw_n=mw_n, fracr=fracr, T=T)


log_level = logging.DEBUG

logging.basicConfig(level=log_level)
test_flash()
