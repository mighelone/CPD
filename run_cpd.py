import cpd
import matplotlib
import matplotlib.pyplot as plt
import logging

matplotlib.use('GTK3Agg')

log_level = logging.INFO

logging.basicConfig(level=log_level)

plt.style.use(['mystyle', 'mystyle-vega'])

if __name__ == '__main__':
    ua = {'C': 74.12, 'H': 4.96, 'O': 13.18, 'N': 1.45, 'S': 0}
    pa = {'FC': 57, 'VM': 43, 'Ash': 0, 'Moist': 0}

    operating_conditions = [[0, 1000], [0.01, 1000], [0.05, 1000]]

    c = cpd.CPD(ultimate_analysis=ua, proximate_analysis=pa,
                pressure=101325, name='Test')
    c.operating_conditions = operating_conditions

    t, y, fract = c._bridge_evolution()

    fig, axes = plt.subplots(
        nrows=2, ncols=1, sharex='col', figsize=(8, 12))
    ax1, ax2 = axes
    l = y[:, 0]
    delta = y[:, 1]
    c = y[:, 2]
    p = l + c
    f = 1 - p
    g1 = 2 * f - delta
    g2 = 2 * (c - c[0])
    ax1.plot(t, l, label='$l$')
    ax1.plot(t, delta / 2, label='$\delta/2$')
    ax1.plot(t, c, label='$c$')
    ax1.plot(t, p, label='$p$')
    ax1.plot(t, g1 / 2, label='$g_1/2$')
    ax1.plot(t, g2 / 2, label='$g_2/2$')
    ax1.set_ylabel('Fraction of bridges')
    ax1.legend()

    labels = ['solid', 'gas', 'tar', 'meta', 'cross']
    for i, l in enumerate(labels):
        ax2.plot(t, fract[:, i], label=l)
    ax2.legend()
    ax2.set_xlim([0, 0.05])
    ax2.set_xlabel('Time, s')
    ax2.set_ylabel('Yield')

    fig.savefig('cpd.png')

    plt.show()
