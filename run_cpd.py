import cpd
import pkp
import matplotlib
import matplotlib.pyplot as plt
import logging
import pandas as pd

matplotlib.use('GTK3Agg')

log_level = logging.WARNING

logging.basicConfig(level=log_level)

plt.style.use(['mystyle', 'mystyle-vega'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


if __name__ == '__main__':
    ua = {'C': 74.12, 'H': 4.96, 'O': 13.18, 'N': 1.45, 'S': 0}
    pa = {'FC': 57, 'VM': 43, 'Ash': 0, 'Moist': 0}

    operating_conditions = [[0, 1000], [0.01, 1000], [0.05, 1000]]

    try:
        results0 = pd.read_csv('CPD_Test.csv', index_col=0)
    except FileNotFoundError:
        c0 = pkp.cpd.CPD(ultimate_analysis=ua, proximate_analysis=pa,
                         pressure=101325, name='Test')
        c0.operating_conditions = operating_conditions
        c0.set_parameters(dt=1e-5, dt_max=1e-4, increment=2)
        results0 = c0.run()

    coal = cpd.CPD(ultimate_analysis=ua, proximate_analysis=pa,
                   pressure=101325, name='Test')
    coal.operating_conditions = operating_conditions

    results = coal.run()

    fig, axes = plt.subplots(
        nrows=2, ncols=1, sharex='col', figsize=(8, 12))
    ax1, ax2 = axes
    results['delta/2'] = 0.5 * results['delta']

    results.plot(x='t', y='l', label='$l$', color=colors[0], ax=ax1)
    results.plot(x='t', y='delta/2', label='$delta/2$', color=colors[1],
                 ax=ax1)
    results.plot(x='t', y='p', label='$c$', color=colors[2], ax=ax1)
    #ax1.plot(t, l, label='$l$', color=colors[0])
    #ax1.plot(t, delta / 2, label='$\delta/2$', color=colors[1])
    #ax1.plot(t, c + l, label='$p$', color=colors[2])

    ax1.plot(results0.index, results0['l'], color=colors[0],
             linestyle='dashed')
    ax1.plot(results0.index, results0['del/2'], color=colors[1],
             linestyle='dashed')
    ax1.plot(results0.index, results0['p'], color=colors[2],
             linestyle='dashed')
    #ax1.plot(t, p, label='$p$')
    #ax1.plot(t, g1 / 2, label='$g_1/2$')
    #ax1.plot(t, g2 / 2, label='$g_2/2$')
    ax1.set_ylabel('Fraction of bridges')
    ax1.legend()

    labels = ['char', 'light_gas', 'tar']
    for i, l in enumerate(labels):
        results.plot(x='t', y=l, label=l, color=colors[i], ax=ax2)
        ax2.plot(results0.index, results0[l],
                 color=colors[i], linestyle='dashed')
    ax2.legend()
    ax2.set_xlim([0, 0.05])
    ax2.set_xlabel('Time, s')
    ax2.set_ylabel('Yield')

    fig.savefig('cpd.png')

    plt.show()
