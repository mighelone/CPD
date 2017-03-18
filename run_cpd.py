import cpd
import pkp
import matplotlib
import matplotlib.pyplot as plt
import logging
import pandas as pd

# matplotlib.use('GTK3Agg')

log_level = logging.WARNING

logging.basicConfig(level=log_level)

plt.style.use(['mystyle', 'mystyle-vega'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

plot = False


def main():
    ua = {'C': 74.12, 'H': 4.96, 'O': 13.18, 'N': 1.45, 'S': 0}
    pa = {'FC': 57, 'VM': 43, 'Ash': 0, 'Moist': 0}

    operating_conditions = [[0, 1000], [0.01, 1000], [0.05, 1000]]

    try:
        results0 = pd.read_csv('CPD_Test.csv')
    except:
        c0 = pkp.cpd.CPD(ultimate_analysis=ua, proximate_analysis=pa,
                         pressure=101325, name='Test')
        c0.operating_conditions = operating_conditions
        c0.set_parameters(dt=1e-5, dt_max=1e-4, increment=2)
        results0 = c0.run().reset_index()

    coal = cpd.CPD(ultimate_analysis=ua, proximate_analysis=pa,
                   pressure=101325, name='Test')
    coal.set_parameters(dt=1e-6, dt_max=1e-4)
    coal.operating_conditions = operating_conditions

    results = coal.run(n_frag=20)
    results['delta/2'] = 0.5 * results['delta']
    labels = ['$l$', '$\delta/2$', '$p$']
    vars_new = ['l', 'delta/2', 'p']
    vars_old = ['l', 'del/2', 'p']

    fig, axes = plt.subplots(
        nrows=2, ncols=1, sharex='col', figsize=(8, 12))
    ax1, ax2 = axes

    for vn, vo, c in zip(vars_new, vars_old, colors):
        results.plot(x='t', y=vn, color=c, ax=ax1)
        results0.plot(x='t', y=vo, color=c, linestyle='dashed', ax=ax1)
    ax1.set_ylabel('Fraction of bridges')
    # add an extra legend
    # http://matplotlib.org/users/legend_guide.html#multiple-legend
    # ax1.add_artist(plt.legend(ax1.lines[0:2], ['Python', 'Original'],
    #                          loc='lower right', frameon=False))
    ax1.legend(ax1.lines[::2], labels, frameon=False,
               loc='upper right')

    # ax2 = plt.subplots(212)
    labels = ['char', 'light_gas', 'tar']
    for i, l in enumerate(labels):
        results.plot(x='t', y=l, label=l, color=colors[i], ax=ax2)
        results0.plot(x='t', y=l, color=colors[i], linestyle='dashed',
                      ax=ax2)
    ax2.legend(ax2.lines[::2], labels, frameon=False, loc='upper right')
    ax2.set_xlim([0, 0.05])
    ax2.set_xlabel('Time, s')
    ax2.set_ylabel('Yield')

    # fig.savefig('cpd.png')

    # plt.show()


if __name__ == '__main__':
    main()
