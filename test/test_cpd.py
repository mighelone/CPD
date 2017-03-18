import cpd
import pkp
import matplotlib
import matplotlib.pyplot as plt
import logging
import pandas as pd
import numpy as np
import pytest
import matplotlib.pyplot as plt

plt.style.use('bmh')


colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

ua = {'C': 74.12, 'H': 4.96, 'O': 13.18, 'N': 1.45, 'S': 0}
pa = {'FC': 57, 'VM': 43, 'Ash': 0, 'Moist': 0}

operating_conditions = [[0, 600], [0.01, 1200], [0.02, 1200]]

pressure = 101325

coals = [
    {
        'name': 'Pittsburgh',
        'nmr': {'mdel': 24, 'mw': 294, 'p0': 0.62, 'sig': 4.5, 'c0': 0}
    },
    {
        'name': '1443-Lignite-ACER',
        'nmr': {'mdel': 36, 'mw': 297, 'p0': 0.59, 'sig': 4.8, 'c0': 0.2}
    },
    {
        'name': '1448-Subbituminous-ACER)',
        'nmr': {'mdel': 37, 'mw': 310, 'p0': 0.54, 'sig': 4.7, 'c0': 0.15}
    },
    {
        'name': 'Illinois-6',
        'nmr': {'mdel': 27, 'mw': 316, 'p0': 0.63, 'sig': 5, 'c0': 0}
    },
    {
        'name': 'Zap',
        'nmr': {'mdel': 40, 'mw': 277, 'p0': 0.63, 'sig': 3.9, 'c0': 0.2}
    }
]


def run_cpd(solver, ua, pa, oc, t=None, **parameters):
    coal = solver(ultimate_analysis=ua, proximate_analysis=pa,
                  pressure=pressure, name='Test')
    coal.set_parameters(**parameters)
    coal.operating_conditions = oc
    res = coal.run(time=t)
    return res.reset_index()


def compare_results(res0, res, variable, atol=1e-2):
    """
    Compare results from the two solutions. The second solution is interpolate
    against the time of the first.
    """
    var_ip = np.interp(res0['t'], res['t'], res[variable])
    try:
        return np.allclose(res0[variable], var_ip, atol=atol)
    except ValueError as e:
        raise(e)


@pytest.fixture(params=coals)
def coal(request):
    return request.param


def test_cpd(coal):
    res_ref = run_cpd(pkp.cpd.CPD, ua, pa, operating_conditions,
                      dt=1e-6, dt_max=1e-4, nmr_parameters=coal['nmr'])
    res_ref['delta/2'] = res_ref['del/2']
    res = run_cpd(cpd.CPD, ua, pa, operating_conditions,
                  dt=1e-6, dt_max=1e-4, nmr_parameters=coal['nmr'])
    res['delta/2'] = 0.5 * res['delta']

    # compare main species
    fig, axes = plt.subplots(
        nrows=2, ncols=1, sharex='col', figsize=(8, 12))
    colors = plt.rcParams['axes.prop_cycle']()
    ax1, ax2 = axes
    ax1.set_title('Coal: {}'.format(coal['name']))
    labels = ('l', 'delta/2', 'p')
    bridge_comparison = []
    for sp in labels:  # , 'tar'):  # , 'tar'):
        print('compare {}'.format(sp))
        color = next(colors)
        res.plot('t', sp, linestyle='solid', ax=ax1, **color)
        res_ref.plot('t', sp, linestyle='dashed', ax=ax1, **color)
        bridge_comparison.append(compare_results(res_ref, res, sp, atol=5e-2))
    ax1.set_ylabel('Fraction of bridges')
    ax1.legend(ax1.lines[::2], labels, frameon=False,
               loc='upper right')

    yield_comparison = []
    colors = plt.rcParams['axes.prop_cycle']()
    labels = ('char', 'tar')
    for sp in labels:  # , 'tar'):  # , 'tar'):
        print('compare {}'.format(sp))
        color = next(colors)
        res.plot('t', sp, linestyle='solid', ax=ax2, **color)
        res_ref.plot('t', sp, linestyle='dashed', ax=ax2, **color)
        yield_comparison.append(compare_results(res_ref, res, sp, atol=5e-2))
    ax2.legend(ax2.lines[::2], labels, frameon=False, loc='upper right')
    ax2.set_xlim([0, 0.02])
    ax2.set_xlabel('Time, s')
    ax2.set_ylabel('Yield')

    fig.savefig('test/cpd-{}.png'.format(coal['name']))

    assert all(bridge_comparison)
    assert all(yield_comparison)
