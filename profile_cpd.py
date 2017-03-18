import cpd
from datetime import datetime

if __name__ == '__main__':
    ua = {'C': 74.12, 'H': 4.96, 'O': 13.18, 'N': 1.45, 'S': 0}
    pa = {'FC': 57, 'VM': 43, 'Ash': 0, 'Moist': 0}

    operating_conditions = [[0, 1000], [0.01, 1000], [0.05, 1000]]

    coal = cpd.CPD(ultimate_analysis=ua, proximate_analysis=pa,
                   pressure=101325, name='Test')
    coal.set_parameters(dt=1e-6, dt_max=1e-4)
    coal.operating_conditions = operating_conditions

    start_time = datetime.now()
    results = coal.run(n_frag=20)
    print('Ex. time = {}'.format(datetime.now() - start_time))
