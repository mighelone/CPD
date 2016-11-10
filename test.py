import cpd
import matplotlib.pyplot as plt

plt.style.use(['mystyle', 'mystyle-vega'])

ua = {'C': 80, 'H': 8, 'O': 12, 'N': 0, 'S': 0}
pa = {'FC': 45.1, 'VM': 50.6, 'Ash': 4.3, 'Moist': 19.0}

operating_conditions = [[0, 300], [0.01, 1200], [0.02, 1200]]

c = cpd.CPD(ultimate_analysis=ua, proximate_analysis=pa,
            pressure=101325, name='Test')
print('c0', c.c0)
print('p0', c.p0)
c.operating_conditions = operating_conditions

t, y = c._bridge_evolution()

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
