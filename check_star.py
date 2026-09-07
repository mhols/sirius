import nextrapol as nx
import matplotlib.pyplot as plt

import sirius_reduction

tharfile = sirius_reduction.DATADIR + '/1164658o.fits.fz'

myext = nx.reduce_star(tharfile, **sirius_reduction.settings_sirius_espadons.get_kwargs())

plt.figure()
plt.title('voie 1')
myext.plot_voie(1, myext.ORDERS)

plt.figure()
plt.title('voie 2')
myext.plot_voie(2, myext.ORDERS)



plt.show()