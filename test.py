import nextrapol as nx
import settings_sirius_espadons
ls /Users/boehm/Desktop/Observations/Sirius_polarbase/Espadons/RAW/02feb10/
N=nx.Night("/Users/boehm/Desktop/Observations/Sirius_polarbase/Espadons/RAW/02feb10/",**settings_sirius_espadons.get_kwargs())
N.star4
N.star4[0].plot_voie(1,[52,53,54,55,56])
plt.show()

