import nextrapol as nx
import settings_sirius_espadons
import os
import matplotlib.pyplot as plt

DATADIR = os.environ['SIRIUS02FEB10'] 


def run1():
    N=nx.Night(DATADIR,**settings_sirius_espadons.get_kwargs())
    N.star4
    #N.star4[0].plot_voie(1,[52,53,54,55,56])
    #plt.show()
    N.stokes

    for o in N.kwargs['ORDERS']:
        I = N.stokes.true_order_number==o

        plt.plot(N.stokes.lambdagrid[I], N.stokes.intensity[I])

    plt.show()



if __name__ == '__main__':

    print (nx.store.keys)

    run1()

