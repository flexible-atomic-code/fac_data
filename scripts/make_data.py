from pfac.atom import atomic_data
from pfac.spm import *
from pfac import const
import os


THIS_DIR = os.path.abspath(os.path.join(__file__, os.pardir))


def make_spectrum(asym, neles, z):
    """ Save the spectral calculation for
    asym: string of element such as 'Li', 'He'...
    """
    dir0 = THIS_DIR + '/../{}/'.format(asym)
    if not os.path.exists(dir0):
        os.mkdir(dir0)

    atomic_data(neles, asym, iprint=1, dir=dir0)

    den = [1.0]  # density
    for nele in neles:
        dir1 = dir0 + '/{0:02d}/'.format(nele)
        if not os.path.exists(dir1):
            os.mkdir(dir1)

        (temp, logt, population) = get_tgrid(z, nele, dt=1)
        print('NION = 1')
        spectrum([nele], temp, den, population,
                 asym, dir0=dir0, dir1=dir1, nion=1)


make_spectrum('He', [1, 2], 2)
make_spectrum('Li', [1, 2, 3], 3)
make_spectrum('Be', [1, 2, 3, 4], 4)
