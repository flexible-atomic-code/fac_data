from pylab import *
from pfac import fac, rfac

clf()
for i in range(3,10):
    d = loadtxt('Cr%02da.es'%i, unpack=1)
    plot(d[0], -d[3], marker='o', label='%s-like'%fac.ATOMICSYMBOL[i])
legend()
xlabel(r'$n_e (10^{24} $cm$^{-3})$')
ylabel(r'redshift (eV)')
