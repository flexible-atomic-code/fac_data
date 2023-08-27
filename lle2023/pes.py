from pylab import *
from pfac import fac, rfac, const

clf()
eb = []
for i in range(3,10):
    d = loadtxt('Cr%02da.es'%i, unpack=1)
    loglog(d[0], -d[3], marker='o', label='%s-like'%fac.ATOMICSYMBOL[i])
    x = loadtxt('Cr%02db.rp'%i, usecols=6)
    eb.append(x)
    for j in range(10):
        x = loadtxt('Cr%02dm0d%02db.rp'%(i,j), usecols=6)
        eb.append(x)
eb = array(eb).reshape((7,11,4))*const.Hartree_eV

legend()
xlabel(r'$\frac{n_e}{10^{24}}$')
ylabel(r'redshift (eV)')
