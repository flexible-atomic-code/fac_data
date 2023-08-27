import os, sys
from pylab import *

z = int(sys.argv[1]) #atomin number
k0 = int(sys.argv[2]) #number of electrons 
k1 = int(sys.argv[3]) #number of electrons, calculation for ions in k0-k1 range
te = float(sys.argv[4]) #temp in eV
nm = int(sys.argv[5])

ds = exp(linspace(log(0.05), log(50), 11)) #electron density grid in 10^24 cm^-3

for i in range(len(ds)):
    for k in range(k0, k1+1):
        c = 'python plasma.py --z=%d --n=%d --d=%g --t=%g --nmax=%d --ce=0 --id=d%02d'%(z,k,ds[i],te,nm,i)
        print(c)
        if i == 0:
            os.system('%s --m=-1'%c)
            os.system('%s --m=-2'%c)
        os.system('%s --m=0'%c)

