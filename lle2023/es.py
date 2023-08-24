from pylab import *
from pfac import rfac
from pfac import fac
import sys

z = int(sys.argv[1])
k0 = int(sys.argv[2])
k1 = int(sys.argv[3])

a = fac.ATOMICSYMBOL[z]

ds = linspace(0.1**(1/3), 10**(1/3), 15)**3.0 #electron density grid in 10^24 cm^-3

for k in range(k0, k1+1):
    d0 = rfac.FLEV('%s%02da.en'%(a,k))
    i2 = []
    for i in range(len(d0.e)):
        if d0.c[i][:3] == '1*1' and d0.v[i]//100 == 2 and d0.nele[i] == k:
            i2.append(i)

    de = zeros(len(ds))
    er = zeros((2,len(ds)))

    r0 = rfac.load_fac('%s%02dcfga.tr'%(a,k))
    w0 = where((r0[2]<i2[0])&(r0[0]>=i2[0])&(r0[0]<=i2[-1]))[0]
    e0 = sum(r0[5][w0]*r0[4][w0])/sum(r0[5][w0])
    er[0,:] = e0
    r0 = rfac.load_fac('%s%02da.tr'%(a,k))
    w0 = where((r0[2]<i2[0])&(r0[0]>=i2[0])&(r0[0]<=i2[-1]))[0]
    e0 = sum(r0[5][w0]*r0[4][w0])/sum(r0[5][w0])
    er[1,:] = e0
    for i in range(len(ds)):
        r = rfac.load_fac('%s%02dm0d%02da.tr'%(a,k,i))
        w = where((r[2]<i2[0])&(r[0]>=i2[0])&(r[0]<=i2[-1]))[0]
        e = sum(r[5][w]*r[4][w])/sum(r[5][w])
        de[i] = e-e0

    savetxt('%s%02da.es'%(a,k), transpose((ds,er[0],er[1],de)))
