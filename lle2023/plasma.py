from pfac.fac import *
from pfac import rfac
from pfac import const
from optparse import OptionParser
from numpy import *
import sys, os, warnings

#ConvertToSFAC('plasma.sf')
#the model used in AA624A74 paper would most resemble these options (for ne=1e22 and te=50 eV)
#python plasma.py 1 0.01 50.0 8 0

p = OptionParser()
p.add_option('-z', '--z', dest='z', type='int',
             default=36, help='atomic number')
p.add_option('-n', '--n', dest='n', type='int',
             default=2, help='number of electrons')
p.add_option('-m', '--m', dest='m', type='int',
             default=-1, help='screening mode')
p.add_option('-d', '--d', dest='d', type='float',
             default=1.0, help='electron density in 1e24')
p.add_option('-t', '--t', dest='t', type='float',
             default=1e3, help='electron temperature in eV')
p.add_option('--zs', dest='zs', type='float',
             default=1, help='effective screen charge in Stewart-Pyatt')
p.add_option('--zp', dest='zp', type='float',
             default=0, help='number of free electrons per ion')
p.add_option('--vxf', dest='vxf', type='int',
             default=0, help='exchange potential for free electrons')
p.add_option('--id', dest='id', type='string',
             default='', help='id for the output file names')
p.add_option('--i0', dest='i0', type='int',
             default=-1, help='lower level n')
p.add_option('--i1', dest='i1', type='int',
             default=-1, help='upper level n')
p.add_option('--nmax', dest='nmax', type='int',
             default=5, help='max n for excitation')
p.add_option('--ce', dest='ce', type='int',
             default=0, help='do collisional excitation')
p.add_option('--slev', dest='slev', type='int',
             default=2, help='level sort mode')
p.add_option('--np', dest='np', type='int',
             default=0, help='openmp threads')
p.add_option('--pwb', dest='pwb', type='int',
             default=0, help='plane-wave born')
p.add_option('--ci', dest='ci', type='int',
             default=0, help='do collisional ionization')
p.add_option('--ng', dest='ng', type='int',
             default=5000, help='radial grid points')
p.add_option('--ws', dest='ws', type='int',
             default=0, help='calculate width/shift')
p.add_option('--bms', dest='bms', type='float',
             default=0, help='collision mass')
p.add_option('--dhm', dest='dhm', type='int',
             default=0, help='debye-huckel mode')
p.add_option('--psp', dest='psp', type='int',
             default=0, help='print SP mode')

(opts, args) = p.parse_args()
print(opts)

m = opts.m
d = opts.d
t = opts.t
z = opts.z
zs = opts.zs
nmax = opts.nmax
n = opts.n
id = opts.id
i0 = opts.i0
i1 = opts.i1
bms = opts.bms
                    
if id == ' ' or id == 'None' or id == 'none':
    id = ''

if opts.zp <= 0:
    zp = z-n
    if zp < 1:
        zp = 1
else:
    zp = opts.zp
    
a = ATOMICSYMBOL[z]
if m >= 0:
    p = '%s%02dm%d%s'%(a,n,m,id)
    p0 = '%s%02d'%(a,n)
elif m == -1:
    p = '%s%02d'%(a,n)
    p0 = '%s%02d'%(a,n)
else:
    p = '%s%02dcfg'%(a,n)
    p0 = '%s%02d'%(a,n)

if opts.ws:
    r = rfac.read_tr(p0+'a.tr')
    e = r[1][0]['Delta E']
    w = r[1][0]['rate']*(1+r[1][0]['upper_2J'])
    x0 = const.hc/(sum(w*e)/sum(w))
    r = rfac.read_tr(p+'a.tr')
    e = r[1][0]['Delta E']
    w = r[1][0]['rate']*(1+r[1][0]['upper_2J'])
    x = const.hc/(sum(w*e)/sum(w))
    de = const.hc*abs(1/x-1/x0)

    if (bms > 0):
        SetBornMass(bms)
        bms *= 1823.0
    else:
        bms = 1.0
    MemENTable(p+'b.en')
    i0 = r[1][0]['lower_index']
    i1 = r[1][0]['upper_index']

    cf = d*1e-4*const.hc/(6*pi)
    wd = sqrt(t/(bms*const.Hartree_eV))*0.529*(4*pi/3*d)**0.333*const.Hartree_eV
    
    ni = len(i0)
    ii = zeros(2*ni,dtype=int32)
    ii[:ni] = i0
    ii[ni:] = i1
    ui = unique(ii)
    nu = len(ui)
    cr0 = zeros(1+max(ui))
    cr = zeros(1+max(ui))
    a = 2.0
    for i in range(nu):
        CERate(p+'b.ce', 'i0.ce', ui[i], -1, t)
        CERate(p+'b.ce', 'i1.ce', -1, ui[i], t)
        os.system('grep "#" i0.ce > i0h.ce')
        os.system('grep "#" i1.ce > i1h.ce')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c0h = transpose(loadtxt('i0h.ce', comments=None, usecols=(1,2,3,4,5)))
            c1h = transpose(loadtxt('i1h.ce', comments=None, usecols=(1,2,3,4,5)))
            c0 = transpose(loadtxt('i0.ce'))
            c1 = transpose(loadtxt('i1.ce'))
        cr[ui[i]] = 0.0
        cr0[ui[i]] = 0.0
        if len(c0) > 0:
            c0 = c0[2]*cf
            cr0[ui[i]] += sum(c0)
            c0r = c0/c0h[4]
            c0 = 0.5*a*c0h[4]*(sqrt(1+4*c0r/a)-1.0)
            cr[ui[i]] += sum(c0)
        if len(c1) > 0:
            c1 = c1[2]*exp(c1h[4]/t)*(1+c1h[1])/(1+c1h[3])*cf
            cr0[ui[i]] += sum(c1)
            c1r = c1/c1h[4]
            c1 = 0.5*a*c1h[4]*(sqrt(1+4*c1r/a)-1.0)
            cr[ui[i]] += sum(c1)

    ui0 = cr0[i0]+cr0[i1]
    wu0 = ui0*const.hc/e**2
    wz = sum(w*wu0)/sum(w)
    wi0 = cr[i0]+cr[i1]
    rid = wi0/wd
    wi = wi0/(1+rid/1.2)
    wr0 = wi0*const.hc/e**2
    wr = wi*const.hc/e**2
    wx0 = sum(w*wr0)/sum(w)
    wx = sum(w*wr)/sum(w)
    print(w/sum(w))
    print(wr)
    print('%d %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E'%(m, d, t, x, wz, wx0, wx, x-x0, 1e8*(1/x-1/x0)))
    exit()

if opts.np > 1:
    InitializeMPI(opts.np)

SetAtom(a)
"""
plasma screening model:
zp: number of free electrons per ion, if 0, the ion's net charge is used
d: electron density in 10^24 cm^-3
t: electron temperature in eV
m: model choice.
   0 -- ion-sphere, if t==0 in this mode, a uniform ion sphere is used.
   1 -- debye-huckel
   2 -- Stewart-Pyatt model.
"""

#SetOption('orbital:relativistic_fermi', 1)
SetOption('orbital:sp_print', 3)
SetOption('orbital:sp_mode', 2)
SetOption('orbital:sp_ofn', 'spn.txt')
SetOption('orbital:fermi_rmf', 'fermi.txt')
SetOption('structure:sort_levs', opts.slev)
if opts.m < -1:
    SetOption('radial:config_energy', 0)
    
SetRadialGrid(opts.ng,-1,-1,-1)
if m >= 0:
    if opts.dhm > 0:
        SetOption('orbital:debye_mode', opts.dhm)
    #print('%d %d %d %g %g'%(m,zp,zs,d,t))
    PlasmaScreen(zp, d, t, m, zs, opts.vxf)
    
nmin = 1
gs=[]
if opts.psp > 0:
    SetOption('radial:print_spm', opts.psp)
if n == 1:
    nmin = 1
    for n0 in range(1,nmax+1):
        Config('g%d'%n0, '%d*1'%n0)        
        gs.append('g%d'%n0)
        if n0 == 1:
            OptimizeRadial(gs[0])
            GetPotential(p+'a.pot')
            if m >= 0:
                SetOrbNMax(20, 0, 10)
    Config('i0', '')
elif n == 2:
    nmin = 1
    for n0 in range(1, nmax+1):
        if n0 == 1:
            Config('g%d'%n0, '1s2')
        else:
            Config('g%d'%n0, '1s1 %d*1'%n0)
        gs.append('g%d'%n0)
        if n0 == 1:
            OptimizeRadial(gs[0])
            GetPotential(p+'a.pot')
            if m >= 0:
                SetOrbNMax(20, 0, 10)
    Config('i0', '1s1')
elif n == 3:
    nmin = 2
    for n0 in range(2, nmax+1):
        Config('g%d'%n0, '1s2 %d*1'%n0)
        gs.append('g%d'%n0)
        Config('e%d'%n0, '1s1 %d*2'%n0)
        gs.append('e%d'%n0)
        if n0 == 2:
            OptimizeRadial(gs[0])
            GetPotential(p+'a.pot')
            if m >= 0:
                SetOrbNMax(20, 0, 10)
    Config('i0', '1s2')
    Config('i1', '1s1 2*1')
elif n < 10:
    nmin = 2
    Config('g2', '1s2 2*%d'%(n-2))
    gs = ['g2']
    Config('e2', '1s1 2*%d'%(n-1))
    gs.append('e2')
    OptimizeRadial(gs[0])
    GetPotential(p+'a.pot')
    if m >= 0:
        SetOrbNMax(20, 0, 10)
    for n0 in range(3, nmax+1):
        Config('g%d'%n0, '1s2 2*%d %d*1'%(n-3,n0))
        gs.append('g%d'%n0)
    for n0 in range(3, nmax+1):
        Config('e%d'%n0, '1s1 2*%d %d*1'%(n-2,n0))
        gs.append('g%d'%n0)
    Config('i0', '1s2 2*%d'%(n-3))
    Config('i1', '1s1 2*%d'%(n-2))
elif n == 10:
    nmin = 2
    Config('g2', '1s2 2*8')
    gs = ['g2']
    OptimizeRadial(gs[0])
    GetPotential(p+'a.pot')
    if m >= 0:
        SetOrbNMax(20, 0, 10)
    for n0 in range(3, nmax+1):
        Config('g%d'%n0, '1s2 2*7 %d*1'%n0)
        gs.append('g%d'%n0)
        Config('e%d'%n0, '1s1 2*8 %d*1'%n0)
        gs.append('e%d'%n0)
    Config('i0', '1s2 2*7')
    Config('i1', '1s1 2*8')
elif n == 11:
    nmin = 3
    for n0 in range(3, nmax+1):
        Config('g%d'%n0, '1s2 2*8 %d*1'%n0)
        gs.append('g%d'%n0)
        if n0 == 3:
            OptimizeRadial(gs[0])
            GetPotential(p+'a.pot')
            if m >= 0:
                SetOrbNMax(20, 0, 10)
    Config('i0', '1s2 2*8')
    Config('i1', '1s2 2*7 3*1')
elif n == 12:
    nmin = 3
    Config('g0', '1s2 2*8 3*2')
    Config('g1', '1s2 2*8 3*1 4*1')
    Config('g1', '1s2 2*8 3*1 5*1')
    Config('i0', '1s2 2*8 3*1')
    Config('i1', '1s2 2*7 3*2')

#ConfigEnergy(0)
#OptimizeRadial(gs[0])
#ConfigEnergy(1)

for g in gs:
    Structure(p+'b.en', [g])
Structure(p+'b.en', ['i0'])
if n > 2:
    Structure(p+'b.en', ['i1'])
MemENTable(p+'b.en')
PrintTable(p+'b.en', p+'a.en')

BasisTable(p+'a.bs')

OrbitalStats(p+'b.rp', nmax)
WaveFuncTable(p+'a.w1s', 1, -1, 0)
WaveFuncTable(p+'a.w2s', 2, -1, 0)
WaveFuncTable(p+'a.w2pm', 2, 1, 0)
WaveFuncTable(p+'a.w2pp', 2, -2, 0)

if i0 < 0 or i1 < 0:
    for n0 in range(nmin, nmax+1):
        for n1 in range(n0, nmax+1):
            TRTable(p+'b.tr', ['g%d'%n0], ['g%d'%n1])
            TRTable(p+'b.tr', ['g%d'%n0], ['e%d'%n1])
            if n0 > nmin:
                TRTable(p+'b.tr', ['g%d'%nmin], ['e%d'%n1])
else:
    TRTable(p+'b.tr', ['g%d'%i0], ['g%d'%i1])
PrintTable(p+'b.tr', p+'a.tr')

if m >= 0:
    if (opts.pwb):
        SetCEBorn(0.0,0)
    if (bms > 0):
        SetBornMass(bms)
        #SetOption('orbital:zcoll', -1)
    if opts.ce:
        for n0 in range(nmin, nmax+1):
            for n1 in range(n0, nmax+1):
                WallTime('CE %d %d'%(n0, n1))
                if i0 >= 0 and i1 >= 0:
                    if (n0 == i0 or n0 == i1 or n1 == i0 or n1 == i1):
                        CETable(p+'b.ce', ['g%d'%n0], ['g%d'%n1])
                else:
                    CETable(p+'b.ce', ['g%d'%n0], ['g%d'%n1])
        PrintTable(p+'b.ce', p+'a.ce')
    if opts.ci:
        CITable(p+'b.ci', ['g%d'%i0], ['i0'])        
        if i1 != i0:
            CITable(p+'b.ci', ['g%d'%i1], ['i0'])
        PrintTable(p+'b.ci', p+'a.ci')
if opts.np > 1:
    FinalizeMPI()
#CloseSFAC()
