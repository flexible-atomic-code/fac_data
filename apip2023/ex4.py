from pfac.fac import *
import sys, os


z = int(sys.argv[1])
k = int(sys.argv[2])

InitializeMPI(4)

a = ATOMICSYMBOL[z]
p = '%s%02d'%(a,k)
pb = p+'b'
pa = p+'a'

WallTime('Beg '+p)

SetAtom(a)

if k == 1:
    Config('g1', '1s1')
    Config('g2', '2*1')
    Config('i1', ' ')
    gcs = ['g1', 'g2']
    dcs = []
    ics = ['i1']
elif k == 2:
    Config('g1', '1s2')
    Config('g2', '1s1 2*1')    
    Config('d2', '2*2')
    Config('i1', '1s1')
    gcs = ['g1','g2']
    dcs = ['d2']
    ics = ['i1']
elif k == 3:
    Config('g1', '1s2 2*1')
    Config('i1', '1s2')
    Config('i2', '1s1 2*1')
    gcs = ['g1']
    dcs = []
    ics = ['i1','i2']
    
ListConfig(pa+'.cfg')

WallTime('OPT')

ConfigEnergy(0)
OptimizeRadial('g1')
ConfigEnergy(1)

WallTime('EN')
Structure(pb+'.en', gcs)
Structure(pb+'.en', ics)
if len(dcs) > 0:
    Structure(pb+'.en', dcs)
MemENTable(pb+'.en')
PrintTable(pb+'.en', pa+'.en')

WallTime('TR')
for i in range(len(gcs)):
    for j in range(i,len(gcs)):
        TRTable(pb+'.tr', [gcs[i]], [gcs[j]])
    for j in range(len(dcs)):
        TRTable(pb+'.tr', [gcs[i]], [dcs[j]])
PrintTable(pb+'.tr', pa+'.tr')
   
WallTime('CE')
for i in range(len(gcs)):
    for j in range(i,len(gcs)):
        CETable(pb+'.ce', [gcs[i]], [gcs[j]])
PrintTable(pb+'.ce', pa+'.ce')

WallTime('CI')
for i in range(len(gcs)):
    for j in range(len(ics)):
        CITable(pb+'.ci', [gcs[i]], [ics[j]])
PrintTable(pb+'.ci', pa+'.ci')

WallTime('RR')
for i in range(len(gcs)):
    for j in range(len(ics)):
        RRTable(pb+'.rr', [gcs[i]], [ics[j]])
PrintTable(pb+'.rr', pa+'.rr')

WallTime('AI')
for i in range(len(dcs)):
    for j in range(len(ics)):
        AITable(pb+'.ai', [dcs[i]], [ics[j]])
PrintTable(pb+'.ai', pa+'.ai')

WallTime('Done')

FinalizeMPI()
