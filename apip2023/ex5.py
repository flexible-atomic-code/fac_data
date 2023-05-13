from pfac.crm import *
from pfac import fac
import sys, os

z = int(sys.argv[1])
k = int(sys.argv[2])
t = float(sys.argv[3])
d = float(sys.argv[4])
ai = float(sys.argv[5])
ar = float(sys.argv[6])

a = fac.ATOMICSYMBOL[z]

p = '%s%02d'%(a,k)
pbi = '%s%02db'%(a,k-1)
pbr = '%s%02db'%(a,k+1)
pb = p+'b'
pa = p+'a'

WallTime('AddIons')
AddIon(k-1, 0.0, pbi)
AddIon(k, 0.0, pb)
AddIon(k+1, 0.0, pbr)

SetBlocks(0)

SetEleDist(0, t, -1, -1)

WallTime('SetRates')
SetTRRates(0)
SetCERates(1)
SetCIRates(1)
SetRRRates(0)
SetAIRates(1)

SetAbund(k-1, ai)
SetAbund(k, 1.0)
SetAbund(k+1, ar)

SetEleDensity(d)

InitBlocks()
LevelPopulation()
DumpRates(pa+'.r0', -1, 0, -1, 1)
DumpRates(pa+'.r1', -1, 1, -1, 1)
DumpRates(pa+'.r3', -1, 3, -1, 1)

SpecTable(pb+'.sp', 0)
PrintTable(pb+'.sp', pa+'.sp')

if os.path.exists(pa+'.ln'):
    os.system('rm '+pa+'.ln')

SelectLines(pb+'.sp', pa+'.ln', 1, 201, 0, 1e5, 0)
SelectLines(pb+'.sp', pa+'.ln', 2, 201, 0, 1e5, 0)
SelectLines(pb+'.sp', pa+'.ln', 2, 20201, 0, 1e5, 0)

WallTime('Done')
