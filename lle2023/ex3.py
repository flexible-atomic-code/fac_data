from pfac.crm import *
import os

AddIon(2, 0.0, 'ex1b')
SetBlocks(-1)

SetEleDist(0, 3e3, -1, -1)

SetTRRates(0)
SetCERates(1)

SetAbund(2, 1.0)
SetEleDensity(1)

InitBlocks()
LevelPopulation()

DumpRates('ex1a.r0', 2, 0, -1, 1)
DumpRates('ex1a.r1', 2, 1, -1, 1)
DumpRates('ex1a.r2', 2, 2, -1, 1)
DumpRates('ex1a.r3', 2, 3, -1, 1)

SpecTable('ex1b.sp', 0)
PrintTable('ex1b.sp', 'ex1a.sp')

if os.path.exists('ex1a.ln'):
    System('rm ex1a.ln')
SelectLines('ex1b.sp', 'ex1a.ln', 2, 0, 6e3, 10e3, 1e-5)
