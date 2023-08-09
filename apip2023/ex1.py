from pfac.fac import *

WallTime('Start')

SetAtom('Fe')

Config('g1', '1s2')
Config('g2', '1s1 2*1')
Config('g3', '1s1 3*1')
ListConfig()
ListConfig('ex1a.cfg')

WallTime('Opt')
ConfigEnergy(0)
OptimizeRadial('g1')
ConfigEnergy(1)

WallTime('EN')
Structure('ex1b.en', ['g1','g2'])
Structure('ex1b.en', ['g3'])

MemENTable('ex1b.en')
PrintTable('ex1b.en', 'ex1a.en')

WallTime('TR')
TRTable('ex1b.tr', ['g1'], ['g2'])
TRTable('ex1b.tr', ['g1'], ['g3'])
TRTable('ex1b.tr', ['g2'], ['g2'])
TRTable('ex1b.tr', ['g2'], ['g3'])
TRTable('ex1b.tr', ['g3'], ['g3'])
PrintTable('ex1b.tr', 'ex1a.tr')

WallTime('CE')
CETable('ex1b.ce', ['g1'], ['g2'])
CETable('ex1b.ce', ['g1'], ['g3'])
CETable('ex1b.ce', ['g2'], ['g2'])
CETable('ex1b.ce', ['g2'], ['g3'])
PrintTable('ex1b.ce', 'ex1a.ce')

WallTime('Done')
