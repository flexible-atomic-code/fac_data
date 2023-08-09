from pfac.fac import *

WallTime('Start')

SetAtom('Fe')

Config('g1', '1s2')
Config('g2', '1s1 2*1')

WallTime('Opt')
OptimizeRadial('g1')

WallTime('Wavefunctions')
WaveFuncTable('w1s.txt', 1, -1, 0)
WaveFuncTable('w2p-.txt', 2, 1, 0)
WaveFuncTable('w3d+.txt', 3, -2, 0)

WaveFuncTable('ws.txt', 0, -1, 2e3)
WaveFuncTable('wp-.txt', 0, 1, 2e3)

WallTime('EN')
Structure('ex2b.en', ['g1','g2'])
MemENTable('ex2b.en')
PrintTable('ex2b.en', 'ex2a.en')

BasisTable('ex2a.bs')
BasisTable('ex2a', 10)

WallTime('Done')
