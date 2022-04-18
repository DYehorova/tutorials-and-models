import LJ

Natoms = 125
density = 0.8442
a = 2.5


delt   = 0.01
tmax   = 100.0
tequil = 2.0
Nstep  = round(tmax/delt)
Nprint = 100

LJ = LJ.LJ(tequil, delt, a, Natoms, Nstep, Nprint, density )

LJ.kernel()
