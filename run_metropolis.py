import metropolis

m       = 1.0
k       = 1.0
xi      = 0.0
temp    = 300
step    = 0.1
Nprint  = 100
Ncycle  = 100000 - 1

ho_calc = metropolis.harm_osc( m, k, step, Nprint, Ncycle, temp, xi )
ho_calc.kernel()
