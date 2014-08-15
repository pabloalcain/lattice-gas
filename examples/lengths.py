#!/usr/bin/env python
from numpy import *
from pylab import *
from latgas import lattice
from latgas.lattice import mdsys, lat, pot

M_t = None
E_t = None
sizes = [100, 50, 30, 20, 10, 5]
for sz in sizes:
  inter = lambda r: -4 + r * 0
  lattice = lat(sz, dim = 2)
  potential = pot(inter, 1.0, 5)
  sys = mdsys(lattice, potential)
  sys.set_mu(8.0)
  E = []
  N = []
  temp = linspace(4, 1.1, 30)
  for i in temp:
    print i
    sys.set_T(i)
    ee, nn = sys.run(1000)
    ee, nn = sys.run(1000)
    E.append(ee)
    N.append(nn)
  E = array(E); N = array(N)
  M = -(2*N - sz**2)/(sz**2)
  
  if E_t == None and M_t == None:
    E_t = E
    M_t = M
  else:
    E_t = vstack([E_t, E])
    M_t = vstack([M_t, M])
