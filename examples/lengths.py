#!/usr/bin/env python
from numpy import *
from pylab import *
from lattice import System, Lattice, Interaction

M_t = None
E_t = None
sizes = [100, 50, 30, 20, 10, 5]
for sz in sizes:
  inter = lambda r: -4
  lat = Lattice(sz, dim = 2)
  pot = Interaction(inter, 1)
  sys = System(lat, pot)
  sys.set_mu(8.0)
  E = []
  N = []
  temp = linspace(4, 1.1, 30)
  for i in temp:
    print i
    sys.set_T(i)
    ee, nn = sys.run(1000 * sz**2)
    ee, nn = sys.run(1000 * sz**2)
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
