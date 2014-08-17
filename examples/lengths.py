#!/usr/bin/env python
from numpy import *
from pylab import *
from latgas import lattice
from latgas.lattice import *

M_t = None
E_t = None
sizes = [1000, 500, 200, 100]
for sz in sizes:
  inter = lambda r: -4
  lattice = Lattice(sz, dim = 2)
  potential = Potential(inter, 1.0, 5)
  sys = System(lattice, potential)
  sys.set_mu(8.0)
  E = []
  N = []
  temp = linspace(2.5, 2.0, 26)
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
for i in M_t:
  plot(temp, i)
legend([str(i) for i in sizes])
