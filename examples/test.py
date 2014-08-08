#!/usr/bin/env python
from numpy import *
from pylab import *
from latgas import lattice 
from latgas.lattice import System, Lattice, Interaction

sz = 30
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
