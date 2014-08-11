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
sys.set_T(20.0)
e, n = sys.run(10 * sz**2)
print "E = {0}, N = {1}".format(e, n)
