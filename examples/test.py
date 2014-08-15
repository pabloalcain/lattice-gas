#!/usr/bin/env python
#from numpy import *
#from pylab import *
from latgas import lattice 
from latgas.lattice import mdsys, lat, pot

sz = 30
inter = lambda r: r * 0 -4
lattice = lat(sz, dim = 2)
potential = pot(inter, 1.0, 5)
sys = mdsys(lattice, potential)
sys.set_mu(8.0)
sys.set_T(20.0)
e, n = sys.run(10)
print "E = {0}, N = {1}".format(e, n)
