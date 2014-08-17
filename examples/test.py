#!/usr/bin/env python
#from numpy import *
#from pylab import *
from latgas import lattice 
from latgas.lattice import mdsys, lat, pot

sz = 30
inter = lambda r: -4
lattice = lat(sz, dim = 2)
potential = pot(inter, 1.0, 5)
sys = mdsys(lattice, potential)
sys.set_mu(8.0)
sys.set_T(20.0)
sys.run(10000)
