#!/usr/bin/env python
#from numpy import *
#from pylab import *
from latgas import lattice 
from latgas.lattice import *

sz = 30
inter = lambda r: -4
lattice = Lattice(sz, dim = 2)
potential = Potential(inter, 1.0, 5)
sys = System(lattice, potential)
sys.set_mu(8.0)
sys.set_T(20.0)
sys.run(10000)
