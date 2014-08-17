#!/usr/bin/env python
"""
test that everything works

This should be superseeded by an eventual unittesting
"""

from latgas import lattice 
from latgas.lattice import Lattice, Potential, System

sz = 30
inter = lambda r: -4
lattice = Lattice(sz, dim = 2)
potential = Potential(inter, 1.0, 5)
sys = System(lattice, potential)
sys.set_mu(8.0)
sys.set_T(20.0)
sys.run(10000)
print "E = {0}, N = {1}".format(sys.E, sys.N)
