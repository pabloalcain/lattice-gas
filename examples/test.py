#!/usr/bin/env python
"""
test that everything works

This should be superseeded by an eventual unittesting
"""

from latgas import lattice 
from latgas.lattice import Lattice, Potential, System

sz = 30
def inter(r):
    if r <= 1.0: return -4
    if r <= 3.0: return 1
lattice = Lattice(sz, dim = 2)
potential = Potential(inter, 3.0, 20)
sys = System(lattice, potential)
sys.lattice.random(1)
sys.set_mu(8.0)
sys.set_T(20.0)
sys.N = sys.tot_population()
sys.E = sys.tot_energy()
sys.run(3000)
print "E = {0}, N = {1}".format(sys.E, sys.N)
