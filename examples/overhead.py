"""
This script is a Q&D way to find the overhead
of calling the "run" function, comparing the
crun method (loops in c) and the run method (loops
in python). We can eventually add an "autocorrelation"
measure to compare times.

"""
#!/usr/bin/env python
from latgas import lattice 
from latgas.lattice import Lattice, Potential, System
from time import time
from pylab import plot, legend, xlabel, ylabel, show
from numpy import array

sz = 30
inter = lambda r: -4
lattice = Lattice(sz, dim = 2)
potential = Potential(inter, 1.0, 5)
sys = System(lattice, potential)
sys.set_mu(8.0)
sys.set_T(20.0)


steps = range(1, 201)
Trun = []
Tpy = []
Tc = []
for i in steps:
    print i
    t1 = time()
    for j in xrange(i):
        sys.whole_lattice()
    t2 = time()
    Tpy.append(1000*(t2-t1))
    t1 = time()
    sys.run(i)
    t2 = time()
    Trun.append(1000*(t2-t1))
    t1 = time()
    sys.crun(i)
    t2 = time()
    Tc.append(1000*(t2-t1))
steps = array(steps)
Trun = array(Trun)
Tpy = array(Tpy)
Tc = array(Tc)
plot(steps, Tpy/steps, steps, Trun/steps, steps, Tc/steps)
legend(["Python stripped", "Python", "c"])
xlabel("Run length")
ylabel("Time per step [ms]")
show()
