from __future__ import division
import numpy as np
import itertools as it
import random as R
import time

class System(object):
    def __init__(self, lattice, interaction, T = 1.0, mu = 0.0):
        """
        Constructor: give the lattice and thermodynamic info
        """
        self.lattice = lattice
        self.interaction = interaction
        self.T = T
        self.mu = mu
        self.E = self.tot_energy()
        self.N = self.tot_population()
    
    def set_T(self, T):
        """
        Encapsulated setter of temperature
        """
        self.T = T

    def set_mu(self, mu):
        """
        Encapsulated setter of chemical potential
        """
        self.mu = mu

    def get_T(self):
        """
        Encapsulated getter of temperature
        """
        return self.T
    
    def get_mu(self):
        """
        Encapsulated getter of chemical potential
        """
        return self.mu

    def tot_energy(self):
        """
        Calculate energy cycling through every position in the lattice
        """
        energy = 0
        for x in range(self.lattice.Lx):
            for y in range(self.lattice.Ly):
                for z in range(self.lattice.Lz):
                    energy = energy + self.energy(x, y, z)/2
        return energy

    def tot_population(self):
        """
        Calculate position cycling through every position in the lattice
        """
        population = 0
        for x in range(self.lattice.Lx):
            for y in range(self.lattice.Ly):
                for z in range(self.lattice.Lz):
                    population = population + self.lattice[x, y, z]
        return population
        
    def energy_if_occupied(self, x, y, z):
        """
        Calculate energy a given point would have *if* that position
        was occupied.
        """

        rmax = self.interaction.rmax
        rmax2 = self.interaction.rmax2
        #sites nearby
        nb = range(-rmax, rmax + 1)
        energy = 0
        d = self.lattice.dim
        s = self.lattice.get_status
        v = self.interaction.V
        
        #sites that we should cycle through
        if d == 3: cyc = it.product(nb, nb, nb)
        elif d == 2: cyc = it.product(nb, nb, [0])
        elif d == 1: cyc = it.product(nb, [0], [0])

        for (i, j, k) in cyc:
            if (i == 0 and j == 0 and k == 0): continue
            xi = x + i
            yj = y + j
            zk = z + k
            if self.lattice.bc == "periodic":
                if xi >= self.lattice.Lx: xi -= self.lattice.Lx
                elif xi < 0: xi += self.lattice.Lx
                
                if yj >= self.lattice.Ly: yj -= self.lattice.Ly
                elif yj < 0: yj += self.lattice.Ly
                
                if zk >= self.lattice.Lz: zk -= self.lattice.Lz
                elif zk < 0: zk += self.lattice.Lz

            elif self.lattice.bc == "free":
                if xi >= self.lattice.Lx or xi < 0: continue
                if yj >= self.lattice.Ly or yj < 0: continue
                if zk >= self.lattice.Lz or zk < 0: continue

            if not self.lattice[xi, yj, zk]: continue
            dist = i**2 + j**2 + k**2
            if dist > rmax2: continue
            energy += v(dist)
        
        return energy

    def energy(self, x, y, z):
        """
        Calculate energy of a given point in the lattice.
        """

        if not self.lattice[x, y, z]:
            return 0
        return self.energy_if_occupied(x, y, z)
    
    def flip(self, x, y, z):
        """
        Method to populate/empty a position. Implements MonteCarlo
        """

        if self.lattice[x, y, z]:
            dn = -1
        else:
            dn = 1
        de = dn * self.energy_if_occupied(x, y, z)
        du = de + self.mu * dn
        if du < 0 or R.random() < np.exp(-du/self.T):
            self.lattice[x, y, z] = not self.lattice[x, y, z]
            self.E+= de
            self.N+= dn
    
    def run(self, Nsteps):
        self.E = self.tot_energy()
        self.N = self.tot_population()
        energy = 0.0
        population = 0.0
        for i in xrange(Nsteps):
            x = int(self.lattice.Lx*R.random())
            y = int(self.lattice.Ly*R.random())
            z = int(self.lattice.Lz*R.random())
            self.flip(x, y, z)
            energy += self.E
            population += self.N
        return energy/Nsteps, population/Nsteps


class Lattice(np.ndarray):
    """
    The lattice is a numpy array with some new methods and attributes.
    The "template" is taken from http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    """
    def __new__(subtype, Lx, buffer=None, offset=0,
                strides=None, order=None, dim=3,
                Ly=None, Lz=None, bc='periodic'):
        """
        Constructor: a square/cubic lattice. Set size and dimension of
        the lattice, as well as the boundary conditions.
        
        Dimensions lower than 3 set "dangling" lengths equal to 1
        """
        if dim > 3:
            raise AttributeError("Dimensions greater than 3 are not supported :( [yet!]")
        if dim < 1:
            raise AttributeError("Dimensions lower than 1 are not supported by reality :( [yet!]")

        if bc != "periodic" and bc != "free":
            raise AttributeError("Only periodic or free boundary conditions are allowed")


        helper = "Dimension is {0}, {1} will not be used"

        if dim < 2:
            if Ly != None: raise Warning(helper.format(dim, "Ly"))
            Ly = 1

        if dim < 3:
            if Lz != None: raise Warning(helper.format(dim, "Lz"))
            Lz = 1
        
        if dim > 1:
            if Ly == None:
                Ly = Lx 
        
        if dim > 2:
            if Lz == None: 
                 Lz = Lx
        
        shape = (Lx, Ly, Lz)
        obj = np.ndarray.__new__(subtype, shape, bool, buffer, offset, strides,
                         order)
        
        obj.Lx = Lx 
        obj.Ly = Ly
        obj.Lz = Lz
        obj.dim = dim
        obj.bc = bc

        return obj

    def  __array_finalize__(self, obj):
        """
        Finalizer: send attributes of obj to class
        """
        if obj is None: return
        self.Lx = getattr(obj, 'Lx', None)
        self.Ly = getattr(obj, 'Ly', None)
        self.Lz = getattr(obj, 'Lz', None)
        self.dim = getattr(obj, 'dim', None)
        self.bc = getattr(obj, 'bc', None)

    def random(self, p = 0.5):
        """
        Fill lattice randomly with probability p.
        """
        for (i, j, k) in it.product(range(self.Lx), range(self.Ly), range(self.Lz)):
            #I'm still surprised this construction works self[...]
            if R.random() < p: 
                self[i, j, k] = True
            else: 
                self[i, j, k] = False
        
        

    def get_status(self, x, y, z):
        if self.bc == "periodic":
            return self[x, y, z]
        
        elif self.bc == "free":
            if x < 0 or x >=self.Lx: return False
            if y < 0 or y >=self.Ly: return False
            if z < 0 or z >=self.Lz: return False
            return self[x, y, z]
        else:
            raise AttributeError("bc = {0} not implemented!".format(self.bc))
        
class Interaction(object):
    """
    Interaction between two lattice sites.
    
    For perfomance issues and typical usage of this code, the
    interaction is defined as a function of r**2

    The rmax is to draw the outer square/cube for each position of
    the lattice. This way, we can only loop among the squares that are
    potentially interacting.

    The potential V is a function that returns interaction given
    distance, for all r less than rmax.
    
    Warning! For r > rmax, V is set to zero.


    """
    def __init__(self, V, rmax):
        """
        Constructor
        """
        self.V = V
        self.rmax = rmax
        self.rmax2 = rmax * rmax
        if type(rmax) != int:
            raise AttributeError("rmax must be an integer, in cell units")
        
    def set_V(self, V):
        """
        Encapsulated setter of the interaction
        """
        self.V = V

    def set_rmax(self, rmax):
        """
        Encapsulated setter of the maximum distance
        """
        self.rmax = rmax
        self.rmax2 = rmax * rmax

    def get_V(self):
        """
        Encapsulated getter of the interaction
        """
        return V

    def get_rmax(self):
        """
        Encapsulated getter of the maximum distance
        """
        return rmax
