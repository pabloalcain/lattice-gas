from __future__ import division
import numpy as np
import itertools as it
import random as R

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

        It does look pretty messed up, but I feel like there is
        somewhere we need to mess it up. It's either here or in the
        flip attempt.


        Something I don't like: the splitting in 3 different yet
        *very* similar methods.
        """

        rmax = self.interaction.rmax
        nearby = range(-rmax, rmax + 1)
        energy = 0

        if self.lattice.dim == 3:
            for (i, j, k) in it.product(nearby, nearby, nearby):
                if (i == 0 and j == 0 and k == 0): continue
                
                if not self.lattice.get_status(x + i, y + j, z + k):
                    continue
                dist = np.sqrt(i**2 + j**2 + k**2)
                if dist > rmax: continue
                energy += self.interaction.V(dist)

        if self.lattice.dim == 2:
            for (i, j) in it.product(nearby, nearby):
                if (i == 0 and j == 0): continue
                
                if not self.lattice.get_status(x + i, y + j, 1):
                    continue
                dist = np.sqrt(i**2 + j**2)
                if dist > rmax: continue
                energy += self.interaction.V(dist)

        if self.lattice.dim == 1:
            for i in nearby:
                if (i == 0): continue
                if not self.lattice.get_status(x + i, 1, 1):
                    continue
                dist = abs(i)
                if dist > rmax: continue
                energy += self.interaction.V(dist)

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
        # Getting advantage of the autocast bool->integer in Python
        dn = - (self.lattice[x, y, z] - 0.5) * 2 
        de = dn * self.energy_if_occupied(x, y, z)
        du = de - self.mu * dn

        if du < 0 or R.random() < np.exp(-du/self.T):
            self.lattice[x, y, z] = not self.lattice[x, y, z]
            self.E+= de
            self.N+= dn
    
    def run(self, Nsteps):
        self.E = self.tot_energy()
        self.N = self.tot_population()
        for i in xrange(Nsteps):
            x = R.randint(0, self.lattice.Lx-1)
            y = R.randint(0, self.lattice.Ly-1)
            z = R.randint(0, self.lattice.Lz-1)
            self.flip(x, y, z)
            if (i*100 % Nsteps) == 0:
                print self.E

    
class Lattice(np.ndarray):
    """
    The lattice is a numpy array with some new methods and attributes.
    The "template" is taken from http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    """
    def __new__(subtype, Lx, buffer=None, offset=0,
                strides=None, order=None, dim=3, Ly=None, Lz=None, bc='periodic'):
        """
        Constructor: a square/cubic lattice. Set size and dimension of
        the lattice, as well as the boundary conditions.
        
        Dimensions lower than 3 set "dangling" lengths equal to 1
        """
        if dim > 3:
            raise AttributeError("Dimensions greater than 3 are not supported :(")
        if dim != 3:
            raise AttributeError("Dimensions other than 3 are not supported :( [yet!]")
        if bc != "periodic" and bc != "free":
            raise AttributeError("Only periodic or free boundary conditions are allowed")


        helper = "Dimension is {0}, {1} will not be used"

        if dim < 2 and Ly!= None:
            Ly = 1
            raise Warning(helper.format(dim, "Ly"))

        if dim < 3 and Lz!= None:
            Lz = 1
            raise Warning(helper.format(dim, "Lz"))
        
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
        Not done yet, since I don't know how to cleanly set these values in the lattice
        """
        for (i, j, k) in it.product(range(self.Lx), range(self.Ly), range(self.Lz)):
            #I'm still surprised this construction works self[...]
            if R.random() < p: 
                self[i, j, k] = True
            else: 
                self[i, j, k] = False
        
        

    def get_status(self, x, y, z):
        if self.bc == "periodic":
            return self[x % self.Lx, y % self.Ly, z % self.Lz]
        
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
    
    The rmax is to draw the outer square/cube for each position of the
    lattice. This way, we can only loop among the squares that are
    potentially interacting.

    The potential V is a function that returns interaction given distance, for all
    r less than rmax.
    
    Warning! For r > rmax, V is set to zero.
    """
    def __init__(self, V, rmax):
        """
        Constructor
        """
        self.V = V
        self.rmax = rmax
        
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
