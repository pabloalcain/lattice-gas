from __future__ import division
import numpy as np
import itertools as it
import random as R
import time
import ctypes as C
from math import ceil

clib = C.CDLL("../liblatgas.so")

        

class lat(C.Structure):
    """
    Basic structure for a lattice.
    """
    _fields_ = [("site", C.POINTER(C.c_bool)),
                ("Lx", C.c_int),
                ("Ly", C.c_int),
                ("Lz", C.c_int),
                ("dim", C.c_int),
                ("periodic", C.c_bool),
                ("free", C.c_bool)]

    def __init__(self, Lx, Ly=None, Lz=None, 
                 dim=3, bc='periodic'):
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

        self.Lx = Lx 
        self.Ly = Ly
        self.Lz = Lz
        
        self.site = ( C.c_bool * ( Lx * Ly * Lz ) ) ()
        self.dim = dim
        self.bc = bc
        
        if bc == "periodic":
            self.periodic = True
            self.free = False
        if bc == "free":
            self.periodic = False
            self.free = True


    def random(self, p = 0.5):
        """
        Fill lattice randomly with probability p.
        """
        for i in range(self.Lx * self.Ly * self.Lz):
            if R.random() < p: 
                self.site[i] = True
            else: 
                self.site[i] = False
        
        

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
        
class pot(C.Structure):
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
    _fields_ = [('rmax', C.c_int),
                ('rcut', C.c_double),
                ('rcore', C.c_double),
                ('invdr', C.c_double),
                ('r', C.POINTER(C.c_double)),
                ('V', C.POINTER(C.c_double))]

    def __init__(self, inter, rcut, npoints, rcore = 0.0):
        """
        Constructor.

        inter is the interaction potential between sites
        """
        self.rcut = rcut
        self.rmax = int(ceil(rcut))
        self.rcore = rcore
        self.invdr = (npoints - 1) / (rcut - rcore)
        self.inter = inter
        self.npoints = npoints
        self.create_tables()

    def create_tables(self):
        """
        Create the lookup tables in r and V pointers
        """
        _r = np.linspace(self.rcore, self.rcut, self.npoints)
        _V = self.inter(_r)
        self.V = (C.c_double * (self.npoints) ) ()
        self.r = (C.c_double * (self.npoints) ) ()
        for i in range(self.npoints):
            self.r[i] = _r[i]
            self.V[i] = _V[i]
        
    def set_V(self, inter):
        """
        Encapsulated setter of the interaction
        """
        self.inter = inter
        self.create_tables(inter)
        
        
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
        return self.inter

    def get_rmax(self):
        """
        Encapsulated getter of the maximum distance
        """
        return self.rmax


class mdsys(C.Structure):
    """
    Structure for the lattice gas system
    """

    _fields_ = [("interaction", pot),
                ("lattice", lat),
                ("T", C.c_double),
                ("mu", C.c_double),
                ("E", C.c_double),
                ("N", C.c_double)]
    
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
        Lx = self.lattice.Lx
        Ly = self.lattice.Ly
        Lz = self.lattice.Lz
        for x in range(Lx):
            for y in range(Ly):
                for z in range(Lz):
                    energy += self.energy(x, y, z)/2
        return energy

    def tot_population(self):
        """
        Calculate position cycling through every position in the lattice
        """
        population = 0
        Lx = self.lattice.Lx
        Ly = self.lattice.Ly
        Lz = self.lattice.Lz
        for i in range(Lx * Ly * Lz):
            population += self.lattice.site[i]
        return population
        
    def energy_if_occupied(self, x, y, z):
        """
        Calculate energy a given point would have *if* that position
        was occupied.
        
        Wrapper to c function.

        """
        return clib.energy_if_occupied(C.byref(self), x, y, z)
        
        
    def energy(self, x, y, z):
        """
        Calculate energy of a given point in the lattice.
        
        Wrapper to c function.
        """
        return clib.energy(C.byref(self), x, y, z)
    
    def flip(self, x, y, z):
        """
        Method to populate/empty a position. Implements MonteCarlo
        """

        clib.flip(C.byref(self), x, y, z)
    
    def run(self, Nsteps):
        """
        Dummy run method, inform mean values
        """
        
        self.E = self.tot_energy()
        self.N = self.tot_population()
        energy = 0.0
        population = 0.0
        for i in xrange(Nsteps):
            self.whole_lattice()
            energy += self.E
            population += self.N
        return energy/Nsteps, population/Nsteps

    def whole_lattice(self):
        """
        Attempt to flip Lx * Ly * Lz positions
        
        Wrapper to c function
        """
        
        clib.whole_lattice(C.byref(self))
