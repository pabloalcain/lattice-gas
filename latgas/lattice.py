"""
Module for lattice gas model. It has inside 3 different classes:

Lattice
Potential
System

that are the basics of the model. System can be linked to different
instances of Lattice and Potential
"""
from __future__ import division
import numpy as np
import random as R
import ctypes as C
from math import ceil

CLIB = C.CDLL("liblatgas.so")

        

class Lattice(C.Structure):
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
            raise AttributeError(("Dimensions greater than 3 are "
                                  "not supported :( [yet!]"))
        if dim < 1:
            raise AttributeError(("Dimensions lower than 1 are not "
                                 "supported by reality :( [yet!]"))

        if bc != "periodic" and bc != "free":
            raise AttributeError(("Only periodic or free boundary "
                                  "conditions are allowed"))

        helper = "Dimension is {0}, {1} will not be used"
        
        if dim == 1:
            if Ly != None:
                raise Warning(helper.format(dim, "Ly"))
            if Lz != None:
                raise Warning(helper.format(dim, "Lz"))
            Ly = 1
            Lz = 1
        
        if dim == 2:
            if Ly == None:
                Ly = Lx
            if Lz != None:
                raise Warning(helper.format(dim, "Lz"))
            Lz = 1

        if dim == 3:
            if Ly == None:
                Ly = Lx
            if Lz == None:
                Lz = Lx

        self.Lx = Lx 
        self.Ly = Ly
        self.Lz = Lz
        
        self.site = ( C.c_bool * ( Lx * Ly * Lz ) ) ()
        self.dim = dim
        self.bc = bc
        
        self.periodic = (bc == "periodic")
        self.free = (bc == "free")
        super(Lattice, self).__init__()


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
        """
        Get status of site. It either belongs to the lattice or we
        need to apply the boundary conditions

        Wrapper to c function
        """
        CLIB.get_status(self, x, y, z)
        
class Potential(C.Structure):
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
        self.rmax2 = self.rmax**2
        self.rcore = rcore
        self.invdr = (npoints - 1) / (rcut - rcore)
        self.inter = inter
        self.npoints = npoints
        self.V = ( C.c_double * (self.npoints) ) ()
        self.r = ( C.c_double * (self.npoints) ) ()
        self.create_tables()
        super(Potential, self).__init__()

    def create_tables(self):
        """
        Create the lookup tables in r and V pointers
        """
        _r = np.linspace(self.rcore, 
                         self.rcut, self.npoints)
        for i in range(self.npoints):
            self.r[i] = _r[i]
            self.V[i] = self.inter(_r[i])
        
    def set_V(self, inter):
        """
        Encapsulated setter of the interaction
        """
        self.inter = inter
        self.create_tables()
        
        
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


class System(C.Structure):
    """
    Structure for the lattice gas system
    """

    _fields_ = [("interaction", Potential),
                ("lattice", Lattice),
                ("T", C.c_double),
                ("mu", C.c_double),
                ("E", C.c_double),
                ("N", C.c_double),
                ("random", C.c_int)]
    
    def __init__(self, lattice, interaction, T = 1.0, mu = 0.0, random = None):
        """
        Constructor: give the lattice and thermodynamic info
        """
        self.lattice = lattice
        self.interaction = interaction
        self.T = T
        self.mu = mu
        self.E = self.tot_energy()
        self.N = self.tot_population()
        if random == None: 
            random = 1000000 * R.random()
        self.random = int(random)
        super(System, self).__init__()
    
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
        return CLIB.energy_if_occupied(C.byref(self), x, y, z)
        
        
    def energy(self, x, y, z):
        """
        Calculate energy of a given point in the lattice.
        
        Wrapper to c function.
        """
        return CLIB.energy(C.byref(self), x, y, z)
    
    def flip(self, x, y, z):
        """
        Method to populate/empty a position. Implements MonteCarlo
        """

        CLIB.flip(C.byref(self), x, y, z)
    
    def run(self, Nsteps):
        """
        Dummy run method, inform mean values
        """
        
        #self.E = self.tot_energy()
        #self.N = self.tot_population()
        energy = 0.0
        population = 0.0
        for _i in xrange(Nsteps):
            self.whole_lattice()
            energy += self.E
            population += self.N
        return energy/Nsteps, population/Nsteps

    def crun(self, Nsteps):
        """
        Wrapper to c function. Useful to find overheads

        Ideal: Overheads << Execution for DT ~ Correlation time
        """
        
        CLIB.run(C.byref(self), Nsteps)
        


    def whole_lattice(self):
        """
        Attempt to flip Lx * Ly * Lz positions
        
        Wrapper to c function
        """
        
        CLIB.whole_lattice(C.byref(self))
