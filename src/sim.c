#include "sim.h"

int main(int argc, char **argv)
{
  System sys;
  int i, nsites;  

  sys.interaction.r = (double *)malloc(5*sizeof(double));
  sys.interaction.V = (double *)malloc(5*sizeof(double));

    for ( i = 0; i < 5; i++) {
      sys.interaction.r[i] = i/4.0;
      sys.interaction.V[i] = -4.0;
    }

  sys.interaction.rcore = 0.0;
  sys.interaction.rcut = 1.0;
  sys.interaction.rmax = 1;
  sys.interaction.invdr = 0.25;

  sys.lattice.Lx = 30;
  sys.lattice.Ly = 30;
  sys.lattice.Lz = 1;
  sys.lattice.dim = 2;
  sys.lattice.periodic = true;
  sys.lattice.free = false;
  nsites = sys.lattice.Lx * sys.lattice.Ly * sys.lattice.Lz;
  sys.lattice.site = (bool *)malloc(nsites*sizeof(bool));

  sys.T = 20.0;
  sys.mu = 8.0;
  sys.E = 0;
  sys.N = 0;


  for ( i = 0; i < nsites; i++) {
    sys.lattice.site[i] = false;
  }
  
  run(&sys, 10000);
  return 0;
}
