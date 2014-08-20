#include "evolution.h"

void run(System *sys, int nsteps)
{
  int i;
  for (i = 0; i < nsteps; i++) {
    whole_lattice(sys);
  }
}

void whole_lattice(System *sys)
{
  if (sys->random != 0) {
    srand(sys->random);
    sys->random = 0;
  }
  int i = 0;
  int N;
  int x, y, z;
  int Lx, Ly, Lz;
  
  Lx = sys->lattice.Lx;
  Ly = sys->lattice.Ly;
  Lz = sys->lattice.Lz;
  
  N = Lx * Ly * Lz;
  /* of course, we don't go through every spin,
     it's just so we do as many montecarlo steps as
     sites in the lattice */
  for (i = 0; i < N; i++) {
    x = rand() % Lx;
    y = rand() % Ly;
    z = rand() % Lz;
    flip(sys, x, y, z);
  }
}

void flip(System *sys, int x, int y, int z)
{
  int dn, de, du;
  int Lx, Ly;
  Lx = sys->lattice.Lx;
  Ly = sys->lattice.Ly;

  if (sys->lattice.site[x + y * Lx + z * Lx * Ly])
    dn = -1;
  else 
    dn = 1;
  de = dn * energy_if_occupied(sys, x, y, z);
  du = de + sys->mu * dn;
  if (du < 0 || 
      rand() < RAND_MAX * exp(-du/sys->T)) {
    sys->lattice.site[x + y * Lx + z * Lx * Ly] = 
      !sys->lattice.site[x + y * Lx + z * Lx * Ly];
    sys->E += de;
    sys->N += dn;
  }
}

double energy_if_occupied(System *sys, int x, int y, int z)
{
  int i, j, k;
  int imax, jmax, kmax;
  double rcut, rcutsq;
  int dist;
  double energy;
  rcut = sys->interaction.rcut;
  imax = rcut;
  jmax = 0;
  kmax = 0;
  energy = 0;
  rcutsq = rcut * rcut;
  if (sys->lattice.dim > 1) 
    jmax = rcut;
  if (sys->lattice.dim > 2) 
    kmax = rcut;
  
#pragma omp parallel for
  for (i = -imax; i <= imax; i++) {
    for (j = -jmax; j <= jmax; j++) {
      for (k = -kmax; k <= kmax; k++) {
        if (i == 0 && j == 0 && k == 0)
	  continue;
	dist = i * i + j * j + k * k;
	if (dist > rcutsq)
	  continue;
        if (!get_status(sys, x + i, y + j, z + k)) 
	  continue;
	dist = sqrt(dist);
        energy += v(sys, dist);
      }
    }
  }
  return energy;
}

double energy(System *sys, int x, int y, int z)
{
  int Lx, Ly;
  Lx = sys->lattice.Lx;
  Ly = sys->lattice.Ly;
  if (!sys->lattice.site[x + y * Lx + z * Lx * Ly])
    return 0;
  return energy_if_occupied(sys, x, y, z);
}
   
bool get_status(System *sys, int x, int y, int z)
{
  int Lx, Ly, Lz;
  Lx = sys->lattice.Lx;
  Ly = sys->lattice.Ly;
  Lz = sys->lattice.Lz;
  if (sys->lattice.periodic) {
    if (x >= Lx)
      x -= Lx;
    else if (x < 0)
      x += Lx;

    if (y >= Ly)
      y -= Ly;
    else if (y < 0)
      y += Ly;

    if (z >= Lz)
      z -= Lz;
    else if (z < 0)
      z += Lz;
  }
  
  if (sys->lattice.free) {
    if (x >= Lx || x < 0)
      return false;
    if (y >= Ly || y < 0)
      return false;
    if (z >= Lz || z < 0)
      return false;
  }
  return sys->lattice.site[x + y * Lx + z * Lx * Ly];
}

double v(System *sys, int dist)
{
  double invdr, rcore;
  int idx;
  
  invdr = sys->interaction.invdr;
  rcore = sys->interaction.rcore;
  
  idx = (int) (dist - rcore) * invdr;
  return sys->interaction.V[idx];
}
