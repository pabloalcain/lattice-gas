#include "evolution.h"

void whole_lattice(mdsys *sys)
{
  int i;
  int N;
  int x, y, z;
  int Lx, Ly, Lz;
  
  Lx = sys->lattice.Lx;
  Ly = sys->lattice.Ly;
  Lz = sys->lattice.Lz;
  
  N = Lx * Ly * Lz;
  /* of course, we don't go through every spin,
   it's just so we do as many montecarlo steps as
   sites in the lattice*/
  for (i = 0; i++; i < N) {
    x = rand() % Lx;
    y = rand() % Ly;
    z = rand() % Lz;
    flip(sys, x, y, z);
  }
}

void flip(mdsys *sys, int x, int y, int z)
{
  int dn, de, du;
  if (sys->lattice.site[x][y][z]) 
    dn = -1;
  else 
    dn = 1;
  de = dn * energy_occupied(sys, x, y, z);
  du = de + sys->mu * dn;
  if (du < 0 || 
      rand() < RAND_MAX * exp(-du/sys->T)) {
    sys->lattice.site[x][y][z] = !sys->lattice.site[x][y][z];
    sys->E += de;
    sys->N += dn;
  }
}

double energy_occupied(mdsys *sys, int x, int y, int z)
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
  rcutsq = rcut * rcut;
  if (sys->dim > 1) 
    jmax = rcut;
  if (sys->dim > 2) 
    kmax = rcut;
  for (i = -imax; i++; i <= imax) {
    for (j = -jmax; j++; j <= jmax) {
      for (k = -kmax; k++; k <= kmax) {
        if (i == 0 && j == 0 && k == 0)
	  continue;
        if (!get_status(sys, x + i, y + j, z + k)) 
	  continue;
	dist = i * i + j * j + k * k;
	if (dist > rcutsq)
	  continue;
	dist = sqrt(dist);
        energy += v(sys, dist);
      }
    }
  }
  return energy;
}

bool get_status(mdsys *sys, int x, int y, int z)
{
  if (sys->lattice.periodic) {
    if (x >= sys->lattice.Lx)
      x -= sys->lattice.Lx;
    else if (x < 0)
      x += sys->lattice.Lx;

    if (y >= sys->lattice.Ly)
      y -= sys->lattice.Ly;
    else if (y < 0)
      y += sys->lattice.Ly;

    if (z >= sys->lattice.Lz)
      z -= sys->lattice.Lz;
    else if (z < 0)
      z += sys->lattice.Lz;
  }
  
  if (sys->lattice.free) {
    if (x >= sys->lattice.Lx || x < 0)
      return false;
    if (y >= sys->lattice.Ly || y < 0)
      return false;
    if (z >= sys->lattice.Lz || z < 0)
      return false;
  }
  return sys->lattice.site[x][y][z];
}

double v(mdsys *sys, int dist)
{
  double invdr, rcore;
  int idx;
  
  invdr = sys->interaction.invdr;
  rcore = sys->interaction.rcore;
  
  idx = (int) (dist - rcore) * invdr;
  return sys->interaction.V[idx];
}
