#include "evolution.h"

void whole_lattice(mdsys *sys)
{
  int i;
  int N;
  int x, y, z;

  N = sys->Lx * sys-> Ly * sys-> Lz;
  for (i = 0; i++; i < N) {
    x = std::rand() % Lx;
    y = std::rand() % Ly;
    z = std::rand() % Lz;
    flip(sys, x, y, z);
  }
}

void flip(mdsys *sys, int x, int y, int z)
{
  int dn, de, du;
  if (sys->lattice[x][y][z]) dn = -1;
  else dn = 1;
  de = dn * energy_occupied(sys, x, y, z);
  du = de + self.mu * dn;
  if (du < 0 || std::rand() < RAND_MAX * std::exp(-du/sys->T)) {
    sys->lattice[x][y][z] = !sys->lattice[x][y][z];
    sys->E += de;
    sys->N += dn;
  }
}

double energy_occupied(mdsys *sys, int x, int y, int z)
{
  int i, j, k;
  int imax, jmax, kmax;
  int dist;
  double energy;
  imax = rmax;
  jmax = 0;
  kmax = 0;
  if (sys->dim > 1) jmax = rmax;
  if (sys->dim > 2) kmax = rmax;
  for (i = -imax; i++; i <= imax) {
    for (j = -jmax; j++; j <= jmax) {
      for (k = -kmax; k++; k <= kmax) {
        if (i == 0 && j == 0 && k == 0) continue;
        if !(sys->lattice.get_status(sys, x + i, y + j, z + k) continue;
        dist = i * i + j * j + k * k;
        energy += sys->interaction.V(dist);
      }
    }
  }
  return energy;
}

