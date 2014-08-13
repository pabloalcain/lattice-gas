#ifndef MDSYS_H
#define MDSYS_H
#include "pot.h"
#include "lat.h"

struct _mdsys {
  pot interaction;
  int Lx, Ly, Lz;
  lat lattice;
  double T, mu;
  double E, N;
  int dim;
};
typedef struct _mdsys mdsys;
#endif
