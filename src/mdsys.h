#ifndef MDSYS_H
#define MDSYS_H
#include "pot.h"
#include "lat.h"

struct _mdsys {
  pot interaction;
  lat lattice;
  double T, mu;
  double E, N;
};
typedef struct _mdsys mdsys;
#endif
