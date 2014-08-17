#ifndef SYS_H
#define SYS_H
#include "pot.h"
#include "lat.h"

struct _sys {
  Potential interaction;
  Lattice lattice;
  double T, mu;
  double E, N;
  int random;
};
typedef struct _sys System;
#endif
