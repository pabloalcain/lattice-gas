#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <stdlib.h>
#include "mdsys.h"

void whole_lattice(mdsys *sys);
void flip(mdsys *sys, int x, int y, int z);
double energy_occupied(mdsys *sys, int x, int y, int z);


#endif
