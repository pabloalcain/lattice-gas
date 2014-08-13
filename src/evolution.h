#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "mdsys.h"

void whole_lattice(mdsys *sys);
void flip(mdsys *sys, int x, int y, int z);
double energy_occupied(mdsys *sys, int x, int y, int z);
bool get_status(mdsys *sys, int x, int y, int z);
double v(mdsys *sys, int dist);

#endif
