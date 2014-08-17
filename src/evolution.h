#ifndef EVOLUTION_H
#define EVOLUTION_H

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "sys.h"

void run(System *sys, int nsteps);
void whole_lattice(System *sys);
void flip(System *sys, int x, int y, int z);
double energy_if_occupied(System *sys, int x, int y, int z);
double energy(System *sys, int x, int y, int z);
bool get_status(System *sys, int x, int y, int z);
double v(System *sys, int dist);

#endif
