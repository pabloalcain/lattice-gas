#ifndef LATTICE_H
#define LATTICE_H

struct _lattice {
  bool*** site;
  int Lx, Ly, Lz;
  bool periodic, free;
};
typedef struct _lattice lattice;
#endif
