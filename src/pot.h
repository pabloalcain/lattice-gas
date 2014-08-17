#ifndef POT_H
#define POT_H

struct _pot {
  int rmax;
  double rcut;
  double rcore;
  double invdr;
  double *r;
  double *V;
};
typedef struct _pot Potential;
#endif
