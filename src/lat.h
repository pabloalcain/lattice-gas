#ifndef LAT_H
#define LAT_H
#include <stdbool.h>

struct _lat {
  bool*** site;
  int Lx, Ly, Lz;
  bool periodic, free;
};
typedef struct _lat lat;
#endif