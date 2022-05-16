#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

typedef struct{
    int id;
    double mass;
    double diam;
    double *pos;
    double *vel;

    double *x,*y,*z;
    double *vx,*vy,*vz;
} particle;

void particle_set(particle *p, int id, double mass, double diam, double pos[3], double vel[3]);

void particle_print(particle b);

void particle_write(particle b, const char* fn);


#endif
