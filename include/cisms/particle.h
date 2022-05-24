#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

typedef struct{
    int id;
    double mass;
    double diam;
    double density;
    double pressure;
    double *pos;
    double *vel;
    double *acc;

} particle;

void particle_init(particle *p);
particle *multi_particle_init(int n_ps);
void multi_particle_free(particle *ps, int n_ps);

void particle_set(particle *p, int id, double mass, double diam, double pos[3], double vel[3]);

void particle_set_zero_acc(particle *p);

void particle_copy(particle *dest, particle *src);

void particle_print(particle b);

void particle_write(particle b, const char* fn);

void particle_write_binary(particle p, const char* fn);

#endif
