#ifndef HYDRO_H_
#define HYDRO_H_

#include "stdio.h"
#include "stdlib.h"
#include "cisms/config.h"
#include "cisms/particle.h"
#include "math.h"

#define PI 3.14159265359
#define POLYTROPIC_INDEX 1.5 // ideal gas
#define PROP_CONST 0.1
#define SMOOTHING_LENGTH 0.05
#define G 6e-3

void hydro_euler(int particle_index, particle *ps, int n_ps);
void gravity(int particle_index, particle *ps, int n_ps);
void calc_density(int particle_index, particle *ps, int n_ps);

double calc_distance(particle *a, particle *b);
double polytrope(double k, double density, double n);
double quintic_spline(double r, double smoothing_length);
double quintic_spline_gradient(int index, double *pos, double smoothing_length);

void velocity_verlet(int particle_index, particle *ps, int n_ps, config conf);

#endif
