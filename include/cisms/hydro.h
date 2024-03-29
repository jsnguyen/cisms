#ifndef HYDRO_H_
#define HYDRO_H_

#include "stdio.h"
#include "stdlib.h"
#include "omp.h"
#include "cisms/config.h"
#include "cisms/particle.h"
#include "math.h"

#define PI 3.14159265359
#define G 6e-3

void hydro_euler(int particle_index, particle *ps, int n_ps, double smooth_len);
void gravity(int particle_index, particle *ps, int n_ps, double max_len, double strength);
void calc_density(int particle_index, particle *ps, int n_ps, double smooth_len, double prop_const, double poly_index);

double calc_distance(particle *a, particle *b);
double polytrope(double k, double density, double n);
double quintic_spline(double r, double smoothing_length);
double quintic_spline_gradient(int index, double *pos, double smoothing_length);

void calc_new_acc(particle *ps, int n_ps, double smooth_len);
void calc_gravity_acc(particle *ps, int n_ps, double smooth_len, double max_len, double strength);
void half_velocity_verlet_position(particle *ps, int n_ps, double td);
void half_velocity_verlet_velocity(particle *ps, int n_ps, double td);

void check_hard_boundaries(int index, particle *ps, int n_ps, double lower_boundary, double upper_boundary);
void check_periodic_boundaries(int index, particle *ps, int n_ps, double lower_boundary, double upper_boundary);

void simple_drag(particle *ps, int n_ps, double drag_coeff);

#endif
