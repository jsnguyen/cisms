#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "time.h"
#include "omp.h"
#include "cisms/config.h"
#include "cisms/particle.h"
#include "cisms/hydro.h"

static inline double get_rand_between(double a, double b) {
    if (a == b) {
        return a;
    } else {
        assert(a<b);
        return a+((double)rand() / (double)RAND_MAX)*(b-a);
    }
}

void rand_pos(particle *ps, int offset, int n_ps, double x_low, double x_high, double y_low, double y_high, double z_low, double z_high){
    for(int i=offset; i<n_ps; i++){
        ps[i].pos[0] = get_rand_between(x_low, x_high);
        ps[i].pos[1] = get_rand_between(y_low, y_high);
        //ps[i].pos[2] = get_rand_between(z_low, z_high);
    }
}

void rand_vel(particle *ps, int offset, int n_ps, double v_low, double v_high){
    for(int i=offset; i<n_ps; i++){
        ps[i].vel[0] = get_rand_between(v_low, v_high);
        ps[i].vel[1] = get_rand_between(v_low, v_high);
        //ps[i].vel[2] = get_rand_between(v_low, v_high);
    }
}


int main(int argc, char* argv[]){

    char *config_filename = "ism.config";

    config *conf;
    conf = config_create();
    config_read(conf,config_filename);
    config_print(conf);

    particle *ps;

    double x_low  = 0.0;
    double x_high = 1.0;
    double y_low  = 0.0;
    double y_high = 1.0;
    double z_low  = 0.0;
    double z_high = 1.0;

    int np = conf->np;
    ps = multi_particle_init(np);

    rand_pos(ps, 0, np, x_low, x_high, y_low, y_high, z_low, z_high);

    double v_low  = -0.01;
    double v_high = 0.01;
    rand_vel(ps, 0, np, v_low, v_high);

    for (int i=0; i<np; i++){
        ps[i].mass = 1.0;
    }

    for(int i=0; i<np; i++){
        calc_density(i, ps, np, conf->smooth_len, conf->prop_const, conf->poly_index);
    }

    FILE *f;
    f = fopen(conf->fn,"w");
    fclose(f);

    for(int i=0; i<np; i++){
        particle_write_binary(ps[i], conf->fn);
    }

    double gravity_strength = 1e-6;

    double start_t = omp_get_wtime();
    double elaps_t;
    double begin_loop_t;
    double cpu_t, loop_t;
    double loop_ravg=0;
     
    printf("Beginning main simulation loop.\n");
    for(int i=0; i < conf->n_iter; i++){

#pragma omp parallel shared(ps, i)
        {

        begin_loop_t = omp_get_wtime();

        calc_new_acc(ps, np, conf->smooth_len);
        calc_gravity_acc(ps, np, conf->smooth_len, 1, gravity_strength);
        half_velocity_verlet_position(ps, np, conf->td);

#pragma omp for
        for(int j=0; j<np; j++){
            calc_density(j, ps, np, conf->smooth_len, conf->prop_const, conf->poly_index);
        }

        calc_new_acc(ps, np, conf->smooth_len);
        calc_gravity_acc(ps, np, conf->smooth_len, 1, gravity_strength);
        half_velocity_verlet_velocity(ps, np, conf->td);

        check_periodic_boundaries(0, ps, np, 0, 1);
        check_periodic_boundaries(1, ps, np, 0, 1);
        //check_periodic_boundaries(2, ps, np, 0, 1);

        simple_drag(ps, np, conf->drag_coeff);

        }

        for(int j=0; j<np; j++){
            particle_write_binary(ps[j], conf->fn);
        }

        elaps_t = omp_get_wtime();
        cpu_t = (double) (elaps_t - start_t);
        loop_t = (double) (elaps_t - begin_loop_t);
        loop_ravg = (loop_t + (i * loop_ravg)) / (i + 1);

        printf("\r[%.2f%%] || elapsed: %.3fs ||  est time comp: %.3fs || avg loop time: %.3fs ||", (double) 100*(i+1)/conf->n_iter, cpu_t, loop_ravg*(conf->n_iter-1-i), loop_ravg);

        fflush(stdout);

    }

    printf("\n");

    multi_particle_free(ps, np);
    config_free(conf);
    
    printf("Done!\n");
    return 0;
}
