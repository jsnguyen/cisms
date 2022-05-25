#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "time.h"
#include "cisms/config.h"
#include "cisms/particle.h"
#include "cisms/hydro.h"

static inline double get_rand_between(double a, double b) {
    assert(a<b);
    return a+((double)rand() / (double)RAND_MAX)*(b-a);
}

void rand_pos_square(particle *ps, int n_ps){
    for(int i=0; i<n_ps; i++){
        ps[i].pos[0] = get_rand_between(-0.9,0.9);
        ps[i].pos[1] = get_rand_between(-0.9,0.9);
        //ps[i].pos[2] = get_rand_between(-0.9,0.9);
    }
}

int main(int argc, char* argv[]){
    if(argc != 2){
        printf("Usage: ./constrained_cube.exe <filename.config>\n");
        return 1;
    }

    config *conf;
    conf = config_create();
    config_read(conf,argv[1]);
    config_print(conf);

    particle *ps;

    ps = multi_particle_init(conf->np);

    rand_pos_square(ps, conf->np);

    for (int i=0; i<conf->np; i++){
        ps[i].mass = 2.0/conf->np;
    }

    for(int i=0; i<conf->np; i++){
        calc_density(i, ps, conf->np, conf->smooth_len, conf->prop_const, conf->poly_index);
    }

    for(int i=0; i<conf->np; i++){
        particle_write(ps[i], conf->fn);
    }

    FILE *f;
    f = fopen(conf->fn,"w");
    fclose(f);

    clock_t start_t = clock();
    clock_t elaps_t;
    clock_t begin_loop_t;
    double cpu_t, loop_t;
    double loop_ravg=0;

    printf("Beginning main simulation loop.\n");

    for(int i=0; i < conf->n_iter; i++){

        begin_loop_t = clock();

        calc_new_acc(ps, conf->np, conf->smooth_len);
        half_velocity_verlet_position(ps, conf->np, conf->td);

        for(int i=0; i<conf->np; i++){
            calc_density(i, ps, conf->np, conf->smooth_len, conf->prop_const, conf->poly_index);
        }

        calc_new_acc(ps, conf->np, conf->smooth_len);
        half_velocity_verlet_velocity(ps, conf->np, conf->td);

        check_hard_boundaries(0, ps, conf->np, -1, 1);
        check_hard_boundaries(1, ps, conf->np, -1, 1);
        check_hard_boundaries(2, ps, conf->np, -1, 1);

        simple_drag(ps, conf->np, conf->drag_coeff);

        for(int i=0; i<conf->np; i++){
            particle_write_binary(ps[i], conf->fn);
        }

        elaps_t = clock();
        cpu_t = (double) (elaps_t - start_t) / CLOCKS_PER_SEC;
        loop_t = (double) (elaps_t - begin_loop_t) / CLOCKS_PER_SEC;
        loop_ravg = (loop_t + (i * loop_ravg)) / (i + 1);

        printf("\r[%.2f%%] || elapsed: %.3fs ||  est time comp: %.3fs || avg loop time: %.3fs ||", (double) 100*(i+1)/conf->n_iter, cpu_t, loop_ravg*(conf->n_iter-1-i), loop_ravg);

        fflush(stdout);
    }

    printf("\n");

    multi_particle_free(ps, conf->np);
    config_free(conf);
    
    printf("Done!\n");
    return 0;
}
