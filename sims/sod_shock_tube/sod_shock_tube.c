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
        ps[i].pos[2] = get_rand_between(z_low, z_high);
    }
}

void uniform_density_pos(particle *ps, int n_ps, int offset, double smooth_len, double number_density, double x_low, double x_high, double y_low, double y_high){

    double dx = x_high - x_low;
    double dy = y_high - y_low;

    int n_x = cbrt(number_density)*dx;
    int n_y = cbrt(number_density)*dy;

    double spacing = 1/cbrt(number_density);

    for(int j=0; j<n_y; j++){
        for(int i=0; i<n_x; i++){
            ps[offset+i+j*n_x].pos[0] = x_low + spacing/2.0 + i*spacing;
            ps[offset+i+j*n_x].pos[1] = y_low + spacing/2.0 + j*spacing;
        }
    }

}

int main(int argc, char* argv[]){

    char *config_filename = "sod_shock_tube.config";

    config *conf;
    conf = config_create();
    config_read(conf,config_filename);
    config_print(conf);

    particle *ps;

    double l_number_density = pow(8,5  );
    double r_number_density = pow(8,5-2);

    double l_x_low  = 0.0;
    double l_x_high = 5.0;
    double l_y_low  = 0.0;
    double l_y_high = 1.0;

    double r_x_low  =  5.0;
    double r_x_high = 10.0;
    double r_y_low  =  0.0;
    double r_y_high =  1.0;

    int n_left  = cbrt(l_number_density) * (l_x_high - l_x_low) * cbrt(l_number_density) * (l_y_high - l_y_low);
    int n_right = cbrt(r_number_density) * (r_x_high - r_x_low) * cbrt(r_number_density) * (r_y_high - r_y_low);

    int np = n_left+n_right;
    ps = multi_particle_init(np);

    uniform_density_pos(ps, np, 0 , conf->smooth_len, l_number_density, l_x_low, l_x_high, l_y_low, l_y_high);
    uniform_density_pos(ps, np, n_left, conf->smooth_len, r_number_density, r_x_low, r_x_high, r_y_low, r_y_high);

    for (int i=0; i<np; i++){
        ps[i].mass = 1.0/np;
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
        half_velocity_verlet_position(ps, np, conf->td);

#pragma omp for
        for(int j=0; j<np; j++){
            calc_density(j, ps, np, conf->smooth_len, conf->prop_const, conf->poly_index);
        }

        calc_new_acc(ps, np, conf->smooth_len);
        half_velocity_verlet_velocity(ps, np, conf->td);

        check_hard_boundaries(0, ps, np, 0, 10);
        check_hard_boundaries(1, ps, np, 0, 1);
        check_hard_boundaries(2, ps, np, 0, 1);

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
