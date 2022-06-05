#include "cisms/hydro.h"

void hydro_euler(int particle_index, particle *ps, int n_ps, double smooth_len){

    double dr[3];
    double new_acc[] = {0,0,0};

    for(int i=0; i<n_ps; i++){

        if(i == particle_index){
            continue;
        }

        dr[0] = ps[particle_index].pos[0] - ps[i].pos[0];
        dr[1] = ps[particle_index].pos[1] - ps[i].pos[1];
        dr[2] = ps[particle_index].pos[2] - ps[i].pos[2];

        new_acc[0] += ps[particle_index].mass * ( ps[particle_index].pressure / (ps[particle_index].density*ps[particle_index].density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(0, dr, smooth_len);
        new_acc[1] += ps[particle_index].mass * ( ps[particle_index].pressure / (ps[particle_index].density*ps[particle_index].density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(1, dr, smooth_len);
        new_acc[2] += ps[particle_index].mass * ( ps[particle_index].pressure / (ps[particle_index].density*ps[particle_index].density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(2, dr, smooth_len);

    }

    ps[particle_index].acc[0] += -new_acc[0];
    ps[particle_index].acc[1] += -new_acc[1];
    ps[particle_index].acc[2] += -new_acc[2];

}

void gravity(int particle_index, particle *ps, int n_ps, double max_len, double strength){

    double dr[3];
    double r;

    for(int i=0; i<n_ps; i++){

        if(i == particle_index){
            continue;
        }

        ps[particle_index].pos[0] = ps[particle_index].pos[0] - floor(ps[particle_index].pos[0] / max_len) * max_len;
        ps[particle_index].pos[1] = ps[particle_index].pos[1] - floor(ps[particle_index].pos[1] / max_len) * max_len;
        ps[particle_index].pos[2] = ps[particle_index].pos[2] - floor(ps[particle_index].pos[2] / max_len) * max_len;

        dr[0] = ps[i].pos[0] - ps[particle_index].pos[0];
        dr[1] = ps[i].pos[1] - ps[particle_index].pos[1];
        dr[2] = ps[i].pos[2] - ps[particle_index].pos[2];

        dr[0] = dr[0] - (int) (dr[0] / max_len) * max_len;
        dr[1] = dr[1] - (int) (dr[1] / max_len) * max_len;
        dr[2] = dr[2] - (int) (dr[2] / max_len) * max_len;

        r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

        if(r < 0.1){
            r = r*10;
        }

        ps[particle_index].acc[0] += dr[0]/r * strength * ps[particle_index].mass * ps[i].mass / (r*r);
        ps[particle_index].acc[1] += dr[1]/r * strength * ps[particle_index].mass * ps[i].mass / (r*r);
        ps[particle_index].acc[2] += dr[2]/r * strength * ps[particle_index].mass * ps[i].mass / (r*r);

    }


}

void calc_density(int particle_index, particle *ps, int n_ps, double smooth_len, double prop_const, double poly_index){

    double r;
    double weight;

    double density = 0;
    for(int i=0; i<n_ps; i++){

        r = calc_distance(&ps[particle_index], &ps[i]);
        weight = quintic_spline(r, smooth_len);

        density += ps[i].mass*weight;

    }

    ps[particle_index].density = density;
    ps[particle_index].pressure = polytrope(prop_const, density, poly_index);

}

double calc_distance(particle *a, particle *b){
    return sqrt(((a->pos[0] - b->pos[0])*(a->pos[0] - b->pos[0]))
              + ((a->pos[1] - b->pos[1])*(a->pos[1] - b->pos[1]))
              + ((a->pos[2] - b->pos[2])*(a->pos[2] - b->pos[2])));
}

double polytrope(double k, double density, double n){
    return k*pow(density, (n+1)/n);
}

double quintic_spline(double r, double smoothing_length){
    double q = r/smoothing_length;

    if ((q >= 0) && (q < 1)) {
        return (7.0/(478.0*PI))*(pow((3-q),5) - 6*pow((2-q),5) + 15*pow((1-q),5));
        
    } else if((q >= 1) && (q < 2)) {
        return (7.0/(478.0*PI))*(pow((3-q),5) - 6*pow((2-q),5));

    } else if((q >= 2) && (q < 3)) {
        return (7.0/(478.0*PI))*pow((3-q),5);

    } else {
        return 0;
    }
}

double quintic_spline_gradient(int index, double *pos, double smoothing_length){
    double q = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2])/smoothing_length;
    double coeff = -5*(pos[index]/smoothing_length)/q;

    if ((q >= 0) && (q < 1)) {
        return (1.0/120.0)*coeff*(pow((3-q),4) - 6*pow((2-q),4) + 15*pow((1-q),4));
        
    } else if((q >= 1) && (q < 2)) {
        return (7.0/(478.0*PI))*coeff*(pow((3-q),4) - 6*pow((2-q),4));

    } else if((q >= 2) && (q < 3)) {
        return (1.0/(120.0*PI))*coeff*pow((3-q),4);

    } else {
        return 0;

    }
}

void calc_new_acc(particle *ps, int n_ps, double smooth_len){

#pragma omp for
    for(int i=0; i<n_ps; i++){

        particle_set_zero_acc(&ps[i]);
        hydro_euler(i, ps, n_ps, smooth_len);
        
    }

}

void calc_gravity_acc(particle *ps, int n_ps, double smooth_len, double max_len, double strength){

#pragma omp for
    for(int i=0; i<n_ps; i++){

        gravity(i, ps, n_ps, max_len, strength);

    }
}

void half_velocity_verlet_position(particle *ps, int n_ps, double td){

#pragma omp for
    for(int i=0; i<n_ps; i++){

        ps[i].pos[0] = ps[i].pos[0] + ps[i].vel[0]*td + 0.5*ps[i].acc[0]*td*td;
        ps[i].pos[1] = ps[i].pos[1] + ps[i].vel[1]*td + 0.5*ps[i].acc[1]*td*td;
        ps[i].pos[2] = ps[i].pos[2] + ps[i].vel[2]*td + 0.5*ps[i].acc[2]*td*td;

    }

}

void half_velocity_verlet_velocity(particle *ps, int n_ps, double td){

#pragma omp for
    for(int i=0; i<n_ps; i++){

        ps[i].vel[0] = ps[i].vel[0] + 0.5*(ps[i].acc[0]+ps[i].acc[0])*td;
        ps[i].vel[1] = ps[i].vel[1] + 0.5*(ps[i].acc[1]+ps[i].acc[1])*td;
        ps[i].vel[2] = ps[i].vel[2] + 0.5*(ps[i].acc[2]+ps[i].acc[2])*td;

    }

}

void check_hard_boundaries(int index, particle *ps, int n_ps, double lower_boundary, double upper_boundary){

#pragma omp for
    for(int i=0; i<n_ps; i++){

        if (ps[i].pos[index] <= lower_boundary){
            ps[i].vel[index] = fabs(ps[i].vel[index]);
        }

        else if (ps[i].pos[index] >= upper_boundary){
            ps[i].vel[index] = -fabs(ps[i].vel[index]);
        }

    }
}

void check_periodic_boundaries(int index, particle *ps, int n_ps, double lower_boundary, double upper_boundary){

    double size = upper_boundary-lower_boundary;
#pragma omp for
    for(int i=0; i<n_ps; i++){

        if (ps[i].pos[index] <= lower_boundary){
            ps[i].pos[index] = ps[i].pos[index] + size;
        }

        else if (ps[i].pos[index] >= upper_boundary){
            ps[i].pos[index] = ps[i].pos[index] - size;
        }

    }
}


void simple_drag(particle *ps, int n_ps, double drag_coeff){

#pragma omp for
    for(int i=0; i<n_ps; i++){

        ps[i].vel[0] -= drag_coeff*ps[i].vel[0];
        ps[i].vel[1] -= drag_coeff*ps[i].vel[1];
        ps[i].vel[2] -= drag_coeff*ps[i].vel[2];

    }

}
