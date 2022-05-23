#include "cisms/hydro.h"

void hydro_euler(int particle_index, particle *ps, int n_ps, double smooth_len){

    particle *p = &ps[particle_index];

    double dr[3];

    double new_acc[] = {0,0,0};

    for(int i=0; i<n_ps; i++){

        if(i == particle_index){
            continue;
        }

        dr[0] = p->pos[0] - ps[i].pos[0];
        dr[1] = p->pos[1] - ps[i].pos[1];
        dr[2] = p->pos[2] - ps[i].pos[2];

        new_acc[0] += p->mass * ( p->pressure / (p->density*p->density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(0, dr, smooth_len);
        new_acc[1] += p->mass * ( p->pressure / (p->density*p->density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(1, dr, smooth_len);
        new_acc[2] += p->mass * ( p->pressure / (p->density*p->density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(2, dr, smooth_len);

    }

    p->acc[0] += -new_acc[0];
    p->acc[1] += -new_acc[1];
    p->acc[2] += -new_acc[2];

}

void gravity(int particle_index, particle *ps, int n_ps){

    particle *p = &ps[particle_index];

    double dr[3];
    double r;

    for(int i=0; i<n_ps; i++){

        if(i == particle_index){
            continue;
        }

        dr[0] = p->pos[0] - ps[i].pos[0];
        dr[1] = p->pos[1] - ps[i].pos[1];
        dr[2] = p->pos[2] - ps[i].pos[2];

        r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

        p->acc[0] += dr[0]/r * G * p->mass * ps[i].mass / (r*r);
        p->acc[1] += dr[1]/r * G * p->mass * ps[i].mass / (r*r);
        p->acc[2] += dr[2]/r * G * p->mass * ps[i].mass / (r*r);

    }


}

void calc_density(int particle_index, particle *ps, int n_ps, double smooth_len, double prop_const, double poly_index){

    particle *p = &ps[particle_index];
    double r;
    double weight;

    double density = 0;
    for(int i=0; i<n_ps; i++){

        r = calc_distance(p, &ps[i]);
        weight = quintic_spline(r, smooth_len);

        density += ps[i].mass*weight;

    }

    p->density = density;
    p->pressure = polytrope(prop_const, density, poly_index);

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
        return (1.0/120.0)*(pow((3-q),5) - 6*pow((2-q),5) + 15*pow((1-q),5));
        
    } else if((q >= 1) && (q < 2)) {
        return (7.0/(478.0*PI))*(pow((3-q),5) - 6*pow((2-q),5));

    } else if((q >= 2) && (q < 3)) {
        return (1/(120.0*PI))*pow((3-q),5);

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

    for(int i=0; i<n_ps; i++){
        particle_set_zero_acc(&ps[i]);
        hydro_euler(i, ps, n_ps, smooth_len);
        //gravity(i, ps, n_ps);
    }

}

void half_velocity_verlet_position(particle *ps, particle *new_ps, int n_ps, double td){

    for(int i=0; i<n_ps; i++){

        particle *p = &ps[i];
        particle *new_p = &new_ps[i];

        new_p->pos[0] = p->pos[0] + p->vel[0]*td + 0.5*p->acc[0]*td*td;
        new_p->pos[1] = p->pos[1] + p->vel[1]*td + 0.5*p->acc[1]*td*td;
        new_p->pos[2] = p->pos[2] + p->vel[2]*td + 0.5*p->acc[2]*td*td;

    }

}

void half_velocity_verlet_velocity(particle *ps, particle *new_ps, int n_ps, double td){

    for(int i=0; i<n_ps; i++){

        particle *p = &ps[i];
        particle *new_p = &new_ps[i];

        new_p->vel[0] = p->vel[0] + 0.5*(p->acc[0]+p->acc[0])*td;
        new_p->vel[1] = p->vel[1] + 0.5*(p->acc[1]+p->acc[1])*td;
        new_p->vel[2] = p->vel[2] + 0.5*(p->acc[2]+p->acc[2])*td;

    }

}

void check_hard_boundaries(particle *ps, int n_ps){
    for(int i=0; i<n_ps; i++){
        particle *p = &ps[i];

        if (p->pos[0] <= -1){
            p->vel[0] = fabs(p->vel[0]);
        }

        else if (p->pos[0] >= 1){
            p->vel[0] = -fabs(p->vel[0]);
        }

        else if (p->pos[1] <= -1){
            p->vel[1] = fabs(p->vel[1]);
        }

        else if (p->pos[1] >= 1){
            p->vel[1] = -fabs(p->vel[1]);
        }

        else if (p->pos[2] <= -1){
            p->vel[2] = fabs(p->vel[2]);
        }

        else if (p->pos[2] >= 1){
            p->vel[2] = -fabs(p->vel[2]);
        }


    }
}

void drag_term(particle *ps, int n_ps, double drag_coeff){

    for(int i=0; i<n_ps; i++){
        particle *p = &ps[i];

        p->vel[0] -= drag_coeff*p->vel[0];
        p->vel[1] -= drag_coeff*p->vel[1];
        p->vel[2] -= drag_coeff*p->vel[2];
    }

}
