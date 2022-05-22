#include "cisms/hydro.h"

void hydro_euler(int particle_index, particle *ps, int n_ps){

    particle *p = &ps[particle_index];

    double dr[3];

    for(int i=0; i<n_ps; i++){

        if(i == particle_index){
            continue;
        }

        dr[0] = p->pos[0] - ps[i].pos[0];
        dr[1] = p->pos[1] - ps[i].pos[1];
        dr[2] = p->pos[2] - ps[i].pos[2];

        p->acc[0] += p->mass * ( p->pressure / (p->density*p->density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(0, dr, SMOOTHING_LENGTH);
        p->acc[1] += p->mass * ( p->pressure / (p->density*p->density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(1, dr, SMOOTHING_LENGTH);
        p->acc[2] += p->mass * ( p->pressure / (p->density*p->density) + ps[i].pressure / (ps[i].density*ps[i].density) ) * quintic_spline_gradient(2, dr, SMOOTHING_LENGTH);

    }

    p->acc[0] = -p->acc[0];
    p->acc[1] = -p->acc[1];
    p->acc[2] = -p->acc[2];

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

void calc_density(int particle_index, particle *ps, int n_ps){

    particle *p = &ps[particle_index];
    double r;
    double weight;

    double density = 0;
    for(int i=0; i<n_ps; i++){

        r = calc_distance(p, &ps[i]);
        weight = quintic_spline(r, SMOOTHING_LENGTH);

        density += ps[i].mass*weight;

    }

    p->density = density;
    p->pressure = polytrope(PROP_CONST, density, POLYTROPIC_INDEX);

}

double calc_distance(particle *a, particle *b){
    return sqrt(((a->pos[0] - b->pos[0])*(a->pos[0] - b->pos[0]))
              + ((a->pos[1] - b->pos[1])*(a->pos[1] - b->pos[1]))
              + ((a->pos[2] - b->pos[2])*(a->pos[2] - b->pos[2])));
}

void velocity_verlet(int particle_index, particle *ps, int n_ps, config conf){

    particle *p = &ps[particle_index];

    for(int i=0; i<conf.np; i++){
        calc_density(i, ps, conf.np);
    }

    for(int i=0; i<conf.np; i++){

        if(i == particle_index){
            continue;
        }

        hydro_euler(i, ps, conf.np);
        gravity(i, ps, conf.np);

    }

    p->pos[0] = p->pos[0] + p->vel[0]*conf.td + 0.5*p->acc[0]*conf.td*conf.td;
    p->pos[1] = p->pos[1] + p->vel[1]*conf.td + 0.5*p->acc[1]*conf.td*conf.td;
    p->pos[2] = p->pos[2] + p->vel[2]*conf.td + 0.5*p->acc[2]*conf.td*conf.td;

    for(int i=0; i<conf.np; i++){
        calc_density(i, ps, conf.np);
    }

    for(int i=0; i<conf.np; i++){

        if(i == particle_index){
            continue;
        }

        hydro_euler(i, ps, conf.np);
        gravity(i, ps, conf.np);

    }

    p->vel[0] = p->vel[0] + 0.5*(p->acc[0]+p->acc[0])*conf.td;
    p->vel[1] = p->vel[1] + 0.5*(p->acc[1]+p->acc[1])*conf.td;
    p->vel[2] = p->vel[2] + 0.5*(p->acc[2]+p->acc[2])*conf.td;

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
