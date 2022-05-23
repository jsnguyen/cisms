#include "cisms/particle.h"

void particle_init(particle *p){
    p = malloc(sizeof(particle));
    p->pos = malloc(3*sizeof(double));
    p->vel = malloc(3*sizeof(double));
    p->acc = malloc(3*sizeof(double));
}

void particle_free(particle *p){

    if(p->pos){
        free(p->pos);
    }

    if(p->vel){
        free(p->pos);
    }

    if(p->acc){
        free(p->pos);
    }

    if(p){
        free(p);
    }
}

particle *multi_particle_init(int n_ps){

    particle *ps;
    ps = malloc(n_ps*sizeof(particle));

    for (int i=0; i<n_ps; i++){
        ps[i].pos = malloc(3*sizeof(double));
        ps[i].vel = malloc(3*sizeof(double));
        ps[i].acc = malloc(3*sizeof(double));

        ps[i].id = 0;
        ps[i].mass = 0;
        ps[i].diam = 0;
        ps[i].density = 0;
        ps[i].pressure = 0;

        ps[i].pos[0] = 0;
        ps[i].pos[1] = 0;
        ps[i].pos[2] = 0;

        ps[i].vel[0] = 0;
        ps[i].vel[1] = 0;
        ps[i].vel[2] = 0;

        ps[i].acc[0] = 0;
        ps[i].acc[1] = 0;
        ps[i].acc[2] = 0;
    }

    return ps;
}

void multi_particle_free(particle *ps, int n_ps){

    for (int i=0; i<n_ps; i++){

        if(ps[i].pos){
            free(ps[i].pos);
        }

        if(ps[i].vel){
            free(ps[i].vel);
        }

        if(ps[i].acc){
            free(ps[i].acc);
        }

    }

    if(ps){
        free(ps);
    }

}


void particle_set(particle *p, int id, double mass, double diam, double pos[3], double vel[3]){
    p->id = id;
    p->mass = mass;
    p->diam = diam;
    memcpy(p->pos, pos, 3*sizeof(double));
    memcpy(p->vel, vel, 3*sizeof(double));
}

void particle_set_zero_acc(particle *p){
    p->acc[0] = 0;
    p->acc[1] = 0;
    p->acc[2] = 0;
}

void particle_copy(particle *dest, particle *src){

    dest->id = src->id;
    dest->mass = src->mass;
    dest->diam = src->diam;
    dest->density = src->density;
    dest->pressure = src->pressure;

    dest->pos[0] = src->pos[0];
    dest->pos[1] = src->pos[1];
    dest->pos[2] = src->pos[2];

    dest->vel[0] = src->vel[0];
    dest->vel[1] = src->vel[1];
    dest->vel[2] = src->vel[2];

    dest->acc[0] = src->acc[0];
    dest->acc[1] = src->acc[1];
    dest->acc[2] = src->acc[2];

}

void particle_print(particle p){
    printf("id     = %i\n", p.id);
    printf("mass   = %f\n", p.mass);
    printf("diam   = %f\n", p.diam);
    printf("dens   = %f\n", p.density);
    printf("press  = %f\n", p.pressure);
    printf("pos    = %.3f,%.3f,%.3f\n", p.pos[0], p.pos[1], p.pos[2]);
    printf("vel    = %.3f,%.3f,%.3f\n", p.vel[0], p.vel[1], p.vel[2]);
}

void particle_write(particle p, const char* fn){
    FILE *f;
    f = fopen(fn,"a");
    fprintf(f,"%f %f %f %e\n",p.pos[0],p.pos[1],p.pos[2], p.density);
    fclose(f);
}
