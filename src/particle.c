#include "cisms/particle.h"


void particle_set(particle *p, int id, double mass, double diam, double pos[3], double vel[3]){
    p->id = id;
    p->mass = mass;
    p->diam = diam;
    memcpy(p->pos, pos, 3*sizeof(double));
    memcpy(p->vel, vel, 3*sizeof(double));
}

void particle_print(particle p){
    printf("id     = %i\n", p.id);
    printf("mass   = %f\n", p.mass);
    printf("diam   = %f\n", p.diam);
    printf("pos    = %.3f,%.3f,%.3f\n", p.pos[0], p.pos[1], p.pos[2]);
    printf("vel    = %.3f,%.3f,%.3f\n", p.vel[0], p.vel[1], p.vel[2]);
}

void particle_write(particle p, const char* fn){
    FILE *f;
    f=fopen(fn,"a");
    fprintf(f,"%f %f %f\n",p.pos[0],p.pos[1],p.pos[2]);
    fclose(f);
}
