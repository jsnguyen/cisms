#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "time.h"
#include "cisms/particle.h"
#include "cisms/config.h"
#define SOL_MASS 1.989e30
#define EAR_MASS 5.972e24
#define AU 1.396e11
#define PARSEC 3.0857e16

int main(int argc, char* argv[]){
    if(argc != 2){
        printf("Usage: ./sim.exe <filename.config>\n");
        return 1;
    }

    config *conf=NULL;
    conf = config_create();
    config_read(conf,argv[1]);
    config_print(conf);

    printf("Beginning main simulation loop.\n");

    clock_t start_t = clock();
    clock_t elaps_t;
    clock_t begin_loop_t;
    double cpu_t, loop_t;
    double loop_ravg=0;

    for(int i=0; i < conf->n_iter; i++){
        begin_loop_t = clock();

        elaps_t = clock();
        cpu_t = (double) (elaps_t - start_t) / CLOCKS_PER_SEC;
        loop_t = (double) (elaps_t - begin_loop_t) / CLOCKS_PER_SEC;
        loop_ravg = (loop_t + (i * loop_ravg)) / (i + 1);

        printf("\r[%.2f%%] || elapsed: %.3fs ||  est time comp: %.3fs || avg loop time: %.3fs ||", (double) 100*(i+1)/conf->n_iter, cpu_t, loop_ravg*(conf->n_iter-1-i), loop_ravg);

        fflush(stdout);
    }

    printf("\n");

    config_destroy(conf);
    
    printf("Success!\n");
    return 0;
}
