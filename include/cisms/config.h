#ifndef CONFIG_H_
#define CONFIG_H_

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"

#define NPARAMS 4
#define BUFFSIZE 255

typedef struct{
    char  *fn;
    int    n_iter;
    int    np;
    double td;
} config;

config* config_create();
void config_destroy(config *conf);
void config_read(config *conf, const char* infile);
void config_print(config *conf);

char *ltrim(char *s);
char *rtrim(char *s);
char *trim(char *s);
#endif
