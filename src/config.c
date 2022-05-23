#include "cisms/config.h"

config* config_create(){
    config *conf;
    conf = malloc(sizeof(config));
    conf->fn = malloc(BUFFSIZE*sizeof(char));
    return conf;
}

void config_free(config *conf){

    if(conf->fn){
        free(conf->fn);
    }

    if(conf){
        free(conf);
    }

}

void config_read(config *conf, const char* infile){
    FILE *fp;
    fp = fopen(infile, "r");

    char line[BUFFSIZE];
    char *num;
    for (int i = 0; i < NPARAMS; i++) {

        if (fgets(line, BUFFSIZE, fp)) {

            num = strchr(line,'=')+1;

            switch (i) {

                case 0:
                    strcpy(conf->fn,trim(num));
                    break;
                case 1:
                    conf->n_iter = (int) atof(num);
                    break;
                case 2:
                    conf->np = (int) atof(num);
                    break;
                case 3:
                    conf->td = atof(num);
                    break;
                case 4:
                    conf->poly_index = atof(num);
                    break;
                case 5:
                    conf->prop_const = atof(num);
                    break;
                case 6:
                    conf->smooth_len = atof(num);
                    break;
                case 7:
                    conf->drag_coeff = atof(num);
                    break;

            }
        }
    }

    fclose(fp);

}
void config_print(config *conf){
    printf("** Parameters **\n");
    printf("fn       = %s\n",conf->fn);
    printf("n_iter   = %i\n",conf->n_iter);
    printf("np       = %i\n",conf->np);
    printf("td       = %f\n",conf->td);
}

char *ltrim(char *s){
    while(isspace(*s)) s++;
    return s;
}

char *rtrim(char *s){
    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}

char *trim(char *s){
    return rtrim(ltrim(s)); 
}
