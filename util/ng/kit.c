#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "kit.h"

//from minimap2
double cputime()
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long peakrss()
{
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
    return r.ru_maxrss * 1024;
#else
    return r.ru_maxrss;
#endif
}

double realtime()
{
    struct timeval tp; 
    gettimeofday(&tp, NULL);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}


void plog(const int level, char *format, ...)
{	
    char *type[5] = {"NOTSET", "DEBUG", "INFO", "WARNING", "ERROR"};
    time_t t;
    struct tm *lt;
    time(&t);
    lt = localtime(&t);

    va_list ap; 
    va_start(ap, format);
    fprintf(stderr, "[%s] %02d-%02d-%02d %02d:%02d:%02d ", type[level], \
    lt->tm_year+1900, lt->tm_mon+1, lt->tm_mday, lt->tm_hour, \
    lt->tm_min, lt->tm_sec);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fprintf(stderr, "\n");
}

void init_sn(ph *s, const int len){
    s->im = len;
    s->i = 0;
    s->n = malloc(s->im * sizeof(UINTL_T));
}

void reallocate_sn(ph *s, const int len){
    s->im += len;
    s->n = realloc(s->n, sizeof(UINTL_T) * s->im);
}

void destroy_sn(ph *s){
    free (s->n);
}

int version(){
    fprintf(stderr, "TYPE_GS: %d, UINT_T: %lu, UINTL_T: %lu, ",\
        GENOME_SIZE, sizeof(UINT_T) * 8, sizeof(UINTL_T) * 8);
    #ifdef RESTRICT
        fprintf(stderr, "RESTRICT_GS: %lu\n", RESTRICT_GS);
    #else
        fprintf(stderr, "RESTRICT_GS: NA\n");
    #endif
    return 1;
}

int quick_select(int *r, int s, int e, int k){
    if (s > e) return e + 1;
    int i, left, pivot;
    while (1){
        left = s;
        pivot = r[e];
        for (i = s; i < e; i++){
            if(r[i] < pivot){
                SWAP(r[i], r[left], int);
                left++;
            }
        }
        SWAP(r[left], r[e], int);

        if (left == k) return pivot;
        if (left < k) s = left + 1;
        else e = left - 1;
    }
}