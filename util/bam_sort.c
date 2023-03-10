/*  bam_sort.c -- sorting and merging.

    Copyright (C) 2008-2016 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>
    Author: Martin Pollard <mp15@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */
#include <time.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/stat.h>
#include "../lib/bsort.h"

#define SORT_DEFAULT_MEGS_PER_THREAD 1024
#define SORT_MIN_MEGS_PER_THREAD 1

static void sort_usage(FILE *fp)
{
    fprintf(fp,
"Usage: bam_sort [options...] [-|in.sam|in.bam]\n"
"Options:\n"
"  -m INT     Set maximum memory per thread; suffix K/M/G recognized [%dM]\n"
"  -o FILE    Write final output to FILE rather than standard output\n"
"  -T PREFIX  Write temporary files to PREFIX.nnnn.bam\n"
"  -@ INT     Number of additional threads to use [0]\n"
"  -i         Write index file\n", SORT_DEFAULT_MEGS_PER_THREAD);
}

static void complain_about_memory_setting(size_t max_mem) {
    char  *suffix = "";
    const size_t nine_k = 9<<10;
    if (max_mem > nine_k) { max_mem >>= 10; suffix = "K"; }
    if (max_mem > nine_k) { max_mem >>= 10; suffix = "M"; }

    fprintf(stderr,
"[bam_sort] -m setting (%zu%s bytes) is less than the minimum required (%dM).\n\n"
"Trying to run with -m too small can lead to the creation of a very large number\n"
"of temporary files.  This may make sort fail due to it exceeding limits on the\n"
"number of files it can have open at the same time.\n\n"
"Please check your -m parameter.  It should be an integer followed by one of the\n"
"letters K (for kilobytes), M (megabytes) or G (gigabytes).  You should ensure it\n"
"is at least the minimum above, and much higher if you are sorting a large file.\n",
            max_mem, suffix, SORT_MIN_MEGS_PER_THREAD);
}

int main(int argc, char *argv[])
{
    size_t max_mem = SORT_DEFAULT_MEGS_PER_THREAD << 20;
    int c, nargs, is_by_qname = 0, ret, o_seen = 0, level = -1;
    char* sort_tag = NULL;
    char *fnout = "-", modeout[12];
    kstring_t tmpprefix = { 0, 0, NULL };
    struct stat st;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    while ((c = getopt(argc, argv, "m:o:T:@:i")) >= 0) {
        switch (c) {
        case 'o': fnout = optarg; o_seen = 1; break;
        case '@': ga.nthreads = atoi(optarg); break;
        case 'm': {
                char *q;
                max_mem = strtol(optarg, &q, 0);
                if (*q == 'k' || *q == 'K') max_mem <<= 10;
                else if (*q == 'm' || *q == 'M') max_mem <<= 20;
                else if (*q == 'g' || *q == 'G') max_mem <<= 30;
                break;
            }
        case 'T': kputs(optarg, &tmpprefix); break;
        case 'i': ga.write_index = 1; break;

        default: sort_usage(stderr); ret = EXIT_FAILURE; goto sort_end;
        }
    }

    nargs = argc - optind;
    if (nargs == 0 && isatty(STDIN_FILENO)) {
        sort_usage(stdout);
        ret = EXIT_SUCCESS;
        goto sort_end;
    }
    else if (nargs >= 2) {
        // If exactly two, user probably tried to specify legacy <out.prefix>
        if (nargs == 2)
            fprintf(stderr, "[bam_sort] Use -T PREFIX / -o FILE to specify temporary and final output files\n");

        sort_usage(stderr);
        ret = EXIT_FAILURE;
        goto sort_end;
    }

    if (max_mem < (SORT_MIN_MEGS_PER_THREAD << 20)) {
        complain_about_memory_setting(max_mem);
        ret = EXIT_FAILURE;
        goto sort_end;
    }

    strcpy(modeout, "wb");
    sam_open_mode(modeout+1, fnout, NULL);
    if (level >= 0) sprintf(strchr(modeout, '\0'), "%d", level < 9? level : 9);

    if (tmpprefix.l == 0) {
        if (strcmp(fnout, "-") != 0) ksprintf(&tmpprefix, "%s.tmp", fnout);
        else kputc('.', &tmpprefix);
    }
    if (stat(tmpprefix.s, &st) == 0 && S_ISDIR(st.st_mode)) {
        unsigned t = ((unsigned) time(NULL)) ^ ((unsigned) clock());
        if (tmpprefix.s[tmpprefix.l-1] != '/') kputc('/', &tmpprefix);
        ksprintf(&tmpprefix, "bam_sort.%d.%u.tmp", (int) getpid(), t % 10000);
    }

    ret = bam_sort_core_ext(is_by_qname, sort_tag, (nargs > 0)? argv[optind] : "-",
                            tmpprefix.s, fnout, modeout, max_mem, ga.nthreads,
                            &ga.in, &ga.out);
    
    int csi = 0;
    int min_shift = 14;
    if (ret >= 0 && ga.write_index && strcmp(fnout, "-") != 0){
        fprintf(stderr, "[S::main] Sort done and indexing...\n");
        ret = sam_index_build3(fnout, 0, csi? min_shift : 0, ga.nthreads);
    }

    if (ret >= 0)
        ret = EXIT_SUCCESS;
    else {
        char dummy[4];
        // If we failed on opening the input file & it has no .bam/.cram/etc
        // extension, the user probably tried legacy -o <infile> <out.prefix>
        if (ret == -2 && o_seen && nargs > 0 && sam_open_mode(dummy, argv[optind], NULL) < 0)
            fprintf(stderr, "[bam_sort] Note the <out.prefix> argument has been replaced by -T/-o options\n");

        ret = EXIT_FAILURE;
    }

sort_end:
    free(tmpprefix.s);
    sam_global_args_free(&ga);

    return ret;
}