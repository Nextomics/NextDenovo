#define _GNU_SOURCE
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <getopt.h>
#include <zlib.h>
#include <libgen.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MIN_SEED_CUTOFF 10000
#define MIN_SEED_DEPTH 20

typedef struct {
	int adjust;
	uint64_t depth;
	uint64_t filter_length;
	uint64_t genome_size;
} opt_t;

typedef struct {
	uint32_t *length;
	uint64_t max_size;
	uint64_t len;
	uint64_t total_bases;
	uint64_t filter_len;
	uint64_t filter_bases;
} reads;

typedef struct {
	uint64_t count[11];
	uint32_t length[11];
} nstat;

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)b - *(int*)a );
}


uint32_t recal_seed_cutoff(opt_t *opt, reads *read, double *seed_depth){
	
	int i;
	int64_t seed_cov_len = 0;
	for (i = 0; i < read->len && read->length[i] >= MIN_SEED_CUTOFF; i ++) seed_cov_len += read->length[i];
	if (seed_cov_len / opt->genome_size < MIN_SEED_DEPTH) 
		while (i < read->len && seed_cov_len < (uint64_t) opt->genome_size * (MIN_SEED_DEPTH + 5))
			seed_cov_len += read->length[i ++];
	
	*seed_depth = (double) seed_cov_len / opt->genome_size;
	return read->length[i - 1];
}


static void out_stat(opt_t *opt, reads *read){
	qsort(read->length, read->len, sizeof(uint32_t), cmpfunc);

	uint64_t bin = read->len/1000;
	if (! bin) bin = 10;
	int bin_step = 1000;
	printf("[Read length histogram ('*' =~ %lu reads)]", bin);
	uint64_t bin_count = 0;
	int bin_index = 1;

	nstat nstat_;
	int nstat_index = 1;
	uint64_t nstat_bases = 0;
	memset(&nstat_, 0, sizeof(nstat));

	uint32_t seed_cutoff = 0;
	double seed_depth = (double) opt->depth;
	int64_t seed_cov_len = (uint64_t) opt->depth * opt->genome_size;
	int i;
	for (i = 0; i < read->len; i ++){
		if (seed_cutoff == 0){
			seed_cov_len -= read->length[i];
			if (seed_cov_len <= read->length[i]) seed_cutoff = read->length[i];
		}

		if (read->length[read->len - 1 - i] < bin_step * bin_index || read->len - 1 - i < bin << 1) {
			bin_count ++;
		}else{
			do {
				printf("\n%7d %7d %10lu  ", bin_step * (bin_index - 1), bin_step * bin_index - 1, bin_count);
				for (int j = 0; j < bin_count/bin; j++){
					printf("*");
				}
				bin_count = 0;
				bin_index ++;
			}
			while (read->length[read->len - 1 - i] > bin_step * bin_index);
			bin_count = 1;
		}
		nstat_bases += read->length[i];
		nstat_.count[nstat_index - 1] ++;
		if (nstat_bases >= nstat_index * 0.1 * read->total_bases){
			nstat_.length[nstat_index - 1] = read->length[i];
			nstat_.count[nstat_index] += nstat_.count[nstat_index - 1];
			nstat_index ++;
		}
	}

	if (opt->adjust && seed_cutoff < MIN_SEED_CUTOFF) seed_cutoff = recal_seed_cutoff(opt, read, &seed_depth);
	else{
		opt->adjust = 0;
		if (!seed_cutoff) {
			seed_depth = (double) read->total_bases/opt->genome_size;
			seed_cutoff = opt->filter_length;
		}
	}
	if (seed_cutoff == opt->filter_length) seed_cutoff ++;

	printf("\n%7d %7d %10lu  ", bin_step * (bin_index - 1), read->length[0], bin_count);
	for (int j = 0; j < bin_count/bin; j++){
		printf("*");
	}

	printf("\n\n[Read length stat]\n");
	printf("%5s %20s %10s\n","Types", "Count (#)", "Length (bp)");
	for (int i = 0 ; i < 9; i ++){
		printf("N%-4d %20lu %7u\n",(i + 1) * 10, nstat_.count[i], nstat_.length[i]);
	}

	printf("\n%-8s %20s %20s %10s\n","Types", "Count (#)", "Bases (bp)", "Depth (X)");
	printf("%-8s %20lu %20lu %10.2f\n","Raw", read->len + read->filter_len, \
		read->total_bases + read->filter_bases, (read->total_bases + read->filter_bases)/(float) opt->genome_size);
	printf("%-8s %20lu %20lu %10.2f\n","Filtered", read->filter_len, read->filter_bases, read->filter_bases/(float) opt->genome_size);
	printf("%-8s %20lu %20lu %10.2f\n","Clean", read->len, read->total_bases, read->total_bases/(float) opt->genome_size);
	
	printf("\n*Suggested seed_cutoff (genome size: %.2fMb, expected seed depth: %lu, real seed depth: %.2f): %u bp\n",\
		(double)opt->genome_size/1000000, opt->depth, seed_depth, seed_cutoff);

	if (seed_cutoff < MIN_SEED_CUTOFF) printf("\033[35m*NOTE:\033[0m The read/seed length is too short, and" 
		" the assembly result is unexpected and please check the assembly quality carefully." 
		" Of course, it's better to sequencing more longer reads and try again.\n");
}

static void count_reads(char *line, opt_t *opt, reads *read)
{
	gzFile fp;
	fp = gzopen(line, "r");
	gzbuffer(fp, 131072);
	if (fp == NULL) {
		fprintf(stderr, "Error! %s does not exist!", line);
		exit(1);
	}
	kseq_t *seq;
	seq = kseq_init(fp);

	while(kseq_read(seq) >= 0) {
		if (seq->seq.l < opt->filter_length ) {
			read->filter_len ++;
			read->filter_bases += seq->seq.l;
		}else{
			read->length[read->len++] = seq->seq.l;
			read->total_bases += seq->seq.l;
			if (read->len >= read->max_size) {
				read->max_size += 1000000;
				read->length = realloc(read->length, read->max_size * sizeof (uint32_t));
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
}

static inline uint64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm')  x *= 1e6;
	else if (*p == 'K' || *p == 'k')  x *= 1e3;
	return (uint64_t)(x + .499);
}

static int usage()
{
	fprintf(stderr, "Usage: seq_stat [options] input.fofn\n");
	fprintf(stderr, "Simple statistics of input files\n\n");

	fprintf(stderr, "Options:\n");
	fprintf(stderr, " -f filter length, default 1kb\n");
	fprintf(stderr, " -g genome size, default: 5Mb\n");
	fprintf(stderr, " -d expected seed depth (30-45), used to be corrected, default: 45\n");
	fprintf(stderr, " -a disable automatic adjustment.\n");
	fprintf(stderr, " -o output file, default: stdout\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int c;
	opt_t opt = {
		.adjust = 1,
		.depth = 45,
		.filter_length = 1000,
		.genome_size = 5000000
	};
	while((c = getopt(argc, argv, "f:g:d:o:a")) != -1) {
		switch (c) {
			case 'f':
				opt.filter_length = mm_parse_num(optarg);
				break;
			case 'g':
				opt.genome_size = mm_parse_num(optarg);
				break;
			case 'd':
				opt.depth = mm_parse_num(optarg);
				break;
			case 'a':
				opt.adjust = 0;
				break;
			case 'o':
				if (freopen(optarg, "wb", stdout) == NULL) {
					fprintf(stderr, "Failed to write output to file %s\n", optarg);
				}else break;
			default:
				return usage();
		}
	}
	if (optind + 1 > argc) return usage();

	char *fofn = strdup(argv[optind]);
	char *path = dirname(fofn);
	FILE *stream = fopen(argv[optind], "r");
	if (stream == NULL) {
		fprintf(stderr, "Failed open input file list!\n");
		return usage();
	}

	reads read;
	memset(&read, 0, sizeof(reads));
	read.max_size = 50000000;
	read.length = malloc(sizeof(uint32_t) * read.max_size);

	char *line = NULL;
	size_t len = 0;
	ssize_t nread;
	while((nread = getline(&line, &len, stream)) != -1) {
		if (nread != 1 && line[0] != '#') {
			line[strlen(line) - 1] = '\0';
			if (line[0] == '/')
				count_reads(line, &opt, &read);
			else {
				char file_path[1024];
				sprintf(file_path, "%s/%s", path, line);
				count_reads(file_path, &opt, &read);
			}
		}
	}
	free(line);
	fclose(stream);
	if (fofn) free (fofn);
	if (read.len) out_stat(&opt, &read);
	free (read.length);
	return 0;
}
