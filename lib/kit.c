#define  _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "bseq.h"

uint64_t calgs(const char *file){
	gzFile fp; 
	fp = gzopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error! %s does not exist!", file); 
		exit(1);
	}
	uint64_t gs = 0;
	kseq_t *seq;
	seq = kseq_init(fp);
	while(kseq_r(seq) >= 0) {
		gs += seq->seq.l;
	}
	kseq_destroy(seq);
	gzclose(fp);
	return gs;
}

///////////////////////////////////
typedef struct record {
	uint32_t *len;
	unsigned int i, im;
} record;

static int cmp_uint(const void * a, const void * b){
	return *(uint32_t *) b - *(uint32_t *) a;
}

void read_data_from_file(record *data, const char *file){
	FILE *fp = fopen(file, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error! %s does not exist!", file);
		exit(1);
	}

	char *line = NULL;
	size_t len = 0;
	ssize_t ret, nread;
	uint32_t name, length;
	uint64_t offset;

	while((nread = getline(&line, &len, fp)) != -1) {
		ret = sscanf(line, "%u\t%lu\t%u", &name, &offset, &length);
		if (ret >= 0){
			data->len[data->i++] = length;
			if (data->i >= data->im) {
				data->im += 100000;
				data->len = realloc(data->len, sizeof(uint32_t) * data->im);
			}
		}
	}
	if (line) free(line);
	fclose (fp);
}

uint32_t cal_minlen_from_idx(char **files, const int file_count, const uint64_t total_len){
	unsigned int i;
	record data = {0};
	data.im = 100000;
	data.len = malloc(sizeof(uint32_t) * data.im);

	for (i = 0; i < file_count; i++){
		read_data_from_file(&data, files[i]);
	}

	qsort(data.len, data.i, sizeof(uint32_t), cmp_uint);
	
	uint64_t total;
	for (i = total = 0; i < data.i && total < total_len; i ++) total += data.len[i];

	total = i >= 1 ? data.len[i - 1] : 0;
	free (data.len);
	return (uint32_t) total;
}

// #include <libgen.h>
// int main(int argc, char *argv[])
// {
// 	// uint64_t gs = calgs(argv[1]);
// 	// printf("genome size: %lu bp\n", gs);

// 	char *path = dirname(strdup(argv[optind]));
// 	FILE *stream = fopen(argv[optind], "r");
// 	if (stream == NULL) {
// 		fprintf(stderr, "Failed open input file list!\n");
// 		exit(1);
// 	}

// 	int file_i = 0, file_m = 1024;
// 	char **files = malloc(file_m * sizeof(char *));

// 	char *line = NULL;
// 	size_t len = 0;
// 	ssize_t nread;
// 	while((nread = getline(&line, &len, stream)) != -1) {
// 		if (nread != 1 && line[0] != '#') {
// 			line[strlen(line) - 1] = '\0';
// 			if (line[0] == '/'){
// 				files[file_i++] = strdup(line);
// 			}else{
// 				char file_path[1024];
// 				sprintf(file_path, "%s/%s", path, line);
// 				files[file_i++] = strdup(file_path);
// 			}
// 			if (file_i >= file_m){
// 				file_m += 1024;
// 				files = realloc(files, file_m * sizeof(char *));
// 			}
// 		}
// 	}
// 	if (line) free(line);
// 	free (path);
// 	fclose(stream);

// 	printf("cal_minlen_from_idx: %u file_i %d\n", cal_minlen_from_idx(files, file_i, 18300000), file_i);
// 	for (file_m = 0; file_m < file_i; file_m ++) free (files[file_m]);
// 	free (files);
// 	return 0;
// }
