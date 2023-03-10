#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "index.h"

void read_index(char *file, index_tag **indexs, uint32_t *index_tag_size, void *f2bit){

	char *line = NULL;
	size_t len = 0;

	FILE *lines = fopen(file, "r");

	if (lines == NULL){
		fprintf(stderr, "Failed open file [%s]!\n", file);
		exit(1);
	}

	uint32_t name;
	uint64_t offset; // may be overflow of max_value
	uint32_t length;
	while (getline(&line, &len, lines) != -1) {

		sscanf(line, "%u\t%lu\t%u", &name, &offset, &length);//%hu\t

		while (*index_tag_size <= name){
			*index_tag_size += 65536;
			*indexs = realloc(*indexs, *index_tag_size * sizeof(index_tag));// need pointer of pointer
		}

		(*indexs)[name].f2bit = f2bit;
		(*indexs)[name].offset = offset;
	}
	free(line);
	fclose(lines);
}


void read_len_index(char *file, index_len_tag **indexs, uint32_t *index_tag_size, FILE *f2bit){

	char *line = NULL;
	size_t len = 0;

	FILE *lines = fopen(file, "r");

	if (lines == NULL){
		fprintf(stderr, "Failed open file [%s]!\n", file);
		exit(1);
	}

	uint32_t name;
	uint64_t offset; // may be overflow of max_value
	uint32_t length;
	while (getline(&line, &len, lines) != -1) {

		sscanf(line, "%u\t%lu\t%u", &name, &offset, &length);//%hu\t

		while (*index_tag_size <= name){
			*index_tag_size += 65536;
			*indexs = realloc(*indexs, *index_tag_size * sizeof(index_len_tag));// need pointer of pointer
		}

		(*indexs)[name].f2bit = f2bit;
		(*indexs)[name].offset = offset;
		(*indexs)[name].length = length;
	}
	free(line);
	fclose(lines);
}

idxs *init_index(FILE *fa){
	FILE *f;
	char idxf[1024];
	char *file = NULL;
	size_t len = 0;
	ssize_t nread;
	idxs *idx = malloc(sizeof(idxs));

	idx->idxm = 65536;
	idx->idx = malloc(idx->idxm * sizeof(index_len_tag));
	idx->fai = 0;
	idx->fam = 1024;
	idx->fa = malloc(idx->fam * sizeof(FILE *));
	while ((nread = getline(&file, &len, fa)) != -1) {
		if (nread != 1 && file[0] != '#'){
			file[nread - 1 ] = '\0';
			f = fopen(file, "r");
			if (f == NULL){
				fprintf(stderr, "Failed open file [%s]\n", file);
				exit(1);
			}

			idx->fa[idx->fai++] = f;
			if (idx->fai >= idx->fam) {
				idx->fam += 1024;
				idx->fa = realloc(idx->fa, idx->fam * sizeof(FILE *));
			}

			nread = sprintf(idxf, "%s.idx", file);
			idxf[nread] = '\0';
			read_len_index(idxf, &idx->idx, &idx->idxm, f);
		}   
	}
	if (file) free(file);
	return idx;
}

void destroy_index(idxs *idx){
	int i;
	for ( i = 0; i < idx->fai; i++){
		fclose(idx->fa[i]);
	}
	free (idx->fa);
	free(idx->idx);
	free (idx);
}

