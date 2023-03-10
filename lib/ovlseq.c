#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include "bseq.h"
#include "ovl.h"
#include "index.h"

typedef struct {
	FILE *ovl;
	FILE **f2bits;
	uint32_t **f2bits_;
	uint32_t *decode_tbl;
	int handle_index;
	index_tag *indexs;
	buffer_t *ovlbuf;
	sbbuf *bitbuf;
} ovls;


static inline void get_2bit_name (char *path, const ssize_t len){
	ssize_t i = len - 1;
	while (path[i] != '/'){
		i --;
	}
	while (++i < len - 4){
		path[i] = path[i + 1];
	}
	path[i++] = '2';
	path[i++] = 'b';
	path[i++] = 'i';
	path[i++] = 't';
	path[i++] = '\0';
}

char *getseq(ovls *ovls_, uint32_t *tmp_arr){
	if (ovls_->f2bits_){
		subbit_(ovls_->bitbuf, (uint32_t *) ovls_->indexs[tmp_arr[4]].f2bit, \
		tmp_arr[1], ovls_->indexs[tmp_arr[4]].offset, tmp_arr[5], tmp_arr[6]);
	}else{
		subbit(ovls_->bitbuf, (FILE *)ovls_->indexs[tmp_arr[4]].f2bit, \
		tmp_arr[1], ovls_->indexs[tmp_arr[4]].offset, tmp_arr[5], tmp_arr[6]);
	}
	return ovls_->bitbuf->seq;
}

ovls *init_ovls(char *idx_fofn, char *ovl_file, int readbit){
	ovls *ovls_;
	ovls_ = malloc(sizeof(ovls));
	ovls_->decode_tbl = init_decode_table(NUM_LIMIT);
	ovls_->handle_index = 0;

	int handle_number = 1024;
	ovls_->f2bits = calloc(handle_number, sizeof(FILE *));
	if (readbit) ovls_->f2bits_ = calloc(handle_number, sizeof(uint32_t *));
	else ovls_->f2bits_ = NULL;

	FILE *idx = fopen(idx_fofn, "r");
	if (idx == NULL){
		fprintf(stderr, "Failed open idx files!\n");
		exit(1);
	}

	char *file = NULL;
	size_t len = 0;
	ssize_t nread;
	void *f = NULL;

	uint64_t fsize;
	uint32_t index_tag_size = 65536;//65536
	ovls_->indexs = calloc(index_tag_size, sizeof(index_tag));
	while ((nread = getline(&file, &len, idx)) != -1) {
		if (nread != 1 && file[0] != '#'){
			file[nread - 1 ] = '\0';
			char *f2bit = calloc(nread + 1, sizeof(char));
			strcpy(f2bit, file);
			get_2bit_name(f2bit, nread - 1);

			ovls_->f2bits[ovls_->handle_index] = fopen(f2bit, "r");
			if (ovls_->f2bits[ovls_->handle_index] == NULL){
				fprintf(stderr, "Failed open file [%s], may be reached the maximum open files [%d].\n", f2bit, ovls_->handle_index);
				exit(1);
			}
			f = (void *)ovls_->f2bits[ovls_->handle_index];

			if (ovls_->f2bits_){
				fseeko(ovls_->f2bits[ovls_->handle_index], 0, SEEK_END);
				fsize = ftello(ovls_->f2bits[ovls_->handle_index]);
				fseeko(ovls_->f2bits[ovls_->handle_index], 2, SEEK_SET);//firt two bytes are lables for bit file
				ovls_->f2bits_[ovls_->handle_index] = malloc(fsize + 2);
				if (ovls_->f2bits_[ovls_->handle_index] == NULL){
					fprintf(stderr, "\033[35m Failed malloc buf for DB file [%s], ", f2bit);
					fsize = 0;
				}else{
					fread(ovls_->f2bits_[ovls_->handle_index], 4, (fsize + 2) >> 2, ovls_->f2bits[ovls_->handle_index]);
					if(!feof(ovls_->f2bits[ovls_->handle_index])){
						fprintf(stderr, "\033[35m Failed read DB file [%s] into buf, ", f2bit);
						fsize = 0;
					}
				}
				if (!fsize){
					fprintf(stderr, "and disable module of reading DB files into buf and continue... \033[0m\n");
					for (ovls_->f2bits_ = 0; len <= ovls_->handle_index; len++){
						if (ovls_->f2bits_[len]) free(ovls_->f2bits_[len]);
					}
					ovls_->f2bits_ = NULL;
				}else f = (void *) ovls_->f2bits_[ovls_->handle_index];
			}

			read_index(file, &ovls_->indexs, &index_tag_size, f);
			free(f2bit);

			ovls_->handle_index ++;
			if (ovls_->handle_index >= handle_number){
				handle_number += 1024;
				ovls_->f2bits = realloc(ovls_->f2bits, handle_number * sizeof(FILE *));
				if (readbit) ovls_->f2bits_ = realloc(ovls_->f2bits_, handle_number * sizeof(uint32_t *));
			}
		}
	}

	if (file){
		free(file);
	}
	fclose(idx);

	ovls_->ovl = fopen(ovl_file, "r");
	if (ovls_->ovl == NULL){
		fprintf(stderr, "Failed open ovl file!\n");
		exit(1);
	}
	ovls_->bitbuf = init_sbbuf(100000);
	ovls_->ovlbuf = init_buffer(DCBUFSIZE);
	return ovls_;
}

void destory_ovls(ovls *ovls_){
	flush_buffer(NULL, ovls_->ovlbuf);
	fclose(ovls_->ovl);
	int i = 0;
	for (; i < ovls_->handle_index; i++){
		fclose(ovls_->f2bits[i]);
		if (ovls_->f2bits_) free(ovls_->f2bits_[i]);
	}
	destroy_sbbbuf(ovls_->bitbuf);
	if (ovls_->f2bits_) free(ovls_->f2bits_);
	// free(ovls_->bases_arr);
	free(ovls_->f2bits);
	free(ovls_->indexs);
	free(ovls_->decode_tbl);
	free (ovls_);
}

int main (int argc, char *argv[]){
	ovls *ovls_ = init_ovls(argv[1], argv[2], 0);
	
	prev_t ids_;
	ids_.prev_qname = ids_.prev_tname = 0;

	uint32_t tmp_arr[8];
	uint8_t init = 0;
	uint32_t last_seed = 0;
	uint32_t t_name, q_name, t_s, t_e;
	while(decode_ovl(ovls_->ovl, ovls_->decode_tbl, &ids_, tmp_arr, ovls_->ovlbuf, 8) >= 0) {
		t_name = tmp_arr[0];
		t_s = tmp_arr[2];
		t_e = tmp_arr[3];
		q_name = tmp_arr[4];
		if (init && t_name != last_seed) fprintf(stdout, "+\n");
		getseq(ovls_, tmp_arr);
		fprintf(stdout, "%u %u %u %s\n", q_name, t_s, t_e, ovls_->bitbuf->seq);
		last_seed = t_name;
		init = 1;
	}
	fprintf(stdout, "+\n");
	destory_ovls(ovls_);
	return 0;
}
