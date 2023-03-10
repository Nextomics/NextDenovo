#include <stdint.h>
#include <stdio.h>

typedef struct {
	void *f2bit;
	uint64_t offset;
} index_tag;

typedef struct {
	FILE *f2bit;
	uint64_t offset;
	uint32_t length;
} index_len_tag;

typedef struct {
	index_len_tag *idx;
	uint32_t idxm;
	FILE **fa;
	int fai;
	int fam;
} idxs;

void read_index(char *file, index_tag **indexs, uint32_t *index_tag_size, void *f2bit);
void read_len_index(char *file, index_len_tag **indexs, uint32_t *index_tag_size, FILE *f2bit);
idxs *init_index(FILE *fa);
void destroy_index(idxs *idx);
