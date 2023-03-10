#ifndef BAM_SORT
#define BAM_SORT

#include "htslib/ksort.h"
#include "htslib/hts_os.h"
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

#define SAM_GLOBAL_ARGS_INIT {{0},{0}}

typedef struct sam_global_args {
    htsFormat in; 
    htsFormat out;
    char *reference;
    int nthreads;
    int write_index;
} sam_global_args;


typedef struct bam1_tag {
    bam1_t *bam_record;
    union {
        const uint8_t *tag;
        uint64_t pos;
    } u;
} bam1_tag;

typedef struct heap1_t{
    int i;
    uint32_t rev;
    uint64_t pos, idx;
    bam1_tag entry;
} heap1_t;

typedef struct {
    int n, init;
    struct heap1_t *heap;
    struct trans_tbl *translation_tbl;
    hts_itr_t **iter;
    samFile **fp;
    bam_hdr_t **hdr;
    char** fn;
    uint64_t idx;
} bam_merge_iter;

int bam_merge_iter_init(int, char*,const char *, const char *, bam_merge_iter *);
int bam_merge_iter_core(bam_merge_iter *);
int bam_merge_iter_destroy(bam_merge_iter *);
int bam_sort_core_ext(int, char*, const char *, const char *, const char *, 
    const char *,size_t, int, const htsFormat *, const htsFormat *);
void sam_global_args_free(sam_global_args *);

#endif
