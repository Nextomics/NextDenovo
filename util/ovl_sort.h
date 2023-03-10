#include "../lib/ovl.h"

#define BIN_OFFSET 6
#define MAX_OVL_COV 150
#define MAX_READ_THREAD 4
#define MAX_SORT_THREAD 10
#define MAX_PATH_LENGTH 1024
#define MAX_BIN_COUNT 31250 //2m
#define MIN_MEM_UNIT 100000000
#define BIN_TOLERANCE_EDGE 50
#define BIN_TOLERANCE_COUNT 5
#define MIN_CONTAINTED_COUNT 2

typedef struct {
	char files[MAX_PATH_LENGTH];
	FILE *fp;
	buffer_t *buf;
	prev_t seed_name;
} file;

typedef struct {
	file *files;
	uint32_t len;
	uint32_t max_size;
} overlap_file;

typedef struct {
	char *outfile;
	overlap *ovls;
	uint8_t sorted;
	uint8_t saved;
	uint32_t len;
	uint32_t max_ovls_len;
	uint32_t max_size;
} overlap_data;

typedef struct {
	FILE *rp;
	uint8_t save_seed;
	char outfile[MAX_PATH_LENGTH];
	overlap_file *input_file;
	overlap_data *input_data;
	threadpool thpool;
	int input_data_len;
	int input_file_index;
	int input_file_step;
} sort_args;

typedef struct {
	int max_thread;
	int input_data_len;
	int max_readfile_in_single_thread;
	uint64_t max_mem_each_file;
} merge_thread_opt;

typedef struct {
	uint8_t pcount;
	uint16_t ovl_i, ovl_m;
	uint16_t *bins;
	overlap *ovls;
	uint32_t qname, qs, qe, qlen, qmaxlen, qcov;
	uint32_t binlen, bincount;
	uint32_t contained, chimera;
	uint32_t max_size, len;
} cov_data;

typedef struct {
	char *idx;
	char *outfile;
	char tempdir[MAX_PATH_LENGTH];
	int thread;
	uint64_t max_mem;
} opt;

static inline uint64_t mm_parse_num(const char *str);
static int get_inputfile_index(sort_args *sort);
static int get_outputfile_index();
static int get_overlapdata_index(sort_args *sort);
static inline int check_chimer(cov_data *cov);
static inline void out_chi_con(cov_data *cov, FILE *fp);
static inline void encode_ovl_filter(FILE *fp, FILE *rp, prev_t *name, overlap *ovl, cov_data *cov, buffer_t *buf);

void init_seeds(char *file);
void sort_ovl(overlap_data *input_data);
void merge_ovl_from_file (sort_args *sort);
void merge_ovl_from_sort(sort_args *sort);
void sort_ovl_file(sort_args *sort);
void print_log(char *format, ...);
void collect_ovl_file(char *file, overlap_file *ovlfiles);
void init_cov(cov_data *cov);
void init_ovlfiles(overlap_file *ovlfiles, char *infile);
void init_mt_opt(merge_thread_opt *mt_opt, opt *opts, int sorted_files_len);
void init_sort(sort_args*, uint8_t, threadpool, overlap_file*, opt*, int, uint64_t, int, int);
void destroy_sort_input(sort_args *sort);
int cmp_ovl(const void *a_, const void *b_);
int load_input_data_from_file(overlap_data *input_data, FILE *fp, buffer_t *buf, prev_t *name);
FILE *init_rp(opt *opts);
