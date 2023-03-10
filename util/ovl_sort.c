#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>
#include "thpool.h"
#include "khash.h"
#include "ovl_sort.h"
#include "../lib/ovl.h"

KHASH_MAP_INIT_INT(id, uint32_t)
static pthread_mutex_t inputlock;
static pthread_mutex_t outputlock;
static pthread_mutex_t sortlock;
static bytes_t *encode_tbl = NULL;
static uint32_t *decode_tbl = NULL;
static int sorted_files_len = 0;
static int min_seed_len = 0;
static int max_bin_cov = 40;
static int max_flank_len = 300;
static int is_hq = 0;
static khash_t(id) *seeds;

void print_log(char *format, ...)
{
	time_t t;
	struct tm *lt;
	time(&t);
	lt = localtime(&t);

	va_list ap;
	va_start(ap, format);
	fprintf(stderr, "[WARNING] %02d-%02d-%02d %02d:%02d:%02d ", \
	lt->tm_year+1900, lt->tm_mon+1, lt->tm_mday, lt->tm_hour, \
	lt->tm_min, lt->tm_sec);
	vfprintf(stderr, format, ap);
	va_end(ap);
	fprintf(stderr, "\n");
}

// void out_cov(cov_data *cov){
// 	if (cov->qname != UINT32_MAX){
// 		printf("qname: %d qlen: %d qcov:%d qsbin:%d qebin:%d\n", cov->qname, cov->qlen/MAX_OVL_COV,
// 				cov->qcov/(cov->qlen/MAX_OVL_COV), cov->qs>> BIN_OFFSET, cov->qe>> BIN_OFFSET);
// 		for (int i = 0; i < cov->len; i++){
// 			printf(" %d", cov->bins[i]);
// 		}
// 		printf("\n");
// 	}
// }

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

static int get_inputfile_index(sort_args *sort){
	// printf("get_inputfile_index..\n");
	pthread_mutex_lock(&inputlock);
	int t = sort->input_file_index += sort->input_file_step;
	pthread_mutex_unlock(&inputlock);
	// printf("get_inputfile_index done %d\n",t);
	return t;
}

static int get_outputfile_index(){
	pthread_mutex_lock(&outputlock);
	int t = ++sorted_files_len;
	pthread_mutex_unlock(&outputlock);
	return t;
}

static int get_overlapdata_index(sort_args *sort){
	// printf("get_overlapdata_index\n" );
	pthread_mutex_lock(&sortlock);
	int i, j;
	i = 1;
	while (i){
		for (j = 0; j < sort->input_data_len; j ++){
			if (sort->input_data[j].sorted == 0) {
				sort->input_data[j].max_ovls_len = 0;
				sort->input_data[j].sorted = 1;
				i = 0;
				break;
			}
		}
		if (i) sleep (1);
	}
	pthread_mutex_unlock(&sortlock);
	// printf("get_overlapdata_index done %d\n",j );
	return j;
}

void init_seeds(char *file)
{
	uint32_t name, length;
	uint64_t offset;
	khint_t k;

	FILE *fp;
	fp = fopen(file, "r");
	if (!fp) {
		print_log("Failed open index file %s!", file); 
		exit(1);
	}

	char *line = NULL;
	size_t len = 0;
	int absent;
	while(getline(&line, &len, fp) != -1) {
		sscanf(line, "%u\t%lu\t%u", &name, &offset, &length);
		k = kh_put(id, seeds, name, &absent);
		if (absent) {
			kh_key(seeds, k) = name;
			kh_val(seeds, k) = length;
			if (length < min_seed_len || !min_seed_len) min_seed_len = length;
		}
	}
	if (line) free(line);
	fclose(fp);
}


void collect_ovl_file(char *file, overlap_file *ovlfiles)
{
	FILE *fp;
	fp = fopen(file, "r");
	if (!fp) {
		print_log("Failed open input overlap file list %s!", file); 
		exit(1);
	}
	
	ssize_t nread;
	size_t len = MAX_PATH_LENGTH;
	char *line = ovlfiles->files[ovlfiles->len].files;
	while((nread = getline(&line, &len, fp)) != -1) {
		if (nread != 1 && ovlfiles->files[ovlfiles->len].files[0] != '#') {
			ovlfiles->files[ovlfiles->len].files[nread - 1] = '\0';
			ovlfiles->len ++;
			line = ovlfiles->files[ovlfiles->len].files;
			if (ovlfiles->len >= ovlfiles->max_size) {
				ovlfiles->max_size += 1000;
				ovlfiles->files = realloc(ovlfiles->files, ovlfiles->max_size * sizeof(file));
			}
		}
	}
	fclose (fp);
}

void init_cov(cov_data *cov){
	memset(cov, -1, sizeof(cov_data));
	cov->max_size = MAX_BIN_COUNT;
	cov->bins = malloc(sizeof(uint16_t) * cov->max_size);
	cov->ovl_m = 5000;
	cov->ovls = malloc(sizeof(overlap) * cov->ovl_m);
	assert(cov->bins && cov->ovls);
}

void init_ovlfiles(overlap_file *ovlfiles, char *infile){
	ovlfiles->len = 0;
	ovlfiles->max_size = 1000;
	ovlfiles->files = malloc(ovlfiles->max_size * sizeof(file));
	collect_ovl_file(infile, ovlfiles);
}

void init_mt_opt(merge_thread_opt *mt_opt, opt *opts, int sorted_files_len){
	int i, j;
	i = j = 2;
	mt_opt->max_thread = mt_opt->max_readfile_in_single_thread = 1;
	for (; mt_opt->max_thread <= opts->thread; mt_opt->max_thread ++){
		for (mt_opt->max_readfile_in_single_thread = (MAX_READ_THREAD + MAX_SORT_THREAD) * 2; \
				mt_opt->max_readfile_in_single_thread > 1; mt_opt->max_readfile_in_single_thread--){
			if (opts->max_mem * opts->thread / (mt_opt->max_readfile_in_single_thread * mt_opt->max_thread) >= MIN_MEM_UNIT && \
				(mt_opt->max_readfile_in_single_thread - 1) * mt_opt->max_thread < sorted_files_len) {
				i = mt_opt->max_thread;
				j = mt_opt->max_readfile_in_single_thread;
				break;
			}
		}
		if (mt_opt->max_readfile_in_single_thread * mt_opt->max_thread >= sorted_files_len) {
			break;
		}
	}

	mt_opt->max_thread = i;
	mt_opt->max_readfile_in_single_thread = j;
	mt_opt->input_data_len = min(mt_opt->max_readfile_in_single_thread * mt_opt->max_thread, sorted_files_len);
	mt_opt->max_mem_each_file = opts->max_mem * opts->thread / mt_opt->input_data_len;
}

FILE *init_rp(opt *opts){
	FILE *rp = NULL;
	if (max_flank_len){
		char outfile[MAX_PATH_LENGTH] = {0};
		assert (sprintf(outfile, "%s.bl", opts->outfile) != -1);
		if ((rp = fopen(outfile, "w")) == NULL) {
			print_log("Failed create rid file %s!", outfile);
			exit(1);
		}
	}
	return rp;
}

void init_sort(sort_args *sort, uint8_t save_seed, threadpool thpool_sort, \
		overlap_file *ovlfiles, opt *opts, int input_data_len, uint64_t max_mem, int input_file_step, int suffix){

	sort->save_seed = save_seed;
	sort->thpool = thpool_sort;
	sort->input_file = ovlfiles;
	sort->input_data_len = input_data_len;
	sort->input_file_step = input_file_step;
	sort->input_file_index = 0;
	memset(sort->outfile, 0, sizeof (sort->outfile));
	assert (sprintf(sort->outfile, "%s/%s.%d.%d.", opts->tempdir, opts->outfile, getpid(), suffix) != -1);
	sort->input_data = calloc(sort->input_data_len, sizeof(overlap_data));
	if (max_mem / sizeof(overlap) >= UINT32_MAX) max_mem = sizeof(overlap) * (UINT32_MAX - 10);
	for (int i = 0; i < sort->input_data_len; i ++){
		// printf("mmalloc: %d %lu %lu %lu\n",i, max_mem,  max_mem / sizeof(overlap), sizeof(overlap));
		sort->input_data[i].outfile = sort->outfile;
		sort->input_data[i].max_size = max_mem / sizeof(overlap);
		sort->input_data[i].ovls = malloc(sort->input_data[i].max_size * sizeof(overlap));
		sort->input_data[i].ovls[0].qe = 0;
		assert (sort->input_data[i].ovls);
	}
}

void destroy_sort_input(sort_args *sort){
	for (int i = 0; i < sort->input_data_len; i ++){
		free (sort->input_data[i].ovls);
	}
	free (sort->input_data);
}

int cmp_ovl(const void *a_, const void *b_){
	overlap *a = (overlap *) a_;
	overlap *b = (overlap *) b_;
	if (a->qname != b->qname){
		return a->qname - b->qname;
	}else if (a->match != b->match){
		return b->match - a->match;
	}else{
		return a->qe - a->qs - b->qe + b->qs;
	}
	// }else if (b->qe - b->qs != a->qe - a->qs){
	// 	return b->qe - b->qs - a->qe + a->qs;
	// }else{
	// 	return a->qs - b->qs;
	// }
}

int load_input_data_from_file(overlap_data *input_data, FILE *fp, buffer_t *buf, prev_t *name){
	uint32_t ovl[8];
	while(decode_ovl(fp, decode_tbl, name, ovl, buf, 8) >= 0) {
		// if (ovl[0] == ovl[4]) continue;
		// if (ovl[3] - ovl[2] < 500 || ovl[6] - ovl[5] < 500) continue;
		// printf("input_data->max_ovls_len: %d %u %d %u\n", input_data->max_ovls_len,ovl[0], ovl[1], input_data->max_size);
		input_data->ovls[input_data->max_ovls_len].qname = ovl[0];
		input_data->ovls[input_data->max_ovls_len].rev = ovl[1];
		input_data->ovls[input_data->max_ovls_len].qs = ovl[2];
		input_data->ovls[input_data->max_ovls_len].qe = ovl[3];
		input_data->ovls[input_data->max_ovls_len].tname = ovl[4];
		input_data->ovls[input_data->max_ovls_len].ts = ovl[5];
		input_data->ovls[input_data->max_ovls_len].te = ovl[6];
		input_data->ovls[input_data->max_ovls_len].match = ovl[7];
		input_data->max_ovls_len ++;
		if (input_data->max_ovls_len >= input_data->max_size){
			input_data->len = 0;
			return 0;
		}
	}
	input_data->len = 0;
	return 1;
}

static inline int check_chimer_hq(cov_data *cov){//TODO: need more test
	int i, j, l, r;

	// for (i=0;i<cov->len;i++){
	// 	fprintf(stderr, "%d(%d) ", cov->bins[i], i);
	// }
	// fprintf(stderr,"\n");
	l = 0;
	while (l < cov->len && cov->bins[l] < 2) l ++;
	r = cov->len;
	while (r > 0 && cov->bins[r - 1] < 2) r --;

	int s, e, flank = 15;
	for (i = l + 1; i < r - 1; i++){
		if (cov->bins[i] <= 1 ) {
			s = i > l + flank ? (i - flank) << BIN_OFFSET : l << BIN_OFFSET;
			e = i + flank < r ? (i + flank) << BIN_OFFSET : r << BIN_OFFSET;
			for (j = 1; j < cov->ovl_i; j ++){
				overlap *ovl = &cov->ovls[j];
				if (ovl->qs < s && ovl->qe > e) break;
			}
			if (j >= cov->ovl_i) {
				// fprintf(stderr, "%u\n", cov->ovls[0].qname);
				return i;
			}
		}
	}
	return 0;
}

static inline int check_chimer(cov_data *cov){
	if (is_hq) return check_chimer_hq(cov);
	int l, r, lable, llable, rlable;
	lable = llable = rlable = 0;
	for (int i = 1; i < cov->len - 1; i++){
		if (cov->bins[i] > 20 && ++llable) {
			if (lable && ++rlable >= 5) break;
		}else{
			l = max(i-5, 0); r = min(i+5, cov->len -1);
			if (llable > 5 && (cov->bins[l] > 20 || cov->bins[r] > 20) && \
				cov->bins[i] <= max(3, min(cov->bins[l], cov->bins[r])/5)){
				lable = i;
			}
		}
	}
	if (rlable < 5) lable = 0;
	return lable;
}

//check hot break ends
static inline int check_chimer2(cov_data *cov){
	int i, j, s, e, t, c;
	overlap *ovl;
	memset(cov->bins, 0, sizeof(uint16_t) * (cov->len/2 + 1));
	for (s = cov->len, c = e = 0, i = 1, j = BIN_OFFSET + 1; i < cov->ovl_i; i ++){
		ovl = &cov->ovls[i];
		if (ovl->qe){
			c ++;
			t = (ovl->qs + 10) >> j;
			if (t < s) s = t;
	 		cov->bins[t] ++;
	 		t = (ovl->qe - 10) >> j;
	 		if (t > e) e = t;
	 		cov->bins[t] ++;
		}
	}

	t = 0;
	if (c > 20){
		int ms, me, m;
		while (s < e && cov->bins[s] < 4) s++;
		while (e > s && cov->bins[e] < 4) e--;
		for (m = 0, ms = cov->bins[s], me = cov->bins[e], i = s; i < e + 1; i ++){
			if (i < s + 5 && cov->bins[i] > ms) ms = cov->bins[i];
			if (i > e - 5 && cov->bins[i] > me) me = cov->bins[i];
			if (cov->bins[i] > cov->bins[m]) m = i;
			// printf("%d\t%d\t%d\n", cov->qname, i, cov->bins[i]);
		}
		
		// printf("%d %d %d %d %d %d %d %d\n", cov->qname,c, m, cov->bins[m],s,ms,e,me);
		float fra = 1;//1.3;
		if (m > s + 5 && m < e - 5 && cov->bins[m] > fra * max(ms, me) && 
				((c > 75 && m > c/5) || (c < 75 && m > c/2))) {//TODO need more test
			t = m << j;
			// printf("%d %d %d %d\n", cov->qname,c, m, cov->bins[m]);
		}
	}
	return t;
}

static inline void out_chi_con(cov_data *cov, FILE *fp){
	if (cov->contained >= MIN_CONTAINTED_COUNT) fprintf(fp, "%u c\n", cov->qname);
	else if (cov->chimera) fprintf(fp, "%u k\n", cov->qname);
}

int cmp_bins(const void *a, const void *b){
	return  ( *(uint16_t*)a < *(uint16_t*)b );
}

void del_repeat_alns(cov_data *cov){//del repeat using break-points
	int i;
	// fprintf(stderr,"#%d %d\n", cov->ovls[0].qname, cov->ovl_i);
	// for (i=0;i<cov->len;i++){
	// 	fprintf(stderr, "%d(%d) ", cov->bins[i], i);
	// }
	// fprintf(stderr,"\n");
	overlap *ovl;
	int offset = 1 + (cov->qlen >> (BIN_OFFSET + 1));
	// int offset = cov->qlen >> (BIN_OFFSET + 1);
	// qsort(cov->bins, (cov->qlen >> (BIN_OFFSET + 1)) + 1, sizeof(uint16_t), cmp_bins);
	// fprintf(stderr, "top5:%d mean:%d offset:%d\n", cov->bins[4], cov->bins[offset/2] + 1, offset);
	int median = 5;//cov->bins[offset/2] > 10 ? cov->bins[offset/2] : 10;
	int flank_len = max_flank_len > 100 ? max_flank_len * 3 : 300;
	// if (cov->bins[0] > median){
		// offset ++;//start + 1
		for (i = 1; i < cov->ovl_i; i ++){
			ovl = &cov->ovls[i];
			if (ovl->qs <= flank_len && ovl->qe + flank_len >= cov->qlen) continue;
			int s = (ovl->qs + 10) >> (BIN_OFFSET + 1);
			int e = (ovl->qe - 10) >> (BIN_OFFSET + 1);
			// fprintf(stderr, "%d %d %d %d %d %d\n",ovl->tname,ovl->qs,ovl->qe,cov->bins[s] ,cov->bins[e + offset],median );
			if (cov->bins[s] >= median && cov->bins[e + offset] >= median) ovl->qe = 0;
		}
	// }
	memset(cov->bins, 0, sizeof(uint16_t) * cov->len);
	for (i = 1; i < cov->ovl_i; i ++){
		ovl = &cov->ovls[i];
		if (!ovl->qe) continue;
		int t, j = (ovl->qs + 10) >> BIN_OFFSET, k = (ovl->qe - 10) >> BIN_OFFSET, mincov = UINT16_MAX;
		for (t = j + 1; t <= k; t ++){
			if (++ cov->bins[t] > UINT16_MAX - 1000) cov->bins[t]--;
			if (cov->bins[t] < mincov) mincov = cov->bins[t];
		}
		if (mincov > 2 * max_bin_cov){
			for (t = j + 1; t <= k; t++) {
				cov->bins[t] --;
			}
			ovl->qe = 0;
		}
	}
	// for (i=0;i<cov->len;i++){
	// 	fprintf(stderr, "%d(%d) ", cov->bins[i], i);
	// }
	// fprintf(stderr,"\n");
}

static void ovl_filter(FILE *fp, FILE *rp, prev_t *name, cov_data *cov, buffer_t *buf){//filter "k" regions
	int i; //, mincov;
	overlap *ovl;
	if (is_hq) del_repeat_alns(cov);
	if (rp){
		int s = 0;
		int e = 0;
		cov->chimera = check_chimer(cov);
		if (cov->chimera || !(cov->contained || is_hq)){
			int j, k, m;
			j = 0;
			if (cov->qcov > cov->qlen * 10){
				for (i = 1; i < cov->len - 1; i ++){
					if (cov->bins[i] < min(4, max_bin_cov/10)){
						if (s == 0) s = i;
						e = i;
					}else if (s) {
						if (cov->chimera && cov->chimera < s && ((!j) || cov->chimera > cov->bins[j - 1])){
							cov->bins[j++] = cov->chimera;
							cov->bins[j++] = cov->chimera;
						}
						cov->bins[j++] = s;
						cov->bins[j++] = e;
						s = e = 0;
					}
				}

				if (s){
					if (cov->chimera && cov->chimera < s && ((!j) || cov->chimera > cov->bins[j - 1])){
						cov->bins[j++] = cov->chimera;
						cov->bins[j++] = cov->chimera;
					}
					cov->bins[j++] = s;
					cov->bins[j++] = e;
				}
				if (cov->chimera && (j == 0 || cov->chimera > cov->bins[j - 1])){
					cov->bins[j++] = cov->chimera;
					cov->bins[j++] = cov->chimera;
				}
			}else if (cov->chimera){
				cov->bins[j++] = cov->chimera;
				cov->bins[j++] = cov->chimera;
			}
			
			// for (i = 0; i < j; i += 2){
			// 	if (i == 0) printf("#%d", cov->qname);
			// 	printf(" (%d %d)",cov->bins[i],cov->bins[i+1]);
			// 	if (i == j - 2) printf("\n");
			// }
			if (j){
				m = j;
				if (cov->bins[0] < 5) m -= 2;
				if (cov->bins[j - 1] > cov->len - 5) m -= 2;

				if (m > 0 ){
					k = 0;
					m = cov->bins[0];
					for (i = 2; i < j; i += 2){
						if (cov->bins[i] - cov->bins[i - 1] > m){
							m = cov->bins[i] - cov->bins[i - 1];
							k = i;
						}
					}

					if (cov->len - cov->bins[i - 1] > m){
						m = cov->len - cov->bins[i - 1];
						s = cov->bins[i - 1];
						e = cov->len;
					}else if (cov->bins[k + 1] > cov->len - 5){
						s = cov->bins[k - 1];
						e = cov->len;
					}else if (k == 0 || cov->bins[k - 2] < 5){
						s = 0;
						e = cov->bins[k];
					}else{
						s = cov->bins[k - 1];
						e = cov->bins[k];
					}
					// printf("$%d %d %d %d %d %d\n", cov->qname, s, e, min_seed_len,m, (min_seed_len >> BIN_OFFSET) * 2 / 3);

					int flank = 5;
					s = s > flank ? (s - flank) << BIN_OFFSET : 0;
					e = (e + flank) << BIN_OFFSET;
					// printf(" %d %d\n",s,e );
					if (m > (min_seed_len >> BIN_OFFSET) * 2 / 3){//TODO CHECK
						cov->chimera = 0;
						for (i = 1; i < cov->ovl_i; i ++){
							ovl = &cov->ovls[i];
							if (ovl->qs < s || ovl->qe > e) ovl->qe = 0;
							// fprintf(stderr, "%u %d %d\n",ovl->qs, s,ovl->qs < s );
						}
					}else cov->chimera = 1;//it's better set 2 or else.
				}else s = e = 0;
			}else assert(s == 0 && e == 0); //TODO remove
		}
		// out_chi_con(cov, rp);
		// printf("%d %f %d %d %d\n",cov->qname, (float)cov->qcov/cov->qlen, s, e, min_seed_len);
		if ((!is_hq) && cov->qcov > cov->qlen * 20 && (!cov->chimera) && cov->contained < MIN_CONTAINTED_COUNT) {
			cov->chimera = check_chimer2(cov);
			if (!e) e = cov->qlen;
			if (cov->chimera <= s + (15 << BIN_OFFSET) || cov->chimera + (15 << BIN_OFFSET) >= e) cov->chimera = 0;
			// if (cov->chimera){
			// 	if (cov->chimera - s > e - cov->chimera) e = cov->chimera + (5 << BIN_OFFSET);
			// 	else s = cov->chimera - (5 << BIN_OFFSET);
			// 	// printf("%d %d %d\n",cov->qname,  s, e);
			// 	if (e - s > min_seed_len * 2 / 3){//TODO CHECK 2/3
			// 		cov->chimera = 0;
			// 		for (i = 1; i < cov->ovl_i; i ++){
			// 			ovl = &cov->ovls[i];
			// 			if (ovl->qs < s || ovl->qe > e) ovl->qe = 0;
			// 		}
			// 	}
			// }
		}
	}
	
	for (i = cov->contained = 0; i < cov->ovl_i; i ++){
		ovl = &cov->ovls[i];
		if (ovl->qe){
			encode_ovl(fp, encode_tbl, name, ovl, buf);
			if (ovl->qname != ovl->tname && ovl->qs <= max_flank_len && ovl->qe + max_flank_len >= cov->qlen) {
				if (!is_hq) cov->contained ++;
				else if (ovl->match >= (ovl->qe - ovl->qs + 1) * 0.9) cov->contained ++;//require check identity for hq reads
			}
		}
	}
	if (rp) out_chi_con(cov, rp);
}

static void flush_cov(FILE *fp, FILE *rp, prev_t *name, cov_data *cov, buffer_t *buf){
	ovl_filter(fp, rp, name, cov, buf);
	flush_buffer(fp, buf);
}

// static inline void encode_ovl_filter(FILE *fp, FILE *rp, prev_t *name, overlap *ovl, cov_data *cov, buffer_t *buf){
// 	int lable = 1;
// 	if (ovl->qname != cov->qname){
// 		if (cov->qname != UINT32_MAX) ovl_filter(fp, rp, name, cov, buf);
// 		cov->qlen = kh_val(seeds, kh_get(id, seeds, ovl->qname));
// 		cov->len = (cov->qlen >> BIN_OFFSET) + 1;
// 		if (cov->len > cov->max_size){
// 			cov->max_size = cov->len;
// 			cov->bins = realloc(cov->bins, cov->max_size * sizeof(uint16_t));
// 		}
// 		memset(cov->bins, 0, cov->len * sizeof(uint16_t));
// 		// fprintf(stderr,"%d %d\n",ovl->qname,cov->len);
// 		cov->qcov = cov->contained = cov->chimera = cov->ovl_i = 0;
// 		cov->qmaxlen = cov->qlen * MAX_OVL_COV;
// 	}else if (cov->qcov > cov->qmaxlen || cov->ovl_i >= MAX_CANDIDATE_OVL){
// 		lable = 0;
// 	}else{
// 		int i;
// 		int j = (ovl->qs + 10) >> BIN_OFFSET;
// 		int k = (ovl->qe - 10) >> BIN_OFFSET;
// 		int mincov = cov->bins[k] + 1;
// 		int maxcov = UINT16_MAX - 10;
// 		// fprintf(stderr, "maxcov:%d\n", maxcov);

// 		for (i = j + 1; i <= k; i++){
// 			// fprintf(stderr,"%d %hhu %d %d %d\n", i, cov->bins[i], ovl->qname,ovl->tname,cov->ovl_i);
// 			if (++cov->bins[i] < mincov) mincov = cov->bins[i];
// 			if (cov->bins[i] > maxcov) cov->bins[i]--;
// 		}
// 		if (mincov >= max_bin_cov){
// 			for (i = j + 1; i <= k; i++) cov->bins[i] --;
// 			lable = 0;
// 		}
// 		// out_cov(cov);
// 	}

// 	if (lable) {
// 		cov->qname = ovl->qname;
// 		cov->qcov += ovl->qe - ovl->qs + 1;
// 		cov->ovls[cov->ovl_i++] = *ovl;		
// 		// for (int i = 0; i <cov->len; i++){
// 		// 	fprintf(stderr,"#%d %hhu %d %d %d\n", i, cov->bins[i], ovl->qname,ovl->tname,cov->ovl_i);
// 		// }
// 		// fprintf(stderr, "1 %u\t%hhu\t%u\t%u\t%u\t%u\t%u\t%u\n", ovl->qname, ovl->rev, ovl->qs, ovl->qe, ovl->tname, ovl->ts, ovl->te, mincov);
// 		// encode_ovl(fp, encode_tbl, name, ovl, buf);
// 		if (ovl->qname != ovl->tname && ovl->qs <= max_flank_len && ovl->qe + max_flank_len >= cov->qlen) cov->contained ++;
// 	}
// }

static inline void encode_ovl_filter_hq(FILE *fp, FILE *rp, prev_t *name, overlap *ovl, cov_data *cov, buffer_t *buf){
	int lable = 1;
	if (ovl->qname != cov->qname){
		if (cov->qname != UINT32_MAX) ovl_filter(fp, rp, name, cov, buf);
		cov->qlen = kh_val(seeds, kh_get(id, seeds, ovl->qname));
		cov->len = (cov->qlen >> BIN_OFFSET) + 2;
		if (cov->len > cov->max_size){
			cov->max_size = cov->len;
			cov->bins = realloc(cov->bins, sizeof(uint16_t) * cov->max_size);
		}
		lable = cov->pcount = 1;
		memset(cov->bins, 0, sizeof(uint16_t) * cov->len);
		cov->qcov = cov->bincount = cov->binlen = cov->contained = cov->chimera = cov->ovl_i = 0;
		cov->qmaxlen = cov->qlen * MAX_OVL_COV * 6;
	}else if (cov->qcov > cov->qmaxlen || cov->ovl_i > UINT16_MAX - 1000){
		lable = 0;
	}else{
		lable = 1;
		int offset = 1 + (cov->qlen >> (BIN_OFFSET + 1));
		int s = (ovl->qs + 10) >> (BIN_OFFSET + 1);
		int e = (ovl->qe - 10) >> (BIN_OFFSET + 1);
		cov->bins[s] ++;
		// cov->bins[s + offset] ++;
		// cov->bins[e] ++;
		cov->bins[e + offset] ++;
	}

	if (lable) {
		cov->qname = ovl->qname; cov->qs = ovl->qs; cov->qe = ovl->qe;
		// fprintf(stderr, "%u\t%hhu\t%u\t%u\t%u\t%u\t%u\t%u\n", ovl->qname, ovl->rev, ovl->qs, ovl->qe, ovl->tname, ovl->ts, ovl->te, ovl->match);
		cov->qcov += cov->qe - cov->qs + 1;
		if (ovl->qname != ovl->tname && ovl->qs <= max_flank_len && ovl->qe + max_flank_len >= cov->qlen) cov->contained ++;
		cov->ovls[cov->ovl_i ++] = *ovl;
		if (cov->ovl_i >= cov->ovl_m) {
			cov->ovl_m += 500;
			cov->ovls = realloc(cov->ovls, sizeof(overlap) * cov->ovl_m);
		}
	}
}


static inline void encode_ovl_filter(FILE *fp, FILE *rp, prev_t *name, overlap *ovl, cov_data *cov, buffer_t *buf){
	if (is_hq) return encode_ovl_filter_hq(fp, rp, name, ovl, cov, buf);
	int lable = 1;
	int mincov = 200;
	if (ovl->qname != cov->qname){
		if (cov->qname != UINT32_MAX) ovl_filter(fp, rp, name, cov, buf);
		cov->qlen = kh_val(seeds, kh_get(id, seeds, ovl->qname));
		cov->len = (cov->qlen >> BIN_OFFSET) + 1;
		if (cov->len > cov->max_size){
			cov->max_size = cov->len;
			cov->bins = realloc(cov->bins, sizeof(uint16_t) * cov->max_size);
		}
		lable = cov->pcount = 1;
		memset(cov->bins, 0, sizeof(uint16_t) * cov->len);
		cov->qcov = cov->bincount = cov->binlen = cov->contained = cov->chimera = cov->ovl_i = 0;
		cov->qmaxlen = cov->qlen * MAX_OVL_COV;
	}else if (cov->qcov > cov->qmaxlen || cov->ovl_i > UINT16_MAX - 1000){
		lable = 0;
	}else{
		int i, m, n;
		m = n = 0;
		int j = (ovl->qs + 10) >> BIN_OFFSET;
		int k = (ovl->qe - 10) >> BIN_OFFSET;
		if ((j > 15 || k < cov->len - 16) && abs(ovl->qs - cov->qs) < BIN_TOLERANCE_EDGE && \
				abs(ovl->qe - cov->qe) < BIN_TOLERANCE_EDGE){
			lable = cov->pcount++ < BIN_TOLERANCE_COUNT ? 2 : 0;
			// printf("qs:%d qe:%d pcount:%d\n", ovl->qs, ovl->qe, cov->pcount);
		}

		if (lable){
			for (i = j + 1; i <= k; i++){
				if (!cov->bins[i]) n ++;
				if (++cov->bins[i] < mincov) mincov = cov->bins[i];
				if (cov->bins[i] > UINT16_MAX - 1000) cov->bins[i]--;
				m += cov->bins[i];
			}
			if ((mincov > max_bin_cov || (float) m/(k-j) > 1.3 * min(max((float) cov->bincount/cov->binlen, 10), max_bin_cov)) && 
					(ovl->qe - ovl->qs <= cov->qlen * 0.8)){
			// if (mincov > max_bin_cov || ((float) m/(k-j) > 1.3 * min(max((float) cov->bincount/cov->binlen, 10), max_bin_cov) && 
			// 		ovl->qe - ovl->qs <= cov->qlen * 0.8)){
				for (i = j + 1; i <= k; i++) {
					cov->bins[i] --;
				}
				lable = 0;
			}else{
				if (lable != 2) cov->pcount = 1;
				lable = 1;
				cov->binlen += n;
				cov->bincount += k - j;
			}
			// out_cov(cov);
		}
	}

	if (lable) {
		cov->qname = ovl->qname; cov->qs = ovl->qs; cov->qe = ovl->qe;
		// printf("%u\t%hhu\t%u\t%u\t%u\t%u\t%u\n", ovl->qname, ovl->rev, ovl->qs, ovl->qe, ovl->tname, ovl->ts, ovl->te);
		// encode_ovl(fp, encode_tbl, name, ovl, buf);
		cov->qcov += cov->qe - cov->qs + 1;
		if (ovl->qname != ovl->tname && ovl->qs <= max_flank_len && ovl->qe + max_flank_len >= cov->qlen) cov->contained ++;
		cov->ovls[cov->ovl_i ++] = *ovl;
		if (cov->ovl_i >= cov->ovl_m) {
			cov->ovl_m += 500;
			cov->ovls = realloc(cov->ovls, sizeof(overlap) * cov->ovl_m);
		}
	}
}

void sort_ovl(overlap_data *input_data){
	qsort(input_data->ovls, input_data->max_ovls_len, sizeof(overlap), cmp_ovl);
	if (input_data->saved) {
		FILE *fp;
		cov_data cov;
		prev_t name = {0, 0};
		char outfile[MAX_PATH_LENGTH] = {0};
		assert (sprintf(outfile, "%stmp%d.ovl", input_data->outfile, get_outputfile_index()) != -1);
		if ((fp = fopen(outfile, "w")) == NULL) {
			print_log("Failed create overlap file %s!", outfile);
			exit(1);
		}

		buffer_t *buf = init_buffer(ECBUFSIZE);
		init_cov(&cov);
		for (int i = 0; i < input_data->max_ovls_len; i ++) {
			encode_ovl_filter(fp, 0, &name, &input_data->ovls[i], &cov, buf);
		}

		flush_cov(fp, 0, &name, &cov, buf);
		fclose(fp);
		free (cov.bins);
		free (cov.ovls);
		input_data->saved = input_data->ovls[0].qe = input_data->sorted = 0;
	}
	// printf("sort done.\n");
	// input_data->sorted = input_data->saved = 0;
}

void merge_ovl_from_file (sort_args *sort){
	FILE *fp, *rp = sort->save_seed ? sort->rp : 0;
	int i, j, k;
	cov_data cov;
	overlap_data *input_data;
	int total_validy = 0;
	prev_t name = {0, 0};
	uint32_t prev_qname = UINT32_MAX;
	uint32_t ovl[] = {UINT32_MAX, 0, 0, 0};
	overlap ovl_temp = {0, 0, 0, 0, 0, 0, 0};
	typedef struct { 
		int filevalidy;
		int datavalidy;
		int infile_index;
		int ovldata_index;
	} overlap_data_index;
	overlap_data_index *ovl_index = malloc(sort->input_file_step * sizeof(overlap_data_index));

	int sort_input_index = get_inputfile_index(sort);
	for (i = sort_input_index - sort->input_file_step; i < sort->input_file->max_size && i < sort_input_index; i++){
		assert(sort->input_file->files[i].fp = fopen(sort->input_file->files[i].files, "r"));
	}
	// printf("merge_ovl_from_file sort_input_index: %d\n", sort_input_index);

	char outfile[MAX_PATH_LENGTH] = {0};
	assert (sprintf(outfile, "%stmp%d.ovl", sort->outfile, get_outputfile_index()) != -1);
	if ((fp = fopen(outfile, "w")) == NULL) {
		print_log("Failed create overlap file %s!", outfile);
		exit(1);
	}
	
	buffer_t *buf = init_buffer(ECBUFSIZE);
	init_cov(&cov);
	for (i = sort_input_index - sort->input_file_step; i < sort->input_file->max_size && i < sort_input_index; i++){
		ovl_index[total_validy].filevalidy = ovl_index[total_validy].datavalidy = 1;
		ovl_index[total_validy].infile_index = i;
		ovl_index[total_validy].ovldata_index = get_overlapdata_index(sort);
		j = load_input_data_from_file(&sort->input_data[ovl_index[total_validy].ovldata_index], sort->input_file->files[i].fp, \
				sort->input_file->files[i].buf, &sort->input_file->files[i].seed_name);

		if (j) {
			ovl_index[total_validy].filevalidy = 0;
			flush_buffer(NULL, sort->input_file->files[i].buf);
			fclose (sort->input_file->files[i].fp);

			// printf("1del file%s\n", sort->input_file->files[i].files);

			assert (! remove (sort->input_file->files[i].files));
		}
		total_validy ++;
	}

	
	k = total_validy;
	while (total_validy){
		// printf("total_validy: %d\n",total_validy );
		for (i = 0; i < k; i++){
			// printf("1 %d\n", i);
			if (ovl_index[i].datavalidy){
				// printf("# %d \n", i);
				input_data = &sort->input_data[ovl_index[i].ovldata_index];
				if ((input_data->ovls[input_data->len].qname == ovl[0] && (input_data->ovls[input_data->len].match > ovl[1] || \
						(input_data->ovls[input_data->len].match == ovl[1] && input_data->ovls[input_data->len].qe - \
						input_data->ovls[input_data->len].qs < ovl[3]))) || input_data->ovls[input_data->len].qname < ovl[0]){
					ovl[0] = input_data->ovls[input_data->len].qname;
					ovl[1] = input_data->ovls[input_data->len].match;
					ovl[2] = i;
					ovl[3] = input_data->ovls[input_data->len].qe - input_data->ovls[input_data->len].qs;
				}
			}
		}

		input_data = &sort->input_data[ovl_index[ovl[2]].ovldata_index];
		if (sort->save_seed && (prev_qname != input_data->ovls[input_data->len].qname)){
			ovl_temp.rev = 0;
			ovl_temp.qname = ovl_temp.tname = input_data->ovls[input_data->len].qname;
			ovl_temp.qe = ovl_temp.te = kh_val(seeds, kh_get(id, seeds, input_data->ovls[input_data->len].qname)) - 1;
			// out_cov(&cov);
			// if (max_flank_len && prev_qname != UINT32_MAX) out_chi_con(&cov, sort->rp);
			encode_ovl_filter(fp, rp, &name, &ovl_temp, &cov, buf);
			prev_qname = input_data->ovls[input_data->len].qname;
		}
		encode_ovl_filter(fp, rp, &name, &input_data->ovls[input_data->len], &cov, buf);
		prev_qname = input_data->ovls[input_data->len].qname;

		ovl[0] = UINT32_MAX;
		input_data->len ++;
		if (input_data->len >= input_data->max_ovls_len){
			input_data->len = input_data->max_ovls_len = 0;
			if (ovl_index[ovl[2]].filevalidy){
				j = load_input_data_from_file(input_data, sort->input_file->files[ovl_index[ovl[2]].infile_index].fp, \
					sort->input_file->files[ovl_index[ovl[2]].infile_index].buf, &sort->input_file->files[ovl_index[ovl[2]].infile_index].seed_name);
				// if (! input_data->ovls[0].qe){
				// 	total_validy --;
				// input_data->sorted = ovl_index[ovl[2]].datavalidy = 0;
				// 	j = 2;
				// }
				if (j) {
					ovl_index[ovl[2]].filevalidy = 0;
					flush_buffer(NULL, sort->input_file->files[ovl_index[ovl[2]].infile_index].buf);
					fclose (sort->input_file->files[ovl_index[ovl[2]].infile_index].fp);
					assert (! remove (sort->input_file->files[ovl_index[ovl[2]].infile_index].files));
				}
				assert(input_data->ovls[0].qe);
			}else{
				total_validy --;
				input_data->sorted = ovl_index[ovl[2]].datavalidy = 0;
			}
		}
	}
	flush_cov(fp, rp, &name, &cov, buf);
	fclose (fp);
	free (ovl_index);
	free (cov.bins);
	free (cov.ovls);
}

void merge_ovl_from_sort(sort_args *sort){
	FILE *fp, *rp = sort->save_seed ? sort->rp : 0;
	cov_data cov;
	overlap_data *input_data;
	prev_t name = {0, 0};
	uint32_t prev_qname = UINT32_MAX;
	uint32_t ovl[] = {UINT32_MAX, 0, 0, 0};
	overlap ovl_temp = {0, 0, 0, 0, 0, 0, 0};
	int total_validy = sort->input_data_len;

	char outfile[MAX_PATH_LENGTH] = {0};
	assert (sprintf(outfile, "%stmp%d.ovl", sort->outfile, get_outputfile_index()) != -1);
	if ((fp = fopen(outfile, "w")) == NULL) {
		print_log("Failed create overlap file %s!", outfile);
		exit(1);
	}
	
	buffer_t *buf = init_buffer(ECBUFSIZE);
	init_cov(&cov);
	while (total_validy){
		for (int i = 0; i < sort->input_data_len; i++){
			if (sort->input_data[i].ovls[0].qe){
				input_data = &sort->input_data[i];
				if ((input_data->ovls[input_data->len].qname == ovl[0] && (input_data->ovls[input_data->len].match > ovl[1] || \
						(input_data->ovls[input_data->len].match == ovl[1] && input_data->ovls[input_data->len].qe - \
						input_data->ovls[input_data->len].qs < ovl[3]))) || input_data->ovls[input_data->len].qname < ovl[0]){
					ovl[0] = input_data->ovls[input_data->len].qname;
					ovl[1] = input_data->ovls[input_data->len].match;
					ovl[2] = i;
					ovl[3] = input_data->ovls[input_data->len].qe - input_data->ovls[input_data->len].qs;
				}
			}
		}
		if (ovl[0] == UINT32_MAX) break;
		input_data = &sort->input_data[ovl[2]];

		if (sort->save_seed && (prev_qname != input_data->ovls[input_data->len].qname)){
			ovl_temp.rev = 0;
			ovl_temp.qname = ovl_temp.tname = input_data->ovls[input_data->len].qname;
			ovl_temp.qe = ovl_temp.te = kh_val(seeds, kh_get(id, seeds, input_data->ovls[input_data->len].qname)) - 1;
			// out_cov(&cov);
			// if (max_flank_len && prev_qname != UINT32_MAX) out_chi_con(&cov, sort->rp);
			encode_ovl_filter(fp, rp, &name, &ovl_temp, &cov, buf);
			prev_qname = input_data->ovls[input_data->len].qname;
		}

		encode_ovl_filter(fp, rp, &name, &input_data->ovls[input_data->len], &cov, buf);
		prev_qname = input_data->ovls[input_data->len].qname;

		ovl[0] = UINT32_MAX;
		input_data->len ++;
		if (input_data->len >= input_data->max_ovls_len){
			total_validy --;
			input_data->ovls[0].qe = 0;
		}
	}
	flush_cov(fp, rp, &name, &cov, buf);
	fclose (fp);
	free (cov.bins);
	free (cov.ovls);
}

void sort_ovl_file(sort_args *sort){
	FILE *fp;
	char *inputfile;
	uint32_t ovl[8];
	prev_t name;
	int key_exist1, key_exist5, input_file_index, input_data_index = -1;
	overlap_data *input_data = NULL;
	buffer_t *buf = init_buffer(DCBUFSIZE);
	khint_t k;

	while (1){
		input_file_index = get_inputfile_index(sort);
		if (input_file_index > sort->input_file->len) break;
		inputfile = sort->input_file->files[input_file_index - 1].files;

		if ((fp = fopen(inputfile, "r")) == NULL) {
			print_log("Failed open overlap file %s!", inputfile);
			exit(1);
		}
		// printf("sort_ovl_filefileindex:%d, file:%s\n", input_file_index, inputfile);
	
		buf->buffer_i = buf->read_n = 0;
		name.prev_qname = name.prev_tname = 0;
		key_exist1 = key_exist5 = 5;//max error record
		while(decode_ovl(fp, decode_tbl, &name, ovl, buf, 8) >= 0) {
			// printf("%u %u %u %u %u %u %u\n", ovl[0], ovl[1], ovl[2], ovl[3],ovl[4],ovl[5],ovl[6]);
			if (ovl[0] == ovl[4]) continue;
			if (ovl[3] - ovl[2] < 500 || ovl[6] - ovl[5] < 500) continue;

			if (input_data_index == -1){
				input_data_index = get_overlapdata_index(sort);
				input_data = &sort->input_data[input_data_index];
			}

			if (input_data->max_ovls_len >= input_data->max_size - 2){
				input_data->saved = 1;
				thpool_add_work(sort->thpool, (void*)sort_ovl, input_data);
				// printf("sort_ovl begin\n");
				input_data_index = get_overlapdata_index(sort);
				input_data = &sort->input_data[input_data_index];
				// printf("%d %d\n", input_data_index, input_data->max_ovls_len);
			}
			//TDDO remove here check
			if (key_exist1 && (k = kh_get(id, seeds, ovl[0])) != kh_end(seeds) && kh_val(seeds, k) >= ovl[3]) {
				input_data->ovls[input_data->max_ovls_len].qname = ovl[0];
				input_data->ovls[input_data->max_ovls_len].rev = ovl[1];
				input_data->ovls[input_data->max_ovls_len].qs = ovl[2];
				input_data->ovls[input_data->max_ovls_len].qe = ovl[3] - 1;//0-based of paf
				input_data->ovls[input_data->max_ovls_len].tname = ovl[4];
				input_data->ovls[input_data->max_ovls_len].ts = ovl[5];
				input_data->ovls[input_data->max_ovls_len].te = ovl[6] - 1;
				input_data->ovls[input_data->max_ovls_len].match = ovl[7];
				input_data->max_ovls_len ++;
			}else if (key_exist1){
				key_exist1 --;
			}

			if (key_exist5 && (k = kh_get(id, seeds, ovl[4])) != kh_end(seeds) && kh_val(seeds, k) >= ovl[6]) { // else if 1 if 2
				input_data->ovls[input_data->max_ovls_len].qname = ovl[4];
				input_data->ovls[input_data->max_ovls_len].rev = ovl[1];
				input_data->ovls[input_data->max_ovls_len].qs = ovl[5];
				input_data->ovls[input_data->max_ovls_len].qe = ovl[6] - 1;
				input_data->ovls[input_data->max_ovls_len].tname = ovl[0];
				input_data->ovls[input_data->max_ovls_len].ts = ovl[2];
				input_data->ovls[input_data->max_ovls_len].te = ovl[3] - 1;
				input_data->ovls[input_data->max_ovls_len].match = ovl[7];
				input_data->max_ovls_len ++;
			}else if (key_exist5){
				key_exist5 --;
			}
		}
		fclose(fp);
	}
	if (input_data_index != -1){
		thpool_add_work(sort->thpool, (void*)sort_ovl, input_data);
	}
	flush_buffer(NULL, buf);
}

static int usage()
{
	fprintf(stderr, "Usage: ovl_sort [options] input.fofn\n");
	fprintf(stderr, "Sort and merge ovl files\n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -i  query idx file [required]\n");
	fprintf(stderr, "  -H  high quality reads, such as HiFi [0]\n");
	fprintf(stderr, "  -m  set max total available buffer size, suffix K/M/G [40G]\n");
	fprintf(stderr, "  -t  number of threads to use [8]\n");
	fprintf(stderr, "  -k  max depth of each overlap [40]\n");
	fprintf(stderr, "  -l  max over hang length to filter [%d]\n", max_flank_len);
	fprintf(stderr, "  -o  output file name [required]\n");
	fprintf(stderr, "  -d  temporary directory [$CWD]\n");
	return 1;
}

int main(int argc, char *argv[])
{	
	opt opts;
	opts.idx = opts.outfile = NULL;
	opts.thread = 8;
	opts.max_mem = 40l << 30;
	getcwd(opts.tempdir, sizeof(opts.tempdir));

	int c;
	while((c = getopt(argc, argv, "i:m:t:k:o:d:l:H")) != -1) {
		switch(c) {
			case 'i': opts.idx = strdup(optarg); break;
			case 't': opts.thread = atoi(optarg); assert (opts.thread >= 2); break;
			case 'k': max_bin_cov = atoi(optarg); break;
			case 'o': opts.outfile = strdup(optarg); break;
			case 'd': strcpy(opts.tempdir, optarg); break;
			case 'l': max_flank_len = atoi(optarg); break;
			case 'H': is_hq = 1; break;
			case 'm': {
				char *q;
				opts.max_mem = strtoll(optarg, &q, 0);
				if (*q == 'k' || *q == 'K') opts.max_mem <<= 10;
				else if (*q == 'm' || *q == 'M') opts.max_mem <<= 20;
				else if (*q == 'g' || *q == 'G') opts.max_mem <<= 30;
				assert (opts.max_mem > 2.2 * MIN_MEM_UNIT);
				break;
			}
			default: return usage();
		}
	}

	if (optind + 1 > argc || !(opts.outfile || opts.idx)) return usage();

	int read_thread = min((opts.thread>>2) + 1, MAX_READ_THREAD);//need test one read thread ?= sort_thread
	int sort_thread = min(opts.thread - read_thread, MAX_SORT_THREAD);
	opts.max_mem = opts.max_mem / 5 / (read_thread + sort_thread) * 4;
	opts.thread = read_thread + sort_thread;
	assert (opts.max_mem > MIN_MEM_UNIT);
	// printf("read_thread:%d sort_thread:%d max_mem:%lu depth:%d\n", read_thread, sort_thread, opts.max_mem, max_bin_cov);

	encode_tbl = init_encode_table(NUM_LIMIT);
	decode_tbl = init_decode_table(NUM_LIMIT);
	seeds = kh_init(id);
	init_seeds(opts.idx);

	threadpool thpool_read = thpool_init(read_thread);
	threadpool thpool_sort = thpool_init(sort_thread);
	assert (!pthread_mutex_init(&inputlock, NULL));
	assert (!pthread_mutex_init(&outputlock, NULL));
	assert (!pthread_mutex_init(&sortlock, NULL));

	overlap_file ovlfiles;
	init_ovlfiles(&ovlfiles, argv[argc - 1]);

	sort_args sort;
	int out_temp_suffix = 1;
	init_sort(&sort, 1, thpool_sort, &ovlfiles, &opts, opts.thread, opts.max_mem, 1, out_temp_suffix++);
	sort.rp = init_rp(&opts);

	for (int i = 0; i < read_thread; i ++){
		thpool_add_work(thpool_read, (void*)sort_ovl_file, &sort);
	}
	// printf("read thpool_add_work done\n");
	
	thpool_wait_all(thpool_read);
	thpool_wait_all(thpool_sort);
	thpool_destroy(thpool_sort);
	thpool_destroy(thpool_read);
	free (ovlfiles.files);

	if (sorted_files_len) sort.save_seed = 0;
	merge_ovl_from_sort(&sort);
	destroy_sort_input(&sort);

	while (sorted_files_len > 1){
		merge_thread_opt mt_opt;
		init_mt_opt(&mt_opt, &opts, sorted_files_len);

		ovlfiles.max_size = sorted_files_len;
		ovlfiles.len = sorted_files_len = 0;
		ovlfiles.files = calloc(ovlfiles.max_size, sizeof(file));

		for (int i = 0; i < ovlfiles.max_size; i++){
			assert (sprintf(ovlfiles.files[i].files, "%stmp%d.ovl", sort.outfile, i + 1) != -1);
			// assert(ovlfiles.files[i].fp = fopen(ovlfiles.files[i].files, "r"));
			ovlfiles.files[i].buf = init_buffer(DCBUFSIZE);
		}

		threadpool thpool_merge = thpool_init(mt_opt.max_thread);
		sort.save_seed = mt_opt.max_thread > 1 ? 0 : 1;
		init_sort(&sort, sort.save_seed, thpool_merge, &ovlfiles, &opts, mt_opt.input_data_len, \
			mt_opt.max_mem_each_file, mt_opt.max_readfile_in_single_thread, out_temp_suffix++);

		for (int i = 0; i < ovlfiles.max_size; i += sort.input_file_step){
			thpool_add_work(thpool_merge, (void*)merge_ovl_from_file, &sort);
			// printf("thpool_merge: %d\n", i);
		}

		thpool_wait_all(thpool_merge);
		thpool_destroy(thpool_merge);
		destroy_sort_input(&sort);
		free(ovlfiles.files);
	}

	char mv[2048] = {0};
	assert (sprintf(mv, "mv %stmp%d.ovl %s", sort.outfile, get_outputfile_index() - 1, opts.outfile) != -1);
	assert(!system(mv));
	if (sort.rp) fclose(sort.rp);
	pthread_mutex_destroy(&inputlock);
	pthread_mutex_destroy(&outputlock);
	pthread_mutex_destroy(&sortlock);
	kh_destroy(id, seeds);
	free(opts.idx); free(opts.outfile);
	free(encode_tbl); free(decode_tbl);
	return 0;
}
