#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h> //getpid
#include <limits.h>

#ifndef LGS_CORRECT
	#define LGS_CORRECT
#endif

#include "ctg_cns.h"
#include "bsort.h"
#include "htslib/bgzf.h"
#include "bseq.h"

#define MAX_PP_PPP 0
#define READS_ONT 1
#define READS_CLR 2
#define READS_HIFI 3
#define READS_RS 4 //pb reads from RS

// #include <time.h>
// static clock_t start;
// static int time_i;
// static inline void time_debug(clock_t t, char *string){
// 	fprintf (stderr, "pid:%d time: %f function:%s seq: %d\n", getpid(), ((double)(t - start))/ CLOCKS_PER_SEC, string, time_i);
// 	fflush(stderr);
// 	start = t;
// }
// static inline void debug(char *string){
// 	fprintf (stderr, "pid:%d log:%s\n", getpid(), string);
// 	fflush(stderr);
// }

// unsigned int lqseq_max_length;

static int READS_TYPE = READS_ONT;
static int GAP_MIN_LEN;
static float GAP_MIN_RATIO1;
static float MAX_CLIP_RATIO;

static uint8_t int_to_base[] = { 
	65, 84, 71, 67, 45, 78, 77 //A, T, G, C, -, N, M
};

static align_tag align_tag_head = {
	.t_pos = -1, 
	.q_base = 0,
	.delta = 0
};

static uint8_t base_to_int[] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//0-15
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//16-31
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//32-47
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//48-63
	4, 0, 4, 3,  4, 4, 4, 2,  4, 4, 4, 4,  4, 6, 5, 4,//64-79
	4, 4, 4, 4,  1, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//80-95
	4, 0, 4, 3,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,//96-111
	4, 4, 4, 4,  1, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4//112-127
};

#define MAX_WAIT_TIME 3600
#define MAX_REMAIN_RAM 10240000000L

void *smalloc(size_t size){
    void *w = NULL, *v = NULL;
    unsigned int t = 0;
    while (t < MAX_WAIT_TIME){
        v = malloc(size);
        w = malloc(MAX_REMAIN_RAM);
        if (v != NULL && w != NULL) {
        	free (w);
        	return v;
        }else{
        	if (v) free (v);
        	if (w) free (w);
        	v = w = NULL;
        }
        t += 10;
        sleep (10);
    }
    return v;
}

void *scalloc(size_t nitems, size_t size){
    void *w = NULL, *v = NULL;
    unsigned int t = 0;
    while (t < MAX_WAIT_TIME){
        v = calloc(nitems, size);
        w = malloc(MAX_REMAIN_RAM);
        if (v != NULL && w != NULL) {
        	free (w);
        	return v;
        }else{
        	if (v) free (v);
        	if (w) free (w);
        	v = w = NULL;
        }
        t += 10;
        sleep (10);
    }
    return v;
}

static void reverse_consensus_base(consensus_base *bases, unsigned int len)
{
	consensus_base tmp;
	consensus_base *p1 = bases;
	consensus_base *p2 = bases + len - 1;
	while (p1 < p2) {
		tmp = *p1;
		*p1++ = *p2;
		*p2-- = tmp;
	}
}

// static void out_align (alignment *aln, char *s){
// 	if (s) printf("%s ", s);
// 	printf("%d %d %d ",aln->aln_t_s, aln->aln_t_e, aln->aln_len);
// 	for (int i = 0; i < aln->aln_len; i++){
// 		printf("%c", aln->t_aln_str[i + aln->shift]);
// 	}
// 	printf("\n");
// 	if (s) printf("%s ", s);
// 	printf("%d %d %d ",aln->aln_t_s, aln->aln_t_e, aln->aln_len);
// 	for (int i = 0; i < aln->aln_len; i++){
// 		printf("%c", aln->q_aln_str[i + aln->shift]);
// 	}
// 	printf("\n");
// }

static inline void get_align_shift(alignment *aln, const int k, const int l){

	int i, j;
	i = 0;
	j = 0;
	while (i < aln->aln_len){
		if (aln->t_aln_str[i] == aln->q_aln_str[i]){
			j ++;
		}else{
			j = 0;
		}

		if (aln->t_aln_str[i] != '-'){
			aln->aln_t_s ++;//true start + 1
		}

		if (l && aln->q_aln_str[i] != '-'){
			aln->aln_q_s ++;
		}

		if (j == k){
			aln->aln_t_s -= k;
			aln->shift = i - k + 1;
			aln->aln_len = aln->aln_len - i + k - 1;
			if (l) aln->aln_q_s -= k;
			break;
		}
		i ++;
	}
	
	if (j == k){

		i = aln->aln_len + i - k;
		j = 0;
		int t = 0;
		while (i >= 0){
			if (aln->t_aln_str[i] == aln->q_aln_str[i]){
				j ++;
			}else{
				j = 0;
			}

			if (aln->t_aln_str[i] != '-'){
				aln->aln_t_e --;//true end - 1
			}

			if (l && aln->q_aln_str[i] != '-'){
				aln->aln_q_e --;
			}

			if (j == k){
				aln->aln_t_e += k;
				aln->aln_len = aln->aln_len - t + k - 1;
				if (l) aln->aln_q_e += k;
				break;
			}
			i --;
			t ++;
		}
	}else{
		aln->aln_len = 0;
	}
}


static void reallocate_cons_trimed_mem(consensus_trimed *cons_trimed, int *cons_trimed_seq_max_len, const int number){
	*cons_trimed_seq_max_len += number;
	cons_trimed->seq = realloc(cons_trimed->seq, *cons_trimed_seq_max_len * sizeof (char));
}


static inline void reallocate_msa_p_d_b_pp_ppp_mem(msa_p_d_b *msa, int re_maxp){
	if (msa->max_size > MAX_PP_PPP){
		msa->max_size += re_maxp;
		msa->pp_ppp = realloc(msa->pp_ppp, sizeof(msa_p_d_b_pp_ppp) * msa->max_size);
	}else{
		msa_p_d_b_pp_ppp * const _ = msa->pp_ppp;
		msa->max_size += re_maxp;
		msa->pp_ppp = malloc(sizeof(msa_p_d_b_pp_ppp) * msa->max_size);
		memcpy(msa->pp_ppp, _, sizeof(msa_p_d_b_pp_ppp) * msa->len);
	}
}

static inline void allocate_msa_mem(msa_p *msa, int len){
	int p, d, b;
	uint64_t max_size;
	int b_ = 6;
	for (p = 0; p < len; p++){
		max_size = msa[p].max_size;
		msa[p].d_b = malloc(max_size * sizeof(msa_p_d_b *) + \
			max_size * b_ * sizeof(msa_p_d_b) + \
			max_size * b_ * MAX_PP_PPP * sizeof(msa_p_d_b_pp_ppp)
		);
		assert (msa[p].d_b);
		
		msa_p_d_b * const _ = (msa_p_d_b *) (msa[p].d_b + max_size);
		msa_p_d_b_pp_ppp * const __ = (msa_p_d_b_pp_ppp *) (_ + max_size * b_);

		for (d = 0; d < max_size; d++){
			msa[p].d_b[d] = _ + d * b_;
			for (b = 0; b < b_; b++){
				msa[p].d_b[d][b].max_size = MAX_PP_PPP;
				msa[p].d_b[d][b].len = 0;
				msa[p].d_b[d][b].pp_ppp = __ + (d * b_ + b) * msa[p].d_b[d][b].max_size;
			}
		}
	}
}

static inline void free_msa(msa_p *msa, int len){
	int b_ = 6;
	for (int p = 0; p < len; p ++){
		for (int d = 0; d < msa[p].max_size; d++){
			for (int b = 0; b < b_; b++){
				if (msa[p].d_b[d][b].len > MAX_PP_PPP) free (msa[p].d_b[d][b].pp_ppp);
			}
		}
		free(msa[p].d_b);
	}
	free (msa);
}


static inline void allocate_lqseq_mem(lqseq *lqseqs, const int len){
	int i, j;
	lqseq *lqseq;
	for (i = 0; i < len; i++){
		lqseq = &lqseqs[i];
		lqseq->sudoseed = NULL;
		lqseq->lqcount = lqseq->len = lqseq->sudoseed_len = 0;
		lqseq->seqs = malloc(LQSEQ_MAX_CAN_COUNT * sizeof(struct seq_) + 
			(LQSEQ_INIT_LEN + 1) * LQSEQ_MAX_CAN_COUNT * sizeof(char));
		assert (lqseq->seqs);
		char * const _ = (char *) (lqseq->seqs + LQSEQ_MAX_CAN_COUNT);
		for (j = 0; j < LQSEQ_MAX_CAN_COUNT; j ++){
			lqseq->seqs[j].len = LQSEQ_INIT_LEN;
			lqseq->seqs[j].seq = _ + j * (LQSEQ_INIT_LEN + 1);
		}
	}
}

static inline void reallocate_lqseq_mem(struct seq_ *lqseq){
	if (lqseq->len > LQSEQ_INIT_LEN){
		lqseq->len += LQSEQ_INIT_LEN;
		lqseq->seq = realloc(lqseq->seq, sizeof(char) * (lqseq->len + 1));
	}else{
		char * const _ = lqseq->seq;
		lqseq->len += LQSEQ_INIT_LEN;
		lqseq->seq = malloc(sizeof(char) * (lqseq->len + 1));
		memcpy(lqseq->seq, _, sizeof(char) * (lqseq->len - LQSEQ_INIT_LEN + 1));
	}
}

static inline void free_lqseq(lqseq *lqseqs, const int len){
	int i, j;
	for (i = 0; i < len; i++){
		for (j = 0; j < LQSEQ_MAX_CAN_COUNT; j ++){
			if (lqseqs[i].seqs[j].len > LQSEQ_INIT_LEN) free(lqseqs[i].seqs[j].seq);
		}
		free (lqseqs[i].seqs);
		if (lqseqs[i].sudoseed) free (lqseqs[i].sudoseed);
	}
	free (lqseqs);
}

static inline int get_align_tag(unsigned int *p, align_tags_t *tags, align_tag *tag){
	uint8_t t = tags->align_tags[*p>>1];
	if (!(*p & 1)) t >>= 4;
	if ((t & 15) == 15) return -1;
	tag->q_base = t & 7;
	if ((*p)++){
		if(t & 8){
			tag->delta++;
		}else{
			tag->delta = 0;
			tag->t_pos ++;
		}
	}else{
		tag->t_pos = tags->aln_t_s;
		tag->delta = 0;
	}
	return 1;
}


static inline void update_msa(msa_p * msa, align_tags_t *tags_list, int aligned_seq_count, int lable){
	unsigned int i, d, p, updated;
	align_tag p1, pp, ppp;
	for (p = 0; p < aligned_seq_count; p++){
		d = 0;
		pp = ppp = align_tag_head;
		while (get_align_tag(&d, &tags_list[p], &p1) > 0){
			if (p1.q_base == 6 || pp.q_base == 6) {
				ppp = pp;
				pp = p1;
				continue;
			}
			msa_p_d_b *p_base = &msa[p1.t_pos].d_b[p1.delta][p1.q_base];
			for (i = 0, updated = 0; i < p_base->len; i++){
				if (p_base->pp_ppp[i].pp.t_pos == pp.t_pos &&
					p_base->pp_ppp[i].pp.delta == pp.delta &&
					p_base->pp_ppp[i].pp.q_base == pp.q_base &&
					p_base->pp_ppp[i].ppp.t_pos == ppp.t_pos &&
					p_base->pp_ppp[i].ppp.delta == ppp.delta && 
					p_base->pp_ppp[i].ppp.q_base == ppp.q_base){
					p_base->pp_ppp[i].link_count ++;
					updated = 1;
					break;
				}
			}
			if (! updated){
				if (i >= p_base->max_size){
					reallocate_msa_p_d_b_pp_ppp_mem(p_base, 10);
				}
				p_base->pp_ppp[p_base->len].pp = pp;
				p_base->pp_ppp[p_base->len].ppp = ppp;
				p_base->pp_ppp[p_base->len].link_count = 1;
				p_base->pp_ppp[p_base->len].score = 0;
				p_base->len ++;
			}
			ppp = pp;
			pp = p1;
		}
		if (lable) free (tags_list[p].align_tags);
	}
	if (lable) free (tags_list);
}	

static int compare_seq_by_kscore(const void * a, const void * b){
	struct seq_ *lqseq1 = (struct seq_ *) a;
	struct seq_ *lqseq2 = (struct seq_ *) b;
	return (lqseq1->kscore < lqseq2->kscore) - (lqseq1->kscore > lqseq2->kscore);
}

static int compare_seq_by_len(const void * a, const void * b){
	struct seq_ *lqseq2 = (struct seq_ *) a;
	struct seq_ *lqseq1 = (struct seq_ *) b;
	return (lqseq1->len < lqseq2->len) - (lqseq1->len > lqseq2->len);
}

static int rev_compare_seq_by_len(const void * b, const void * a){
	struct seq_ *lqseq2 = (struct seq_ *) a;
	struct seq_ *lqseq1 = (struct seq_ *) b;
	return (lqseq1->len < lqseq2->len) - (lqseq1->len > lqseq2->len);
}

static consensus_data *error_seed(int len){
	consensus_data *consensus = malloc(sizeof(consensus_data));
	consensus->len = len;
	consensus->cns_bases = malloc(consensus->len * sizeof(consensus_base));
	return consensus;
}

static void find_ref_lqseq(lqseq *lq){
	int j;
	struct seq_ tmp;
	for (j = 1; j < lq->len; j ++){
		if (!lq->seqs[j].order){
			tmp = lq->seqs[0];
			lq->seqs[0] = lq->seqs[j];
			lq->seqs[j] = tmp;
			break;
		}
	}
}

static void count_kmers(lqseq *lq, uint16_t *kmers, int c, int l){
	int j, k, s, index;
	char *seq;
	uint16_t kmer;
	memset(kmers, 0, sizeof(uint16_t) * KMER_LEN_COUNT);
	for (j = 0; j < min(lq->len, c); j ++){
		if (lq->seqs[j].len < KMER_LEN) continue;
		// s = l ? lq->seqs[j].len - KMER_RANGE : 0;
		s = l && lq->seqs[j].len > KMER_RANGE ? lq->seqs[j].len - KMER_RANGE : 0;
		for (kmer = k = 0; k < min(lq->seqs[j].len, KMER_RANGE) - KMER_LEN; k ++){
			if (k){
				kmer = kmer << 2 | base_to_int[(unsigned char) lq->seqs[j].seq[s + k + KMER_LEN - 1]];
			}else{
				seq = lq->seqs[j].seq;
				for (index = 0; index < KMER_LEN; index ++){
					kmer = kmer << 2 | base_to_int[(unsigned char) seq[s + k + index]];
				}
			}
			kmers[kmer] ++;
		}
	}
}

static void count_kscore(lqseq *lq, uint16_t *kmers, int l){
	int j, k, s, index;
	char *seq;
	uint16_t kmer;
	for (j = 0; j < lq->len; j ++){
		lq->seqs[j].kscore = 0;
		if (lq->seqs[j].len < KMER_LEN) continue;
		// s = l ? lq->seqs[j].len - KMER_RANGE : 0;
		s = l && lq->seqs[j].len > KMER_RANGE ? lq->seqs[j].len - KMER_RANGE : 0;
		for (kmer = k = 0; k < min(lq->seqs[j].len, KMER_RANGE) - KMER_LEN; k ++){
			if (k){
				kmer = kmer << 2 | base_to_int[(unsigned char) lq->seqs[j].seq[s + k + KMER_LEN - 1]];
			}else{
				seq = lq->seqs[j].seq;
				for (index = 0; index < KMER_LEN; index ++){
					kmer = kmer << 2 | base_to_int[(unsigned char) seq[s + k + index]];
				}
			}
			lq->seqs[j].kscore += kmers[kmer];
		}
	}
}

static void save_kscore(lqseq *lq, uint16_t *score){
	int j;
	for (j = 0; j < lq->len; j ++){
		score[lq->seqs[j].order] = lq->seqs[j].kscore;
	}
}

static void merge_kscore(lqseq *lq, uint16_t *score){
	int j;
	for (j = 0; j < lq->len; j ++){
		lq->seqs[j].kscore += score[lq->seqs[j].order];
	}
}

// static int generate_lqseqs_from_tags1(lqseq *lqseqs, const int lqseqs_count, align_tags_t *tags_list, 
// 	const unsigned int aligned_seq_count){

// 	unsigned int i, j, k, p, step;
// 	unsigned int index, max_aln_lqseq_len, max_aln_length;

// 	align_tags_t *tags;
// 	align_tag tag = {0, 0, 0};
// 	max_aln_length = 0;
// 	step = 1;

// 	lqseq *lqseq;
// 	for (i = lqseqs_count - 1; i < lqseqs_count; i--){
// 		max_aln_lqseq_len = 0;
// 		lqseq = &lqseqs[i];
// 		lqseq->len = 0;
// 		lqseq->sudoseed = NULL;
// 		lqseq->seqs = malloc(LQSEQ_MAX_CAN_COUNT * sizeof(struct seq_));
// 		for (j = 0; j < LQSEQ_MAX_CAN_COUNT; j ++){
// 			lqseq->seqs[j].len = 500;
// 			lqseq->seqs[j].seq = malloc(lqseq->seqs[j].len + 1);
// 		}
// 		for (k = 0, j = step; j < aligned_seq_count; j++){ //discard seed
// 			tags = &tags_list[j];
// 			if (tags->aln_t_s <= lqseq->start && tags->aln_t_e >= lqseq->end){
// 				if (!k++ && j > 1) step = j - 1;// need update if step startswith 0
// 				index = p = 0;
// 				while (get_align_tag(&p, tags, &tag) > 0 && tag.t_pos <= lqseq->end){
// 					if (tag.t_pos >= lqseq->start && tag.q_base != 4){
// 						lqseq->seqs[lqseq->len].seq[index++] = int_to_base[tag.q_base];
// 						if (index >= lqseq->seqs[lqseq->len].len){
// 							lqseq->seqs[lqseq->len].len += 500;
// 							lqseq->seqs[lqseq->len].seq = realloc(lqseq->seqs[lqseq->len].seq, 
// 								lqseq->seqs[lqseq->len].len);
// 						}
// 					}
// 				}
// 				if (index > lqseq->end - lqseq->start + 1){
// 					lqseq->seqs[lqseq->len].seq[index] = '\0';
// 					lqseq->seqs[lqseq->len].len = index;
// 					lqseq->seqs[lqseq->len].order = lqseq->len;
// 					lqseq->len ++;
// 					if (index > max_aln_lqseq_len) max_aln_lqseq_len = index;
// 				}
// 				if (lqseqs[i].len >= LQSEQ_MAX_CAN_COUNT) break;
// 			}
// 		}
// 		if (lqseq->len <= 4) {
// 			lqseq->len = 0;
// 		}else{
// 			qsort(lqseq->seqs, lqseq->len, sizeof(struct seq_), compare_seq_by_len);
// 			k = lqseq->len/2;
// 			while (lqseq->seqs[lqseq->len - 1].len > 2 * lqseq->seqs[k].len || 
// 				lqseq->seqs[lqseq->len - 1].len >= 1.4 * lqseq->seqs[lqseq->len - 2].len) lqseq->len--;
// 			uint16_t kmers[KMER_LEN_COUNT];
// 			// count_kmers(lqseq, kmers, KMER_MAX_SEQ, 0);
// 			count_kmers(lqseq, kmers, LQSEQ_MAX_CAN_COUNT, 0);
// 			count_kscore(lqseq, kmers, 0);
// 			unsigned int klastscore, kmaxscore;
// 			unsigned int kmaxlen = lqseq->seqs[0].len;
// 			if (kmaxlen > 100){
// 				uint16_t score[LQSEQ_MAX_CAN_COUNT];
// 				save_kscore(lqseq, score);
// 				// count_kmers(lqseq, kmers, KMER_MAX_SEQ, 1);
// 				count_kmers(lqseq, kmers, LQSEQ_MAX_CAN_COUNT, 1);
// 				count_kscore(lqseq, kmers, 1);
// 				merge_kscore(lqseq, score);
// 			}

// 			qsort(lqseq->seqs, lqseq->len, sizeof(struct seq_), compare_seq_by_kscore);
// 			kmaxlen = lqseq->seqs[0].len;
// 			klastscore = kmaxscore = lqseq->seqs[0].kscore;

// 			for (k = j = 0; j < lqseq->len; j ++){
// 				if (lqseq->seqs[j].kscore * 10 < kmaxscore || j >= LQSEQ_MAX_COUNT || 
// 					lqseq->seqs[j].kscore * 2 < klastscore) break;
// 				klastscore = lqseq->seqs[j].kscore;
// 				if (j < KMER_MAX_SEQ && lqseq->seqs[j].kscore > kmaxscore * 0.8 && lqseq->seqs[j].len > kmaxlen){
// 					kmaxlen = lqseq->seqs[j].len;
// 					k = j;
// 				}
// 			}

// 			lqseq->indexs = 0;
// 			lqseq->indexe = kmaxlen > LQSEQ_MAX_REV_LEN && j > 6 ? 5 : j - 1;
// 			if (lqseq->indexe - lqseq->indexs <= 3) {
// 				lqseq->len = 0;
// 				continue;
// 			}

// 			if (lqseqs[i].seqs[0].len < 3000){
// 				j = lqseqs[i].indexs;
// 				k = j + 6 < lqseqs[i].indexe ? 6 : lqseqs[i].indexe - j + 1;
// 			}else{
// 				j = lqseqs[i].indexs;
// 				k = j + 2 < lqseqs[i].indexe ? 2 : lqseqs[i].indexe - j + 1;
// 			}
// 			lqseqs[i].sudoseed = poa_to_consensus(&lqseqs[i].seqs[j], k);
// 			// lqseqs[i].sudoseed = strdup(lqseqs[i].seqs[k].seq);
// 			lqseqs[i].sudoseed_len = strlen(lqseqs[i].sudoseed);

// 			if (max_aln_lqseq_len + lqseqs[i].sudoseed_len > max_aln_length) {
// 				max_aln_length = max_aln_lqseq_len + lqseqs[i].sudoseed_len;
// 			}
// 		}
// 		// time_debug(clock(), "lq produced... ");
// 	}

// 	// for (i = 0 ; i < lqseqs_count; i++){
// 	// 	if (lqseqs[i].len){
// 	// 		for (k = 0; k < lqseqs[i].len; k++){
// 	// 			j = k >= lqseqs[i].indexs && k <= lqseqs[i].indexe ? 1 : 0;
// 	// 			printf("%d %d %d %d s:%d e:%d order:%d kscore:%d len:%d %s\n",lqseqs_count, i, k, j, lqseqs[i].start, 
// 	// 				lqseqs[i].end, lqseqs[i].seqs[k].order, lqseqs[i].seqs[k].kscore, strlen(lqseqs[i].seqs[k].seq), lqseqs[i].seqs[k].seq);
// 	// 		}
// 	// 		printf("seed%d: %d %s\n",i, lqseqs[i].sudoseed_len, lqseqs[i].sudoseed);
// 	// 	}
// 	// }
// 	// exit(1);
// 	for (i = 0; i < aligned_seq_count; i++){
// 		free (tags_list[i].align_tags);
// 	}
// 	free (tags_list);
// 	return max_aln_length;
// }


static void generate_lqseqs_from_cluster(lqseq *lqseq, gap_cluster *clu){
	uint32_t i, s, index;
	for (i = 0; i < clu->i_m && lqseq->len < LQSEQ_MAX_CAN_COUNT; i++){
		if (clu->gap[i]->l != 2) continue;
		for (index = 0, s = clu->gap[i]->gap.s; s < clu->gap[i]->gap.e; s++){
			lqseq->seqs[lqseq->len].seq[index++] = seq_nt16_str[bam_seqi(clu->gap[i]->dseq, s)];
			if (index > lqseq->seqs[lqseq->len].len) reallocate_lqseq_mem(&lqseq->seqs[lqseq->len]);
		}
		lqseq->seqs[lqseq->len].seq[index] = '\0';
		lqseq->seqs[lqseq->len].len = index;
		lqseq->seqs[lqseq->len].order = lqseq->len;
		if (index > lqseq->lqcount) lqseq->lqcount = index;
		lqseq->len ++;
	}
	assert (lqseq->start == clu->r.s);
}

static void reverse_lqseq(lqseq *lqseq){
	struct seq_ tmp;
	struct seq_ *p1 = lqseq->seqs;
	struct seq_ *p2 = lqseq->seqs + lqseq->len - 1;
	while (p1 < p2) {
		tmp = *p1;
		*p1++ = *p2;
		*p2-- = tmp;
      }
}

static void remove_short_lqseq(lqseq *lqseq){
	int k;
	qsort(lqseq->seqs, lqseq->len, sizeof(struct seq_), rev_compare_seq_by_len);
	k = lqseq->len/4;
	while (lqseq->len > k && (lqseq->seqs[lqseq->len - 1].len < lqseq->seqs[k].len/2 || 
		lqseq->seqs[lqseq->len - 1].len * 1.4 < lqseq->seqs[lqseq->len - 2].len)) lqseq->len--;
	if (k == lqseq->len) lqseq->len = 0;
	if (lqseq->len > LQSEQ_MAX_COUNT) lqseq->len = LQSEQ_MAX_COUNT;
	reverse_lqseq(lqseq);
	// for (i = 0; i < lqseq->len; i++){
	// 	printf("%d s:%d e:%d order:%d len:%d %s\n", i, lqseq->start, 
	// 		lqseq->end, lqseq->seqs[i].order, strlen(lqseq->seqs[i].seq), lqseq->seqs[i].seq);
	// }
}


static int generate_lqseqs_from_tags_kmer(lqseq *lqseqs, const int lqseqs_count, align_tags_t *tags_list, \
	const unsigned int aligned_seq_count, const gap_clusters *clusters){

	allocate_lqseq_mem(lqseqs, lqseqs_count);
	align_tags_t *tags;
	align_tag tag = {0, 0, 0};
	int align_tags_len, align_tags_max = 1000000;
	align_tag *align_tags = malloc(sizeof(align_tag) * align_tags_max);
	
	lqseq *lqseq;
	unsigned int i, p;
	int j, k, index, s = lqseqs_count - 1;
	for (i = 1; i < aligned_seq_count; i ++){
		tags = &tags_list[i];
		while (s >= 0 && (lqseqs[s].start < tags->aln_t_s || lqseqs[s].len >= LQSEQ_MAX_CAN_COUNT)) s--;
		for (j = s; j >= 0 && lqseqs[j].end <= tags->aln_t_e; j--);
		if (j == s) continue;

		align_tags_len = p = 0;
		while (get_align_tag(&p, tags, &tag) > 0){
			align_tags[align_tags_len++] = tag;
			if (tag.t_pos <= lqseqs[j + 1].end){
				if (align_tags_len >= align_tags_max){
					align_tags_max += 1000;
					align_tags = realloc(align_tags, sizeof(align_tag) * align_tags_max);
				}
			}else break;
		}

		for (k = s; k > j; k--){
			index = 0;
			lqseq = &lqseqs[k];
			if (lqseq->len >= LQSEQ_MAX_CAN_COUNT) continue;
			for (p = lqseq->start - tags->aln_t_s; p < align_tags_len && align_tags[p].t_pos <= lqseq->end; p++){
				if (align_tags[p].t_pos >= lqseq->start && align_tags[p].q_base != 4){
					lqseq->seqs[lqseq->len].seq[index++] = int_to_base[align_tags[p].q_base];
					if (index > lqseq->seqs[lqseq->len].len) reallocate_lqseq_mem(&lqseq->seqs[lqseq->len]);
				}
			}
			if (index){
				lqseq->seqs[lqseq->len].seq[index] = '\0';
				lqseq->seqs[lqseq->len].len = index;
				lqseq->seqs[lqseq->len].order = lqseq->len;
				if (index > lqseq->lqcount) lqseq->lqcount = index;
				lqseq->len ++;
			}else lqseq->sudoseed_len ++;
		}
	}
	free (align_tags);

	int clusters_i = clusters->i - 1, max_aln_length = 0;
	for (i = 0; i < lqseqs_count; i++){
		lqseq = &lqseqs[i];
		// if (i) printf("%d %d %d %d %d\n",i, lqseq->start, lqseq->end, lqseq->l,lqseqs[i-1].end - lqseq->start);
		// for (k = 0; k < lqseq->len; k++){
		// 	printf("#%d/%d %d s:%d e:%d order:%d lable:%d len:%d %s\n",i,lqseqs_count, k, lqseq->start, 
		// 		lqseq->end, lqseq->seqs[k].order, lqseq->l, strlen(lqseq->seqs[k].seq), lqseq->seqs[k].seq);
		// }
		if (lqseq->l == 1) {
			while (clusters_i >= 0 && !clusters->clusters[clusters_i].i_m) clusters_i--;
			generate_lqseqs_from_cluster(lqseq, &clusters->clusters[clusters_i--]);
			// if (lqseq->len < clusters_i/2){
			// 	lqseq->len = lqseq->lqcount = lqseq->sudoseed_len = 0;
			// 	generate_lqseqs_from_cluster(lqseq, &clusters->clusters[clusters_i]);
			// }
			// clusters_i--;
		}
		if (!lqseq->len) continue;
		int8_t used[LQSEQ_MAX_CAN_COUNT] = {0};
		for (j = s = 0; j < lqseq->len; j ++){
			lqseq->seqs[j].kscore = 1;
			if (used[j]) continue;
			for (k = j + 1; k < lqseq->len; k ++){
				if (strcmp(lqseq->seqs[j].seq, lqseq->seqs[k].seq) == 0){
					used[k] = 1;
					lqseq->seqs[j].kscore ++;
				}
			}
			if (lqseq->seqs[j].kscore > lqseq->seqs[s].kscore || (
				lqseq->seqs[j].kscore == lqseq->seqs[s].kscore && lqseq->seqs[j].len > lqseq->seqs[s].len)) s = j;
		}

		if ((lqseq->seqs[s].kscore > lqseq->len/3 || lqseq->seqs[s].len < 10 || lqseq->len <= 4) && 
			(lqseq->seqs[s].kscore != 1 || (lqseq->len != 3 && lqseq->len != 4))){
			lqseq->len = -2;//rm for debug
			lqseq->l = 4;
			lqseq->sudoseed = strdup(lqseq->seqs[s].seq);
			lqseq->sudoseed_len = lqseq->seqs[s].len;
		}else{
			if (lqseq->len > 4){
				qsort(lqseq->seqs, lqseq->len, sizeof(struct seq_), compare_seq_by_len);
				k = lqseq->len/2;
				while (lqseq->len > k && (lqseq->seqs[lqseq->len - 1].len > 2 * lqseq->seqs[k].len || 
					lqseq->seqs[lqseq->len - 1].len >= 1.4 * lqseq->seqs[lqseq->len - 2].len)) lqseq->len--;
				if (k == lqseq->len){
					lqseq->len = 0;
					continue;
				}

				j = 0;
				k = lqseq->len/2;
				if (lqseq->seqs[j].len < lqseq->seqs[k].len/2){
					reverse_lqseq(lqseq);
					while (lqseq->seqs[lqseq->len - 1].len < lqseq->seqs[k].len/2) lqseq->len --;
					if (k == lqseq->len){
						lqseq->len = 0;
						continue;
					}
				}
			}

			uint16_t kmers[KMER_LEN_COUNT];
			// count_kmers(lqseq, kmers, KMER_MAX_SEQ, 0);
			count_kmers(lqseq, kmers, LQSEQ_MAX_CAN_COUNT, 0);
			count_kscore(lqseq, kmers, 0);
			unsigned int klastscore, kmaxscore;
			unsigned int kmaxlen = lqseq->seqs[0].len;
			if (kmaxlen > 100){
				uint16_t score[LQSEQ_MAX_CAN_COUNT];
				save_kscore(lqseq, score);
				// count_kmers(lqseq, kmers, KMER_MAX_SEQ, 1);
				count_kmers(lqseq, kmers, LQSEQ_MAX_CAN_COUNT, 1);
				count_kscore(lqseq, kmers, 1);
				merge_kscore(lqseq, score);
			}

			qsort(lqseq->seqs, lqseq->len, sizeof(struct seq_), compare_seq_by_kscore);
			kmaxlen = lqseq->seqs[0].len;
			klastscore = kmaxscore = lqseq->seqs[0].kscore;

			for (k = j = 0; j < lqseq->len; j ++){
				if (lqseq->seqs[j].kscore * 10 < kmaxscore || j >= LQSEQ_MAX_COUNT || 
					lqseq->seqs[j].kscore * 2 < klastscore) break;
				klastscore = lqseq->seqs[j].kscore;
				if (j < KMER_MAX_SEQ && lqseq->seqs[j].kscore > kmaxscore * 0.8 && lqseq->seqs[j].len > kmaxlen){
					kmaxlen = lqseq->seqs[j].len;
					k = j;
				}
			}

			lqseq->indexs = 0;
			lqseq->indexe = kmaxlen > LQSEQ_MAX_REV_LEN && j > 6 ? 5 : j - 1;
			if (lqseq->indexe - lqseq->indexs <= 1 || 
					(lqseq->seqs[0].len > 20000 && lqseq->len < LQSEQ_MAX_CAN_COUNT/3)) {
				lqseq->len = 0;
				continue;
			}

			if (lqseq->seqs[0].len < 3000){
				j = lqseq->indexs;
				k = j + 6 < lqseq->indexe ? 6 : lqseq->indexe - j + 1;
			}else{
				j = lqseq->indexs;
				k = j + 2 < lqseq->indexe ? 2 : lqseq->indexe - j + 1;
			}
			if (lqseq->seqs[0].len < 20000) lqseq->sudoseed = poa_to_consensus(&lqseq->seqs[j], k);
			else lqseq->sudoseed = strdup(lqseqs[i].seqs[0].seq);
			lqseq->sudoseed_len = strlen(lqseq->sudoseed);
		}
		if (lqseq->lqcount + lqseq->sudoseed_len > max_aln_length) {
			max_aln_length = lqseq->lqcount + lqseq->sudoseed_len;
		}
		// time_debug(clock(), "lq produced... ");
	}

	// rm lqseq->len = -2; for debug
	// for (i = 0 ; i < lqseqs_count; i++){
	// 	if (lqseqs[i].len){
	// 		for (k = 0; k < lqseqs[i].len; k++){
	// 			j = k >= lqseqs[i].indexs && k <= lqseqs[i].indexe ? 1 : 0;
	// 			printf("%d %d %d %d s:%d e:%d order:%d lable: %d kscore:%d len:%d %s\n",lqseqs_count, i, k, j, lqseqs[i].start, 
	// 				lqseqs[i].end, lqseqs[i].seqs[k].order, lqseqs[i].l, lqseqs[i].seqs[k].kscore, strlen(lqseqs[i].seqs[k].seq), lqseqs[i].seqs[k].seq);
	// 		}
	// 		printf("seed%d: %d %s\n",i, lqseqs[i].sudoseed_len, lqseqs[i].sudoseed);
	// 	}else{
	// 		printf("%d %d len:%d s:%d e:%d\n",lqseqs_count, i, lqseqs[i].len, lqseqs[i].start,lqseqs[i].end);
	// 	}
	// }
	// exit(1);
	for (i = 0; i < aligned_seq_count; i++){
		free (tags_list[i].align_tags);
	}
	free (tags_list);
	return max_aln_length;
}

static int generate_lqseqs_from_tags(lqseq *lqseqs, const int lqseqs_count, align_tags_t *tags_list, \
	const unsigned int aligned_seq_count, const gap_clusters *clusters){

	allocate_lqseq_mem(lqseqs, lqseqs_count);
	align_tags_t *tags;
	align_tag tag = {0, 0, 0};
	int align_tags_len, align_tags_max = 1000000;
	align_tag *align_tags = malloc(sizeof(align_tag) * align_tags_max);
	
	lqseq *lqseq;
	unsigned int i, p;
	int j, k, index, s = lqseqs_count - 1;
	for (i = 1; i < aligned_seq_count; i ++){
		tags = &tags_list[i];
		while (s >= 0 && (lqseqs[s].start < tags->aln_t_s || lqseqs[s].len >= LQSEQ_MAX_CAN_COUNT)) s--;
		for (j = s; j >= 0 && lqseqs[j].end <= tags->aln_t_e; j--);
		if (j == s) continue;

		align_tags_len = p = 0;
		while (get_align_tag(&p, tags, &tag) > 0){
			align_tags[align_tags_len++] = tag;
			if (tag.t_pos <= lqseqs[j + 1].end){
				if (align_tags_len >= align_tags_max){
					align_tags_max += 1000;
					align_tags = realloc(align_tags, sizeof(align_tag) * align_tags_max);
				}
			}else break;
		}

		for (k = s; k > j; k--){
			index = 0;
			lqseq = &lqseqs[k];
			if (lqseq->len >= LQSEQ_MAX_CAN_COUNT) continue;
			for (p = lqseq->start - tags->aln_t_s; p < align_tags_len && align_tags[p].t_pos <= lqseq->end; p++){
				if (align_tags[p].t_pos >= lqseq->start && align_tags[p].q_base != 4){
					lqseq->seqs[lqseq->len].seq[index++] = int_to_base[align_tags[p].q_base];
					if (index > lqseq->seqs[lqseq->len].len) reallocate_lqseq_mem(&lqseq->seqs[lqseq->len]);
				}
			}
			if ((lqseq->l && index) || index > lqseq->end - lqseq->start + 1){
				lqseq->seqs[lqseq->len].seq[index] = '\0';
				lqseq->seqs[lqseq->len].len = index;
				lqseq->seqs[lqseq->len].order = lqseq->len;
				if (index > lqseq->lqcount) lqseq->lqcount = index;
				lqseq->len ++;
			}else lqseq->sudoseed_len ++;
		}
	}
	free (align_tags);

	int clusters_i = clusters->i - 1, max_aln_length = 0;
	for (i = 0; i < lqseqs_count; i++){
		lqseq = &lqseqs[i];
		// if (i) printf("%d %d %d %d %d\n",i, lqseq->start, lqseq->end, lqseq->l,lqseqs[i-1].end - lqseq->start);
		// for (k = 0; k < lqseq->len; k++){
		// 	printf("#%d/%d %d s:%d e:%d order:%d lable:%d len:%d %s\n",i,lqseqs_count, k, lqseq->start, 
		// 		lqseq->end, lqseq->seqs[k].order, lqseq->l, strlen(lqseq->seqs[k].seq), lqseq->seqs[k].seq);
		// }
		if (lqseq->l == 1) {
			while (clusters_i >= 0 && !clusters->clusters[clusters_i].i_m) clusters_i--;
			generate_lqseqs_from_cluster(lqseq, &clusters->clusters[clusters_i--]);
			// if (lqseq->len < clusters_i/2){
			// 	lqseq->len = lqseq->lqcount = lqseq->sudoseed_len = 0;
			// 	generate_lqseqs_from_cluster(lqseq, &clusters->clusters[clusters_i]);
			// }
			// clusters_i--;
		}else if (lqseq->l > 1 && lqseq->len > 4) remove_short_lqseq(lqseq);
		if (lqseq->len <= 4 || lqseq->len < lqseq->sudoseed_len * 0.5) {
			lqseq->len = 0;
		}else{
			qsort(lqseq->seqs, lqseq->len, sizeof(struct seq_), compare_seq_by_len);
			k = lqseq->len/2;
			while (lqseq->len > k && (lqseq->seqs[lqseq->len - 1].len > 2 * lqseq->seqs[k].len || 
				lqseq->seqs[lqseq->len - 1].len >= 1.4 * lqseq->seqs[lqseq->len - 2].len)) lqseq->len--;
			if (k == lqseq->len){
				lqseq->len = 0;
				continue;
			}

			j = 0;
			k = lqseq->len/2;
			if (lqseq->seqs[j].len < lqseq->seqs[k].len/2){
				reverse_lqseq(lqseq);
				while (lqseq->seqs[lqseq->len - 1].len < lqseq->seqs[k].len/2) lqseq->len --;
				if (k == lqseq->len){
					lqseq->len = 0;
					continue;
				}
			}

			uint16_t kmers[KMER_LEN_COUNT];
			// count_kmers(lqseq, kmers, KMER_MAX_SEQ, 0);
			count_kmers(lqseq, kmers, LQSEQ_MAX_CAN_COUNT, 0);
			count_kscore(lqseq, kmers, 0);
			unsigned int klastscore, kmaxscore;
			unsigned int kmaxlen = lqseq->seqs[0].len;
			if (kmaxlen > 100){
				uint16_t score[LQSEQ_MAX_CAN_COUNT];
				save_kscore(lqseq, score);
				// count_kmers(lqseq, kmers, KMER_MAX_SEQ, 1);
				count_kmers(lqseq, kmers, LQSEQ_MAX_CAN_COUNT, 1);
				count_kscore(lqseq, kmers, 1);
				merge_kscore(lqseq, score);
			}

			qsort(lqseq->seqs, lqseq->len, sizeof(struct seq_), compare_seq_by_kscore);
			kmaxlen = lqseq->seqs[0].len;
			klastscore = kmaxscore = lqseq->seqs[0].kscore;

			for (k = j = 0; j < lqseq->len; j ++){
				if (lqseq->seqs[j].kscore * 10 < kmaxscore || j >= LQSEQ_MAX_COUNT || 
					lqseq->seqs[j].kscore * 2 < klastscore) break;
				klastscore = lqseq->seqs[j].kscore;
				if (j < KMER_MAX_SEQ && lqseq->seqs[j].kscore > kmaxscore * 0.8 && lqseq->seqs[j].len > kmaxlen){
					kmaxlen = lqseq->seqs[j].len;
					k = j;
				}
			}

			lqseq->indexs = 0;
			lqseq->indexe = kmaxlen > LQSEQ_MAX_REV_LEN && j > 6 ? 5 : j - 1;
			if (lqseq->indexe - lqseq->indexs <= 3 || 
					(lqseq->seqs[0].len > 20000 && lqseq->len < LQSEQ_MAX_CAN_COUNT/3)) {
				lqseq->len = 0;
				continue;
			}

			if (lqseq->seqs[0].len < 3000){
				j = lqseq->indexs;
				k = j + 6 < lqseq->indexe ? 6 : lqseq->indexe - j + 1;
			}else{
				j = lqseq->indexs;
				k = j + 2 < lqseq->indexe ? 2 : lqseq->indexe - j + 1;
			}
			if (lqseq->seqs[0].len < 20000) lqseq->sudoseed = poa_to_consensus(&lqseq->seqs[j], k);
			else lqseq->sudoseed = strdup(lqseqs[i].seqs[0].seq);
			lqseq->sudoseed_len = strlen(lqseq->sudoseed);

			if (lqseq->lqcount + lqseq->sudoseed_len > max_aln_length) {
				max_aln_length = lqseq->lqcount + lqseq->sudoseed_len;
			}
		}
		// time_debug(clock(), "lq produced... ");
	}
	// for (i = 0 ; i < lqseqs_count; i++){
	// 	if (lqseqs[i].len){
	// 		for (k = 0; k < lqseqs[i].len; k++){
	// 			j = k >= lqseqs[i].indexs && k <= lqseqs[i].indexe ? 1 : 0;
	// 			printf("%d %d %d %d s:%d e:%d order:%d lable: %d kscore:%d len:%d %s\n",lqseqs_count, i, k, j, lqseqs[i].start, 
	// 				lqseqs[i].end, lqseqs[i].seqs[k].order, lqseqs[i].l, lqseqs[i].seqs[k].kscore, strlen(lqseqs[i].seqs[k].seq), lqseqs[i].seqs[k].seq);
	// 		}
	// 		printf("seed%d: %d %s\n",i, lqseqs[i].sudoseed_len, lqseqs[i].sudoseed);
	// 	}else{
	// 		printf("%d %d len:%d s:%d e:%d\n",lqseqs_count, i, lqseqs[i].len, lqseqs[i].start,lqseqs[i].end);
	// 	}
	// }
	// exit(1);
	for (i = 0; i < aligned_seq_count; i++){
		free (tags_list[i].align_tags);
	}
	free (tags_list);
	return max_aln_length;
}

static consensus_trimed *get_lqseqs_from_align_tags(align_tags_t *tags_list, msa_p *msa,\
		int aligned_seq_count, int len){
	int p, d, b, m, n;
	int64_t p_pp_score, p_pp_score_, pp_ppp_score = -10;
	// int64_t p_base_score = -10;
	align_tag global_best_p = {0, 0, -1};
	
	// msa_p_d_b_pp_ppps pp_ppps;
	// pp_ppps.len = pp_ppps.max_size = 0;
	allocate_msa_mem(msa, len);
	update_msa(msa, tags_list, aligned_seq_count, 1);
	msa_p_d_b_pp_ppp *pp_ppp_m, *pp_ppp_n;

	if (READS_TYPE == READS_HIFI){
		for (p = 0; p < len; p++){
			for ( d = 0; d < msa[p].max_size; d++ ){
				for (b = 0; b < 6; b++){
					msa_p_d_b *p_base = &msa[p].d_b[d][b];
					p_base->max_size = 0;
					p_pp_score_ = p_pp_score = INT64_MIN;
					for (m = 0; m < p_base->len; m++){
						pp_ppp_m = &p_base->pp_ppp[m];
						if (pp_ppp_m->pp.t_pos == -1){
							pp_ppp_m->score = 10 * pp_ppp_m->link_count - 4 * msa[p].coverage;
						}else{

							msa_p_d_b *pp_base = &msa[pp_ppp_m->pp.t_pos].d_b[pp_ppp_m->pp.delta][pp_ppp_m->pp.q_base];

							for (n = 0; n < pp_base->len; n++){
								pp_ppp_n = &pp_base->pp_ppp[n];
								if (pp_ppp_n->pp.t_pos == pp_ppp_m->ppp.t_pos &&
									pp_ppp_n->pp.delta == pp_ppp_m->ppp.delta &&
									pp_ppp_n->pp.q_base == pp_ppp_m->ppp.q_base){
									pp_ppp_score = pp_ppp_n->score + 10 * pp_ppp_m->link_count - 4 * msa[p].coverage;
									if (pp_ppp_score > pp_ppp_m->score) {
										pp_ppp_m->score = pp_ppp_score;
										p_pp_score_ = pp_ppp_n->score;
									}
									if (pp_ppp_n->score > p_pp_score || 
										(pp_ppp_n->score == p_pp_score && pp_ppp_m->pp.q_base != 4)){
										p_base->max_size = m;
										p_pp_score = pp_ppp_n->score;
									}
								}
							}
						}
						if (pp_ppp_m->score > p_base->pp_ppp[p_base->max_size].score || \
							(pp_ppp_m->score == p_base->pp_ppp[p_base->max_size].score && pp_ppp_m->pp.q_base != 4)){
							p_base->max_size = m;
							p_pp_score = p_pp_score_;
						}
					}
					global_best_p.t_pos = p;
					global_best_p.delta = d;
					global_best_p.q_base = b;
				}
			}
		}
	}else{
		for (p = 0; p < len; p++){
			for ( d = 0; d < msa[p].max_size; d++ ){
				for (b = 0; b < 6; b++){
					msa_p_d_b *p_base = &msa[p].d_b[d][b];
					p_base->max_size = 0;
					p_pp_score_ = p_pp_score = INT64_MIN;
					// p_base_score = INT64_MIN;
					// p_base->best_score = -10;
					// p_base->best_pp.t_pos = -1;
					for (m = 0; m < p_base->len; m++){
						pp_ppp_m = &p_base->pp_ppp[m];
						if (pp_ppp_m->pp.t_pos == -1){
							pp_ppp_m->score = 10 * pp_ppp_m->link_count - 2 * msa[p].coverage;
						}else{

							msa_p_d_b *pp_base = &msa[pp_ppp_m->pp.t_pos].d_b[pp_ppp_m->pp.delta][pp_ppp_m->pp.q_base];

							for (n = 0; n < pp_base->len; n++){
								pp_ppp_n = &pp_base->pp_ppp[n];
								if (pp_ppp_n->pp.t_pos == pp_ppp_m->ppp.t_pos &&
									pp_ppp_n->pp.delta == pp_ppp_m->ppp.delta &&
									pp_ppp_n->pp.q_base == pp_ppp_m->ppp.q_base){
									pp_ppp_score = pp_ppp_n->score + 10 * pp_ppp_m->link_count - 2 * msa[p].coverage;
									if (pp_ppp_score > pp_ppp_m->score) {
										pp_ppp_m->score = pp_ppp_score;
										p_pp_score_ = pp_ppp_n->score;
									}
									if (pp_ppp_m->link_count > p_base->pp_ppp[p_base->max_size].link_count/2 &&
										pp_ppp_n->score > p_pp_score && (pp_ppp_m->pp.q_base == 4 || pp_ppp_m->pp.q_base == b
										|| pp_ppp_m->ppp.q_base == b || pp_ppp_m->pp.q_base == pp_ppp_m->ppp.q_base)){
										p_base->max_size = m;
										p_pp_score = pp_ppp_n->score;
									}
									// if (pp_ppp_n->score > p_pp_score || //test2.log6
									// 	(pp_ppp_n->score == p_pp_score && pp_ppp_m->pp.q_base != 4)){
									// 	p_base->max_size = m;
									// 	p_pp_score = pp_ppp_n->score;
									// }
									// if (pp_ppp_n->score > p_pp_score){//test2.log7
									// 	p_base->max_size = m;
									// 	p_pp_score = pp_ppp_n->score;
									// }
								}
							}
						}
						if (pp_ppp_m->score > p_base->pp_ppp[p_base->max_size].score || \
							(pp_ppp_m->score == p_base->pp_ppp[p_base->max_size].score && pp_ppp_m->pp.q_base != 4)){
							p_base->max_size = m;
							p_pp_score = p_pp_score_;
						// if (pp_ppp_m->score > p_base->best_score || (pp_ppp_m->score == p_base->best_score && pp_ppp_m->pp.q_base != 4)){// find the bese score for a given p	
							// p_base->best_score = pp_ppp_m->score;
							// p_base->best_pp = pp_ppp_m->pp;
							// p_base->best_link_count = pp_ppp_m->link_count;
						}
	// printf("%d %d %c depth:%d %d %d %c %d %d %c %ld %d %ld\n",p,d,int_to_base[b],msa[p].coverage, 
	// 	pp_ppp_m->pp.t_pos,pp_ppp_m->pp.delta,int_to_base[pp_ppp_m->pp.q_base],pp_ppp_m->ppp.t_pos,pp_ppp_m->ppp.delta,
	// 	int_to_base[pp_ppp_m->ppp.q_base],pp_ppp_m->score, pp_ppp_m->link_count, p_base->pp_ppp[p_base->max_size].score);
					}
					global_best_p.t_pos = p;
					global_best_p.delta = d;
					global_best_p.q_base = b;
				}
			}
		}
	}
	// free (pp_ppps.pp_ppp);

	consensus_trimed *cons_trimed;
	cons_trimed = malloc(sizeof(consensus_trimed));
	int cons_trimed_seq_max_len = len * 2 + 1;
	cons_trimed->len = 0;
	cons_trimed->seq = malloc(cons_trimed_seq_max_len * sizeof(char));
	
	msa_p_d_b *global_best_d_b;
	uint16_t global_best_quality;
	while (1){
		if (global_best_p.q_base != 4){
			if (cons_trimed->len >= cons_trimed_seq_max_len){
				cons_trimed_seq_max_len += 500;
				cons_trimed->seq = realloc(cons_trimed->seq, cons_trimed_seq_max_len * sizeof(char));
			}
			global_best_d_b = &msa[global_best_p.t_pos].d_b[global_best_p.delta][global_best_p.q_base];
			
			// printf("#%d %d %d %c global_best_d_b->max_size:%d linkcount:%d cov:%d\n", cons_trimed->len,global_best_p.t_pos,
			// 	global_best_p.delta,int_to_base[global_best_p.q_base],global_best_d_b->max_size,
			// 	global_best_d_b->pp_ppp[global_best_d_b->max_size].link_count,
			// 	msa[global_best_p.t_pos].coverage);

			global_best_quality = global_best_d_b->pp_ppp[global_best_d_b->max_size].link_count;
			cons_trimed->seq[cons_trimed->len ++] = 
		 		// msa[global_best_p.t_pos].d_b[global_best_p.delta][global_best_p.q_base].best_link_count * 5
		 		global_best_quality * 5
		 		> msa[global_best_p.t_pos].coverage || int_to_base[global_best_p.q_base] == 'N' ? \
		 		int_to_base[global_best_p.q_base] : tolower(int_to_base[global_best_p.q_base]);
		}

		global_best_d_b = &msa[global_best_p.t_pos].d_b[global_best_p.delta][global_best_p.q_base];
		global_best_p = global_best_d_b->pp_ppp[global_best_d_b->max_size].pp;


		// global_best_p = msa[global_best_p.t_pos].d_b[global_best_p.delta][global_best_p.q_base].best_pp;
		// printf("%d %d %c\n", global_best_p.t_pos, global_best_p.delta, int_to_base[global_best_p.q_base]);
		if (global_best_p.t_pos == -1){
			break;
		}
	}
	free_msa(msa, len);
	cons_trimed->seq[cons_trimed->len] = '\0';
	// printf("%s\n", cons_trimed->seq);

	// int t = -1;
	// for (int i = 0; i < cons_trimed->len; i++){
	// 	if (cons_trimed->seq[i] == 'N') {t++; printf("\n%d ", t);continue;}
	// 	printf("%c", cons_trimed->seq[i]);
	// }
	// printf("\n");

	return cons_trimed;
}

static consensus_data *update_consensus_trimed(lqseq *lqseqs, const int lqseqs_count, \
		consensus_data *consensus, const int seed_len){
	
	consensus_data *consensus_ = calloc(1, sizeof(consensus_data));
	consensus_->max_size = seed_len/5 * 6 + 1;
	consensus_->cns_bases = malloc(consensus_->max_size * sizeof(consensus_base));
	unsigned int p = 0;
	unsigned int i = consensus->lstrip, j;
	int update, lqseqs_index;
	update = 1;
	lqseqs_index = lqseqs_count - 1;
	
	while (i < consensus->len - consensus->rstrip){
		p = consensus->cns_bases[i].pos;
		
		if (lqseqs_index >= 0 && ((lqseqs[lqseqs_index].len <= 0 && lqseqs[lqseqs_index].len != -2) || p > lqseqs[lqseqs_index].end)) {
			lqseqs_index--;
			update = 1;
		}
		if (lqseqs_index >= 0 && (lqseqs[lqseqs_index].len > 0 || lqseqs[lqseqs_index].len == -2) && p >= lqseqs[lqseqs_index].start && 
			p <= lqseqs[lqseqs_index].end){
			if (update){
				// printf("\n%d %d ", lqseqs_index, lqseqs[lqseqs_index].sudoseed_len);
				for (j = 0; j < lqseqs[lqseqs_index].sudoseed_len; j++){
					consensus_->cns_bases[consensus_->len].base = lqseqs[lqseqs_index].sudoseed[j];
					consensus_->cns_bases[consensus_->len++].pos = lqseqs[lqseqs_index].start;
					// printf("%c", lqseqs[lqseqs_index].sudoseed[j]);
					if (consensus_->len >= consensus_->max_size){
						consensus_->max_size += 100000;
						consensus_->cns_bases = realloc(consensus_->cns_bases, consensus_->max_size * sizeof(consensus_base));
					}
				}
				update = 0;
			}
		}else{
			consensus_->cns_bases[consensus_->len++] = consensus->cns_bases[i];
			update = 1;
			if (consensus_->len >= consensus_->max_size){
				consensus_->max_size += 100000;
				consensus_->cns_bases = realloc(consensus_->cns_bases, consensus_->max_size * sizeof(consensus_base));
			}
		}
		i ++;
	}
	
	return consensus_;
}

static inline void get_align_tags(alignment *aln, align_tags_t *tags, msa_p *msa){
	//A--TC--GATCATGC--
	//ATT-CGG-ATCATGCAT
	
	int p, i, len = (aln->aln_len + 1) >> 1;
	tags->aln_t_s = tags->aln_t_e = aln->aln_t_s;
	tags->align_tags = calloc(len + 1, sizeof(uint8_t));

	uint8_t b, l = 0;
	uint16_t delta = 0;
	tags->aln_t_e --;
	for (p = 0, i = aln->shift; i < aln->aln_len + aln->shift; i++, p++){
		// printf("%d %c\n", i, aln->q_aln_str[i]);
		b = base_to_int[(unsigned char)aln->q_aln_str[i]];
		if (aln->t_aln_str[i] == '-') {
			b |= 8;
			delta ++;
		}else {
			tags->aln_t_e ++;
			l = delta = 0;
		}
		if (!(p & 1)) b <<= 4;
		tags->align_tags[p>>1] |= b;
		if (delta == 0 && aln->q_aln_str[i] != 'M'){
			msa[tags->aln_t_e].coverage ++;
		}

		if (delta >= msa[tags->aln_t_e].max_size){
			msa[tags->aln_t_e].max_size = delta + 1;
		}

		if (delta >= GAP_MIN_LEN && !l) {
			msa[tags->aln_t_e].l_ins ++;
			l = 1;
		}

		if (delta == 0 && aln->q_aln_str[i] == '-'){
			msa[tags->aln_t_e].l_del ++;
		}
	}
	
	if ((p - 1) & 1) tags->align_tags[p>>1] |= 255;
	else tags->align_tags[p>>1] |= 15;
}

static void fill_aln_with_seed(alignment *aln, int seed_len){
	char _[seed_len + 1];
	memset(_, 'M', seed_len * sizeof(char));
	_[seed_len] = '\0';
	
	strcpy(&aln->t_aln_str[aln->aln_len], _);
	strcpy(&aln->q_aln_str[aln->aln_len], _);
	aln->aln_len += seed_len;
}

static void fill_aln_with_lqseq(alignment *aln, char *seed, int seed_len, char *lqseq, int lqseq_len){
	if (lqseq_len > seed_len){
		char _[lqseq_len + 1];
		strcpy(_, seed);
		while (seed_len < lqseq_len) _[seed_len++] = '-';
		_[seed_len] = '\0';
		strcpy(&aln->t_aln_str[aln->aln_len], _);
		strcpy(&aln->q_aln_str[aln->aln_len], lqseq);
	}else{
		char _[seed_len + 1];
		strcpy(_, lqseq);
		while (lqseq_len < seed_len) _[lqseq_len++] = '-';
		_[lqseq_len] = '\0';
		strcpy(&aln->t_aln_str[aln->aln_len], seed);
		strcpy(&aln->q_aln_str[aln->aln_len], _);
	}
	aln->aln_len += lqseq_len;
}

static consensus_trimed *generate_consensus_trimed(lqseq *lqseqs, const int lqseqs_count,\
		const int max_aln_length){

	int seq_count = LQSEQ_MAX_COUNT;
	alignment aln, aln_link;
	align_tags_t *tags_list;
	int seed_len, query_len;
	tags_list = malloc(seq_count * sizeof(align_tags_t));

	msa_p * msa = NULL;

	int max_mem_len = max_aln_length + 2;
	uint64_t max_mem_d = (int) max_mem_len * 0.4 + 1;
	aln.t_aln_str = malloc(max_mem_len * sizeof(char));
	aln.q_aln_str = malloc(max_mem_len * sizeof(char));
	int aligned_linkseq_max = 65536 * 2 + 1;
	int aligned_linkseq_len = 0;
	int delta;
	aln_link.shift = 0;
	aln_link.t_aln_str = malloc(aligned_linkseq_max * sizeof(char));
	aln_link.q_aln_str = malloc(aligned_linkseq_max * sizeof(char));

	int *V = malloc( max_mem_d * 2 * sizeof(int));
	uint8_t **D = malloc( max_mem_d * sizeof(uint8_t *) + max_mem_d * (max_mem_d + 1)/2 * sizeof(uint8_t));
	assert (V && D);
	uint8_t * const _ = (uint8_t *) (D + max_mem_d);
	for (uint64_t d = 0; d < max_mem_d; d ++ ) {
		D[d] = _ + d * (d + 1)/2;
	}

	int i, j;
	for (i = 0; i < lqseqs_count; i++) lqseqs[i].lqcount = 0;

	for (i = 0; i < seq_count; i++){
		aligned_linkseq_len = 0;
		aln_link.aln_t_s = 0;
		aln_link.aln_t_e = 0;
		aln_link.aln_len = 0;
		for (j = lqseqs_count - 1; j >= 0; j--){
			if (lqseqs[j].len <= 0) continue;

			seed_len = lqseqs[j].sudoseed_len;
			aligned_linkseq_len += seed_len + 1;
			aln_link.aln_t_e += seed_len + 1;
			aln_link.t_aln_str[aln_link.aln_len] = 'N';
			aln_link.q_aln_str[aln_link.aln_len] = 'N';
			aln_link.aln_len ++;

			query_len = (i + lqseqs[j].indexs) > lqseqs[j].indexe ? seed_len : lqseqs[j].seqs[i + lqseqs[j].indexs].len;
			int mlen = query_len > seed_len ? query_len : seed_len;
			if (aln_link.aln_len + mlen * 2 + 1 >= aligned_linkseq_max){
				aligned_linkseq_max += query_len + mlen * 2 + 1;
				aln_link.t_aln_str = realloc(aln_link.t_aln_str, aligned_linkseq_max * sizeof(char));
				aln_link.q_aln_str = realloc(aln_link.q_aln_str, aligned_linkseq_max * sizeof(char));
			}

			if (i + lqseqs[j].indexs > lqseqs[j].indexe) lqseqs[j].lqcount = 0;
			if (i + lqseqs[j].indexs > lqseqs[j].indexe || (i && (query_len < seed_len * 0.5 || query_len > seed_len * 1.3))){
				if (lqseqs[j].lqcount ++ < lqseqs[j].indexe - lqseqs[j].indexs) fill_aln_with_seed(&aln_link, seed_len);
				else fill_aln_with_lqseq(&aln_link, lqseqs[j].sudoseed, seed_len, lqseqs[j].seqs[lqseqs[j].indexs].seq, 
						lqseqs[j].seqs[lqseqs[j].indexs].len);
			}else{
				aln.aln_t_s = 0;
				aln.aln_t_e = seed_len;
				aln.aln_len = 0;
				aln.shift = 0;

				memset(V, 0, max_mem_d * 2 * sizeof(int));
				align(lqseqs[j].seqs[i + lqseqs[j].indexs].seq, query_len,
							lqseqs[j].sudoseed, seed_len, &aln, V, D);
				// align_nd(lqseqs[j].seqs[i + lqseqs[j].indexs].seq, query_len,
				// lqseqs[j].sudoseed, seed_len, &aln);
				// aln.aln_t_len = seed_len;
				// aln.aln_q_len = lqseqs[j].seqs[i + lqseqs[j].indexs].len;
				// printf("%d %d %d %d\n", i, j, query_len, seed_len);
				// printf("#%d %d %d %s %s\n", i, j, query_len, lqseqs[j].seqs[i + lqseqs[j].indexs].seq, lqseqs[j].sudoseed);
				// out_align(&aln, 0);

				if (aln.aln_len > 2){

					strcpy(&aln_link.t_aln_str[aln_link.aln_len], aln.t_aln_str);
					strcpy(&aln_link.q_aln_str[aln_link.aln_len], aln.q_aln_str);
					aln_link.aln_len += aln.aln_len;

					while (aln.aln_t_len < seed_len){
						aln_link.t_aln_str[aln_link.aln_len] = lqseqs[j].sudoseed[aln.aln_t_len++];
						aln_link.q_aln_str[aln_link.aln_len] = '-';
						aln_link.aln_len ++;
					}
					delta = 0;
					while (aln.aln_q_len < lqseqs[j].seqs[i + lqseqs[j].indexs].len && delta++ < 250){
						aln_link.q_aln_str[aln_link.aln_len] = lqseqs[j].seqs[i + lqseqs[j].indexs].seq[aln.aln_q_len++];
						aln_link.t_aln_str[aln_link.aln_len] = '-';
						aln_link.aln_len ++;
					}

					assert (aln_link.aln_len <= aligned_linkseq_max);
				}else{
					if (lqseqs[j].lqcount ++ < lqseqs[j].indexe - lqseqs[j].indexs) fill_aln_with_seed(&aln_link, seed_len);
					else fill_aln_with_lqseq(&aln_link, lqseqs[j].sudoseed, seed_len, lqseqs[j].seqs[lqseqs[j].indexs].seq, 
							lqseqs[j].seqs[lqseqs[j].indexs].len);
				}
			}
		}

		aligned_linkseq_len ++;
		aln_link.aln_t_e ++;
		aln_link.t_aln_str[aln_link.aln_len] = 'N';
		aln_link.q_aln_str[aln_link.aln_len] = 'N';
		aln_link.aln_len ++;
		aln_link.t_aln_str[aln_link.aln_len] = '\0';
		aln_link.q_aln_str[aln_link.aln_len] = '\0';

		if (! msa) msa = calloc(aligned_linkseq_len, sizeof(msa_p));
		// out_align(&aln_link, 0);
		get_align_tags(&aln_link, &tags_list[i], msa);

	}
	// exit(1);
	free(aln.t_aln_str);
	free(aln.q_aln_str);
	free(aln_link.t_aln_str);
	free(aln_link.q_aln_str);
	free(V);
	free(D);
	consensus_trimed *cons_trimed = get_lqseqs_from_align_tags( tags_list, msa, seq_count, aligned_linkseq_len);
	return cons_trimed;
}

static uint32_t get_lqseq_minlen(lqseq *lqseq){
	int i;
	uint32_t t = lqseq->seqs[0].len;
	for ( i = 1; i < lqseq->len; i ++){
		if (lqseq->seqs[i].len < t) t  = lqseq->seqs[i].len;
	}
	return t;
}

static consensus_data *iterate_generate_consensus_trimed(lqseq *lqseqs, const int lqseqs_count,\
	consensus_data *consensus, const int len, int max_aln_length, const int iterate){
	
	int i, j, psed_len, max_psed_len, max_sed_len, max_dif_len;
	unsigned int k;
	psed_len = max_psed_len = max_sed_len = max_dif_len = 0;
	consensus_trimed *cons_trimed;
	for (i = 1; i <= iterate; i++){
		max_aln_length += max_dif_len;
		j = lqseqs_count;
		psed_len = max_sed_len = max_dif_len = 0;
		cons_trimed = generate_consensus_trimed(lqseqs, lqseqs_count, max_aln_length);
		// printf("lqseqs_count %d\n%s\n",lqseqs_count, cons_trimed->seq);
		for (k = cons_trimed->len; k; k --){
			if (cons_trimed->seq[k - 1] != 'N'){
				if (cons_trimed->seq[k - 1] < 'a'){
					// printf("%d %d\n",j,k );
					lqseqs[j].sudoseed[lqseqs[j].sudoseed_len++] = cons_trimed->seq[k - 1];
				}else{
					lqseqs[j].sudoseed[lqseqs[j].sudoseed_len++] = toupper(cons_trimed->seq[k - 1]);
					lqseqs[j].lqcount ++;
				}
				if (lqseqs[j].sudoseed_len >= max_psed_len) {
					max_psed_len += 1000;
					lqseqs[j].sudoseed = realloc(lqseqs[j].sudoseed, max_psed_len);
				}
			}else{
				if (j != lqseqs_count) {
					lqseqs[j].sudoseed[lqseqs[j].sudoseed_len] = '\0';
					if (lqseqs[j].sudoseed_len > psed_len + max_dif_len) max_dif_len = lqseqs[j].sudoseed_len - psed_len;
					if (((lqseqs[j].sudoseed_len <= lqseqs[j].end - lqseqs[j].start + 1 || 
						lqseqs[j].lqcount > lqseqs[j].sudoseed_len * 4/5) && !lqseqs[j].l)||(
						lqseqs[j].l && lqseqs[j].sudoseed_len * 1.3 < get_lqseq_minlen(&lqseqs[j]))) lqseqs[j].len = -1;
					// printf("%d #%d %d %d %d %d %s\n",j, k, i, lqseqs[j].len, lqseqs[j].lqcount, 
					// 	lqseqs[j].sudoseed_len, lqseqs[j].sudoseed);
				}
				j --;
				while (j >= 0 && lqseqs[j].len <= 0) j--;
				if (j < 0) continue;
				max_psed_len = psed_len = lqseqs[j].sudoseed_len;
				lqseqs[j].sudoseed_len = lqseqs[j].lqcount = 0;
				// if (sed_len > max_psed_len) max_psed_len = sed_len;
			}
		}
		free_consensus_trimed(cons_trimed);
	}
	consensus_data *consensus_= update_consensus_trimed(lqseqs, lqseqs_count, consensus, len);
	return consensus_;
}

static consensus_data *generate_cns_from_best_score_fast( msa_p *msa, align_tag *global_best_p,\
		int aligned_seq_count, int len, int min_cov, float min_error_corrected_ratio){
	consensus_data *consensus = calloc(1, sizeof(consensus_data));
	consensus->max_size = len/5 * 6 + 1;
	consensus->cns_bases = malloc(consensus->max_size * sizeof(consensus_base));

	msa_p_d_b *global_best_d_b;
	while (1){
		if (global_best_p->q_base != 4){

// global_best_d_b = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base];
// int global_best_quality = global_best_d_b->pp_ppp[global_best_d_b->max_size].link_count;
// int pqv = 100 * global_best_quality/msa[global_best_p->t_pos].coverage;
// printf("%d %d %c %d %d %d\n",global_best_p->t_pos,  consensus->len, 
// int_to_base[global_best_p->q_base], pqv, msa[global_best_p->t_pos].l_ins,
// msa[global_best_p->t_pos].coverage);

			consensus->cns_bases[consensus->len].base = msa[global_best_p->t_pos].coverage > min_cov ? 
				int_to_base[global_best_p->q_base] : tolower(int_to_base[global_best_p->q_base]);
			consensus->cns_bases[consensus->len++].pos = global_best_p->t_pos;
			if (consensus->len >= consensus->max_size){
				consensus->max_size += 100000;
				consensus->cns_bases = realloc(consensus->cns_bases, consensus->max_size * sizeof(consensus_base));
			}
		}
		global_best_d_b = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base];
		global_best_p = &global_best_d_b->pp_ppp[global_best_d_b->max_size].pp;
		if (global_best_p->t_pos == -1){
			break;
		}
	}
	free_msa(msa, len);
	// reverse_str(cons_trimed->seq, cons_trimed->len);
	// cons_trimed->seq[cons_trimed->len] = '\0';
	return consensus;
}

typedef struct {
	pos gap;
	int l;
} del;

static int get_lqseqs_from_cluster(gap_cluster *clu, lqseq *lqseqs, int index){
	if (clu->i_m){
		if (index >= 0){
			while (index > 0 && lqseqs[index].start <= clu->r.e) index --;
			// printf("%d %u\n",index,  lqseqs[index].start);
			if (lqseqs[index].start > clu->r.e) index ++;
			lqseqs[index].start = clu->r.s;
			lqseqs[index].end = clu->r.e;
			lqseqs[index].l = 1;
			// printf("gap:%u %u %u\n",index,lqseqs[index].start,lqseqs[index].end  );
		}else{
			lqseqs[++index].start = clu->r.s;
			lqseqs[index].end = clu->r.e;
			lqseqs[index].l = 1;
		}
	}
	return index;
}

static int get_lqseqs_from_dels(del *del, lqseq *lqseqs, int index){
	if (index >= 0){
		int s;
		s = del->gap.s < lqseqs[index].start ? del->gap.s : lqseqs[index].start;
		while (index > 0 && lqseqs[index].start <= del->gap.e && !lqseqs[index].l) index --;
		if (lqseqs[index].start > del->gap.e) lqseqs[++index].end = 0;
		else if (lqseqs[index].l) return index;
		lqseqs[index].start = s;
		lqseqs[index].end = del->gap.e > lqseqs[index].end ? del->gap.e : lqseqs[index].end;
		lqseqs[index].l = del->l;
	}else{
		lqseqs[++index].start = del->gap.s;
		lqseqs[index].end = del->gap.e;
		lqseqs[index].l = del->l;
	}
	return index;
}

static int cal_del_pos(msa_p *msa, int s, int e){
	int i, validy;
	for (validy = 0, i = s; i <= e; i++){
		if (msa[i].l_del > msa[i].coverage * 0.6) validy++;
	}
	return validy;
}

static int get_l_del_regions(msa_p *msa, consensus_data *consensus, del **dels_){
	int i, s, e, p, ps, pe, l;
	int index = ps = pe = 0;
	int index_m = 500;
	del *dels = malloc(sizeof(del) * index_m);
	for (i = 1; i < consensus->len; i++){
		if (msa[consensus->cns_bases[i].pos].l_del < 
			msa[consensus->cns_bases[i].pos].coverage * DEL_MIN_DEPTH_RATIO && 
			consensus->cns_bases[i].pos < consensus->cns_bases[i - 1].pos + DEL_MIN_LEN) continue;
		if (i >= ps && i <= pe) continue;
		s = i - 1;
		while (s > 0 && msa[consensus->cns_bases[s].pos].l_del > 
			msa[consensus->cns_bases[s].pos].coverage * DEL_MIN_DEPTH_RATIO) s--;
		e = i + 1;
		while (e < consensus->len - 1 && msa[consensus->cns_bases[e].pos].l_del >
			msa[consensus->cns_bases[e].pos].coverage * DEL_MIN_DEPTH_RATIO) e++;
		if (consensus->cns_bases[e].pos - consensus->cns_bases[s].pos < 10) continue;
		// printf("%d %d\n",s, e );
		p = cal_del_pos(msa, consensus->cns_bases[s].pos, consensus->cns_bases[e].pos);
		l = consensus->cns_bases[e].pos - consensus->cns_bases[s].pos + 1;
		if ((READS_TYPE == READS_CLR || READS_TYPE == READS_RS) && p < l * 0.05) continue;

		l = p > l / 3 ? 2 : 3;
		ps = s;
		pe = e;

		for (p = 0, s = i - LQSEQ_MIN_LEN/2; s > 0; s --){
			if (consensus->cns_bases[s].qv >= HQ_MIN_QV && msa[consensus->cns_bases[s].pos].l_del < 
				msa[consensus->cns_bases[s].pos].coverage * DEL_MIN_DEPTH_RATIO) p ++;
			else p = 0;
			if (p >= HQSEQ_MIN_LEN && 
				base_to_int[(unsigned char) consensus->cns_bases[s].base] != 
					base_to_int[(unsigned char) consensus->cns_bases[s - 1].base] &&
				msa[consensus->cns_bases[s].pos].l_ins <= 0) break;
		}

		for (p = 0, e = i + LQSEQ_MIN_LEN/2; e < consensus->len - 1; e ++){
			if (consensus->cns_bases[e].qv >= HQ_MIN_QV && msa[consensus->cns_bases[e].pos].l_del < 
				msa[consensus->cns_bases[e].pos].coverage * DEL_MIN_DEPTH_RATIO) p ++;
			else p = 0;
			if (p >= HQSEQ_MIN_LEN && 
				base_to_int[(unsigned char) consensus->cns_bases[e].base] != 
					base_to_int[(unsigned char) consensus->cns_bases[e + 1].base] &&
				msa[consensus->cns_bases[e].pos].l_ins <= 0) break;
		}

		s = s >= 0 ? consensus->cns_bases[s].pos : consensus->cns_bases[0].pos;
		e = e < consensus->len - 1 ? consensus->cns_bases[e].pos : consensus->cns_bases[consensus->len - 1].pos;
		if (e - s < DEL_MIN_LEN) continue;  //TODO CHECK
		if (index == 0 || s > dels[index - 1].gap.e){
			dels[index].gap.s = s;
			dels[index].gap.e = e;
			dels[index].l = l;//TODO CHECK
			if (++index >= index_m) {
				index_m += 100;
				dels = realloc(dels, index_m * sizeof(del));
			}
		}else dels[index - 1].gap.e = e;

	}
	*dels_ = dels;
	// for (i = 0; i < index; i ++){
	// 	printf("1:%d %d %d %d\n",i, dels[i].gap.s,dels[i].gap.e,dels[i].l);
	// }
	// exit(1);
	return index;
}

static lqseq *get_lqseqs_from_gap(msa_p *msa, consensus_data *consensus, int *lqseqs_count, 
		gap_clusters *clusters){//TODO OPTIMIZE SPEED

	int i, s, e, p, l_ins, index, lqseqs_max_size = 200;
	lqseq *lqseqs = malloc(lqseqs_max_size * sizeof(lqseq));
	s = e = index = lqseqs[0].start = lqseqs[0].end = 0;
	int clusters_i = clusters->i;
	// int index_m = 500;
	del *dels = NULL;
	int dels_i = get_l_del_regions(msa, consensus, &dels);

	for (i = consensus->len - 1; i >= 0; i--){
		if (msa[consensus->cns_bases[i].pos].l_ins < 
			msa[consensus->cns_bases[i].pos].coverage * GAP_MIN_RATIO1) continue;
		if (msa[consensus->cns_bases[i].pos].l_ins < 
				msa[consensus->cns_bases[i].pos].coverage * GAP_MIN_RATIO2){
			s = (int) consensus->cns_bases[i].pos - GAP_FLANK_LEN;
			e = consensus->cns_bases[i].pos + GAP_FLANK_LEN;
			l_ins = msa[consensus->cns_bases[i].pos].l_ins;

			for (p = i - 1; p >= 0 && consensus->cns_bases[p].pos >= s; p --){
				if (consensus->cns_bases[p].pos != consensus->cns_bases[p + 1].pos) {
					l_ins += msa[consensus->cns_bases[p].pos].l_ins;
				}
			}
			for (p = i + 1; p < consensus->len && consensus->cns_bases[p].pos <= e; p++){
				if (consensus->cns_bases[p].pos != consensus->cns_bases[p - 1].pos){
					l_ins += msa[consensus->cns_bases[p].pos].l_ins;
				}
			}
			// printf("%d %d %c %d %d %d %d\n",consensus->cns_bases[i].pos, i, 
			// consensus->cns_bases[i].base, consensus->cns_bases[i].qv, msa[consensus->cns_bases[i].pos].l_ins,
			// msa[consensus->cns_bases[i].pos].coverage, l_ins);
			if (l_ins < msa[consensus->cns_bases[i].pos].coverage * GAP_MIN_RATIO3) continue;
		}
		
		// if (reads_type == 2 && msa[consensus->cns_bases[i].pos].l_ins < 
		// 		msa[consensus->cns_bases[i].pos].coverage * 0.6) continue;
		// printf("l_ins:%d %d %c %d %d %d %f\n",consensus->cns_bases[i].pos, i, 
		// 	consensus->cns_bases[i].base, consensus->cns_bases[i].qv, msa[consensus->cns_bases[i].pos].l_ins,
		// 	msa[consensus->cns_bases[i].pos].coverage, (float)msa[consensus->cns_bases[i].pos].l_ins/ msa[consensus->cns_bases[i].pos].coverage);

		for (p = 0, s = i - LQSEQ_MIN_LEN/2; s > 0; s --){
			if (consensus->cns_bases[s].qv >= HQ_MIN_QV) p ++;
			else p = 0;
			if (p >= HQSEQ_MIN_LEN && 
				base_to_int[(unsigned char) consensus->cns_bases[s].base] != 
					base_to_int[(unsigned char) consensus->cns_bases[s - 1].base] &&
				msa[consensus->cns_bases[s].pos].l_ins <= 0) break;
		}

		for (p = 0, e = i + LQSEQ_MIN_LEN/2; e < consensus->len - 1; e ++){
			if (consensus->cns_bases[e].qv >= HQ_MIN_QV) p ++;
			else p = 0;
			if (p >= HQSEQ_MIN_LEN && 
				base_to_int[(unsigned char) consensus->cns_bases[e].base] != 
					base_to_int[(unsigned char) consensus->cns_bases[e + 1].base] &&
				msa[consensus->cns_bases[e].pos].l_ins <= 0) break;
		}

		s = s >= 0 ? consensus->cns_bases[s].pos : consensus->cns_bases[0].pos;
		e = e < consensus->len - 1 ? consensus->cns_bases[e].pos : consensus->cns_bases[consensus->len - 1].pos;
		// printf("%u %u\n",s,e );
		if (index == 0 || e + GAP_BETWEEN_LEN < lqseqs[index - 1].start){
			while (dels_i && e < dels[dels_i - 1].gap.s){
				index = get_lqseqs_from_dels(&dels[dels_i - 1], lqseqs, index - 1);
				dels_i --;
				if (++index >= lqseqs_max_size) {
					lqseqs_max_size += 100;
					lqseqs = realloc(lqseqs, lqseqs_max_size * sizeof(lqseq));
				}
			}

			while (clusters_i > 0 && e < clusters->clusters[clusters_i - 1].r.s){
				index = get_lqseqs_from_cluster(&clusters->clusters[clusters_i - 1], lqseqs, index - 1);
				clusters_i --;
				while (clusters_i > 0 && !clusters->clusters[clusters_i - 1].i_m) clusters_i--;
				if (++index >= lqseqs_max_size) {
					lqseqs_max_size += 100;
					lqseqs = realloc(lqseqs, lqseqs_max_size * sizeof(lqseq));
				}
			}

			lqseqs[index].start = s;
			lqseqs[index].end = e;
			lqseqs[index].l = 0;
			if (++index >= lqseqs_max_size) {
				lqseqs_max_size += 100;
				lqseqs = realloc(lqseqs, lqseqs_max_size * sizeof(lqseq));
			}
		}else lqseqs[index - 1].start = s;
	}
	if (dels) free (dels);
	*lqseqs_count = index;
	return lqseqs;
}

static consensus_data *generate_cns_from_best_score_lq( msa_p *msa, align_tags_t *tags_list, align_tag *global_best_p,\
		int aligned_seq_count, int len, int min_cov, float min_error_corrected_ratio, gap_clusters *clusters){
	
	int lq_min_length = 2;
	int p, qv, lq, lq_s, lq_e;
	lq_s = lq_e = -1;
	p = qv = lq = 0;

	int lqseqs_max_size = 200;
	int lqseqs_index = 0;
	lqseq *lqseqs = malloc(lqseqs_max_size * sizeof(lqseq));

	msa_p_d_b *global_best_d_b;
	consensus_data consensus;
	memset(&consensus, 0, sizeof(consensus));
	consensus.max_size = len/5 * 6 + 1;
	consensus.cns_bases = malloc( consensus.max_size * sizeof(consensus_base) );

	int DAG_MIN_QV = 80;
	int clusters_i = clusters->i;
	while (1){
		if (global_best_p->q_base != 4){

			if (p >= consensus.max_size){
				consensus.max_size += 5000;
				consensus.cns_bases = realloc(consensus.cns_bases, consensus.max_size * sizeof(consensus_base));
			}
			consensus.cns_bases[p].pos = global_best_p->t_pos;
			global_best_d_b = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base];
			qv = 100 * global_best_d_b->pp_ppp[global_best_d_b->max_size].link_count / msa[global_best_p->t_pos].coverage;
			// printf("p %d %d qn:%d depth:%d\n",p, global_best_p->t_pos, qv, msa[global_best_p->t_pos].coverage);
			if (msa[global_best_p->t_pos].coverage < 4){
				lq = 0;
				lq_s = -1;
			}else if (qv < DAG_MIN_QV){
				if (lq_s == -1) lq_s = p;
				lq_e = p;
				lq = 1;
			}else if (lq && p - lq_e > 2 * lq_min_length && consensus.cns_bases[p].pos != consensus.cns_bases[p - 1].pos){
				lq_e = p - lq_min_length - 1;
				lq_s = lq_s > lq_min_length ? lq_s - lq_min_length : 1;
				// if (lqseqs_index >= 1) printf("%d %d %d %d %d %d\n",lq_s,lq_e,consensus.cns_bases[lq_s].pos,consensus.cns_bases[lq_e].pos,lqseqs[lqseqs_index - 1].start,lqseqs[lqseqs_index - 1].end );
				if (lqseqs_index >= 1 && consensus.cns_bases[lq_s].pos >= lqseqs[lqseqs_index - 1].start){
					lqseqs[lqseqs_index - 1].start = consensus.cns_bases[lq_e].pos;
				}else{
					while (clusters_i > 0 && consensus.cns_bases[lq_s].pos < clusters->clusters[clusters_i - 1].r.s){
						lqseqs_index = get_lqseqs_from_cluster(&clusters->clusters[clusters_i - 1], lqseqs, lqseqs_index - 1);
						clusters_i --;
						while (clusters_i > 0 && !clusters->clusters[clusters_i - 1].i_m) clusters_i--;
						// printf("clusters_i %d\n", clusters_i);
						if (++lqseqs_index >= lqseqs_max_size){
							lqseqs_max_size += 100;
							lqseqs = realloc(lqseqs, lqseqs_max_size * sizeof(lqseq));
						}
					}

					lqseqs[lqseqs_index].end = consensus.cns_bases[lq_s].pos;
					lqseqs[lqseqs_index].start = consensus.cns_bases[lq_e].pos;
					lqseqs[lqseqs_index].l = 4;
					if (++lqseqs_index >= lqseqs_max_size) {
						lqseqs_max_size += 100;
						lqseqs = realloc(lqseqs, lqseqs_max_size * sizeof(lqseq));
					}
				}
				lq = 0;
				lq_s = -1;
			}

			if (msa[global_best_p->t_pos].coverage > min_cov && qv > DAG_MIN_QV){
				consensus.cns_bases[p].base = int_to_base[global_best_p->q_base];
			}else{
				consensus.cns_bases[p].base = tolower(int_to_base[global_best_p->q_base]); 
			}
			p ++;
		}

		global_best_d_b = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base];
		global_best_p = &global_best_d_b->pp_ppp[global_best_d_b->max_size].pp;

		if (global_best_p->t_pos == -1){
			break;
		}
	}
	consensus.len = p;
	free_msa(msa, len);
	reverse_consensus_base(consensus.cns_bases, consensus.len);

	// for (p = 0; p < lqseqs_index; p ++){
	// 	printf("lqseq %d: %d %d %d\n",p, lqseqs[p].start,lqseqs[p].end,
	// 		lqseqs[p].end-lqseqs[p].start+1 );
	// }
	// exit(1);

	int max_aln_length = generate_lqseqs_from_tags_kmer(lqseqs, lqseqs_index, tags_list, aligned_seq_count, clusters);
	consensus_data *consensus_ = iterate_generate_consensus_trimed(lqseqs, lqseqs_index, &consensus, len, max_aln_length, 2);
	if (lqseqs) free_lqseq(lqseqs, lqseqs_index);
	free (consensus.cns_bases);
	// free (consensus);
	return consensus_;
}

static consensus_data *generate_cns_from_best_score( msa_p *msa, align_tags_t *tags_list, align_tag *global_best_p,\
		int aligned_seq_count, int len, int min_cov, float min_error_corrected_ratio, gap_clusters *clusters){
	msa_p_d_b *global_best_d_b;
	consensus_data consensus;
	memset(&consensus, 0, sizeof(consensus));
	consensus.max_size = len/5 * 6 + 1;
	consensus.cns_bases = malloc( consensus.max_size * sizeof(consensus_base) );
	// time_debug(clock(), "generate_cns_from_best_score start.......");

	while (1){
		if (global_best_p->q_base != 4){
			global_best_d_b = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base];
			consensus.cns_bases[consensus.len].qv = 100 * global_best_d_b->pp_ppp[global_best_d_b->max_size].link_count / 
				msa[global_best_p->t_pos].coverage;

			if (msa[global_best_p->t_pos].coverage > min_cov && consensus.cns_bases[consensus.len].qv > LQBASE_MIN_QV){
				consensus.cns_bases[consensus.len].base = int_to_base[global_best_p->q_base];
			}else consensus.cns_bases[consensus.len].base = tolower(int_to_base[global_best_p->q_base]);
			consensus.cns_bases[consensus.len++].pos = global_best_p->t_pos;
			if (consensus.len >= consensus.max_size){
				consensus.max_size += 100000;
				consensus.cns_bases = realloc(consensus.cns_bases, consensus.max_size * sizeof(consensus_base));
			}
		}
		global_best_d_b = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base];
		global_best_p = &global_best_d_b->pp_ppp[global_best_d_b->max_size].pp;
		if (global_best_p->t_pos == -1){
			break;
		}
	}

	// time_debug(clock(), "generate_cns_from_best_score done...");
	reverse_consensus_base(consensus.cns_bases, consensus.len);
	
	int lqseqs_count = 0;
	lqseq *lqseqs = get_lqseqs_from_gap(msa, &consensus, &lqseqs_count, clusters);
	free_msa(msa, len);
	// time_debug(clock(), "get_lqseqs_from_gap done...");
	int max_aln_length = generate_lqseqs_from_tags(lqseqs, lqseqs_count, tags_list, aligned_seq_count, clusters);
	// time_debug(clock(), "generate_lqseqs_from_tags done...");
	consensus_data *consensus_ = iterate_generate_consensus_trimed(lqseqs, lqseqs_count, &consensus, len, max_aln_length, 2);
	if (lqseqs) free_lqseq(lqseqs, lqseqs_count);
	// time_debug(clock(), "iterate_generate_consensus_trimed done...");
	free (consensus.cns_bases);
	// free (consensus);
	return consensus_;
}

static consensus_data *get_cns_from_align_tags( align_tags_t *tags_list, msa_p *msa, int aligned_seq_count, \
		int len, int min_cov, float min_error_corrected_ratio, int fast, gap_clusters *clusters){
	int p, d, b, m, n;
	// time_debug(clock(), "malloc msa begin...");
	allocate_msa_mem(msa, len);
	// time_debug(clock(), "update msa loop begin...");
	update_msa(msa, tags_list, aligned_seq_count, 0);
	// time_debug(clock(), "update msa done...");

	int64_t p_pp_score, p_pp_score_, pp_ppp_score = -10;
	int64_t global_best_score = INT64_MIN;
	align_tag global_best_p = {0, 0, -1};

	msa_p_d_b_pp_ppp *pp_ppp_m, *pp_ppp_n;
	if (READS_TYPE == READS_RS){
		for (p = 0; p < len; p++){
			for ( d = 0; d < msa[p].max_size; d++ ){
				for (b = 0; b < 6; b++){
					msa_p_d_b *p_base = &msa[p].d_b[d][b];
					p_base->max_size = 0;//index of best pp
					p_pp_score_ = p_pp_score = INT64_MIN;
					for (m = 0; m < p_base->len; m++){
						pp_ppp_m = &p_base->pp_ppp[m];
						if (pp_ppp_m->pp.t_pos == -1){
							pp_ppp_m->score = 10 * pp_ppp_m->link_count - 3 * msa[p].coverage;
						}else{

							msa_p_d_b *pp_base = &msa[pp_ppp_m->pp.t_pos].d_b[pp_ppp_m->pp.delta][pp_ppp_m->pp.q_base];

							for (n = 0; n < pp_base->len; n++){
								pp_ppp_n = &pp_base->pp_ppp[n];

								if (pp_ppp_n->pp.t_pos == pp_ppp_m->ppp.t_pos &&
									pp_ppp_n->pp.delta == pp_ppp_m->ppp.delta &&
									pp_ppp_n->pp.q_base == pp_ppp_m->ppp.q_base){
									pp_ppp_score = pp_ppp_n->score + 10 * pp_ppp_m->link_count - 3 * msa[p].coverage;
									if (pp_ppp_score > pp_ppp_m->score) {
										pp_ppp_m->score = pp_ppp_score;
										p_pp_score_ = pp_ppp_n->score;
									}
								}
							}
						}
						if (pp_ppp_m->score >= p_base->pp_ppp[p_base->max_size].score){
							p_base->max_size = m;
							p_pp_score = p_pp_score_;
						}
					}
					if (m && p == len - 1 && p_base->pp_ppp[p_base->max_size].score >= global_best_score){
						global_best_p.t_pos = p;
						global_best_p.delta = d;
						global_best_p.q_base = b;
						if (p_base->pp_ppp[p_base->max_size].score > global_best_score){
							global_best_score = p_base->pp_ppp[p_base->max_size].score;
						}
					}
				}
			}
		}
	}else if (READS_TYPE == READS_CLR){
		for (p = 0; p < len; p++){
			for ( d = 0; d < msa[p].max_size; d++ ){
				for (b = 0; b < 6; b++){
					msa_p_d_b *p_base = &msa[p].d_b[d][b];
					p_base->max_size = 0;//index of best pp
					p_pp_score_ = p_pp_score = INT64_MIN;
					// int tmp = 0;
					// for (m = 0; m < p_base->len; m++){//TODO FIX AND REMOVE
					// 	n = p_base->pp_ppp[m].link_count;
					// 	if (n > tmp) tmp = n;
					// }
					for (m = 0; m < p_base->len; m++){
						pp_ppp_m = &p_base->pp_ppp[m];
						if (pp_ppp_m->pp.t_pos == -1){
							pp_ppp_m->score = 10 * pp_ppp_m->link_count - 3 * msa[p].coverage;
						}else{

							msa_p_d_b *pp_base = &msa[pp_ppp_m->pp.t_pos].d_b[pp_ppp_m->pp.delta][pp_ppp_m->pp.q_base];

							for (n = 0; n < pp_base->len; n++){
								pp_ppp_n = &pp_base->pp_ppp[n];

								if (pp_ppp_n->pp.t_pos == pp_ppp_m->ppp.t_pos &&
									pp_ppp_n->pp.delta == pp_ppp_m->ppp.delta &&
									pp_ppp_n->pp.q_base == pp_ppp_m->ppp.q_base){
									pp_ppp_score = pp_ppp_n->score + 10 * pp_ppp_m->link_count - 3 * msa[p].coverage;
									if (pp_ppp_score > pp_ppp_m->score) {
										pp_ppp_m->score = pp_ppp_score;
										p_pp_score_ = pp_ppp_n->score;
									}
									if (pp_ppp_n->score > p_pp_score ||
										(pp_ppp_n->score == p_pp_score && pp_ppp_m->pp.q_base != 4)){
									// if (pp_ppp_n->score > p_pp_score && (pp_ppp_m->pp.q_base == 4 || pp_ppp_m->pp.q_base == b)){
										p_base->max_size = m;
										p_pp_score = pp_ppp_n->score;
									}
								}
							}
						}
						if (pp_ppp_m->score > p_base->pp_ppp[p_base->max_size].score || \
							(pp_ppp_m->score == p_base->pp_ppp[p_base->max_size].score && pp_ppp_m->pp.q_base != 4)){
							p_base->max_size = m;
							p_pp_score = p_pp_score_;
						// if (pp_ppp_m->score > p_base->best_score || (pp_ppp_m->score == p_base->best_score && 
							// pp_ppp_m->pp.q_base != 4)){// find the bese score for a given p
							// p_base->best_score = pp_ppp_m->score;
							// p_base->best_pp = pp_ppp_m->pp;
							// p_base->best_link_count = pp_ppp_m->link_count;
						}
	// printf("%d %d %c %d %d %c %d %d %c %ld %d %ld\n",p,d,int_to_base[b], pp_ppp_m->pp.t_pos,pp_ppp_m->pp.delta,
	// 	int_to_base[pp_ppp_m->pp.q_base],pp_ppp_m->ppp.t_pos,pp_ppp_m->ppp.delta,
	// 	int_to_base[pp_ppp_m->ppp.q_base],pp_ppp_m->score, pp_ppp_m->link_count, p_base->pp_ppp[p_base->max_size].score);
					}
					if (m && p == len - 1 && p_base->pp_ppp[p_base->max_size].score >= global_best_score){
						global_best_p.t_pos = p;
						global_best_p.delta = d;
						global_best_p.q_base = b;
						if (p_base->pp_ppp[p_base->max_size].score > global_best_score){
							global_best_score = p_base->pp_ppp[p_base->max_size].score;
						}
					}
				}
			}
		}
	}else if (READS_TYPE == READS_HIFI){
		for (p = 0; p < len; p++){
			for ( d = 0; d < msa[p].max_size; d++ ){
				for (b = 0; b < 6; b++){
					msa_p_d_b *p_base = &msa[p].d_b[d][b];
					p_base->max_size = 0;//index of best pp
					p_pp_score_ = p_pp_score = INT64_MIN;
					// int tmp = 0;
					// for (m = 0; m < p_base->len; m++){//TODO FIX AND REMOVE
					// 	n = p_base->pp_ppp[m].link_count;
					// 	if (n > tmp) tmp = n;
					// }
					for (m = 0; m < p_base->len; m++){
						pp_ppp_m = &p_base->pp_ppp[m];
						if (pp_ppp_m->pp.t_pos == -1){
							pp_ppp_m->score = 10 * pp_ppp_m->link_count - 4 * msa[p].coverage;
						}else{

							msa_p_d_b *pp_base = &msa[pp_ppp_m->pp.t_pos].d_b[pp_ppp_m->pp.delta][pp_ppp_m->pp.q_base];

							for (n = 0; n < pp_base->len; n++){
								pp_ppp_n = &pp_base->pp_ppp[n];

								if (pp_ppp_n->pp.t_pos == pp_ppp_m->ppp.t_pos &&
									pp_ppp_n->pp.delta == pp_ppp_m->ppp.delta &&
									pp_ppp_n->pp.q_base == pp_ppp_m->ppp.q_base){
									pp_ppp_score = pp_ppp_n->score + 10 * pp_ppp_m->link_count - 4 * msa[p].coverage;
									if (pp_ppp_score > pp_ppp_m->score) {
										pp_ppp_m->score = pp_ppp_score;
										p_pp_score_ = pp_ppp_n->score;
									}
									if (pp_ppp_n->score > p_pp_score ||
										(pp_ppp_n->score == p_pp_score && pp_ppp_m->pp.q_base != 4)){
									// if (pp_ppp_n->score > p_pp_score && (pp_ppp_m->pp.q_base == 4 || pp_ppp_m->pp.q_base == b)){
										p_base->max_size = m;
										p_pp_score = pp_ppp_n->score;
									}
								}
							}
						}
						if (pp_ppp_m->score > p_base->pp_ppp[p_base->max_size].score || \
							(pp_ppp_m->score == p_base->pp_ppp[p_base->max_size].score && pp_ppp_m->pp.q_base != 4)){
							p_base->max_size = m;
							p_pp_score = p_pp_score_;
						// if (pp_ppp_m->score > p_base->best_score || (pp_ppp_m->score == p_base->best_score && 
							// pp_ppp_m->pp.q_base != 4)){// find the bese score for a given p
							// p_base->best_score = pp_ppp_m->score;
							// p_base->best_pp = pp_ppp_m->pp;
							// p_base->best_link_count = pp_ppp_m->link_count;
						}
	// printf("%d %d %c %d %d %c %d %d %c %ld %d %ld\n",p,d,int_to_base[b], pp_ppp_m->pp.t_pos,pp_ppp_m->pp.delta,
	// 	int_to_base[pp_ppp_m->pp.q_base],pp_ppp_m->ppp.t_pos,pp_ppp_m->ppp.delta,
	// 	int_to_base[pp_ppp_m->ppp.q_base],pp_ppp_m->score, pp_ppp_m->link_count, p_base->pp_ppp[p_base->max_size].score);
					}
					if (m && p == len - 1 && p_base->pp_ppp[p_base->max_size].score >= global_best_score){
						global_best_p.t_pos = p;
						global_best_p.delta = d;
						global_best_p.q_base = b;
						if (p_base->pp_ppp[p_base->max_size].score > global_best_score){
							global_best_score = p_base->pp_ppp[p_base->max_size].score;
						}
					}
				}
			}
		}
	}else{
		for (p = 0; p < len; p++){
			for ( d = 0; d < msa[p].max_size; d++ ){
				for (b = 0; b < 6; b++){
					msa_p_d_b *p_base = &msa[p].d_b[d][b];
					p_base->max_size = 0;//index of best pp
					p_pp_score_ = p_pp_score = INT64_MIN;
					int tmp = 0;
					for (m = 0; m < p_base->len; m++){//TODO FIX AND REMOVE
						n = p_base->pp_ppp[m].link_count;
						if (n > tmp) tmp = n;
					}
					for (m = 0; m < p_base->len; m++){
						pp_ppp_m = &p_base->pp_ppp[m];
						if (pp_ppp_m->pp.t_pos == -1){
							pp_ppp_m->score = 10 * pp_ppp_m->link_count - 3 * msa[p].coverage;
						}else{

							msa_p_d_b *pp_base = &msa[pp_ppp_m->pp.t_pos].d_b[pp_ppp_m->pp.delta][pp_ppp_m->pp.q_base];

							for (n = 0; n < pp_base->len; n++){
								pp_ppp_n = &pp_base->pp_ppp[n];

								if (pp_ppp_n->pp.t_pos == pp_ppp_m->ppp.t_pos &&
									pp_ppp_n->pp.delta == pp_ppp_m->ppp.delta &&
									pp_ppp_n->pp.q_base == pp_ppp_m->ppp.q_base){
									pp_ppp_score = pp_ppp_n->score + 10 * pp_ppp_m->link_count - 3 * msa[p].coverage;
									if (pp_ppp_score > pp_ppp_m->score) {
										pp_ppp_m->score = pp_ppp_score;
										p_pp_score_ = pp_ppp_n->score;
									}
		
									if (((pp_ppp_m->ppp.delta > 1 || pp_ppp_m->pp.delta > 0) && 
										(pp_ppp_m->link_count > msa[p].coverage * 0.2 || pp_ppp_m->link_count > tmp/2))|| 
										(pp_ppp_m->link_count > p_base->pp_ppp[p_base->max_size].link_count/2 && 
										pp_ppp_n->score > p_pp_score && (pp_ppp_m->pp.q_base == 4 || pp_ppp_m->pp.q_base == b
										|| pp_ppp_m->ppp.q_base == b || pp_ppp_m->pp.q_base == pp_ppp_m->ppp.q_base))){
										p_base->max_size = m;
										p_pp_score = pp_ppp_n->score;
									}
								}
							}
						}
						if (pp_ppp_m->score > p_base->pp_ppp[p_base->max_size].score || \
							(pp_ppp_m->score == p_base->pp_ppp[p_base->max_size].score && pp_ppp_m->pp.q_base != 4)){
							p_base->max_size = m;
							p_pp_score = p_pp_score_;
						}
					}
					if (m && p == len - 1 && p_base->pp_ppp[p_base->max_size].score >= global_best_score){
						global_best_p.t_pos = p;
						global_best_p.delta = d;
						global_best_p.q_base = b;
						if (p_base->pp_ppp[p_base->max_size].score > global_best_score){
							global_best_score = p_base->pp_ppp[p_base->max_size].score;
						}
					}
				}
			}
		}
	}
	// printf("%d %d %d\n",global_best_p.t_pos,global_best_p.delta,global_best_p.q_base );
	// free (pp_ppps.pp_ppp);
	// time_debug(clock(), "update score chain done...");

	consensus_data *consensus;
	if (fast){
		for (p = 0; p < aligned_seq_count; p++) free (tags_list[p].align_tags);
		free (tags_list);
		consensus = generate_cns_from_best_score_fast(msa, &global_best_p, \
			aligned_seq_count, len, min_cov, min_error_corrected_ratio);
	}else if (READS_TYPE == READS_HIFI){
		consensus = generate_cns_from_best_score_lq(msa, tags_list, &global_best_p, \
			aligned_seq_count, len, min_cov, min_error_corrected_ratio, clusters);
	}else{
		consensus = generate_cns_from_best_score(msa, tags_list, &global_best_p, \
			aligned_seq_count, len, min_cov, min_error_corrected_ratio, clusters);
	}
	return consensus;
}

void free_consensus_trimed(consensus_trimed *consensus_trimed){
	free(consensus_trimed->seq);
	free (consensus_trimed);
}

void free_consensus_trimed_data(consensus_trimed_data *cons_trimed_data){
	int i;
	for (i = 0; i < cons_trimed_data->i_m; i++) free (cons_trimed_data->data[i].seq);
	free (cons_trimed_data->data);
	free (cons_trimed_data);
}

void set_satags(uint8_t *s, satags *sas){
	satag *sa;
	sas->i = 0;
	char *token = strtok(bam_aux2Z(s), ",");
	while (1){
		if (!token) break;
		if (sas->i >= sas->i_m) {
			sas->i_m += 10;
			if (sas->i_m - 10) sas->sa = realloc(sas->sa, sas->i_m * sizeof(satag));
			else sas->sa = malloc(sas->i_m * sizeof(satag));
		}
		sa = &sas->sa[sas->i++];
		sa->rname = token;
		token = strtok(NULL, ",");
		sa->pos = atoll(token) - 1;
		token = strtok(NULL, ",");
		sa->strand = *token == '+' ? 0 : 1;
		sa->cigar = strtok(NULL, ",");
		strtok(NULL, ",");strtok(NULL, ";");
		token = strtok(NULL, ",");
	}
}

refs_ *refs_init(int l){
	refs_ *refs = (refs_ *) malloc(sizeof(refs_));
	refs->i = 0;
	refs->i_m = l;
	refs->ref = calloc(refs->i_m, sizeof(ref_));
	return refs;
}

void refs_realloc(refs_ *refs, int l){
	refs->i_m += l;
	refs->ref = realloc(refs->ref, refs->i_m * sizeof(ref_));
	memset(refs->ref + refs->i_m - l, 0, l * sizeof(ref_));
}

void refs_destroy(refs_ *refs){
	uint32_t i;
	for (i = 0; i < refs->i; i ++){
		if (refs->ref[i].n) free(refs->ref[i].n);
		if (refs->ref[i].s) free(refs->ref[i].s);
		if (refs->ref[i].qv) free(refs->ref[i].qv);
	}
	free (refs->ref);
	free (refs);
}

int binary_search(char *target, char **array, int n) {
	int low = 0, high = n, middle = 0;
	while (low < high) {
		middle = (low + high)/2;
		if (strcmp(target, array[middle]) == 0){
			return 1;
		} else if (strcmp(target, array[middle]) < 0) {
			high = middle;
		} else if (strcmp(target, array[middle]) > 0) {
			low = middle + 1;
		}
	}
	return 0;
}

static int str_cmp(const void* a, const void* b){
	return strcmp(*(char**)a, *(char**)b);
}


static int startswith(const char *pre, const char *str)
{
    size_t lenpre = strlen(pre),
           lenstr = strlen(str);
    return lenstr < lenpre ? 0 : memcmp(pre, str, lenpre) == 0;
}

void set_ref_qv(char *s, ref_ *r){
	char *token, *qv = NULL, sep[2] = " ";

	r->qv_l = 0;
	if (s){
		token = strtok(s, sep);
		while (token != NULL){
			if (startswith("node", token)) r->qv_l = atoi(token + 7);
			if (startswith("qv", token)) qv = token + 5;
			token = strtok(NULL, sep);
		}
	}
	// printf("%d %s\n", r->qv_l, qv);
	uint64_t t, i;
	ref_qv *q;
	if (r->qv_l && qv){
		i = 0;
		sep[0] = ':';
		r->qv = malloc(r->qv_l * sizeof(ref_qv));
		token = strtok(qv, sep);
		while (token != NULL){
			q = &r->qv[i++];
			t = strtoull(token, NULL, 16);
			q->p = t >> 32;
			q->ide = t >> 20 & 0x3ff;
			q->ort = t >> 10 & 0x3ff;
			q->irt = t & 0x3ff;
			// printf("%lu %u %u %u %u\n",t, q->p, q->ide, q->ort, q->irt);
			token = strtok(NULL, sep);
		}
	}else{
		r->qv_l = 0;
		r->qv = NULL;
	}
}

refs_ *read_ref(char *reff, char **accept_names, int accept_names_len){
	gzFile fp = gzopen(reff, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error! %s does not exist!", reff);
		exit(1);
	}

	if (accept_names_len) qsort(accept_names, accept_names_len, sizeof(char *), str_cmp);

	ref_ *r = NULL;
	refs_ *refs = refs_init(1000);
	kseq_t *seq = kseq_init(fp);
	while(kseq_read(seq) >= 0) {
		r = &refs->ref[refs->i];
		if (accept_names_len && !binary_search(seq->name.s, accept_names, accept_names_len)) continue;
		r->n = strdup(seq->name.s);
		r->length = seq->seq.l;
		r->s = malloc(sizeof(uint32_t) * (r->length/16 + 1));
		seq2bit1(r->s, r->length, seq->seq.s);
		set_ref_qv(seq->comment.s, r);
		if (++refs->i >= refs->i_m) refs_realloc(refs, 100);
	}
	kseq_destroy(seq);
	gzclose(fp);
	// exit(1);
	return refs;
}

void out_ref(refs_ *refs){
	uint32_t i;
	ref_ *r = NULL;
	uint32_t l_m = 100000;
	char *s = malloc(l_m);
	for (i = 0; i < refs->i; i++){
		r = &refs->ref[i];
		if (r->length > l_m){
			l_m = r->length + 1;
			s = realloc(s, l_m);
		}
		bit2seq1(r->s, r->length, s);
		printf(">%s %u\n%s\n", r->n, r->length, s);
	}
	free (s);
}

uint32_t cigarint2ul(const uint32_t *cigar, const uint32_t n_cigar, const int end){
	int index = end ? n_cigar - 1 : 0;
	int32_t len = bam_cigar_oplen(cigar[index]);
	uint8_t curcigar = bam_cigar_op(cigar[index]);
	if (curcigar == BAM_CSOFT_CLIP || curcigar == BAM_CHARD_CLIP) return len;
	else return 0;
}

inline double cal_clip_len(bam1_t *read){
	uint32_t *cigar = bam_get_cigar(read);
	int32_t len = bam_cigar_oplen(cigar[0]), cliplen = 0;
	uint8_t curcigar = bam_cigar_op(cigar[0]);
	if (curcigar == BAM_CSOFT_CLIP || curcigar == BAM_CHARD_CLIP) cliplen += len;
	
	len = bam_cigar_oplen(cigar[read->core.n_cigar - 1]);
	curcigar = bam_cigar_op(cigar[read->core.n_cigar - 1]);
	if (curcigar == BAM_CSOFT_CLIP || curcigar == BAM_CHARD_CLIP) cliplen += len;
	
	return read->core.l_qseq > 0 ? cliplen / (double) read->core.l_qseq : 0;
}

int32_t cal_l_qseq_from_cigar(const uint32_t *cigar, const uint32_t n_cigar){
	uint32_t i, n, c, rlen = 0;
	for (i = 0; i < n_cigar; i++){
		n = bam_cigar_oplen(cigar[i]);
		c = bam_cigar_op(cigar[i]);
		switch (c){
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
			case BAM_CMATCH:
			case BAM_CEQUAL:
			case BAM_CDIFF:
			case BAM_CINS:
				rlen += n;
				break;
		}
	}
	return rlen;
}

int32_t cal_l_qseq(bam1_t *read){
	uint32_t *cigar = bam_get_cigar(read);
	uint8_t curcigar = bam_cigar_op(cigar[0]);
	if (!read->core.l_qseq) return cal_l_qseq_from_cigar(cigar, read->core.n_cigar);
	if (curcigar == BAM_CSOFT_CLIP) return read->core.l_qseq;
	int32_t len = bam_cigar_oplen(cigar[0]), rlen = read->core.l_qseq;
	if (curcigar == BAM_CHARD_CLIP) rlen += len;
	
	len = bam_cigar_oplen(cigar[read->core.n_cigar - 1]);
	curcigar = bam_cigar_op(cigar[read->core.n_cigar - 1]);
	if (curcigar == BAM_CHARD_CLIP) rlen += len;
	return rlen;
}

uint32_t cigarstr2ul(const char *s, const int end){
	if (end){
		int index = 0;
		while (*(s + 1) != '\0'){
			if (*s >= '0' &&  *s <= '9') index ++;
			else index = 0;
			s ++;
		}
		s -= index;
	}

	uint32_t result = 0;
	while (*s >= '0' &&  *s <= '9') {
		result = result * 10 + *s - '0';
		s ++;
	}
	if (*s != 'H' && *s != 'S') result = 0;
	return result;
}

int32_t cigarstr2rlen(const char *s){
	uint32_t rlen, clen;
	rlen = clen = 0;
	while (*s != '\0'){
		if (*s >= '0' &&  *s <= '9'){
			clen = clen * 10 + *s - '0';
		}else{
			if(*s == 'M' || *s == 'D') rlen += clen;
			clen = 0;
		}
		s ++;
	}
	return rlen;
}

int32_t bam2aln(alignment *aln, const char *rfseq, const uint8_t *rdseqi, const uint32_t *cigar, \
		const uint32_t n_cigar){
	uint32_t i, n, c, rdi = 0, rfi = aln->aln_t_s;
	for (i = 0; i < n_cigar; i++){
		n = bam_cigar_oplen(cigar[i]);
		c = bam_cigar_op(cigar[i]);
		switch (c){
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				rdi += n;
				break;
			case BAM_CREF_SKIP:
				aln->aln_t_s += n;
				rfi += n; 
				break;
			case BAM_CMATCH:
			case BAM_CINS:
			case BAM_CDEL:
				while (aln->aln_len + n >= aln->max_aln_len) {
					aln->max_aln_len += 10000;
					aln->t_aln_str = realloc(aln->t_aln_str, aln->max_aln_len * sizeof(char));
					aln->q_aln_str = realloc(aln->q_aln_str, aln->max_aln_len * sizeof(char));
				}
				switch (c){
					case BAM_CMATCH:
						while (n--){
							aln->t_aln_str[aln->aln_len] = rfseq[rfi++];
							aln->q_aln_str[aln->aln_len++] = seq_nt16_str[bam_seqi(rdseqi,rdi)];
							// printf("M%d %d %c\n",rdi,aln->aln_len, aln->q_aln_str[aln->aln_len-1]);
							rdi ++;
						}
						break;
					case BAM_CINS:
						while (n--){
							aln->t_aln_str[aln->aln_len] = '-';
							aln->q_aln_str[aln->aln_len++] = seq_nt16_str[bam_seqi(rdseqi,rdi)];
							rdi ++;
						}
						break;
					case BAM_CDEL:
						while (n--){
							aln->t_aln_str[aln->aln_len] = rfseq[rfi++];
							aln->q_aln_str[aln->aln_len++] = '-';
						}
						break;
				}
				break;
			default:
				fprintf(stderr, "Error, unexpected cigar: %d%c", n, BAM_CIGAR_STR[c]);
				return 1;
		}
	}
	return rfi;
}


#define TEM_CLIP_RATIO 0.1
#define MAX_GAP_LEN 30000
#define mabs(x, y) ((x) > (y) ? (x) - (y) : (y) - (x))

void check_indel(gap *g, const int32_t rlen, const pos *rfp1, const pos *rdp1, const pos *rfp2, const pos *rdp2){
	int l = 0;
	int32_t mclen = rlen * TEM_CLIP_RATIO;
	if (rfp1->s > rfp2->s){
		l = 1;
		const pos *t = rfp1; rfp1 = rfp2; rfp2 = t;
		t = rdp1; rdp1 = rdp2; rdp2 = t;
	}
	// printf("SA1: %u %u %u %u %u %u %u %u %u\n",rlen, rfp1->s, rfp1->e, rdp1->s, rdp1->e, rfp2->s, rfp2->e, rdp2->s, rdp2->e);
	if (rfp2->e > rfp1->e && rdp2->e > rdp1->e && \
		rdp1->s < mclen && rdp2->e > rlen - mclen && \
		mabs(rfp2->s, rfp1->e) < MAX_GAP_LEN && mabs(rdp2->s, rdp1->e) < MAX_GAP_LEN && rfp1->s != rfp2->s){//TODO check
		uint32_t score = rdp1->s + rlen - rdp2->e + mabs(rfp2->s, rfp1->e) + mabs(rdp2->s, rdp1->e);
		if (score < g->score || !g->score){
			// g->s1 = rfp1->s;
			// g->e1 = rfp1->e;
			// g->s2 = rfp2->s;
			g->score = score;
			g->ds = l ? rdp1->s : rdp2->s;
			g->fs = l ? rfp1->s : rfp2->s;
			if (rfp1->e < rfp2->s){
				g->gap.s = rfp1->e;
				g->gap.e = rfp2->s;
			}else{
				g->gap.s = rfp2->s;
				g->gap.e = rfp1->e;
			}
		}
	}
}

static int compare_gap_by_media(const void * a, const void * b){
	gap_ **g1 = (gap_ **) a;
	gap_ **g2 = (gap_ **) b;
	uint32_t p1 = (*g1)->gap.s + (*g1)->gap.e;
	uint32_t p2 = (*g2)->gap.s + (*g2)->gap.e;
	return (p1 > p2) - (p1 < p2);
}

static int compare_gap_by_start(const void * a, const void * b){
	gap_ *g1 = (gap_ *) a;
	gap_ *g2 = (gap_ *) b;
	if (g1->gap.s != g2->gap.s) return (g1->gap.s > g2->gap.s) - (g1->gap.s < g2->gap.s);
	else return (g1->gap.e > g2->gap.e) - (g1->gap.e < g2->gap.e);
}

static void cal_gap_cluster_median(gap_cluster *clu){
	int32_t i, j, count_m, count_t, count_t_diff;
	uint32_t s, e, median, median_t, count_mc, offset = 10;
	uint64_t count_m_diff;
	while (offset <= 100){
		clu->median = count_m = count_mc = count_m_diff = 0;
		for (i = 0; i < clu->i_m; i ++){	
			median = (clu->gap[i]->gap.s + clu->gap[i]->gap.e)/2;
			if (median == clu->median) continue;
			s = median > offset ? median - offset : 0;
			e = median + offset;
			
			for (count_t = count_t_diff = 0, j = i - 1; j >= 0; j--){
				median_t = (clu->gap[j]->gap.s + clu->gap[j]->gap.e)/2;
				if (median_t >= s) {
					count_t ++;
					count_t_diff += mabs(median_t, median);
				}else break;
			}
			for (j = i + 1; j < clu->i_m; j++){
				median_t = (clu->gap[j]->gap.s + clu->gap[j]->gap.e)/2;
				if (median_t <= e) {
					count_t ++;
					count_t_diff += mabs(median_t, median);
				}else break;
			}

			if (count_t > count_m || (count_t == count_m && count_m_diff > count_t_diff)) {
				count_m = count_t;
				count_mc = median;
				count_m_diff = count_t_diff;
			}
			// printf("%d %d %d %d %d %d %d %d\n",i,clu->gap[i]->gap.s,clu->gap[i]->gap.e,median,count_t, count_m, count_mc, offset);
		}
		if (count_m >= max(3, clu->i_m/6)){
			clu->median = count_mc;
			break;
		}
		offset += 10;
	}
	if (offset > 100) clu->median = (clu->gap[clu->i_m/2]->gap.s + clu->gap[clu->i_m/2]->gap.e)/2;
}

int update_gap_cluster(gaps *gs, gap_clusters *clusters, uint16_t *ref_ds, const int w, const int d, const int32_t ref_s){
	if (d < 10) return 0;
	if (!clusters->i_m) {
		clusters->i_m = 10;
		clusters->clusters = malloc(clusters->i_m * sizeof(gap_cluster));
	}

	int32_t i, j, t, e, p;
	int md = d * CLUSTER_MIN_DEPTH_RATIO;
	gap_cluster *clu;

	qsort(gs->gap, gs->i, sizeof(gap_), compare_gap_by_start);
	// for (i = 0; i < gs->i; i++){
	// 	p = (gs->gap[i].gap.s + gs->gap[i].gap.e)/2 - ref_s;
	// 	printf("#gaps: %d %d %d %d %d %d %d\n",i, gs->gap[i].gap.s,  gs->gap[i].gap.e, p, ref_ds[p/INS_WIN_STEP], ref_s, w);
	// }
	// exit(1);

	for (clusters->i = i = 0; i < (int32_t)gs->i - md; i++){
		p = (gs->gap[i].gap.s + gs->gap[i].gap.e)/2 - ref_s;
		if (p < w || ref_ds[p/INS_WIN_STEP] >= d/2) continue;
		e = gs->gap[i].gap.e;
		clu = &clusters->clusters[clusters->i];
		clu->i_m = 0;
		t = 1;
		for (j = i + 1; j < gs->i && gs->gap[j].gap.s <= e; j ++){
			if (ref_ds[((gs->gap[j].gap.s + gs->gap[j].gap.e)/2 - ref_s)/INS_WIN_STEP] >= d/2) continue;
			t ++;
			if (gs->gap[j].gap.e > e) e = gs->gap[j].gap.e;
			if (clu->i_m < LQSEQ_MAX_CAN_COUNT<<1) clu->gap[clu->i_m++] = &gs->gap[j];
		}
		// printf("%d %d %d %d %d\n",i,j, gs->gap[i].gap.s,  gs->gap[i].gap.e, t);
		i = j - 1;
		if (clu->i_m > md && ref_ds[p/INS_WIN_STEP] < t) {
			clusters->i++;
			if (clusters->i >= clusters->i_m){
				clusters->i_m += 10;
				clusters->clusters = realloc(clusters->clusters, clusters->i_m * sizeof(gap_cluster));
			}
		}
	}

	for (i = t = 0; i < clusters->i; i ++){
		clu = &clusters->clusters[i];
		t += clu->i_m;
		qsort(clu->gap, clu->i_m, sizeof(gap_ *), compare_gap_by_media);
		//here we may need to calculate two or more medians because one cluster may include two or more sub-clusters.
		cal_gap_cluster_median(clu);
	}

	// for (i = 0; i < clusters->i; i ++){
	// 	clu = &clusters->clusters[i];
	// 	for (j = 0; j < clu->i_m; j ++){
	// 		printf("gap_clusters: %d %d s:%d e:%d media:%d depth:%d %d %d\n",i, j, clu->gap[j]->gap.s,
	// 			clu->gap[j]->gap.e,(clu->gap[j]->gap.s + clu->gap[j]->gap.e)/2, 
	// 			ref_ds[(clu->gap[j]->gap.s + clu->gap[j]->gap.e-ref_s-ref_s)/2/INS_WIN_STEP],clu->median,clu->gap[j]->l);
	// 	}
	// }
	// exit(1);
	return t;
}

void free_gap_clusters(gap_clusters *clusters){
	free (clusters->clusters);
}

void clean_gap_clusters_lable(gap_clusters *clusters){
	uint32_t i, j;
	gap_cluster *clu;
	for (i = 0; i < clusters->i; i ++){
		clu = &clusters->clusters[i];
		for (j = 0; j < clu->i_m; j ++) clu->gap[j]->l = 0;
	}
}

void update_gap_info(gaps *gs, alignment *aln, gap *g, uint8_t *s, uint32_t s_l, uint32_t p_id){
	uint32_t i, m = gs->i_m;
	if (gs->i >= gs->i_m) {
		gs->i_m += 1000;
		if (m) gs->gap = realloc(gs->gap, gs->i_m * sizeof(gap_));
		else gs->gap = malloc(gs->i_m * sizeof(gap_));
		for (i = m; i < gs->i_m; i ++) {
			gs->gap[i].dl_m = 50000;
			gs->gap[i].dseq = malloc(gs->gap[i].dl_m * sizeof(char));
		}
	}
	gap_ *g_ = &gs->gap[gs->i++];
	g_->gap = g->gap;
	g_->p_id = p_id;
	g_->p_s = aln->aln_q_s;
	g_->s_id = g->fs;
	g_->s_s = g->ds;
	g_->l = 0;

	i = (s_l + 1)/2;
	if (i > g_->dl_m){
		g_->dl_m = i;
		g_->dseq = realloc(g_->dseq, g_->dl_m * sizeof(uint8_t));
	}
	memcpy(g_->dseq, s, i * sizeof(uint8_t));
}

void free_gap_info(gaps *gs){
	uint32_t i;
	for (i = 0; i < gs->i_m; i ++) free (gs->gap[i].dseq);
	free (gs->gap);
}

void update_sup_alns(sup_alns *aln, uint32_t fs, uint32_t ds, uint32_t l, uint32_t *cigar){
	uint32_t i, m = aln->i_m;
	if (aln->i >= aln->i_m) {
		aln->i_m += 1000;
		if (m) aln->aln = realloc(aln->aln, aln->i_m * sizeof(sup_aln));
		else aln->aln = malloc(aln->i_m * sizeof(sup_aln));
		for (i = m; i < aln->i_m; i ++) {
			aln->aln[i].i_m = 20000;
			aln->aln[i].cigar = malloc(aln->aln[i].i_m * sizeof(uint32_t));
		}
	}

	sup_aln *aln_ = &aln->aln[aln->i++];
	aln_->fs = fs;
	aln_->ds = ds;
	aln_->i = l;
	if (aln_->i > aln_->i_m){
		aln_->i_m = aln_->i;
		aln_->cigar = realloc(aln_->cigar, aln_->i_m * sizeof(uint32_t));
	}
	memcpy(aln_->cigar, cigar, l * sizeof(uint32_t));
}

void free_sup_alns(sup_alns *aln){
	uint32_t i;
	for (i = 0; i < aln->i_m; i ++) free (aln->aln[i].cigar);
	free (aln->aln);
}

static int32_t find_low_depth_edge(const uint16_t *r, int32_t s, const int l, const int d, const int lable){
	int md = d * INS_MIN_DEPTH_RATIO * 2;
	if (lable) while (s > 1 && r[s] <= md) s--;
	else while (s < l && r[s] <= md) s++;
	return s;
}

void update_ld_regs(ld_regs *regs, const uint16_t *r, const int32_t l, const int w, const int d, const int32_t s){
	if (!regs->i_m) {
		regs->i_m = 50;
		regs->reg = malloc(regs->i_m * sizeof(pos));
	}
	regs->i = 0;

	int init_data, has_data;
	init_data = has_data = 0;
	int32_t i, t, md = d * INS_MIN_DEPTH_RATIO;
	for (i = 0; i < l; i ++){
		if (r[i] <= md){
			if (!init_data) {
				t = find_low_depth_edge(r, i, l, d, 1);
				// regs->reg[regs->i].s = t > 1 ? t * INS_WIN_STEP + w : 0;
				regs->reg[regs->i].s = t > 1 ? t * INS_WIN_STEP: 0;
				t = find_low_depth_edge(r, i, l, d, 0);
				// regs->reg[regs->i].e = (t - 1) * INS_WIN_STEP;
				regs->reg[regs->i].e = (t - 1) * INS_WIN_STEP + w;
				i = t;
				init_data = 1;
			}else{
				t = find_low_depth_edge(r, i, l, d, 1);
				t = t * INS_WIN_STEP;
				if (t > regs->reg[regs->i].e + INS_WIN_DIV/2 * w) {
					regs->reg[++regs->i].s = t;
					if (regs->i >= regs->i_m - 1){
						regs->i_m += 50;
						regs->reg = realloc(regs->reg, regs->i_m * sizeof(pos));
					}
				}
				t = find_low_depth_edge(r, i, l, d, 0);
				regs->reg[regs->i].e = (t - 1) * INS_WIN_STEP + w;
				i = t;
			}
			if (regs->reg[regs->i].s > regs->reg[regs->i].e) SWAP(regs->reg[regs->i].s, regs->reg[regs->i].e, uint32_t);
		}
	}

	if (init_data) regs->i ++;

	// for (i = 0; i < regs->i; i ++){
	// 	printf("ld_regs1: %d s:%d(%d) e:%d(%d) media:%d len:%d depth:%d bin_len:%d mean_depth:%d\n",i,
	// 		regs->reg[i].s,regs->reg[i].s+s,regs->reg[i].e, regs->reg[i].e+s,
	// 		(regs->reg[i].s + regs->reg[i].e)/2, l * INS_WIN_STEP,
	// 		r[(regs->reg[i].s + regs->reg[i].e)/2/INS_WIN_STEP],
	// 		 w, d);
	// }
	// exit(1);
}

static int compare_pos_by_start(const void * a, const void * b){
	pos *g1 = (pos *) a;
	pos *g2 = (pos *) b;
	return g1->s == g2->s ? g1->e - g2->e : g1->s - g2->s;
}

void update_ld_regs_with_refqv(ld_regs *regs, const uint16_t *r, ref_ *ref, const int32_t w, const int32_t s_t,
		const int32_t e_t, const int32_t d_t, const uint32_t ide_t, const uint32_t ort_t, const uint32_t irt_t){
	
	int32_t i, p, s, e, l, t = 0;
	for (i = 0; i < ref->qv_l && ref->qv[i].p < e_t; i ++){
		if (ref->qv[i].p < s_t) continue;
		// printf("%d %d %d %d %d %d %d\n",ref->qv[i].p, ref->qv[i].ide,ref->qv[i].ort,ref->qv[i].irt, ide_t, ort_t,  irt_t );
		if (ref->qv[i].ide < ide_t && ref->qv[i].ort < ort_t && ref->qv[i].irt < irt_t){
			s = ref->qv[i].p > w * 2 + s_t ? (ref->qv[i].p - w * 2 - s_t) / INS_WIN_STEP : 0;
			e = ref->qv[i].p + w * 2 < e_t ?  (ref->qv[i].p + w * 2 - s_t) / INS_WIN_STEP : (e_t - s_t)/INS_WIN_STEP;
			for (l = 0, p = s; p <= e && !l; p++) if (r[p] <= d_t) l = 1;
			if (l){
				t ++;
				regs->reg[regs->i].s = ref->qv[i].p - s_t;
				regs->reg[regs->i].e = ref->qv[i].p + 1 - s_t;
				if (++ regs->i >= regs->i_m){
					regs->i_m += 50;
					regs->reg = realloc(regs->reg, regs->i_m * sizeof(pos));
				}
			}
		}
	}

	if (t){
		qsort (regs->reg, regs->i, sizeof(pos), compare_pos_by_start);
		for (i = 1; i < regs->i; i ++){
			if (regs->reg[i].s < regs->reg[i - 1].e + INS_WIN_DIV/2 * w){
				regs->reg[i].s = regs->reg[i - 1].s;
				if (regs->reg[i].e < regs->reg[i - 1].e) regs->reg[i].e = regs->reg[i - 1].e;
				regs->reg[i - 1].s = regs->reg[i - 1].e = 0;
			}
		}
		// for (i = 0; i < regs->i; i ++){
		// 	printf("ld_regs2: %d s:%d(%d) e:%d(%d) media:%d len:%d depth:%d bin_len:%d lq_depth:%d\n",i,
		// 		regs->reg[i].s,regs->reg[i].s+s_t,regs->reg[i].e, regs->reg[i].e+s_t,
		// 		(regs->reg[i].s + regs->reg[i].e)/2, l * INS_WIN_STEP,
		// 		r[(regs->reg[i].s + regs->reg[i].e)/2/INS_WIN_STEP],
		// 		w, d_t);
		// }
		// exit(1);
	}
}

void free_ld_regs(ld_regs *regs){
	free (regs->reg);
}

static int cal_win_len(const int w, const int s, const uint64_t l){
	int b = l;
	if (l > w){
		int n = (int)((float)(l - s)/(w - s) + 0.999);
		b = (int)((float)(l + (n - 1) * s)/n + 0.999);
	}
	return b;
}

static void clip_aln(alignment *aln, const int32_t s, const int32_t e, const int l){
	int32_t s_ = 0, e_ = aln->aln_len - 1;
	while (aln->aln_t_s < s){
		if (aln->t_aln_str[s_++] != '-') aln->aln_t_s ++;
		if (l && aln->q_aln_str[s_ - 1] != '-') aln->aln_q_s ++;
	}
	while (aln->t_aln_str[s_] == '-') s_++;

	while (aln->aln_t_e > e){
		if (aln->t_aln_str[e_--] != '-') aln->aln_t_e --;
		if (l && aln->q_aln_str[e_ + 1] != '-') aln->aln_q_e --;
	}
	if (e_ > s_ + 500){
		aln->aln_len = e_ - s_ + 1;
		memmove(aln->q_aln_str, aln->q_aln_str + s_, aln->aln_len);
		memmove(aln->t_aln_str, aln->t_aln_str + s_, aln->aln_len);
	}else aln->aln_len = 10;
}

sup_aln *find_sup_alns(sup_alns *sup_aln, const uint32_t fs, const uint32_t ds){
	uint32_t i;
	for (i = 0; i < sup_aln->i; i ++){
		if (sup_aln->aln[i].fs == fs && sup_aln->aln[i].ds == ds) return &sup_aln->aln[i];
	}
	assert(i != sup_aln->i);
	return NULL;
}

uint32_t update_align_tags(gap_clusters *clusters, sup_alns *sup_alns, 
	align_tags_t *tags_list, uint32_t seq_count, char *rfseq, alignment *aln,
	const int32_t ref_s, const int32_t ref_e, msa_p *msa){
		
	// clean_gap_clusters_lable(clusters);
	uint32_t i, j, s, e, median, lqseq_count, offset;
	gap_cluster *clu;
	gap_ *gap;
	sup_aln *s_aln;

	for (i = 0; i < clusters->i; i ++){
		clu = &clusters->clusters[i];
		for (lqseq_count = 0, offset = 20; lqseq_count < LQSEQ_MAX_CAN_COUNT && 
				lqseq_count < clu->i_m * 0.8 && offset < 300; offset += 20){
			s = clu->median > offset ? clu->median - offset : 0;
			e = clu->median + offset;
			for (j = 0; j < clu->i_m; j ++){
				if (clu->gap[j]->l) continue;
				median = (clu->gap[j]->gap.s + clu->gap[j]->gap.e)/2;
				// printf("i:%d j:%d offset:%d median:%d s:%d e:%d median_:%d\n",i,j,offset, clu->median,s,e,median);
				if (median < s || median > e) continue;
				gap = clu->gap[j];
				s_aln = find_sup_alns(sup_alns, gap->s_id, gap->s_s);
				aln->aln_len = aln->shift = 0;
				aln->aln_t_s = s_aln->fs;
				aln->aln_t_e = bam2aln(aln, rfseq, gap->dseq, s_aln->cigar, s_aln->i);
				aln->aln_q_s = s_aln->ds;
				aln->aln_q_e = s_aln->ds;
				// out_align(aln, 0);
				if (aln->aln_t_s < ref_s || aln->aln_t_e > ref_e) clip_aln(aln, ref_s, ref_e, 1);
				get_align_shift(aln, 8, 1);
				if (aln->aln_t_s > aln->aln_t_e - 500) continue;
				aln->aln_t_s -= ref_s;
				aln->aln_t_e -= ref_s;
				get_align_tags(aln, &tags_list[seq_count++], msa);
				gap->l = offset/20;
				gap->s_id = seq_count - 1;
				gap->s_s = aln->aln_q_s;
				lqseq_count ++;
			}
		}
	}
	// for (i = 0; i < clusters->i; i ++){
	// 	clu = &clusters->clusters[i];
	// 	for (j = 0; j < clu->i_m; j ++){
	// 		printf("gap_clusters:update_align_tags: %d %d s:%d e:%d (s+e)/2:%d median:%d, lable:%d\n",i, j, 
	// 			clu->gap[j]->gap.s, clu->gap[j]->gap.e,(clu->gap[j]->gap.s + clu->gap[j]->gap.e)/2, 
	// 			clu->median, clu->gap[j]->l);
	// 	}
	// }
	return seq_count;
}

static int cal_valid_gap(gap_cluster *clu){
	int i, j;
	for (i = j = 0; j < clu->i_m; j ++){
		if (clu->gap[j]->l) i++;
	}
	return i;
}

void generate_gapseqs(gap_clusters *clusters, align_tags_t *tags_list, const int32_t s_){
	uint32_t i, j, s, e, lqseq_count, lqseq_pcount, 
		lqseq_rmcount, offset, offset_step;
	gap_cluster *clu;
	align_tag tag = {0, 0, 0};
	align_tags_t *f_tag, *l_tag;
	for (i = 0; i < clusters->i; i ++){
		clu = &clusters->clusters[i];
		offset = 10;
		lqseq_rmcount = clu->r.s = clu->r.e = 0;
		while (1){
			for (lqseq_pcount = lqseq_count = 0; offset < 30000 && 
					lqseq_pcount < clu->i_m - lqseq_rmcount && 
					(lqseq_count >= lqseq_pcount || lqseq_pcount < clu->i_m/2); offset += 10){
				s = clu->median > offset ? clu->median - offset - s_: 0;
				e = clu->median + offset - s_;
				lqseq_pcount = lqseq_count;
				// printf("median:%d offset:%d s:%d e:%d lqseq_rmcount:%d %d %d %d\n",clu->median, offset, s + s_,e + s_,lqseq_rmcount,clu->r.s+s_, clu->r.e+s_,s_ );
				for (lqseq_rmcount = lqseq_count = j = 0; j < clu->i_m; j ++){
					if (!clu->gap[j]->l) {lqseq_rmcount ++; continue;}
					f_tag = &tags_list[clu->gap[j]->p_id];
					l_tag = &tags_list[clu->gap[j]->s_id];
					if (f_tag->aln_t_s > l_tag->aln_t_s){
						SWAP(clu->gap[j]->p_id, clu->gap[j]->s_id, uint32_t);
						SWAP(clu->gap[j]->p_s, clu->gap[j]->s_s, uint32_t);
						SWAP(f_tag, l_tag, align_tags_t*);
					}
					if (f_tag->aln_t_s < s && f_tag->aln_t_e > s
						&& l_tag->aln_t_s < e && l_tag->aln_t_e > e && s < l_tag->aln_t_s &&
						e > f_tag->aln_t_e) lqseq_count ++;
			// printf("gap_clusters: %d %d s:%d e:%d media:%d lable:%d tag: %d %d %d %d\n",i, j, 
			// 		clu->gap[j]->gap.s, clu->gap[j]->gap.e,(clu->gap[j]->gap.s + clu->gap[j]->gap.e)/2, 
			// 		clu->gap[j]->l, f_tag->aln_t_s, f_tag->aln_t_e, 
			// 		l_tag->aln_t_s, l_tag->aln_t_e);
				}
				if (lqseq_count > lqseq_pcount) {clu->r.s = s;clu->r.e = e;}
			}
			// printf("%d %d %d %d %d %d %d\n",i, clu->i_m, lqseq_count, lqseq_pcount, lqseq_rmcount, clu->r.s+s_, clu->r.e+s_);
			for (offset_step = UINT32_MAX, lqseq_count = j = 0; j < clu->i_m; j ++){
				if (!clu->gap[j]->l) continue;
				f_tag = &tags_list[clu->gap[j]->p_id];
				l_tag = &tags_list[clu->gap[j]->s_id];
				if (f_tag->aln_t_s > clu->r.s || f_tag->aln_t_e < clu->r.s || l_tag->aln_t_s > clu->r.e 
						|| l_tag->aln_t_e < clu->r.e) {
					clu->gap[j]->l = 1;
					continue;
				}
				s = 0;
				e = clu->gap[j]->p_s - 1;
				while (get_align_tag(&s, f_tag, &tag) > 0){
					if (tag.q_base != 4) e++;
					if (tag.t_pos == clu->r.s) break;
				}
				clu->gap[j]->gap.s = e;

				s = 0;
				e = clu->gap[j]->s_s;
				while (get_align_tag(&s, l_tag, &tag) > 0){
					if (tag.t_pos == clu->r.e + 1) break;
					if (tag.q_base != 4) e++;
				}
				clu->gap[j]->gap.e = e;//TODO CHECK
				if (clu->gap[j]->gap.e > clu->gap[j]->gap.s + 10) {
					lqseq_count ++;
					clu->gap[j]->l = 2;
				}else clu->gap[j]->l = 1;
				if (mabs(clu->gap[j]->gap.s, clu->gap[j]->gap.e) < offset_step){
					offset_step = mabs(clu->gap[j]->gap.s, clu->gap[j]->gap.e);
				}
			}
			// printf("i:%d ###lqseq_count: %d lqseq_pcount: %d offset_step: %d\n",i, lqseq_count, lqseq_pcount, offset_step);
			if (lqseq_count >= lqseq_pcount/2 || lqseq_count >= 10) break;
			else offset += offset_step/2 + 20;//TODO CHECK 
		}
	}

	for (i = 0; i < clusters->i; i ++){
		clu = &clusters->clusters[i];
		if (!clu->i_m) continue;
		if (i < clusters->i - 1 && clu->r.e + 500 >= clusters->clusters[i + 1].r.s){
			if (cal_valid_gap(&clusters->clusters[i + 1]) > cal_valid_gap(clu)) {
				clu->i_m = 0;
				continue;
			}else clusters->clusters[i + 1].i_m = 0;
		}
		// printf("gap_cluster1:%d s:%d(%d) e:%d(%d) lqseq_count:%d\n",i,clu->r.s,clu->r.s+s_,
		// 	clu->r.e,clu->r.e+s_, clu->i_m);
		// for (j = 0; j < clu->i_m; j ++){
		// 	if (clu->gap[j]->l != 2) continue;
		// 	printf("gap_clusters1: %d refs:%d refe:%d %d s:%d e:%d lable:%d seq:",i, clu->r.s, clu->r.e,
		// 		j, clu->gap[j]->gap.s, clu->gap[j]->gap.e, clu->gap[j]->l);
		// 	for (s = clu->gap[j]->gap.s; s < clu->gap[j]->gap.e; s++){
		// 		printf("%c", seq_nt16_str[bam_seqi(clu->gap[j]->dseq,s)]);
		// 	}
		// 	printf("\n");
		// }
	}
	// exit(1);
}


void update_split_p(ld_regs *split_ps, gap_clusters *clusters, ld_regs *regs, const int32_t s, const int32_t l, ref_ *ref){
	if (!split_ps->i_m) {
		split_ps->i_m = 50;
		split_ps->reg = malloc(split_ps->i_m * sizeof(pos));
	}

	int i, j, split;
	int ENDING_FLANK = 1000;
	gap_cluster *clu;
	pos *reg;
	for (i = j = 0; i < regs->i; i ++){
		reg = &regs->reg[i];
		if (reg->s < ENDING_FLANK || reg->e + ENDING_FLANK > l) continue;
		// printf("#%u %u %u %d\n",reg->s+s, reg->e+s,l,j );
		j = j > 1 ? j - 1 : 0;
		for (split = 1; j < clusters->i && split; j ++){
			clu = &clusters->clusters[j];
			// printf("%u %u %u %u %d %d\n", reg->s+s, reg->e+s,clu->r.s+s,clu->r.e+s, j, clusters->i);
			if (clu->r.s > reg->e) break;
			else if ((reg->s <= clu->r.s && clu->r.s <= reg->e) || (reg->s <= clu->r.e && clu->r.e <= reg->e) ||
			(clu->r.s <= reg->s && reg->s <= clu->r.e) || (clu->r.s <= reg->e && reg->e <= clu->r.e)) split = 0;
		}
		if (split) {
			if (!split_ps->i || reg->s + s > split_ps->reg[split_ps->i - 1].e + 10000){
				split_ps->reg[split_ps->i].s = reg->s + s;
				split_ps->reg[split_ps->i].e = reg->e + s;
				if (++split_ps->i >= split_ps->i_m) {
					split_ps->i_m += 50;
					split_ps->reg = realloc(split_ps->reg, split_ps->i_m * sizeof(pos));
				}
			}else split_ps->reg[split_ps->i - 1].e = reg->e + s;
		}
        // printf("ld_regs: %d s:%d(%d) e:%d(%d) media:%d len:%d split:%d\n",i,
        //     regs->reg[i].s,regs->reg[i].s+s,regs->reg[i].e, regs->reg[i].e+s,
        //     (regs->reg[i].s + regs->reg[i].e)/2 + s, l * INS_WIN_STEP, split);
    }

    int p;
    uint32_t sco;
    for (i = 0; i < split_ps->i; i++){
    	reg = &split_ps->reg[i];
    	for (sco = j = p = 0; j < ref->qv_l && ref->qv[j].p <= reg->e; j ++){
    		if (ref->qv[j].p >= reg->s){
    			if (sco == 0 || ref->qv[j].ide + ref->qv[j].ort + ref->qv[j].irt < sco){
    				sco = ref->qv[j].ide + ref->qv[j].ort + ref->qv[j].irt;
    				p = j;
    			}
    		}
    	}
    	if (sco && sco < 2900) reg->s = reg->e = ref->qv[p].p;
    	// printf("split_ps: %d start: %d end: %d\n", i, reg->s, reg->e);
    }
}

static consensus_trimed_data *link_consensus_fast(consensuss_data *consensus_t, const int len, const int k){
	int i, j, l, p, s = consensus_t->s/2;
	consensus_data *consensus = NULL, *consensusnext = NULL;
	for (i = consensus_t->i - 1; i; i --){
		consensus = consensus_t->consensus[i];
		consensusnext = consensus_t->consensus[i - 1];
		consensus->rstrip = consensusnext->lstrip = s;
		while (consensus->cns_bases[consensus->len - consensus->rstrip].pos <
			consensus->cns_bases[consensus->len - 1].pos + s) consensus->rstrip++;
		while (consensus->cns_bases[consensus->len - consensus->rstrip].pos >
			consensus->cns_bases[consensus->len - 1].pos + s) consensus->rstrip--;
		while (consensusnext->cns_bases[consensusnext->lstrip].pos <
			consensusnext->cns_bases[0].pos - s) consensusnext->lstrip--;
		while (consensusnext->cns_bases[consensusnext->lstrip].pos >
			consensusnext->cns_bases[0].pos - s) consensusnext->lstrip++;

		l = 0;
		p = consensus->uncorrected_len - consensusnext->uncorrected_len;
		while (l < k){
			j = consensusnext->cns_bases[consensusnext->lstrip].pos - 
				consensus->cns_bases[consensus->len - consensus->rstrip].pos;
			// printf("%d %d %d %d %d %d %d %d\n",j,p,consensus->len - consensus->rstrip,
			// 	consensus->cns_bases[consensus->len - consensus->rstrip].pos,
			// 	consensus->cns_bases[consensus->len - consensus->rstrip].pos+consensus->uncorrected_len,
			// 	consensusnext->lstrip,consensusnext->cns_bases[consensusnext->lstrip].pos,
			// 	consensusnext->cns_bases[consensusnext->lstrip].pos + 
			// 	consensusnext->uncorrected_len);
			if (j == p && consensus->cns_bases[consensus->len - consensus->rstrip].base == 
				consensusnext->cns_bases[consensusnext->lstrip].base) {
				l++;
				consensusnext->lstrip --;
				consensus->rstrip ++;
			}else{
				l = 0;
				if (j >= p) consensusnext->lstrip ++;
				else consensusnext->lstrip --;
			}
		}
	}

	if (consensus_t->i > 1){
		assert (l == k);
		consensus->rstrip -= k;
		consensusnext->lstrip += k;
	}

	consensus_trimed_data *cons_trimed_data = calloc(1, sizeof(consensus_trimed_data));
	cons_trimed_data->i_m = 1;
	cons_trimed_data->data = calloc(1, sizeof(consensus_trimed));
	consensus_trimed *cons_trimed = &cons_trimed_data->data[0];
	int cons_trimed_seq_max_len = len/5 * 6 + 1;
	cons_trimed->seq = malloc(cons_trimed_seq_max_len * sizeof(char));
	for (i = 0; i < consensus_t->i; i ++){
		consensus = consensus_t->consensus[i];
		for (j = consensus->len - consensus->rstrip - 1; j >= (int) consensus->lstrip; j--){
			cons_trimed->seq[cons_trimed->len++] = consensus->cns_bases[j].base;
			if (cons_trimed->len >= cons_trimed_seq_max_len){
				cons_trimed_seq_max_len += 100000;
				cons_trimed->seq = realloc(cons_trimed->seq, cons_trimed_seq_max_len * sizeof(char));
			}
		}
		free (consensus->cns_bases);
		free (consensus);
	}
	cons_trimed->seq[cons_trimed->len] = '\0';
	return cons_trimed_data;
}

static consensus_trimed_data *link_consensus(consensuss_data *consensus_t, ld_regs *split_ps, const int len, const int k, const int split){
	int i, j, l, p, s = consensus_t->s/2;
	consensus_data *consensus = NULL, *consensusnext = NULL;
	for (i = 0; i < consensus_t->i - 1; i ++){
		consensus = consensus_t->consensus[i];
		consensusnext = consensus_t->consensus[i + 1];
		consensus->rstrip = consensusnext->lstrip = s;
		while (consensus->cns_bases[consensus->len - consensus->rstrip].pos <
			consensus->cns_bases[consensus->len - 1].pos - s) consensus->rstrip--;
		while (consensus->cns_bases[consensus->len - consensus->rstrip].pos >
			consensus->cns_bases[consensus->len - 1].pos - s) consensus->rstrip++;
		while (consensusnext->cns_bases[consensusnext->lstrip].pos <
			consensusnext->cns_bases[0].pos + s) consensusnext->lstrip++;
		while (consensusnext->cns_bases[consensusnext->lstrip].pos >
			consensusnext->cns_bases[0].pos + s) consensusnext->lstrip--;

		l = 0;
		p = consensusnext->uncorrected_len - consensus->uncorrected_len;
		while (l < k){
			j = consensus->cns_bases[consensus->len - consensus->rstrip].pos - 
				consensusnext->cns_bases[consensusnext->lstrip].pos;
				// printf("%d %d %d %d %d %d %d %d %d\n",i, j,p,consensus->len - consensus->rstrip,
				// consensus->cns_bases[consensus->len - consensus->rstrip].pos,
				// consensus->cns_bases[consensus->len - consensus->rstrip].pos+consensus->uncorrected_len,
				// consensusnext->lstrip,consensusnext->cns_bases[consensusnext->lstrip].pos,
				// consensusnext->cns_bases[consensusnext->lstrip].pos + 
				// consensusnext->uncorrected_len);
			if (j == p && consensus->cns_bases[consensus->len - consensus->rstrip].base == 
				consensusnext->cns_bases[consensusnext->lstrip].base) {
				l++;
				consensusnext->lstrip --;
				consensus->rstrip ++;
			}else{
				l = 0;
				if (j > p) consensusnext->lstrip ++;
				else if (j < p) consensusnext->lstrip --;
				else{
					int d = consensus->cns_bases[consensus->len - consensus->rstrip].pos + consensus->uncorrected_len - 1;
					while (consensus->cns_bases[consensus->len - consensus->rstrip].pos + consensus->uncorrected_len > d) consensus->rstrip ++;
					while (consensusnext->cns_bases[consensusnext->lstrip].pos + consensusnext->uncorrected_len > d) consensusnext->lstrip --;
				}
			}
		}
	}

	if (consensus_t->i > 1){
		assert (l == k);
		consensus->rstrip -= k;
		consensusnext->lstrip += k;
	}

	// for (i = 0; i < split_ps->i; i ++){
	// 	printf("split_reg: s:%d e:%d media:%d\n", split_ps->reg[i].s, split_ps->reg[i].e, (split_ps->reg[i].s+split_ps->reg[i].e)/2);
	// }

	consensus_trimed_data *cons_trimed_data = malloc(sizeof(consensus_trimed_data));
	cons_trimed_data->i_m = split ? split_ps->i + 1 : 1;
	cons_trimed_data->data = calloc(cons_trimed_data->i_m, sizeof(consensus_trimed));

	int cons_trimed_index = 0, cons_trimed_seq_max_len = len/5 * 6 + 1;
	consensus_trimed *cons_trimed = &cons_trimed_data->data[cons_trimed_index ++];
	cons_trimed->seq = malloc(cons_trimed_seq_max_len * sizeof(char));
	// consensus_trimed *cons_trimed = calloc(1, sizeof(consensus_trimed));
	// int cons_trimed_seq_max_len = len/5 * 6 + 1;
	// cons_trimed->seq = malloc(cons_trimed_seq_max_len * sizeof(char));

	// int N_len = 10;
	l = 0;
	s = l < split_ps->i ? (split_ps->reg[l].s + split_ps->reg[l].e) / 2 : -1;
	l ++;
	for (i = 0; i < consensus_t->i; i ++){
		consensus = consensus_t->consensus[i];
		p = consensus_t->consensus[i]->uncorrected_len;
		for (j = consensus->lstrip; j < consensus->len - consensus->rstrip; j++){
			if (split && consensus->cns_bases[j].pos + p >= s && j >= 1 && consensus->cns_bases[j - 1].pos + p < s ){
				if (split == 1 && cons_trimed->len){
					cons_trimed->seq[cons_trimed->len] = '\0';
					cons_trimed_seq_max_len = len/5 * 6 + 1;
					cons_trimed = &cons_trimed_data->data[cons_trimed_index++];
					cons_trimed->seq = malloc(cons_trimed_seq_max_len * sizeof(char));
				}else if (split == 2) cons_trimed->seq[cons_trimed->len++] = 'N';
				while (consensus->cns_bases[j].pos + p == s) j++;
			}else cons_trimed->seq[cons_trimed->len++] = consensus->cns_bases[j].base;
			
			if (consensus->cns_bases[j].pos + p > s && l < split_ps->i){
				s = (split_ps->reg[l].s + split_ps->reg[l].e) / 2;
				l ++;
			}
			// printf("old: %d new: %d %c\n", consensus->cns_bases[j].pos + p, 
			// 	cons_trimed->len - 1, cons_trimed->seq[cons_trimed->len-1]);

			if (cons_trimed->len >= cons_trimed_seq_max_len - 1){
				cons_trimed_seq_max_len += 100000;
				cons_trimed->seq = realloc(cons_trimed->seq, cons_trimed_seq_max_len * sizeof(char));
			}
		}
		free (consensus->cns_bases);
		free (consensus);
	}
	cons_trimed_data->i_m = cons_trimed_index;
	cons_trimed->seq[cons_trimed->len] = '\0';
	return cons_trimed_data;
}

int cal_rreads_w(pos *rs, int l){
	uint32_t pivot;
	int i, left, s = 0, e = l - 1, k = l/2;
	
	while (1){
		left = s;
		pivot = rs[e].e - rs[e].s;
		for (i = s; i < e; i ++){
			if(rs[i].e - rs[i].s < pivot){
				SWAP(rs[left], rs[i], pos);
				left ++;
			}
		}
		SWAP(rs[left], rs[e], pos);

		if (left == k) {
			pivot = (pivot + 1)/INS_WIN_DIV;
			return pivot > INS_WIN_MIN_SIZE ? pivot : INS_WIN_MIN_SIZE;
		}
		if (left < k) s = left + 1;
		else e = left - 1;
	}
}


int quick_select(int *r, uint32_t s, uint32_t e, uint32_t k){ 
    if (s > e) return e + 1;
    uint32_t i, left, pivot;
    while (1){
        left = s;
        pivot = r[e];
        // printf("%d %d %d\n",s,e,r[s]);
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

int cal_ref_ide(ref_qv *qv, uint32_t l){
	if (l == 0 || qv == NULL) return 0;
	uint32_t i = 0;
	int *t = malloc(sizeof(int) * l);
	for (i = 0; i < l; i ++) t[i] = qv[i].ide;
	int media = quick_select(t, 0, l - 1, l/2);
	free (t);
	return media;
}
int cal_ref_d_ave(const uint16_t *r, const int32_t l, const int clip){
	
	int32_t i;
	uint64_t j = 1, t = 150, h = 0;

	while (j && t/j > h/3){
		h = t/j * 3;
		t = j = 0;
		for (i = clip; i < l - clip; i += 10){
			if (r[i] && r[i] < h) {//TODO CHECK
				t += r[i];
				j ++;
			}
		}
	}
	return j ? t/j : 0;
}

int cal_ref_d(const uint16_t *r, const int32_t l){
	uint32_t i, j, e;
	int media, ignore5 = l > 20000 ? 10000 : l > 200 ? 100 : 20, ignore3 = 0;
	while (!r[ignore5 ++]);
	while (!r[l - 1 - ignore3 ++]);
	int *t = malloc(sizeof(int) * (l - ignore5 - ignore3));
	for (j = e = 0, i = ignore5; i < l - ignore3; i ++, j++) {
		t[j] = r[i];
		if (t[j] < 4) e ++;
	}
	if (!j) media = 0;
	else if (l > 50000 && (double) e / j > 0.2) media = cal_ref_d_ave(r, l, ignore5);
	else media = quick_select(t, 0, j - 1, j/2);
	free (t);
	return media;
}

void inline update_ref_d(uint16_t *r, int w, pos *p, int32_t s){
	uint32_t s_ = p->s > s ? p->s - s : 0;
	uint32_t e_ = p->e - s;
	if (e_ - s_ + 1 >= w * 3){
		s_ = (s_ + w)/INS_WIN_STEP;
		e_ = (e_ - 2 * w)/INS_WIN_STEP;
		while (s_ <= e_) r[s_++] ++;
	}
}

void inline update_ref_ds(uint16_t *r, int w, pos *p, int l, int32_t s){
	int i;
	for (i = 0; i < l; i ++) update_ref_d(r, w, &p[i], s);
}

void out_ref_ds(uint16_t *r, int l, int w, int s){
	int i = 0;
	for (i = 0; i < l; i ++){
		printf("%d %d %d %d\n",i, i * INS_WIN_STEP + s, i * INS_WIN_STEP + w + s, r[i]);
	}
}

typedef struct ctg_cns_cfg {
	pos *rreads;
	uint16_t *ref_ds;
	int reads_type;
	int split;
	float ide_t;
	uint32_t ort_t;
	uint32_t irt_t;
	alignment aln;
	consensuss_data consensus_t;
	sup_alns sup_aln;
	ld_regs ld_cluster;
	gap_clusters g_clusters;
	ld_regs split_ps;
	satags sas;
	gaps gs;
} ctg_cns_cfg;

ctg_cns_cfg *ctg_cns_init(int consensus_w, int reads_type, int split, float ide_t, float ort_t, float irt_t){
	ctg_cns_cfg *cfg = calloc(1, sizeof(ctg_cns_cfg));

	cfg->reads_type = reads_type;
	cfg->split = split;
	cfg->ide_t = ide_t;
	cfg->ort_t = 1000 * ort_t;
	cfg->irt_t = 1000 * irt_t;

	cfg->aln.max_aln_len = 100000;
	cfg->aln.t_aln_str = malloc(cfg->aln.max_aln_len * sizeof(char));
	cfg->aln.q_aln_str = malloc(cfg->aln.max_aln_len * sizeof(char));
	
	cfg->consensus_t.i_m = 5;
	cfg->consensus_t.s = 1000000;
	if (!consensus_w) consensus_w = 40000000;
	else assert (consensus_w > cfg->consensus_t.s * 4);
	cfg->consensus_t.w = consensus_w;
	cfg->consensus_t.consensus = malloc(cfg->consensus_t.i_m * sizeof(consensus_data *));

	cfg->rreads = malloc(INS_RADOM_COUNT * sizeof(pos));
	int w = cfg->consensus_t.w > INS_RADOM_LEN ? cfg->consensus_t.w : INS_RADOM_LEN;
	cfg->ref_ds = malloc((w/INS_WIN_STEP + 200000) * sizeof(uint16_t));

	//TODO init bam_list

	return cfg;
}

void ctg_cns_destroy(ctg_cns_cfg *cfg){
	free (cfg->aln.t_aln_str);
	free (cfg->aln.q_aln_str);
	free (cfg->consensus_t.consensus);
	free (cfg->rreads);
	free (cfg->ref_ds);
    free_gap_info(&cfg->gs);
    free_sup_alns (&cfg->sup_aln);
    free_ld_regs(&cfg->ld_cluster);
    free_ld_regs(&cfg->split_ps);
    free_gap_clusters(&cfg->g_clusters);
    if (cfg->sas.i_m) free (cfg->sas.sa);
	free (cfg);
}

consensus_trimed_data *ctg_cns_core(ctg_cns_cfg *cfg, ref_ *ref, char *bam_list){

	bam_merge_iter bam_iter;
	bam1_t *brecord;

	int32_t j, l, p, s, e, l_qseq;
	uint8_t strand, *satag_;
	satag *sa;
	satags *sas = &cfg->sas;
	pos rfp1, rdp1, rfp2, rdp2;
	uint32_t *cigar;

	uint8_t *rdseqi;
	alignment aln_, *aln = &cfg->aln;

	msa_p *msa;
	align_tags_t *tags_list;
	uint32_t seq_count, seq_count_m;
	// lqseq_max_length = UINT16_MAX;//TODO remove DAG_MAX_LENGTH; DO NOT CHANGE

	consensuss_data *consensus_t =& cfg->consensus_t;

	int brk_g, rreads_i, rreads_w, ref_d, ref_ide = 0;
	sup_alns *sup_aln = &cfg->sup_aln;
	ld_regs *ld_cluster = &cfg->ld_cluster;
	pos *rreads = cfg->rreads;
	uint16_t *ref_ds = cfg->ref_ds;

	gap g;
	gaps *gs = &cfg->gs;
	gap_clusters *g_clusters = &cfg->g_clusters;

	ld_regs *split_ps= &cfg->split_ps;

	int rege, fast = 0, fra_map = 0, total_map = 0;
	
	READS_TYPE = cfg->reads_type;
	if (READS_TYPE != READS_ONT){
		GAP_MIN_LEN = 5;
		GAP_MIN_RATIO1 = 0.3;
	}else{
		GAP_MIN_LEN = 3;
		GAP_MIN_RATIO1 = 0.01;
	}
	MAX_CLIP_RATIO = READS_TYPE == READS_HIFI ? 0.1 : 0.7;
	
	assert (ref->length < INT_MAX);
	char *rfseq = malloc(ref->length + 1);//TODO remove
	bit2seq1(ref->s, ref->length, rfseq);
	int32_t b = cal_win_len(consensus_t->w, consensus_t->s, ref->length);

	brk_g = ref->length > INS_MIN_CHECK_LEN && !fast ? 1 : 0;
	consensus_t->i = s = e = ref_d = rreads_i = rreads_w = split_ps->i = 0;

	char reg[1024];
	if (brk_g) ref_ide = cal_ref_ide(ref->qv, ref->qv_l);
	while (e < ref->length){
		e = s + b > ref->length ? ref->length : s + b;//not include
		l = e - s;
		aln_.shift = aln_.aln_t_s = 0;
		aln_.q_aln_str = aln_.t_aln_str = rfseq + s;
		aln_.aln_q_len = aln_.aln_t_len = l;
		aln_.aln_t_e = aln_.aln_len = l;

		g_clusters->i = ld_cluster->i = gs->i = sup_aln->i = seq_count = p = 0;
		seq_count_m = max(l / 1000, 2);
		msa = calloc(l + 1, sizeof(msa_p));
		tags_list = malloc(seq_count_m * sizeof(align_tags_t));
		assert (msa && tags_list);
		get_align_tags(&aln_, &tags_list[seq_count++], msa);

		if (brk_g) memset(ref_ds, 0, l/INS_WIN_STEP * sizeof(uint16_t));

		rege = s == 0 ? (e > INS_RADOM_LEN ? e : INS_RADOM_LEN) : e;
		sprintf(reg, "%s:%d-%d", ref->n, s, rege);
		bam_merge_iter_init(0, NULL, bam_list, reg, &bam_iter);// TODO FIX HEAD READ
		while (bam_merge_iter_core (&bam_iter) > 0){
			brecord = bam_iter.heap->entry.bam_record;
			assert (brecord->core.pos >= p);
			p = brecord->core.pos;
			if (p >= e) rege = 0;
			cigar = bam_get_cigar(brecord);
			l_qseq = cal_l_qseq(brecord);
			rfp1.s = brecord->core.pos;
			rfp1.e = bam_endpos(brecord);//not include
			rdp1.s = cigarint2ul(cigar, brecord->core.n_cigar, 0);
			rdp1.e = l_qseq - cigarint2ul(cigar, brecord->core.n_cigar, 1);////not include
			// printf("reads:%s %u %u %u %u %u\n",bam_get_qname(brecord), l_qseq, rfp1.s, rfp1.e, rdp1.s,rdp1.e );
			satag_ = bam_aux_get(brecord, "SA");
			g.score = 0;
			if (satag_){
				set_satags(satag_, sas);
				strand = brecord->core.flag & 16 ? 1: 0;
				for (j = 0; j < sas->i; j ++){
					sa = &sas->sa[j];
					// if (bam_name2id(bhead, sa->rname) == tid && sa->strand == strand){
					if (strcmp(sa->rname, ref->n) == 0 && sa->strand == strand){
						rfp2.s = sa->pos;
						rfp2.e = sa->pos + cigarstr2rlen(sa->cigar);
						rdp2.s = cigarstr2ul(sa->cigar, 0);
						rdp2.e = l_qseq - cigarstr2ul(sa->cigar, 1);
						// printf("%u %u %u %u\n",rdp2.s,rdp2.e,l_qseq, cigarstr2ul(sa->cigar, 1));
						check_indel(&g, l_qseq, &rfp1, &rdp1, &rfp2, &rdp2);
					}
				}
				// if (g.score) printf("gap_reads:%s %u %u %u %d %d %d\n",bam_get_qname(brecord), g.gap.s, g.gap.e, g.score, rege, s, e);
			}
			// printf("0:%s %d %d\n", bam_get_qname(brecord), rfp1.s, g.score);
			if (rege && brecord->core.flag & 0xD04 && brk_g && g.score)
				update_sup_alns(sup_aln, rfp1.s, rdp1.s, brecord->core.n_cigar, cigar);
			if (brecord->core.flag & 0xD04) continue;

			// printf("%s %d %d %d %f\n", bam_get_qname(brecord), rdp1.s, rdp1.e, l_qseq, (rdp1.e - rdp1.s)/(double)l_qseq);
			total_map ++;
			if ((rdp1.e - rdp1.s) / (double) l_qseq < 0.7) fra_map ++;
			if ((!g.score) && (rdp1.e - rdp1.s) / (double) l_qseq <= MAX_CLIP_RATIO) continue;

			if (brk_g){
				if (!rreads_w) {
					rreads[rreads_i++] = rfp1;
					if (rreads_i >= INS_RADOM_COUNT) {
						rreads_w = cal_rreads_w(rreads, rreads_i);
						update_ref_ds(ref_ds, rreads_w, rreads, rreads_i, s);
					}
				}else update_ref_d(ref_ds, rreads_w, &rfp1, s);
			}
			if (!rege) continue;

			aln->aln_t_s = rfp1.s;
			aln->aln_t_e = rfp1.e;//not include
			aln->aln_q_s = rdp1.s;
			aln->aln_q_e = rdp1.e;//not include
			aln->aln_len = aln->shift = 0;
			rdseqi = bam_get_seq(brecord);
			l = bam2aln(aln, rfseq, rdseqi, cigar, brecord->core.n_cigar);
			if (l != aln->aln_t_e) {
				fprintf(stderr, "bamaln error, %s\n", ref->n); 
				exit(1);
			}
			if (aln->aln_t_s < s || aln->aln_t_e > e) clip_aln(aln, s, e, g.score);
			get_align_shift(aln, 8, g.score);
			if (aln->aln_t_s > aln->aln_t_e - 500) continue;
			aln->aln_t_s -= s;
			aln->aln_t_e -= s;
			// out_align(aln, bam_get_qname(brecord));
			if ((msa[aln->aln_t_s].coverage > 3000 && msa[aln->aln_t_e].coverage > 3000) || \
				(msa[aln->aln_t_s].coverage > 500 && msa[aln->aln_t_e].coverage > 500 && rdp1.e - rdp1.s < l_qseq * 0.9)) continue;

			get_align_tags(aln, &tags_list[seq_count++], msa);
			if (seq_count >= seq_count_m){
				seq_count_m += 1000;
				tags_list = realloc(tags_list, seq_count_m * sizeof(align_tags_t));
			}
			if (brk_g && g.score && g.gap.s >= s && g.gap.e <= e) update_gap_info(gs, aln, &g, rdseqi, l_qseq, seq_count - 1);
				// printf("2: %s %d\n", bam_get_qname(brecord),rfp1.s);}
		}
		// hts_itr_destroy(iter);
		bam_merge_iter_destroy(&bam_iter);
		if (seq_count < 150 || rreads_i < 150 || sup_aln->i == 0) brk_g = 0;
		// time_debug(clock(), "read bam done...");
		if (brk_g){
			if (!rreads_w) {
				rreads_w = cal_rreads_w(rreads, rreads_i);
				// printf("update_ref_ds init\n");
				update_ref_ds(ref_ds, rreads_w, rreads, rreads_i, s);
			}
			// out_ref_ds(ref_ds, (e - s)/INS_WIN_STEP, rreads_w, s);exit(1);
			// printf("update_ref_ds done\n");
			if (!ref_d) ref_d = cal_ref_d(ref_ds, (e - s)/INS_WIN_STEP);
			// printf("rreads_w:%d ref_d:%d\n",rreads_w, ref_d);
			update_ld_regs(ld_cluster, ref_ds, (e - s)/INS_WIN_STEP, rreads_w, ref_d, s);
			if (ref_ide) update_ld_regs_with_refqv(ld_cluster, ref_ds, ref, rreads_w * INS_WIN_DIV, s, e,
				ref_d * INS_MIN_DEPTH_RATIO_REFQV, ref_ide * cfg->ide_t, cfg->ort_t, cfg->irt_t);
			// printf("update_ld_regs done\n");
			l = update_gap_cluster(gs, g_clusters, ref_ds, rreads_w, ref_d, s);
			// printf("update_gap_cluster done\n");
			if (l + seq_count >= seq_count_m){
				seq_count_m = l + seq_count + 1;
				tags_list = realloc(tags_list, seq_count_m * sizeof(align_tags_t));
			}
			seq_count = update_align_tags(g_clusters, sup_aln, tags_list, seq_count, rfseq, aln, s, e, msa);
			generate_gapseqs(g_clusters, tags_list, s);
			if (ref_d > 15) update_split_p(split_ps, g_clusters, ld_cluster, s, e - s, ref);
		}

		// fprintf(stderr, "gap_aln:%d sup_aln:%d depth_cluster:%d gap_cluster:%d split_count:%d bin_len:%d median_depth:%d fra_map:%.3f\n", 
			// gs->i, sup_aln->i, ld_cluster->i, g_clusters->i, split_ps->i, rreads_w, ref_d, (double) fra_map / (total_map + 1));
		// exit(1);
		consensus_t->consensus[consensus_t->i] = get_cns_from_align_tags(tags_list, msa, 
			seq_count, e - s, 4, 0, fast, g_clusters);
		consensus_t->consensus[consensus_t->i]->uncorrected_len = s;
		if (++consensus_t->i >= consensus_t->i_m) {
			consensus_t->i_m += 5;
			consensus_t->consensus = realloc(consensus_t->consensus, consensus_t->i_m * sizeof(consensus_data *));
		}
		s = e - consensus_t->s;
	}
	free (rfseq);
	
	double fra = (double) fra_map / (total_map + 1);
	if (READS_TYPE == READS_HIFI && fra > 0.1 ){
		fprintf(stderr, "Warning, Too many (%.3f%%) fragment mappings in %s, please polish the genome with other reads first, or"
			" adjust the mapping parameters to tolerate more errors, such as use asm20/map-pb instead of asm5 for minimap2,"
			" continue anyway...\n",
			(double) fra_map * 100 / (total_map + 1), ref->n);
		// if (total_map > 100 && fra > 0.5){
		// 	free(msa);
		// 	free (rfseq);
		// 	for (p = 0; p < seq_count; p++) free (tags_list[p].align_tags);
		// 	free (tags_list);
		// 	consensus_trimed_data *cons_trimed_data = calloc(1, sizeof(consensus_trimed_data));
		// 	cons_trimed_data->i_m = 1;
		// 	cons_trimed_data->data = calloc(cons_trimed_data->i_m, sizeof(consensus_trimed));
		// 	consensus_trimed *cons_trimed = &cons_trimed_data->data[0];
		// 	cons_trimed->seq = calloc(2, sizeof(char));
		// 	cons_trimed->len = cons_trimed->identity = 1;
		// 	return cons_trimed_data;
		// }
	}

	consensus_trimed_data *cons_trimed_data;
	if (fast) cons_trimed_data = link_consensus_fast(consensus_t, ref->length, 50);
	else cons_trimed_data = link_consensus(consensus_t, split_ps, ref->length, 50, cfg->split);
	return cons_trimed_data;
}

int main(int argc, char *argv[]){
	// start = clock();
	// time_i = 0;

	char *reff = argv[1];
	char *bam_list = argv[2];
	refs_ *refs = read_ref(reff, NULL, 0);

	int i, j;
	ctg_cns_cfg *cfg = ctg_cns_init(5000000, READS_ONT, 0, 0.8, 0.8, 0.8);
	//ctg_cns_cfg *cfg = ctg_cns_init(5000000, 3, 0, 0.8, 0.8, 0.8);
	for (i = 0; i < refs->i; i++){
		assert (refs->ref[i].length < INT_MAX);
		consensus_trimed_data *cons_trimed_data = ctg_cns_core(cfg, &refs->ref[i], bam_list);
		if (cons_trimed_data->i_m > 1){
			for (j = 0; j < cons_trimed_data->i_m; j++){
				printf(">%s_%d_lgs %d %f\n%s\n", refs->ref[i].n, j, cons_trimed_data->data[j].len, 
					cons_trimed_data->data[j].identity, cons_trimed_data->data[j].seq);
			}
		}else{
			printf(">%s_lgs %d %f\n%s\n", refs->ref[i].n, cons_trimed_data->data[0].len, 
				cons_trimed_data->data[0].identity, cons_trimed_data->data[0].seq);
		}
		free_consensus_trimed_data(cons_trimed_data);
	}

	ctg_cns_destroy(cfg);
	refs_destroy(refs);
}
