#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h> //getpid
#include "nextcorrect.h"

#ifndef SELF_CORRECT
	#define SELF_CORRECT
#endif


#define MAX_PP_PPP 0
#define READS_ONT 1
#define READS_CLR 2
#define READS_HIFI 3

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

unsigned int lqseq_max_length;
static int READS_TYPE = READS_ONT;

static uint8_t int_to_base[] = { 
	65, 84, 71, 67, 45, 78, 77
};

static align_tag align_tag_head = {
	.t_pos = -1,  
	.q_base = 0, 
	.delta = 0
};

static uint8_t base_to_int[] = {
	// A 0 T 1 G 2 C 3
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//0-15
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//16-31
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//32-47
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//48-63
	4, 0, 4, 3,  4, 4, 4, 2,  4, 4, 4, 4,  4, 6, 5, 4,//64-79
	4, 4, 4, 4,  1, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,//80-95
	4, 0, 4, 3,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,//96-111
	4, 4, 4, 4,  1, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4//112-127
};

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

static void string_lower(char *p) {
	for ( ; *p; ++p) *p = tolower(*p);
}

// static void out_align (alignment *aln){
// 	printf("%d %d %d ", aln->aln_t_s, aln->aln_t_e, aln->aln_len);
// 	for (int i = 0; i < aln->aln_len; i++){
// 		printf("%c", aln->t_aln_str[i + aln->shift]);
// 	}
// 	printf("\n");
// 	printf("%d %d %d ", aln->aln_t_s, aln->aln_t_e, aln->aln_len);
// 	for (int i = 0; i < aln->aln_len; i++){
// 		printf("%c", aln->q_aln_str[i + aln->shift]);
// 	}
// 	printf("\n");
// }

static float get_align_ide(alignment *aln){
	int m = 0;
	for (int i = 0; i < aln->aln_len; i++){
		if (aln->q_aln_str[i] == aln->t_aln_str[i]) m ++;
	}

	return aln->aln_len == 0 ? 0 : (float) m / aln->aln_len;
}

static inline void get_align_shift(alignment *aln, int k){

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

		if (j == k){
			aln->aln_t_s -= k;
			aln->shift = i - k + 1;
			aln->aln_len = aln->aln_len - i + k - 1;
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

			if (j == k){
				aln->aln_t_e += k;
				aln->aln_len = aln->aln_len - t + k - 1;
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

	for (p = 0; p < len; p++){
		max_size = msa[p].max_size;
		msa[p].d_b = calloc(1, max_size * sizeof(msa_p_d_b *) + \
			max_size * 6 * sizeof(msa_p_d_b) + \
			max_size * 6 * MAX_PP_PPP * sizeof(msa_p_d_b_pp_ppp)
		);

		msa_p_d_b * const _ = (msa_p_d_b *) (msa[p].d_b + max_size);
		msa_p_d_b_pp_ppp * const __ = (msa_p_d_b_pp_ppp *) (_ + max_size * 6);

		for (d = 0; d < max_size; d++){
			msa[p].d_b[d] = _ + d * 6;
			for (b = 0; b < 6; b++){
				msa[p].d_b[d][b].max_size = MAX_PP_PPP;
				msa[p].d_b[d][b].len = 0;
				msa[p].d_b[d][b].pp_ppp = __ + (d * 6 + b) * msa[p].d_b[d][b].max_size;
			}
		}
	}
}

static inline void free_msa(msa_p *msa, int len){
	for (int p = 0; p < len; p ++){
		for (int d = 0; d < msa[p].max_size; d++){
			for (int b = 0; b < 6; b++){
				if (msa[p].d_b[d][b].max_size > MAX_PP_PPP) free (msa[p].d_b[d][b].pp_ppp);
			}
		}
		free(msa[p].d_b);
	}
	free (msa);
}

static inline void update_msa(msa_p * msa, align_tags_t *tags_list, int aligned_seq_count, int lable){
	int i, d, p, updated;
	align_tag *p1, *pp, *ppp;
	for (p = 0; p < aligned_seq_count; p++){
		for (d = 0; d < tags_list[p].len; d++){
			p1 = &tags_list[p].align_tags[d];
			pp = d > 0 ? &tags_list[p].align_tags[d - 1] : &align_tag_head;
			ppp = d > 1 ? &tags_list[p].align_tags[d - 2] : &align_tag_head;

			if (p1->q_base == 6 || pp->q_base == 6) continue;
			msa_p_d_b *p_base = &msa[p1->t_pos].d_b[p1->delta][p1->q_base];
			for (i = 0, updated = 0; i < p_base->len; i++){
				if (p_base->pp_ppp[i].pp.t_pos == pp->t_pos &&
					p_base->pp_ppp[i].pp.delta == pp->delta &&
					p_base->pp_ppp[i].pp.q_base == pp->q_base &&
					p_base->pp_ppp[i].ppp.t_pos == ppp->t_pos &&
					p_base->pp_ppp[i].ppp.delta == ppp->delta && 
					p_base->pp_ppp[i].ppp.q_base == ppp->q_base){
					p_base->pp_ppp[i].link_count ++;
					updated = 1;
					break;
				}
			}
			if (! updated){
				if (i >= p_base->max_size){
					reallocate_msa_p_d_b_pp_ppp_mem(p_base, 5);
				}
				p_base->pp_ppp[p_base->len].pp = *pp;
				p_base->pp_ppp[p_base->len].ppp = *ppp;
				p_base->pp_ppp[p_base->len].link_count = 1;
				p_base->pp_ppp[p_base->len].score = 0;
				// p_base->p_base[i] = pp_ppps->len++;	
				p_base->len ++;
			}
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


static consensus_trimed *error_seed(int len){
	consensus_trimed *cons_trimed = malloc(sizeof(consensus_trimed));
	cons_trimed->len = len;
	cons_trimed->seq = malloc(cons_trimed->len * sizeof(char));
	return cons_trimed;
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

static int compare_seq_by_len(const void * a, const void * b){
	struct seq_ *lqseq2 = (struct seq_ *) a;
	struct seq_ *lqseq1 = (struct seq_ *) b;
	return (lqseq1->len < lqseq2->len) - (lqseq1->len > lqseq2->len);
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

static int generate_lqseqs_from_tags(lqseq *lqseqs, const int lqseqs_count, align_tags_t *tags_list, \
	const unsigned int aligned_seq_count, const int split){

	int i, j, k, lable, large_seq, start, end, index, max_aln_lqseq_len;

	int max_aln_length = 0;
	// int aln_m, aln_l;

	align_tags_t *tags;
	for (i = 0; i < lqseqs_count; i++){
		max_aln_lqseq_len = large_seq = 0;
		start = lqseqs[i].start;
		end = lqseqs[i].end;
		lqseqs[i].len = 0;
		// lqseqs[i].seqs = malloc(aligned_seq_count * sizeof(struct seq_));
		lqseqs[i].seqs = malloc(LQSEQ_MAX_CAN_COUNT * sizeof(struct seq_));
		lqseqs[i].sudoseed = NULL;
		for (j = 0; j < aligned_seq_count; j++){ // select seed reads as template
			// aln_m = aln_l = 0;
			index = 0;
			tags = &tags_list[j];
			if (tags->align_tags[0].t_pos <= start && tags->align_tags[tags->len - 1].t_pos >= end){
				lable = 0;
				for (k = start - tags->align_tags[0].t_pos; k < tags->len && tags->align_tags[k].t_pos <= end; k++){
					if (tags->align_tags[k].t_pos >= start && tags->align_tags[k].q_base != 4){
						lqseqs[i].seqs[lqseqs[i].len].seq[index++] = int_to_base[tags->align_tags[k].q_base];
						if (index > lqseq_max_length - 1) {
							large_seq++;
							lable = 1;
							break;
						}
					}
					// if (tags->align_tags[k].t_pos >= start){
					// 	aln_l ++;
					// 	if (tags->align_tags[k].q_base != 4 && tags->align_tags[k].delta == 0) aln_m ++;
					// }
				}
				if (index && ! lable){
					lqseqs[i].seqs[lqseqs[i].len].seq[index] = '\0';
					lqseqs[i].seqs[lqseqs[i].len].len = index;
					lqseqs[i].seqs[lqseqs[i].len].order = lqseqs[i].len;
					// lqseqs[i].seqs[lqseqs[i].len].ide = (float) aln_m/aln_l;
					lqseqs[i].len ++;
					if (index > max_aln_lqseq_len) max_aln_lqseq_len = index;
				}
				// if (lqseqs[i].len >= 60) break;
				if (lqseqs[i].len >= LQSEQ_MAX_CAN_COUNT) break;
			}
		}

		if ((float) large_seq/(lqseqs[i].len + large_seq) > 1.0/3 || lqseqs[i].len <= 4 || (split && lqseqs[i].len < 10)) {
			lqseqs[i].len = 0;
		}else{

			uint16_t kmers[KMER_LEN_COUNT];
			count_kmers(&lqseqs[i], kmers, 1, 0);
			count_kscore(&lqseqs[i], kmers, 0);
			qsort(lqseqs[i].seqs, lqseqs[i].len, sizeof(struct seq_), compare_seq_by_kscore);
			count_kmers(&lqseqs[i], kmers, KMER_MAX_SEQ, 0);
			count_kscore(&lqseqs[i], kmers, 0);
			unsigned int klastscore, kmaxscore = lqseqs[i].seqs[0].kscore;
			unsigned int kmaxlen = lqseqs[i].seqs[0].len, kminlen;
			if (kmaxlen > 500 || (kmaxlen > 200 && kmaxscore < 200)){//LQSEQ_MAX_REV_LEN){
				uint16_t score[LQSEQ_MAX_CAN_COUNT];
				if (lqseqs[i].seqs[0].order) find_ref_lqseq(&lqseqs[i]);
				save_kscore(&lqseqs[i], score);
				count_kmers(&lqseqs[i], kmers, 1, 1);
				count_kscore(&lqseqs[i], kmers, 1);
				qsort(lqseqs[i].seqs, lqseqs[i].len, sizeof(struct seq_), compare_seq_by_kscore);
				count_kmers(&lqseqs[i], kmers, KMER_MAX_SEQ, 1);
				count_kscore(&lqseqs[i], kmers, 1);
				merge_kscore(&lqseqs[i], score);
			}

			qsort(lqseqs[i].seqs, lqseqs[i].len, sizeof(struct seq_), compare_seq_by_kscore);
			kminlen = kmaxlen = lqseqs[i].seqs[0].len;
			klastscore = kmaxscore = lqseqs[i].seqs[0].kscore;

			for (k = j = 0; j < lqseqs[i].len; j ++){
				if (lqseqs[i].seqs[j].kscore * 10 < kmaxscore || j >= LQSEQ_MAX_COUNT || lqseqs[i].seqs[j].kscore * 2 < klastscore || 
					(j > 4 && kmaxlen > 200 && lqseqs[i].seqs[j].kscore < kmaxscore * 0.6 && lqseqs[i].seqs[j].len < kminlen * 0.8)) break;
				klastscore = lqseqs[i].seqs[j].kscore;
				if (j < KMER_MAX_SEQ && lqseqs[i].seqs[j].kscore > kmaxscore * 0.8){
					if (lqseqs[i].seqs[j].len > kmaxlen) kmaxlen = lqseqs[i].seqs[j].len;
					else if (lqseqs[i].seqs[j].len < kminlen) kminlen = lqseqs[i].seqs[j].len;
				}
			}

			lqseqs[i].indexs = 0;
			lqseqs[i].indexe = kmaxlen > LQSEQ_MAX_REV_LEN && j > 6 ? 5 : j - 1;
			if (lqseqs[i].indexe - lqseqs[i].indexs <= 3) {
				lqseqs[i].len = 0;
				continue;
			}
			// if (lqseqs[i].seqs[0].kscore < (min(KMER_RANGE, kmaxlen) - KMER_LEN) * KMER_MAX_SEQ/2){
				if (lqseqs[i].seqs[0].len < 3000){
					j = lqseqs[i].indexs;
					k = j + 6 < lqseqs[i].indexe ? 6 : lqseqs[i].indexe - j + 1;
				}else{
					j = lqseqs[i].indexs;
					k = j + 2 < lqseqs[i].indexe ? 2 : lqseqs[i].indexe - j + 1;
				}
				lqseqs[i].sudoseed = poa_to_consensus(&lqseqs[i].seqs[j], k);
			// }else{
			// 	lqseqs[i].sudoseed = strdup(lqseqs[i].seqs[k].seq);
			// }
			lqseqs[i].sudoseed_len = strlen(lqseqs[i].sudoseed);

			if (lqseqs[i].sudoseed_len > 500){
				int kmaxlen, kminlen;
				k = kmaxlen = kminlen = lqseqs[i].seqs[lqseqs[i].indexs].len;
				for (j = lqseqs[i].indexs + 1; j <= lqseqs[i].indexe && j <= lqseqs[i].indexs + 4; j ++){
					k += lqseqs[i].seqs[j].len;
					if (lqseqs[i].seqs[j].len > kmaxlen) kmaxlen = lqseqs[i].seqs[j].len;
					else if (lqseqs[i].seqs[j].len < kminlen) kminlen = lqseqs[i].seqs[j].len;
				}
				k = kmaxlen != kminlen ?  (k - kmaxlen - kminlen)/ (j - lqseqs[i].indexs - 2) : k / (j - lqseqs[i].indexs);
				if (lqseqs[i].sudoseed_len > k + k / 10){
					// printf("i:%d %d %d %d %f\n",i,lqseqs[i].sudoseed_len,lqseqs[i].seqs[j].len,k,(float)k/lqseqs[i].sudoseed_len);
					for (kminlen = k, k = lqseqs[i].indexs; k < j; k ++){
						if (lqseqs[i].seqs[k].len != kmaxlen && lqseqs[i].seqs[k].len >= kminlen) break;
					}
					free (lqseqs[i].sudoseed);
					if (j == k) for (k = 0; k < lqseqs[i].len && lqseqs[i].seqs[k].order; k ++);
					lqseqs[i].sudoseed = strdup(lqseqs[i].seqs[k].seq);
					lqseqs[i].sudoseed_len = lqseqs[i].seqs[k].len;
				}
			}

			if (max_aln_lqseq_len + lqseqs[i].sudoseed_len > max_aln_length) {
				max_aln_length = max_aln_lqseq_len + lqseqs[i].sudoseed_len;
			}
		}
	}

	// for (i = 0 ; i< lqseqs_count; i++){
	// 	if (lqseqs[i].len){
	// 		for (int k = 0; k < lqseqs[i].len; k++){
	// 			j = k >= lqseqs[i].indexs && k <= lqseqs[i].indexe ? 1 : 0;
	// 			printf("%d %d %d %d s:%d e:%d order:%d kscore:%d len:%d %s\n",lqseqs_count, i, k, j, lqseqs[i].start, 
	// 				lqseqs[i].end, lqseqs[i].seqs[k].order, lqseqs[i].seqs[k].kscore, strlen(lqseqs[i].seqs[k].seq), lqseqs[i].seqs[k].seq);
	// 		}
	// 		printf("seed%d: %d %s\n",i, lqseqs[i].sudoseed_len, lqseqs[i].sudoseed);
	// 	}
	// }
	// exit(1);
	
	for (i = 0; i < aligned_seq_count; i++){
		free (tags_list[i].align_tags);
	}
	free (tags_list);
	return max_aln_length;
	
	// time_debug(clock(), "lq produced...");
}

static int remove_differ_len_lqseq(lqseq *lqseq){
	int s, j, k = lqseq->end - lqseq->start + 1;
	int offset = min(max(30, k/10), k/3);
	int8_t dif_len[LQSEQ_MAX_CAN_COUNT] = {0};
	for (j = s = 0; j < lqseq->len; j ++){
		if (lqseq->seqs[j].len + offset >= k && lqseq->seqs[j].len <= k + offset) {
			s ++;
		}else{
			dif_len[j] = 1;
		}
	}

	if (s != lqseq->len && (s >= lqseq->len/2 || (s >= lqseq->len/3 && s >= 3))){
		k = lqseq->len;
		for (j = 0; j < lqseq->len && j < k; j ++){ 
			if (dif_len[j]) {
				for (k --; k > j; k--){
					if (!dif_len[k]) {
						SWAP(lqseq->seqs[j], lqseq->seqs[k], struct seq_);
						break;
					}
				}
			}
		}
		lqseq->len = k;
	}
	return s;
}

typedef struct {
	uint16_t d;// diff
	uint16_t s;// same
	int del;
} phs;
static void remove_differ_phase_lqseq1(lqseq *lqseq, phs *phase_sco){
	
	int i, j, k;
	int8_t dif[LQSEQ_MAX_CAN_COUNT] = {0};
	for (j = 0; j < lqseq->len; j ++){
		i = lqseq->seqs[j].order;
		if (phase_sco[i].s < phase_sco[i].d) dif[j] = 1;//NB:can not <=
		else if (phase_sco[i].d >= 10) dif[j] = 1;//NB:for >=3 kmer-types such as polyploid or mult-repeats
	}

	for (k = lqseq->len, j = 0; j < lqseq->len && j < k; j ++){
		if (dif[j]) {
			for (k --; k > j; k--){
				if (!dif[k]) {
					SWAP(lqseq->seqs[j], lqseq->seqs[k], struct seq_);
					break;
				}
			}
		}
	}
	lqseq->len = k;
}


static void mark_del_lqseq(lqseq *lqseq, phs *phase_sco){
	
	int i, j, k;
	for (k = 0, j = 1; j < lqseq->len; j ++){//some diff sites are wrong
		i = lqseq->seqs[j].order;
		if (phase_sco[i].s >= 3 && !phase_sco[i].d) k ++;
	}
	if (k >= 2){
		for (j = 0; j < lqseq->len; j ++){
			i = lqseq->seqs[j].order;
			if (phase_sco[i].d) phase_sco[i].del = 1;
		}
	}else{
		for (j = 0; j < lqseq->len; j ++){
			i = lqseq->seqs[j].order;
			if (phase_sco[i].s < phase_sco[i].d || phase_sco[i].d >= 3) phase_sco[i].del = 1;
		}
	}
}

static void remove_differ_phase_lqseq(lqseq *lqseq, phs *phase_sco){
	
	int i, j, k;
	int8_t dif[LQSEQ_MAX_CAN_COUNT] = {0};
	for (j = 0; j < lqseq->len; j ++){
		i = lqseq->seqs[j].order;
		if (phase_sco[i].del) dif[j] = 1;
	}

	for (k = lqseq->len, j = 0; j < lqseq->len && j < k; j ++){
		if (dif[j]) {
			for (k --; k > j; k--){
				if (!dif[k]) {
					SWAP(lqseq->seqs[j], lqseq->seqs[k], struct seq_);
					break;
				}
			}
		}
	}
	lqseq->len = k;
}

static int select_most_lqseq(lqseq *lqseq, int len){
	int j, s, k;
	int8_t used[LQSEQ_MAX_CAN_COUNT] = {0};
	for (j = s = 0; j < min(lqseq->len, len); j ++){
		lqseq->seqs[j].kscore = 1;
		if (used[j]) continue;
		for (k = j + 1; k < lqseq->len; k ++){
			if (lqseq->seqs[j].len == lqseq->seqs[k].len && strcmp(lqseq->seqs[j].seq, lqseq->seqs[k].seq) == 0){
				used[k] = 1;
				lqseq->seqs[j].kscore ++;
			}
		}
	
		if (lqseq->seqs[j].kscore > lqseq->seqs[s].kscore || (lqseq->seqs[s].order > lqseq->seqs[j].order &&
			lqseq->seqs[j].kscore == lqseq->seqs[s].kscore)) s = j;
			// lqseq->seqs[j].kscore == lqseq->seqs[s].kscore && lqseq->seqs[j].len > lqseq->seqs[s].len)) s = j;
	}
	return s;
}

static void select_most2_lqseq(lqseq *lqseq, int len, int *m1, int *m2){
	//*m1 == *m2 means only one lqseq type
	int j, k;
	int8_t used[LQSEQ_MAX_CAN_COUNT] = {0};
	for (*m1 = *m2 = j = 0; j < min(lqseq->len, len); j ++){
		lqseq->seqs[j].kscore = 1;
		if (used[j]) continue;
		for (k = j + 1; k < lqseq->len; k ++){
			if (lqseq->seqs[j].len == lqseq->seqs[k].len && strcmp(lqseq->seqs[j].seq, lqseq->seqs[k].seq) == 0){
				used[k] = 1;
				lqseq->seqs[j].kscore ++;
			}
		}
		if (lqseq->seqs[j].kscore > lqseq->seqs[*m1].kscore || (lqseq->seqs[j].kscore == lqseq->seqs[*m1].kscore && 
			lqseq->seqs[j].order < lqseq->seqs[*m1].order)){
			*m2 = *m1;
			*m1 = j;
		}else if (*m2 == *m1 || lqseq->seqs[j].kscore > lqseq->seqs[*m2].kscore){
			*m2 = j;
		}
	}
}

//For kscore already be counted
static void select_most2_lqseq_with_kscore(lqseq *lqseq, int len, int *m1, int *m2){
	//*m1 == *m2 means only one lqseq type
	int j;
	for (*m1 = *m2 = j = 0; j < min(lqseq->len, len); j ++){
		if (lqseq->seqs[j].kscore > lqseq->seqs[*m1].kscore || (lqseq->seqs[j].kscore == lqseq->seqs[*m1].kscore && 
			lqseq->seqs[j].order < lqseq->seqs[*m1].order)){
			*m2 = *m1;
			*m1 = j;
		}else if (*m2 == *m1 || lqseq->seqs[j].kscore > lqseq->seqs[*m2].kscore){
			*m2 = j;
		}
	}
}

// static int select_longest_lqseq(lqseq *lqseq, int len){
// 	int j, l = 0;
// 	for (j = 1; j < min(lqseq->len, len); j ++){
// 		if (lqseq->seqs[j].len > lqseq->seqs[l].len) l = j;
// 	}
// 	return l;
// }

static void _set_s_e(char *seq, int seq_len, int *s, int *e){
	*s = 0;
	while (*s + 1 < seq_len && seq[*s] == seq[*s + 1]) (*s) ++;
	*e = seq_len - 1;
	while (*e > 0 && seq[*e - 1] == seq[*e]) (*e) --;
}
static int homo_end_compress_is_same(char *seq1, int seq1_len, char *seq2, int seq2_len){
	//remove heter: TCAAAAA TCAAAA 
	int s1_s, s1_e, s2_s, s2_e;
	_set_s_e(seq1, seq1_len, &s1_s, &s1_e);
	// printf("%d %d\n", s1_s, s1_e);
	_set_s_e(seq2, seq2_len, &s2_s, &s2_e);
	if (s1_e <= s1_s && s2_e <= s2_s) return 1;
	if (s1_e - s1_s != s2_e - s2_s) return 0;
	int i;
	for (i = 0; i <= s1_e - s1_s; i ++){
		if (seq1[i + s1_s] != seq2[i + s2_s]) return 0;
	}
	return 1;
}

static int prefixhomo_compress_is_same(char *seq1, int seq1_len, char *seq2, int seq2_len){

	int i = 0, j = 0;

	while (i < seq1_len && j < seq2_len){
		// printf("seq1:%s seq2:%s %d %d %c %c %d %d\n",seq1, seq2, i, j, seq1[i], seq2[j],pi,pj);
		if (seq1[i] != seq2[j]) return 0;
		while (i + 1 < seq1_len && seq1[i] == seq1[i + 1]) i ++;
		while (j + 1 < seq2_len && seq2[j] == seq2[j + 1]) j ++;

		i ++;
		j ++;
	}

	return 1;
}

static int trim_endssr_is_same(char *seq1, int seq1_len, char *seq2, int seq2_len){
	//remove heter: TCTTTG TCTTTGTG;CATG CATGTG
	if (seq1_len < seq2_len){
		char *tmp = seq1;
		int tmp_len = seq1_len;
		seq1 = seq2;
		seq1_len = seq2_len;
		seq2 = tmp;
		seq2_len = tmp_len;
	}

	int i;
	for (i = 0; i < seq2_len; i ++){ 
		if (seq1[i] != seq2[i]) return 0;
	}

	int j;
	for (j = seq1_len - 1; j >= i; j--){
		if (seq1[j] != seq2[seq2_len - (seq1_len - j)]) return 0;
	}
	return 1;
}


static int generate_lqseqs_from_tags_kmer(lqseq *lqseqs, const int lqseqs_count, align_tags_t *tags_list, \
	const unsigned int aligned_seq_count, const int split){

	int i, j, k, s, lable, start, end, index, max_aln_lqseq_len;
	int max_aln_length = max_aln_lqseq_len = 0;

	lqseq *lqseq;
	align_tags_t *tags;
	assert (aligned_seq_count < 65536);
	
	for (i = 0; i < lqseqs_count; i++){
		lqseq = &lqseqs[i];
		start = lqseq->start;
		end = lqseq->end;
		lqseq->len = 0;
		lqseq->seqs = malloc(LQSEQ_MAX_CAN_COUNT * sizeof(struct seq_));
		lqseq->sudoseed = NULL;
		for (j = 0; j < aligned_seq_count; j++){ // select seed reads as template
			index = 0;
			tags = &tags_list[j];
			if (tags->align_tags[0].t_pos <= start && tags->align_tags[tags->len - 1].t_pos >= end){
				lable = 0;
				for (k = start - tags->align_tags[0].t_pos; k < tags->len && tags->align_tags[k].t_pos <= end; k++){
					if (tags->align_tags[k].t_pos >= start && tags->align_tags[k].q_base != 4){
						lqseq->seqs[lqseq->len].seq[index++] = int_to_base[tags->align_tags[k].q_base];
						if ((j && index > lqseq_max_length - 1) || index > DAG_MAX_LENGTH - 1) {
							lable = 1;
							break;
						}
					}
				}
				if (index && ! lable){
					lqseq->seqs[lqseq->len].seq[index] = '\0';
					lqseq->seqs[lqseq->len].len = index;
					lqseq->seqs[lqseq->len].order = j;//seq source
					lqseq->seqs[lqseq->len].kscore = 0;
					lqseq->len ++;
					// phase_sco[j].t ++;
					if (index > max_aln_lqseq_len) max_aln_lqseq_len = index;
					// if (lqseq->len == 1 || strcmp(lqseq->seqs[0].seq, lqseq->seqs[lqseq->len - 1].seq) == 0) phase_sco[j].c ++;
				}
				if (lqseqs[i].len >= LQSEQ_MAX_CAN_COUNT) break;
				// if (lqseq->len >= 20) break;
			}
		}
	}

	phs *phase_sco = calloc(aligned_seq_count, sizeof(phs));//phased score of each reads
	int has_heter = 0;
	for (i = 0; i < lqseqs_count; i++){//select heter. sites, first.
		lqseq = &lqseqs[i];
		if (!lqseq->len) continue;
		select_most2_lqseq(lqseq, lqseq->len, &s, &k);
		if (s != k && lqseq->seqs[k].kscore >= 3 && lqseq->seqs[s].len == lqseq->seqs[k].len){//only consider SNPs
			if (s == 0 || k == 0){
				int heter = s == 0 ? k : s;
				for (j = 0; j < lqseq->len; j ++){
					index = lqseq->seqs[j].order;
					if (lqseq->seqs[0].len == lqseq->seqs[j].len && strcmp(lqseq->seqs[0].seq, lqseq->seqs[j].seq) == 0) phase_sco[index].s ++;
					else if (lqseq->seqs[heter].len == lqseq->seqs[j].len && strcmp(lqseq->seqs[heter].seq, lqseq->seqs[j].seq) == 0) phase_sco[index].d ++;
				}
			}
			lqseq->indexs = 1;// heter flag
		}else lqseq->indexs = 0;// homo flag
		
		if (!has_heter && (lqseq->indexs == 1 || (s != k && lqseq->seqs[k].kscore >= 5 && 
			lqseq->seqs[s].kscore + lqseq->seqs[k].kscore >= lqseq->len * 0.8 && 
			!prefixhomo_compress_is_same(lqseq->seqs[s].seq, lqseq->seqs[s].len, lqseq->seqs[k].seq, lqseq->seqs[k].len)))){//candidate heter. site.
			has_heter = 1;
		}
	}

	if (has_heter && !phase_sco[0].s){//no SNPs were detected, but has_heter
		for (i = 0; i < lqseqs_count; i++){
			lqseq = &lqseqs[i];
			if (!lqseq->len) continue;
			select_most2_lqseq_with_kscore(lqseq, lqseq->len, &s, &k);
			if (s != k && lqseq->seqs[k].kscore >= 5 && (lqseq->seqs[s].kscore + lqseq->seqs[k].kscore) >= lqseq->len * 0.8 && 
				(lqseq->seqs[s].len >= lqseq->seqs[k].len + 5 || lqseq->seqs[k].len >= lqseq->seqs[s].len + 5 || 
				!prefixhomo_compress_is_same(lqseq->seqs[s].seq, lqseq->seqs[s].len, lqseq->seqs[k].seq, lqseq->seqs[k].len))){//candidate heter. site.
				
				int s_, k_;
				if (s == 0) {
					s_ = 1;
					k_ = 0;
				}else if (k == 0){
					s_ = 0;
					k_ = 1;
				}else{
					s_ = (homo_end_compress_is_same(lqseq->seqs[s].seq, lqseq->seqs[s].len, lqseq->seqs[0].seq, lqseq->seqs[0].len) 
						|| trim_endssr_is_same(lqseq->seqs[s].seq, lqseq->seqs[s].len, lqseq->seqs[0].seq, lqseq->seqs[0].len)
						|| prefixhomo_compress_is_same(lqseq->seqs[s].seq, lqseq->seqs[s].len, lqseq->seqs[0].seq, lqseq->seqs[0].len));
					k_ = (homo_end_compress_is_same(lqseq->seqs[k].seq, lqseq->seqs[k].len, lqseq->seqs[0].seq, lqseq->seqs[0].len) 
						|| trim_endssr_is_same(lqseq->seqs[k].seq, lqseq->seqs[k].len, lqseq->seqs[0].seq, lqseq->seqs[0].len) ||
						prefixhomo_compress_is_same(lqseq->seqs[k].seq, lqseq->seqs[k].len, lqseq->seqs[0].seq, lqseq->seqs[0].len));
				}

				int same, heter;
				if (s_ && !k_) {
					same = s;
					heter = k;
				}else if (k_ && !s_){
					same = k;
					heter = s;
				}else continue;

				for (j = 0; j < lqseq->len; j ++){
					index = lqseq->seqs[j].order;
					if (lqseq->seqs[same].len == lqseq->seqs[j].len && strcmp(lqseq->seqs[same].seq, lqseq->seqs[j].seq) == 0) phase_sco[index].s ++;
					else if (lqseq->seqs[heter].len == lqseq->seqs[j].len && strcmp(lqseq->seqs[heter].seq, lqseq->seqs[j].seq) == 0) phase_sco[index].d ++;
				}
				lqseq->indexs = 2;// heter flag
			}else lqseq->indexs = 0;// homo flag
		}
	}
	
	for (i = 0; i < lqseqs_count; i++){
		lqseq = &lqseqs[i];
		if (!lqseq->len) continue;
		mark_del_lqseq(lqseq, phase_sco);
	}

	for (i = 0; i < lqseqs_count; i++){
		lqseq = &lqseqs[i];
		if (!lqseq->len) continue;
		// for (j = 0; j < lqseq->len; j ++){
		// 	printf("p1#%d/%d heter:%d %d s:%d e:%d order:%d score:%d s/d:%d/%d len:%d %s\n",i,lqseqs_count, 
		// 	lqseq->indexs, j, lqseq->start, lqseq->end, lqseq->seqs[j].order,lqseq->seqs[j].kscore, 
		// 	phase_sco[lqseq->seqs[j].order].s, phase_sco[lqseq->seqs[j].order].d, 
		// 	strlen(lqseq->seqs[j].seq), lqseq->seqs[j].seq);
		// }printf("\n");
		remove_differ_phase_lqseq(lqseq, phase_sco);
	}

	for (i = 0; i < lqseqs_count; i++){
		lqseq = &lqseqs[i];
		if (!lqseq->len) continue;
		// s = select_most_lqseq(lqseq, lqseq->len);
		select_most2_lqseq(lqseq, lqseq->len, &s, &k);
		index = lqseq->seqs[s].order;
		if (lqseq->indexs && s != k && s != 0 && lqseq->seqs[k].kscore >= 3 && phase_sco[index].s >= phase_sco[index].d + 3){//NB:for >=3 kmer-types such as polyploid or mult-repeats
			int sphase_sco, kphase_sco;
			for (sphase_sco = kphase_sco = 0, j = 1; j < lqseq->len; j ++){
				index = lqseq->seqs[j].order;
				if (phase_sco[index].d >= 3) continue;
				if (lqseq->seqs[s].len == lqseq->seqs[j].len && strcmp(lqseq->seqs[s].seq, lqseq->seqs[j].seq) == 0) 
					sphase_sco += phase_sco[index].s - phase_sco[index].d;
				else if (lqseq->seqs[k].len == lqseq->seqs[j].len && strcmp(lqseq->seqs[k].seq, lqseq->seqs[j].seq) == 0) 
					kphase_sco += phase_sco[index].s - phase_sco[index].d;
			}
			// printf("%d %d %d %d\n",s,k, sphase_sco, kphase_sco);
			if (sphase_sco < kphase_sco) s = k;
		}else if (lqseq->seqs[0].len > 50 && lqseq->seqs[s].kscore < lqseq->len/3 && lqseq->seqs[s].kscore < 3){
			int sl = remove_differ_len_lqseq(lqseq);
			if (sl <= 3) {//NB: if this region includes lq seqs with large SD., just do not correct.
				s = 0;
				lqseq->seqs[s].kscore = 65534;
			}
		}

		// for (k = 0; k < lqseq->len; k++){
		// 	printf("f#%d/%d %d s:%d e:%d order:%d score:%d s/d:%d/%d len:%d %s\n",i,lqseqs_count, k, lqseq->start, 
		// 		lqseq->end, lqseq->seqs[k].order, lqseq->seqs[k].kscore, phase_sco[lqseq->seqs[k].order].s, 
		// 		phase_sco[lqseq->seqs[k].order].d, strlen(lqseq->seqs[k].seq), lqseq->seqs[k].seq);
		// }
		// printf("most:%d score:%d\n\n",s, lqseq->seqs[s].kscore);
		if (lqseq->seqs[s].kscore > 2 || lqseq->seqs[s].kscore >= lqseq->len/2){//TODO check?
			lqseq->sudoseed = strdup(lqseq->seqs[s].seq);
			lqseq->sudoseed_len = lqseq->seqs[s].len;
			if (lqseq->seqs[s].kscore < lqseq->len /2) string_lower(lqseq->sudoseed);
			lqseq->len = -2;//rm for debug
		}else{
			remove_differ_len_lqseq(lqseq);
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
			count_kmers(lqseq, kmers, LQSEQ_MAX_CAN_COUNT, 0);
			count_kscore(lqseq, kmers, 0);
			unsigned int klastscore, kmaxscore;
			unsigned int kmaxlen = lqseq->seqs[0].len;
			if (kmaxlen > 100){
				uint16_t score[LQSEQ_MAX_CAN_COUNT];
				save_kscore(lqseq, score);
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
		if (max_aln_lqseq_len + lqseq->sudoseed_len > max_aln_length) {
			max_aln_length = max_aln_lqseq_len + lqseq->sudoseed_len;
		}
	}

	// rm lqseq->len = -2; for debug
	// for (i = 0 ; i < lqseqs_count; i++){
	// 	if (lqseqs[i].len){
	// 		for (k = 0; k < lqseqs[i].len; k++){
	// 			j = k >= lqseqs[i].indexs && k <= lqseqs[i].indexe ? 1 : 0;
	// 			printf("%d %d %d %d s:%d e:%d order:%d kscore:%d len:%d %s\n",lqseqs_count, i, k, j, lqseqs[i].start, 
	// 				lqseqs[i].end, lqseqs[i].seqs[k].order,  lqseqs[i].seqs[k].kscore, strlen(lqseqs[i].seqs[k].seq), lqseqs[i].seqs[k].seq);
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
	free (phase_sco);

	return max_aln_length;
}

static int generate_lqseqs_from_tags_kmer1(lqseq *lqseqs, const int lqseqs_count, align_tags_t *tags_list, \
	const unsigned int aligned_seq_count, const int split){

	int i, j, k, s, lable, start, end, index, max_aln_lqseq_len;
	int max_aln_length = max_aln_lqseq_len = 0;

	align_tags_t *tags;
	for (i = 0; i < lqseqs_count; i++){
		start = lqseqs[i].start;
		end = lqseqs[i].end;
		lqseqs[i].len = 0;
		lqseqs[i].seqs = malloc(LQSEQ_MAX_CAN_COUNT * sizeof(struct seq_));
		lqseqs[i].sudoseed = NULL;
		for (j = 0; j < aligned_seq_count; j++){ // select seed reads as template
			index = 0;
			tags = &tags_list[j];
			if (tags->align_tags[0].t_pos <= start && tags->align_tags[tags->len - 1].t_pos >= end){
				lable = 0;
				for (k = start - tags->align_tags[0].t_pos; k < tags->len && tags->align_tags[k].t_pos <= end; k++){
					if (tags->align_tags[k].t_pos >= start && tags->align_tags[k].q_base != 4){
						lqseqs[i].seqs[lqseqs[i].len].seq[index++] = int_to_base[tags->align_tags[k].q_base];
						if (index > lqseq_max_length - 1) {
							lable = 1;
							break;
						}
					}
				}
				if (index && ! lable){
					lqseqs[i].seqs[lqseqs[i].len].seq[index] = '\0';
					lqseqs[i].seqs[lqseqs[i].len].len = index;
					lqseqs[i].seqs[lqseqs[i].len].order = lqseqs[i].len;
					lqseqs[i].len ++;
					if (index > max_aln_lqseq_len) max_aln_lqseq_len = index;
				}
				// if (lqseqs[i].len >= LQSEQ_MAX_CAN_COUNT) break;
				if (lqseqs[i].len >= 15) break;
			}
		}
	}
	
	int max_aln_len = 2000;
	int max_mem_len = max_aln_len + 2;
	uint64_t max_mem_d = (int) max_mem_len * 0.2 + 1;
	alignment aln;
	aln.t_aln_str = malloc(max_mem_len * sizeof(char));
	aln.q_aln_str = malloc(max_mem_len * sizeof(char));
	int *V = malloc(max_mem_d * 2 * sizeof(int));
	uint8_t **D = NULL;
	for (int i = 0; i < 3 && ! D; i++){
		D = malloc( max_mem_d * sizeof(uint8_t *) + max_mem_d * (max_mem_d + 1)/2 * sizeof(uint8_t));
	}
	if (D != NULL){
		uint8_t * const _ = (uint8_t *) (D + max_mem_d);
		for (uint64_t d = 0; d < max_mem_d; d ++ ) {
			D[d] = _ + d * (d + 1)/2;
		}
	}

	lqseq *lqseq;
	for (i = 0; i < lqseqs_count; i++){
		lqseq = &lqseqs[i];
		// if (i) printf("%d %d %d  %d\n",i, lqseq->start, lqseq->end, lqseqs[i-1].end - lqseq->start);
		// for (k = 0; k < lqseq->len; k++){
		// 	printf("#%d/%d %d s:%d e:%d order:%d len:%d %s\n",i,lqseqs_count, k, lqseq->start, 
		// 		lqseq->end, lqseq->seqs[k].order,  strlen(lqseq->seqs[k].seq), lqseq->seqs[k].seq);
		// }
		remove_differ_len_lqseq(lqseq);
		// printf("after remove_differ_len_lqseq merge:\n");
		// for (k = 0; k < lqseq->len; k++){
		// 	printf("#%d/%d %d s:%d e:%d order:%d len:%d %s\n",i,lqseqs_count, k, lqseq->start, 
		// 		lqseq->end, lqseq->seqs[k].order,  strlen(lqseq->seqs[k].seq), lqseq->seqs[k].seq);
		// }

		int select_most = 0;
		if (lqseq->len <= 2){
			select_most = 1;
			s = 0;
			lqseq->seqs[s].kscore = 0;
		}else{
			s = select_most_lqseq(lqseq, lqseq->len);
			for (j = k = 0; j < lqseq->len; j ++){//get second max 2.
				if (j != s && (lqseq->seqs[j].kscore > lqseq->seqs[k].kscore || (lqseq->seqs[j].kscore == \
						lqseq->seqs[k].kscore && lqseq->seqs[j].len > lqseq->seqs[k].len))) k = j;
			}

			if (lqseq->seqs[s].kscore >= lqseq->len * 0.8){
				select_most = 2;
			}else if (lqseq->seqs[s].kscore + lqseq->seqs[k].kscore >= lqseq->len * 0.8){
				select_most = 3;
				if (s != 0){
					if (lqseq->seqs[k].kscore >= lqseq->seqs[s].kscore/2){
						if (k == 0) s = 0;
						else if (D != NULL && lqseq->seqs[0].len < max_aln_len && lqseq->seqs[s].len < max_aln_len && \
							lqseq->seqs[k].len < max_aln_len){
							aln.aln_t_s = aln.aln_len = aln.shift = 0;
							memset(V, 0, max_mem_d * 2 * sizeof(int));
							align_hq(lqseq->seqs[0].seq, lqseq->seqs[0].len, \
								lqseq->seqs[s].seq, lqseq->seqs[s].len, &aln, V, D);
							float s_ide = get_align_ide(&aln);
							aln.aln_t_s = aln.aln_len = aln.shift = 0;
							memset(V, 0, max_mem_d * 2 * sizeof(int));
							align_hq(lqseq->seqs[0].seq, lqseq->seqs[0].len, \
								lqseq->seqs[k].seq, lqseq->seqs[k].len, &aln, V, D);
							float k_ide = get_align_ide(&aln);
							s = k_ide > s_ide ? k : s;
						}
					}
				}
			}else if (k < 3 || lqseq->seqs[s].kscore < lqseq->len/2){ // for high similarity repeat
				for (j = 1, k = 0; j < min(lqseq->len, 3); j ++){
					if (mabs(lqseq->seqs[j].len, lqseq->seqs[0].len) < mabs(lqseq->seqs[j].len, lqseq->seqs[s].len)) k ++;
				}
				if (k >= min(lqseq->len, 3) - 1){
					s = lqseq->seqs[s].len;
					for (j = 1; j < lqseq->len; j ++){
						if (mabs(lqseq->seqs[j].len, lqseq->seqs[0].len) > mabs(lqseq->seqs[j].len, s)){
							for (k = j + 1; k < lqseq->len; k++){
								if (mabs(lqseq->seqs[k].len, lqseq->seqs[0].len) < mabs(lqseq->seqs[k].len, s)) {
									SWAP(lqseq->seqs[j], lqseq->seqs[k], struct seq_);
									break;
								}
							}
							if (k == lqseq->len) break;
						}
					}
					lqseq->len = j;
				}
				// printf("after re-remove_differ_len_lqseq2 merge:\n");
				// for (k = 0; k < lqseq->len; k++){
				// 	printf("#%d/%d %d s:%d e:%d order:%d len:%d %s\n",i,lqseqs_count, k, lqseq->start, 
				// 		lqseq->end, lqseq->seqs[k].order,  strlen(lqseq->seqs[k].seq), lqseq->seqs[k].seq);
				// }
			}else select_most = 4;
		}
		// printf("select order:%d count:%d/%d select_most:%d\n", s, lqseq->seqs[s].kscore,lqseq->len, select_most);
		if (select_most){
			lqseq->sudoseed = strdup(lqseq->seqs[s].seq);
			lqseq->sudoseed_len = lqseq->seqs[s].len;
			if (lqseq->seqs[s].kscore < lqseq->len /2) string_lower(lqseq->sudoseed);
			lqseq->len = -2;//rm for debug
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
			count_kmers(lqseq, kmers, LQSEQ_MAX_CAN_COUNT, 0);
			count_kscore(lqseq, kmers, 0);
			unsigned int klastscore, kmaxscore;
			unsigned int kmaxlen = lqseq->seqs[0].len;
			if (kmaxlen > 100){
				uint16_t score[LQSEQ_MAX_CAN_COUNT];
				save_kscore(lqseq, score);
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
		if (max_aln_lqseq_len + lqseq->sudoseed_len > max_aln_length) {
			max_aln_length = max_aln_lqseq_len + lqseq->sudoseed_len;
		}
	}

	// rm lqseq->len = -2; for debug
	// for (i = 0 ; i < lqseqs_count; i++){
	// 	if (lqseqs[i].len){
	// 		for (k = 0; k < lqseqs[i].len; k++){
	// 			j = k >= lqseqs[i].indexs && k <= lqseqs[i].indexe ? 1 : 0;
	// 			printf("%d %d %d %d s:%d e:%d order:%d kscore:%d len:%d %s\n",lqseqs_count, i, k, j, lqseqs[i].start, 
	// 				lqseqs[i].end, lqseqs[i].seqs[k].order,  lqseqs[i].seqs[k].kscore, strlen(lqseqs[i].seqs[k].seq), lqseqs[i].seqs[k].seq);
	// 		}
	// 		printf("seed%d: %d %s\n",i, lqseqs[i].sudoseed_len, lqseqs[i].sudoseed);
	// 	}else{
	// 		printf("%d %d len:%d s:%d e:%d\n",lqseqs_count, i, lqseqs[i].len, lqseqs[i].start,lqseqs[i].end);
	// 	}
	// }
	// exit(1);
	free(aln.t_aln_str);
	free(aln.q_aln_str);
	free(V);
	if (D) free(D);
	for (i = 0; i < aligned_seq_count; i++){
		free (tags_list[i].align_tags);
	}
	free (tags_list);

	return max_aln_length;
}

static consensus_trimed *get_lqseqs_from_align_tags(align_tags_t *tags_list, msa_p *msa,\
		int aligned_seq_count, int len){
	int p, d, b, m, n;
	int64_t pp_ppp_score = -10;
	align_tag global_best_p = {0, 0, -1};
	
	// msa_p_d_b_pp_ppps pp_ppps;
	// pp_ppps.len = pp_ppps.max_size = 0;
	allocate_msa_mem(msa, len);
	update_msa(msa, tags_list, aligned_seq_count, 1);
	msa_p_d_b_pp_ppp *pp_ppp_m, *pp_ppp_n;

	int factor = READS_TYPE == READS_HIFI ? 4 : 2;
	for (p = 0; p < len; p++){
		for ( d = 0; d < msa[p].max_size; d++ ){
			for (b = 0; b < 6; b++){
				msa_p_d_b *p_base = &msa[p].d_b[d][b];
				p_base->best_score = -10;
				p_base->best_pp.t_pos = -1;
				for (m = 0; m < p_base->len; m++){
					pp_ppp_m = &p_base->pp_ppp[m];
					if (pp_ppp_m->pp.t_pos == -1){
						pp_ppp_m->score = 10 * pp_ppp_m->link_count - factor * msa[p].coverage;
					}else{

						msa_p_d_b *pp_base = &msa[pp_ppp_m->pp.t_pos].d_b[pp_ppp_m->pp.delta][pp_ppp_m->pp.q_base];

						for (n = 0; n < pp_base->len; n++){
							pp_ppp_n = &pp_base->pp_ppp[n];
							if (pp_ppp_n->pp.t_pos == pp_ppp_m->ppp.t_pos &&
								pp_ppp_n->pp.delta == pp_ppp_m->ppp.delta &&
								pp_ppp_n->pp.q_base == pp_ppp_m->ppp.q_base){
								pp_ppp_score = pp_ppp_n->score + 10 * pp_ppp_m->link_count - factor * msa[p].coverage;
								if (pp_ppp_score > pp_ppp_m->score) pp_ppp_m->score = pp_ppp_score;
							}
						}
					}

					if (pp_ppp_m->score > p_base->best_score || (pp_ppp_m->score == p_base->best_score && pp_ppp_m->pp.q_base != 4)){// find the bese score for a given p
						p_base->best_score = pp_ppp_m->score;
						p_base->best_pp = pp_ppp_m->pp;
						// printf("%d %d %d %c\n", p, p_base->best_pp.t_pos, p_base->best_pp.delta, int_to_base[p_base->best_pp.q_base]);
						p_base->best_link_count = pp_ppp_m->link_count;
					}
				}
				global_best_p.t_pos = p;
				global_best_p.delta = d;
				global_best_p.q_base = b;
			}
		}
	}
	// free (pp_ppps.pp_ppp);

	consensus_trimed *cons_trimed;
	cons_trimed = malloc(sizeof(consensus_trimed));
	int cons_trimed_seq_max_len = len * 2 + 1;
	cons_trimed->len = 0;
	cons_trimed->seq = malloc(cons_trimed_seq_max_len * sizeof(char));

	int min_qv_factor = READS_TYPE == READS_HIFI ? 2: 5;
	while (1){
		if (global_best_p.q_base != 4){
			if (cons_trimed->len >= cons_trimed_seq_max_len){
				cons_trimed_seq_max_len += 500;
				cons_trimed->seq = realloc(cons_trimed->seq, cons_trimed_seq_max_len * sizeof(char));
			}
			cons_trimed->seq[cons_trimed->len ++] = 
				msa[global_best_p.t_pos].d_b[global_best_p.delta][global_best_p.q_base].best_link_count * min_qv_factor \
				> msa[global_best_p.t_pos].coverage || int_to_base[global_best_p.q_base] == 'N' ? \
				int_to_base[global_best_p.q_base] : tolower(int_to_base[global_best_p.q_base]);
		}
		global_best_p = msa[global_best_p.t_pos].d_b[global_best_p.delta][global_best_p.q_base].best_pp;
		// printf("%d %d %c\n", global_best_p.t_pos, global_best_p.delta, int_to_base[global_best_p.q_base]);
		if (global_best_p.t_pos == -1){
			break;
		}
	}
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

static int update_lqreg(lqreg lq[], char *seq, unsigned int p, int lq_i, unsigned int *hq_m, unsigned int *lq_m){
	if (seq[p] >= 'a') {
		if (!lq[lq_i].lqlen) lq[lq_i].start = p;
		if ((*lq_m)++ > 2) *hq_m = 0;
		lq[lq_i].end = p;
		lq[lq_i].lqlen ++;
		lq[lq_i].lq_total_len ++;
	}else{
		if (lq[lq_i].lqlen && lq[lq_i].start == 0) {// TODO check, trim start
			lq_i++;
			*hq_m = 0;
		}else if (*hq_m + lq[lq_i].start > lq[lq_i].end || (*hq_m)++ > LQREG_MAX_GAP){
			//here lq[lq_i].end - lq[lq_i].start + 1 != lq[lq_i].lqlen
			if (lq[lq_i].end > lq[lq_i].start + LQREG_MAX_LEN) lq_i ++;
			else lq[lq_i].lqlen = lq[lq_i].end = 0;
			*hq_m = 0;
		}else if (*hq_m >= lq[lq_i].lqlen){
			lq[lq_i].lqlen = lq[lq_i].end = 0;
			*hq_m = 0;
		}
		*lq_m = 0;
	}
	return lq_i;
}

static consensus_trimed *update_consensus_trimed(lqseq *lqseqs, const int lqseqs_count, \
		consensus_data *consensus, const int seed_len){
	
	consensus_trimed *cons_trimed_;
	cons_trimed_ = malloc(sizeof(consensus_trimed));
	int cons_trimed_seq_max_len = seed_len * 2 + 1;
	cons_trimed_->len = 0;
	cons_trimed_->seq = malloc(cons_trimed_seq_max_len * sizeof(char));
	unsigned int p = 0;
	unsigned int i = consensus->lstrip, j;
	int update, lqseqs_index;
	update = 1;
	lqseqs_index = lqseqs_count - 1;

	lqreg lq[LQREG_MAX_COUNT];
	memset(lq, 0, sizeof(lq));
	unsigned int lq_m, hq_m;
	int lq_i = lq_m = hq_m = 0;
	
	while (i < consensus->len - consensus->rstrip){
		p = consensus->cns_bases[i].pos;

		if (lqseqs_index >= 0 && ((lqseqs[lqseqs_index].len <= 0 && lqseqs[lqseqs_index].len != -2) || p > lqseqs[lqseqs_index].end)) {
			lqseqs_index--;
			update = 1;
		}
		if (lqseqs_index >= 0 && (lqseqs[lqseqs_index].len > 0 || lqseqs[lqseqs_index].len == -2) && p >= lqseqs[lqseqs_index].start && p <= lqseqs[lqseqs_index].end){
			if (update){
				// printf("index:%d seed:%s\n", lqseqs_index, lqseqs[lqseqs_index].sudoseed);
				for (j = 0; j < lqseqs[lqseqs_index].sudoseed_len; j++){
					cons_trimed_->seq[cons_trimed_->len++] = lqseqs[lqseqs_index].sudoseed[j];
					// cons_trimed_->seq[cons_trimed_->len++] = 'N';
					if (cons_trimed_->len >= cons_trimed_seq_max_len) reallocate_cons_trimed_mem(cons_trimed_, &cons_trimed_seq_max_len, 500);
					lq_i = update_lqreg(lq, cons_trimed_->seq, cons_trimed_->len - 1, lq_i, &hq_m, &lq_m);
					if (lq_i >= LQREG_MAX_COUNT) break;
				}
				update = 0;
			}
		}else{
			cons_trimed_->seq[cons_trimed_->len++] = consensus->cns_bases[i].base;
			update = 1;
			if (cons_trimed_->len >= cons_trimed_seq_max_len) reallocate_cons_trimed_mem(cons_trimed_, &cons_trimed_seq_max_len, 500);
			lq_i = update_lqreg(lq, cons_trimed_->seq, cons_trimed_->len - 1, lq_i, &hq_m, &lq_m);
			if (lq_i >= LQREG_MAX_COUNT) break;
		}
		i ++;
	}
	if (lq[lq_i].end == cons_trimed_->len - 1) lq_i ++;//TODO check, trim end

	cons_trimed_->seq[cons_trimed_->len] = '\0';
	
	for (i = 0; i < lqseqs_count; i++) {
		free (lqseqs[i].seqs);
		if (lqseqs[i].sudoseed) free (lqseqs[i].sudoseed);
	}
	if (READS_TYPE == READS_HIFI){
		int tolerate_len = 0;
		unsigned int lq_total_len = 0;
		i = 0;
		while (i < cons_trimed_->len && cons_trimed_->seq[i] >= 'a') i ++;
		if (i < tolerate_len) i = 0;
		j = 0;
		while (j < cons_trimed_->len && cons_trimed_->seq[cons_trimed_->len - 1 -j] >= 'a') j ++;
		if (j < tolerate_len) j = 0;
		if (i + j < cons_trimed_->len){
			if (i > 0 || j > 0){
				cons_trimed_->len = cons_trimed_->len - i - j;
				memmove(cons_trimed_->seq, cons_trimed_->seq + i, cons_trimed_->len);
				cons_trimed_->seq[cons_trimed_->len] = '\0';
			}

			int remove_lq_total_len = i + j;
			for (i = 0; i < lq_i; i ++) lq_total_len += lq[i].lq_total_len;
			cons_trimed_->identity = 1 - (float) (lq_total_len - remove_lq_total_len) / cons_trimed_->len;
		}else{
			cons_trimed_->len = 2;
			cons_trimed_->identity = 0;
		}
		cons_trimed_->seq[cons_trimed_->len] = '\0';
		return cons_trimed_;
	}else if (lq_i){
		lq_m = 0;// start
		hq_m = lq[0].start;// end + 1
		p = lq[0].start;// length
		unsigned int lq_total_len = lq[0].lq_total_len - lq[0].lqlen;
		for (i = 1; i < LQREG_MAX_COUNT && lq[i].end; i ++){
			if (lq[i].start - lq[i - 1].end > p){
				lq_m = lq[i - 1].end + 1;
				hq_m = lq[i].start;
				lq_total_len = lq[i].lq_total_len - lq[i].lqlen;
				p = lq[i].start - lq[i - 1].end;
			}
		}

		if (i < LQREG_MAX_COUNT && cons_trimed_->len - lq[i - 1].end > p) {
			lq_m = lq[i - 1].end + 1;
			hq_m = cons_trimed_->len;
			lq_total_len = lq[i].lq_total_len;
		}
		cons_trimed_->len = hq_m - lq_m;
		if (lq_m) memmove(cons_trimed_->seq, cons_trimed_->seq + lq_m, cons_trimed_->len);
		cons_trimed_->identity = 1 - (float) lq_total_len / cons_trimed_->len;
		cons_trimed_->seq[cons_trimed_->len] = '\0';
		return cons_trimed_;
	}else{

		if (cons_trimed_->seq[0] >= 'a'){
			i = 0;
			while (cons_trimed_->seq[i] >= 'a'){i ++;}
			cons_trimed_->len -= i;
			memmove(cons_trimed_->seq, cons_trimed_->seq + i, cons_trimed_->len);
			cons_trimed_->seq[cons_trimed_->len] = '\0';
			lq[0].lq_total_len -= i;
		}
		cons_trimed_->identity = 1 - (float) lq[0].lq_total_len/cons_trimed_->len;
		return cons_trimed_;
	}
}


static inline void get_align_tags(alignment *aln, align_tags_t *tags, msa_p *msa){
	//A--TC--GATCATGC--
	//ATT-CGG-ATCATGCAT

	int i = aln->shift;
	int p = 0;
	uint16_t delta = 0;

	tags->len = aln->aln_len;
	tags->aln_t_s = aln->aln_t_s;
	tags->align_tags = malloc(tags->len * sizeof(align_tag));

	aln->aln_t_s --;

	while (i < aln->aln_len + aln->shift){

		if (aln->t_aln_str[i] != '-'){
			aln->aln_t_s ++;
			delta = 0;
		}

		// printf("get_align_tags: %d %d %d %d\n", i, aln->aln_len, aln->shift, aln->aln_t_s);
		p = i - aln->shift;
		tags->align_tags[p].t_pos = aln->aln_t_s;
		tags->align_tags[p].delta = delta++;
		tags->align_tags[p].q_base = base_to_int[(unsigned char)aln->q_aln_str[i]];

		if (tags->align_tags[p].delta == 0 && aln->q_aln_str[i] != 'M'){
			msa[tags->align_tags[p].t_pos].coverage ++;
		}
		if (tags->align_tags[p].delta >= msa[tags->align_tags[p].t_pos].max_size){
			msa[tags->align_tags[p].t_pos].max_size = tags->align_tags[p].delta + 1;
		}
		// if (p > 1){
		// 	tags->align_tags[p].pp = tags->align_tags[p - 1].p;
		// 	tags->align_tags[p].ppp = tags->align_tags[p - 2].p;
		// }else if (p == 0){
		// 	tags->align_tags[p].pp.t_pos = -1;
		// 	tags->align_tags[p].pp.delta = 0;
		// 	tags->align_tags[p].pp.q_base = 0;
		// 	tags->align_tags[p].ppp.t_pos = -1;
		// 	tags->align_tags[p].ppp.delta = 0;
		// 	tags->align_tags[p].ppp.q_base = 0;
		// }else{
		// 	tags->align_tags[p].pp = tags->align_tags[p - 1].p;
		// 	tags->align_tags[p].ppp.t_pos = -1;
		// 	tags->align_tags[p].ppp.delta = 0;
		// 	tags->align_tags[p].ppp.q_base = 0;
		// }
		i ++ ;
	}
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
	int aligned_linkseq_max = lqseq_max_length * 2 + 1;
	int aligned_linkseq_len = 0;
	int delta;
	aln_link.shift = 0;
	aln_link.t_aln_str = malloc(aligned_linkseq_max * sizeof(char));
	aln_link.q_aln_str = malloc(aligned_linkseq_max * sizeof(char));

	int *V = malloc( max_mem_d * 2 * sizeof(int));
	uint8_t **D = malloc( max_mem_d * sizeof(uint8_t *) + max_mem_d * (max_mem_d + 1)/2 * sizeof(uint8_t));
	uint8_t * const _ = (uint8_t *) (D + max_mem_d);
	for (uint64_t d = 0; d < max_mem_d; d ++ ) {
		D[d] = _ + d * (d + 1)/2;
	}

	void (*aln_fun)(char *, int, char *, int, alignment *, int *, uint8_t **);
	aln_fun = READS_TYPE == READS_HIFI ? align_hq : align;

	for (int i = 0; i < seq_count; i++){
		aligned_linkseq_len = 0;
		aln_link.aln_t_s = 0;
		aln_link.aln_t_e = 0;
		aln_link.aln_len = 0;

		for (int j = lqseqs_count - 1; j >= 0; j--){
			if (lqseqs[j].len <= 0) continue;

			seed_len = lqseqs[j].sudoseed_len;

			aligned_linkseq_len += seed_len + 1;
			aln_link.aln_t_e += seed_len + 1;
			aln_link.t_aln_str[aln_link.aln_len] = 'N';
			aln_link.q_aln_str[aln_link.aln_len] = 'N';
			aln_link.aln_len ++;

			query_len = (i + lqseqs[j].indexs) > lqseqs[j].indexe ? seed_len : lqseqs[j].seqs[i + lqseqs[j].indexs].len;
			if (aln_link.aln_len + query_len + seed_len + 1 >= aligned_linkseq_max){
				aligned_linkseq_max += query_len + seed_len + 1;
				aln_link.t_aln_str = realloc(aln_link.t_aln_str, aligned_linkseq_max * sizeof(char));
				aln_link.q_aln_str = realloc(aln_link.q_aln_str, aligned_linkseq_max * sizeof(char));
			}

			if (i + lqseqs[j].indexs > lqseqs[j].indexe || (i && (query_len < seed_len * 0.5 || query_len > seed_len * 1.3))){
				char _[seed_len + 1];
				memset(_, 'M', seed_len * sizeof(char));
				_[seed_len] = '\0';
				
				strcpy(&aln_link.t_aln_str[aln_link.aln_len], _);
				strcpy(&aln_link.q_aln_str[aln_link.aln_len], _);
				aln_link.aln_len += seed_len;
			}else{
				aln.aln_t_s = 0;
				aln.aln_t_e = seed_len;
				aln.aln_len = 0;
				aln.shift = 0;

				memset(V, 0, max_mem_d * 2 * sizeof(int));
				aln_fun(lqseqs[j].seqs[i + lqseqs[j].indexs].seq, query_len,\
							lqseqs[j].sudoseed, seed_len, &aln, V, D);

				// if (lqseqs[j].seqs[i + lqseqs[j].indexs].order == 0 && seed_len > 1000) {
					// printf("%d %d %d %d\n", i, j, query_len, seed_len);
					// printf("#%d %d %d %s\n", i, j, query_len, lqseqs[j].seqs[i + lqseqs[j].indexs].seq);
					// out_align(&aln);}
				
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
					char _[seed_len + 1];
					memset(_, 'M', seed_len * sizeof(char));
					_[seed_len] = '\0';
					strcpy(&aln_link.t_aln_str[aln_link.aln_len], _);
					strcpy(&aln_link.q_aln_str[aln_link.aln_len], _);
					aln_link.aln_len += seed_len;
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
		// out_align(&aln_link);
		get_align_tags(&aln_link, &tags_list[i], msa);

	}
	free(aln.t_aln_str);
	free(aln.q_aln_str);
	free(aln_link.t_aln_str);
	free(aln_link.q_aln_str);
	free(V);
	free(D);
	consensus_trimed *cons_trimed = get_lqseqs_from_align_tags( tags_list, msa, seq_count, aligned_linkseq_len);	
	free_msa(msa, aligned_linkseq_len);
	return cons_trimed;
}

static consensus_trimed *iterate_generate_consensus_trimed(lqseq *lqseqs, const int lqseqs_count,\
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
		// printf("cons_trimed->seq: %s\n", cons_trimed->seq);
		for (k = cons_trimed->len; k; k --){
			if (cons_trimed->seq[k - 1] != 'N'){
				if (cons_trimed->seq[k - 1] < 'a'){
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
					if (lqseqs[j].lqcount > lqseqs[j].sudoseed_len * 4/5) lqseqs[j].len = -1;
					// printf("%d #%d %d %d %s\n",k, j, i, lqseqs[j].sudoseed_len, lqseqs[j].sudoseed);
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
	cons_trimed = update_consensus_trimed(lqseqs, lqseqs_count, consensus, len);
	return cons_trimed;
}

static consensus_trimed *generate_cns_from_best_score_fast( msa_p *msa, align_tags_t *tags_list, align_tag *global_best_p,\
		int aligned_seq_count, int len, int min_cov, float min_error_corrected_ratio){

	int cons_trimed_seq_max_len = len * 2 + 1;
	consensus_trimed *cons_trimed = calloc(1, sizeof(consensus_trimed));
	cons_trimed->seq = malloc(cons_trimed_seq_max_len * sizeof(char));

	lqreg lq[LQREG_MAX_COUNT];
	memset(lq, 0, sizeof(lq));
	int lq_i = 0;

	while (1){
		if (global_best_p->q_base != 4){
			if (msa[global_best_p->t_pos].coverage > min_cov){
				cons_trimed->seq[cons_trimed->len++] = int_to_base[global_best_p->q_base];
				if (lq[lq_i].end >= lq[lq_i].start + 50 || !lq_i) {
					if (++lq_i >= LQREG_MAX_COUNT) break;
				}else lq[lq_i].end = 0;
			}else{
				cons_trimed->seq[cons_trimed->len++] = tolower(int_to_base[global_best_p->q_base]); // lower-case for the base with low coverage (maybe uncorrected).
				if (!lq[lq_i].end){
					lq[lq_i].start = cons_trimed->len - 1;
					lq[lq_i].lqlen = 0;
				}
				lq[lq_i].end = cons_trimed->len - 1;
				lq[lq_i].lq_total_len ++;
				lq[lq_i].lqlen ++;
			}
			if (cons_trimed->len >= cons_trimed_seq_max_len){
				cons_trimed_seq_max_len += 500;
				cons_trimed->seq = realloc(cons_trimed->seq, cons_trimed_seq_max_len * sizeof(char));
			}
		}
		global_best_p = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base].best_pp;
		if (global_best_p->t_pos == -1){
			break;
		}
	}

	int i, lq_m = 0, hq_m = lq[0].start, l = hq_m;
	unsigned int lq_total_len = lq[0].lq_total_len - lq[0].lqlen;
	for (i = 1; i < LQREG_MAX_COUNT && lq[i].end; i ++){
		if (lq[i].start - lq[i - 1].end > l){
			lq_m = lq[i - 1].end + 1;
			hq_m = lq[i].start;
			lq_total_len = lq[i].lq_total_len - lq[i].lqlen;
			l = hq_m - lq_m;
		}
	}

	if (i < LQREG_MAX_COUNT && cons_trimed->len - lq[i - 1].end > l) {
		lq_m = lq[i - 1].end + 1;
		hq_m = cons_trimed->len;
		lq_total_len = lq[i].lq_total_len;
	}

	cons_trimed->len = hq_m - lq_m;
	cons_trimed->identity = 1 - (float) lq_total_len / cons_trimed->len;
	if (lq_m) memmove(cons_trimed->seq, cons_trimed->seq + lq_m, cons_trimed->len);
	reverse_str(cons_trimed->seq, cons_trimed->len);
	cons_trimed->seq[cons_trimed->len] = '\0';

	for (int p = 0; p < aligned_seq_count; p++){
		free (tags_list[p].align_tags);
	}
	free (tags_list);
	return cons_trimed;
}

static consensus_trimed *generate_cns_from_best_score_kmer( msa_p *msa, align_tags_t *tags_list, align_tag *global_best_p,\
		int aligned_seq_count, int len, int min_cov, float min_error_corrected_ratio, int split){
	
	int lq_min_length = 2;
	int p, qv, lq, lq_s, lq_e;
	lq_s = lq_e = -1;
	p = qv = lq = 0;

	int lqseqs_max_size = 200;
	int lqseqs_index = 0;
	int lqseq_total_length = 0;
	lqseq *lqseqs = malloc(lqseqs_max_size * sizeof(lqseq));

	msa_p_d_b *global_best_d_b;
	consensus_data consensus;
	memset(&consensus, 0, sizeof(consensus));
	consensus.max_size = len/5 * 6 + 1;
	consensus.cns_bases = malloc( consensus.max_size * sizeof(consensus_base) );

	int dag_min_qv = 80;
	align_tags_t *ref_tag = &tags_list[0];
	while (1){
		if (global_best_p->q_base != 4){

			if (p >= consensus.max_size){
				consensus.max_size += 5000;
				consensus.cns_bases = realloc(consensus.cns_bases, consensus.max_size * sizeof(consensus_base));
			}
			consensus.cns_bases[p].pos = global_best_p->t_pos;
			global_best_d_b = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base];
			qv = 100 * global_best_d_b->best_link_count / msa[global_best_p->t_pos].coverage;
			// printf("p %d %d qn:%d depth:%d\n",p, global_best_p->t_pos, qv, msa[global_best_p->t_pos].coverage);
			if (msa[global_best_p->t_pos].coverage < 4){
				lq = 0;
				lq_s = -1;
				lqseq_total_length ++;
			}else if (qv < dag_min_qv || global_best_p->q_base != ref_tag->align_tags[global_best_p->t_pos].q_base){
			// }else if (qv < dag_min_qv || (global_best_p->q_base != ref_tag->align_tags[global_best_p->t_pos].q_base && msa[global_best_p->t_pos].coverage >= global_best_d_b->best_link_count + 2)){

				if (lq_s == -1) lq_s = p;
				lq_e = p;
				lq = 1;
				lqseq_total_length ++;
			}else if (lq && p - lq_e > 2 * lq_min_length && consensus.cns_bases[p].pos != consensus.cns_bases[p - 1].pos){
				lq_e = p - lq_min_length - 1;
				lq_s = lq_s > lq_min_length ? lq_s - lq_min_length : 1;
				// if (lqseqs_index >= 1) printf("%d %d %d %d %d %d\n",lq_s,lq_e,consensus.cns_bases[lq_s].pos,consensus.cns_bases[lq_e].pos,lqseqs[lqseqs_index - 1].start,lqseqs[lqseqs_index - 1].end );
				if (lqseqs_index >= 1 && consensus.cns_bases[lq_s].pos >= lqseqs[lqseqs_index - 1].start){
					lqseqs[lqseqs_index - 1].start = consensus.cns_bases[lq_e].pos;
				}else{			
					lqseqs[lqseqs_index].end = consensus.cns_bases[lq_s].pos;
					lqseqs[lqseqs_index].start = consensus.cns_bases[lq_e].pos;
					if (++lqseqs_index >= lqseqs_max_size) {
						lqseqs_max_size += 100;
						lqseqs = realloc(lqseqs, lqseqs_max_size * sizeof(lqseq));
					}
				}
				lq = 0;
				lq_s = -1;
			}

			// if (msa[global_best_p->t_pos].coverage > min_cov && qv > dag_min_qv){
			if (msa[global_best_p->t_pos].coverage > min_cov){
				consensus.cns_bases[p].base = int_to_base[global_best_p->q_base];
			}else{
				consensus.cns_bases[p].base = tolower(int_to_base[global_best_p->q_base]); 
			}
			p ++;
		}

		global_best_p = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base].best_pp;
		if (global_best_p->t_pos == -1){
			break;
		}
	}

	consensus.len = p;//include lstrip & rstrip
	consensus_trimed *cons_trimed;
	if (consensus.len > 2 && lqseq_total_length < consensus.len * DAG_MAX_RATIO && consensus.uncorrected_len -\
		consensus.lstrip - consensus.rstrip < (consensus.len - consensus.lstrip - consensus.rstrip) * (1 - min_error_corrected_ratio)){
		// time_debug(clock(), "lq start...");
		reverse_consensus_base(consensus.cns_bases, consensus.len);
		int max_aln_length = generate_lqseqs_from_tags_kmer(lqseqs, lqseqs_index, tags_list, aligned_seq_count, split);
		cons_trimed = iterate_generate_consensus_trimed(lqseqs, lqseqs_index, &consensus, len, max_aln_length, 2);
		// time_debug(clock(), "lq end...");
	}else{
		for (p = 0; p < aligned_seq_count; p++){
			free (tags_list[p].align_tags);
		}
		free (tags_list);
		cons_trimed = error_seed(2);
	}

	free (lqseqs);
	free (consensus.cns_bases);
	// free (consensus);
	return cons_trimed;
}

static consensus_trimed *generate_cns_from_best_score( msa_p *msa, align_tags_t *tags_list, align_tag *global_best_p,\
		int aligned_seq_count, int len, int min_cov, float min_error_corrected_ratio, int split){
	int p = 0;
	int lable = 1;
	// int lq_min_qv = 40;
	int lq_min_length = 8;
	int qv, pqv, hq, lq, lq_l, lq_s, lq_e;
	lq_s = lq_e = -1;
	hq = qv = lq_l = lq = 0;
	int lqseqs_max_size = 200;
	int lqseqs_index = 0;
	int lqseq_total_length = 0;
	int global_best_quality = 0;

	lqseq *lqseqs = malloc(lqseqs_max_size * sizeof(lqseq));

	consensus_data consensus;
	memset(&consensus, 0, sizeof(consensus_data));
	// consensus = calloc(1, sizeof(consensus_data));
	consensus.max_size = len * 2 + 1;
	consensus.cns_bases = malloc( consensus.max_size * sizeof(consensus_base) );
	// time_debug(clock(), "lq prodecing.......");
	while (1){
		if (global_best_p->q_base != 4){

			if (p >= consensus.max_size){
				consensus.max_size += 500;
				consensus.cns_bases = realloc(consensus.cns_bases, consensus.max_size * sizeof(consensus_base));
			}
			consensus.cns_bases[p].pos = global_best_p->t_pos;
			global_best_quality = msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base].best_link_count;
			pqv = 100 * global_best_quality/msa[global_best_p->t_pos].coverage;

			if (pqv > DAG_MIN_QV){
				hq++;
			}else{
				hq = 0;
				lqseq_total_length ++;
			}
			
			if (hq > lq_min_length/2 && lq_e - lq_s < lq_min_length/2){
				qv = lq_l = lq = 0;
				lq_s = -1;
			}

			if ((qv + pqv)/(lq_l + 1) < DAG_MIN_QV){
				if (lq_s == -1) lq_s = p;
				lq_e = p;
				lq = 1;
				lq_l ++;
				qv += pqv;
			}else if (lq && p - lq_e > 2 * lq_min_length && consensus.cns_bases[p].pos != consensus.cns_bases[p - 1].pos){
				if (lq_e - lq_s + 1 > lq_min_length && lq_e - lq_s + 1 < lqseq_max_length){
					lq_e = p - lq_min_length - 1;
					lq_s = lq_s > lq_min_length ? lq_s - lq_min_length : 1;

					lqseqs[lqseqs_index].end = consensus.cns_bases[lq_s].pos;
					lqseqs[lqseqs_index].start = consensus.cns_bases[lq_e].pos;
					if (lqseqs_index >= 1 && lqseqs[lqseqs_index].end == lqseqs[lqseqs_index -1].start){
						while (lqseqs[lqseqs_index].end == lqseqs[lqseqs_index - 1].start && lq_s < p - 4){
							lqseqs[lqseqs_index].end = consensus.cns_bases[++lq_s].pos;
						}
					}

					if (++lqseqs_index >= lqseqs_max_size) {
						lqseqs_max_size += 100;
						lqseqs = realloc(lqseqs, lqseqs_max_size * sizeof(lqseq));
					}
				// }else if (lq_e - lq_s + 1 >= lqseq_max_length){
				// 	p = 2;//discard corrected seed with length of low quality region >=lqseq_max_length bp.
				// 	break;
				}
				qv = lq_l = lq = 0;
				lq_s = -1;
			}else if (lq && consensus.cns_bases[p].pos != consensus.cns_bases[p - 1].pos){
				qv = lq_l = 0;
			}

			if (msa[global_best_p->t_pos].coverage > min_cov && pqv > LQBASE_MIN_QV){
			// if (msa[global_best_p->t_pos].coverage > min_cov){
				consensus.cns_bases[p].base = int_to_base[global_best_p->q_base];
				lable = 0;
				consensus.lstrip = 0;
			}else{
				consensus.cns_bases[p].base = tolower(int_to_base[global_best_p->q_base]); 
				// lower-case for the base with low coverage or quality (maybe uncorrected).
				consensus.uncorrected_len ++;
				consensus.lstrip ++;
				if (lable) consensus.rstrip ++;
			}
			p ++;
		}
		
		global_best_p = &msa[global_best_p->t_pos].d_b[global_best_p->delta][global_best_p->q_base].best_pp;
		if (global_best_p->t_pos == -1){
			break;
		}
	}

	consensus.len = p;//include lstrip & rstrip
	consensus_trimed *cons_trimed;
	if (consensus.len > 2 && lqseq_total_length < consensus.len * DAG_MAX_RATIO && consensus.uncorrected_len -\
		consensus.lstrip - consensus.rstrip < (consensus.len - consensus.lstrip - consensus.rstrip) * (1 - min_error_corrected_ratio)){
		// time_debug(clock(), "lq start...");
		reverse_consensus_base(consensus.cns_bases, consensus.len);
		int max_aln_length = generate_lqseqs_from_tags(lqseqs, lqseqs_index, tags_list, aligned_seq_count, split);

		cons_trimed = iterate_generate_consensus_trimed(lqseqs, lqseqs_index, &consensus, len, max_aln_length, 2);
		// time_debug(clock(), "lq end...");
	}else{
		for (p = 0; p < aligned_seq_count; p++){
			free (tags_list[p].align_tags);
		}
		free (tags_list);
		cons_trimed = error_seed(2);
	}

	free (lqseqs);
	free (consensus.cns_bases);
	// free (consensus);
	return cons_trimed;
}

int get_terminal_ssr(int kmers[256], const int ssr_range, const int ssr_len, const char *seq, const int s){

	int i, k;
	uint8_t kmer;
	memset(kmers, 0, sizeof(int) * 256);
	for (i = kmer = 0; i < ssr_range; i ++){
		if (i){
			kmer = kmer << 2 | base_to_int[(unsigned char) seq[s + i + ssr_len - 1]];
		}else{
			for (k = 0; k < ssr_len; k ++){
				kmer = kmer << 2 | base_to_int[(unsigned char) seq[s + k]];
			}
		}
		kmers[kmer] ++;
	}

	for (i = k = 0; i < 256; i++){
		if (kmers[i] > k){
			k = kmers[i];
			kmer = i;
		}
	}
	// printf("%d %d\n",kmer,k );
	return kmer;
}

static int clip_terminal_ssr(const char *seq, int seq_len, int ssr_len, int kmer, int dire){
	
	int i, p, p1, p2, k, gap = 20;
	uint8_t kmer_t;
	if (dire){//3->5
		for (i = kmer_t = 0; i < 8; i += 2){
			kmer_t = kmer_t << 2 | (kmer >> i & 3);
		}
		// printf("%d %d\n",kmer_t, kmer );
		seq_len --;
		for (kmer = kmer_t, i = p = p1 = p2 = kmer_t = 0; i < seq_len - ssr_len; i ++){
			if (i){
				kmer_t = kmer_t << 2 | base_to_int[(unsigned char) seq[seq_len - i - ssr_len + 1]];
			}else{
				for (k = 0; k < ssr_len; k ++){
					kmer_t = kmer_t << 2 | base_to_int[(unsigned char) seq[seq_len - k]];
				}
			}
			// printf("%d %d %d %d %d\n",seq_len-i,kmer_t,kmer,p,i);
			if (kmer_t != kmer){
				// if (i - p > gap && seq[seq_len - i] < 'a') break;
				if (i - p > gap){
					if (!p1) p1 = p;
					else if (p2){
						if (i - p2 < 100){ p = p1; break;}
						else p1 = p2 = 0;
					}
				}
			}else {
				p = i;
				if (p1 && p2 == 0) p2 = p;
			}
		}
		p = p > 100 ? p + ssr_len : 0;
	}else{//5->3
		for (i = p = p1 = p2 = kmer_t = 0; i < seq_len - ssr_len; i ++){
			if (i){
				kmer_t = kmer_t << 2 | base_to_int[(unsigned char) seq[i + ssr_len - 1]];
			}else{
				for (k = 0; k < ssr_len; k ++){
					kmer_t = kmer_t << 2 | base_to_int[(unsigned char) seq[k]];
				}
			}
			// printf("%d %d %d\n",i,kmer_t,kmer);
			// if (kmer_t != kmer){
			// 	if ( i - p > gap && seq[i] < 'a') break;
			// }else p = i;

			if (kmer_t != kmer){
				if (i - p > gap){
					if (!p1) p1 = p;
					else if (p2){
						if (i - p2 < 100){ p = p1; break;}
						else p1 = p2 = 0;
					}
				}
			}else {
				p = i;
				if (p1 && p2 == 0) p2 = p;
			}
		}
		p = p > 100 ? p + ssr_len : 0;
	}
	return p;
}

static void trim_terminal_ssr(consensus_trimed *cons_trimed){
	int kmer, clip_s, clip_e, ssr_range = 24, ssr_len = 4, kmers[256] = {0};
	clip_s = clip_e = 0;
	kmer = get_terminal_ssr(kmers, ssr_range, ssr_len, cons_trimed->seq, 0);
	if (kmers[kmer] >= 4){
		clip_s = clip_terminal_ssr(cons_trimed->seq, cons_trimed->len, ssr_len, kmer, 0);
		while (clip_s < cons_trimed->len && cons_trimed->seq[clip_s] >= 'a') clip_s ++;
		// if (clip_s){
		// 	printf("s1: %d %d %s\n",kmer, kmers[kmer], cons_trimed->seq);
		// 	printf("s2: %d %d %s\n\n", kmer, clip_s, cons_trimed->seq + clip_s);
		// }
	}
	kmer = get_terminal_ssr(kmers, ssr_range, ssr_len, cons_trimed->seq, cons_trimed->len - ssr_range - ssr_len + 1);
	if (kmers[kmer] >= 4){
		clip_e = clip_terminal_ssr(cons_trimed->seq, cons_trimed->len, ssr_len, kmer, 1);
		while (clip_e < cons_trimed->len && cons_trimed->seq[cons_trimed->len - clip_e - 1] >= 'a') clip_e ++;
		// if (clip_e){
		// 	printf("e1: %d %d %s\n",kmer, kmers[kmer], cons_trimed->seq);
		// 	cons_trimed->seq[cons_trimed->len - clip_e] = '\0';
		// 	printf("e2: %d %d %s\n", kmer, clip_e, cons_trimed->seq);
		// }
	}
	if (clip_s + clip_e < cons_trimed->len - 10){
		if (clip_e) cons_trimed->seq[cons_trimed->len - clip_e] = '\0';
		cons_trimed->len -= clip_s + clip_e;
		if (clip_s) memmove(cons_trimed->seq, cons_trimed->seq + clip_s, cons_trimed->len + 1);
	}else cons_trimed->len = 4;
	// printf("%d\n",strlen(cons_trimed->seq) );
}

static consensus_trimed *get_cns_from_align_tags( align_tags_t *tags_list, msa_p *msa, int aligned_seq_count, \
		int len, int min_cov, float min_error_corrected_ratio, int split, int fast){
	int p, d, b, m, n;
	// time_debug(clock(), "malloc msa begin...");

	allocate_msa_mem(msa, len);
	// time_debug(clock(), "update msa loop begin...");
	update_msa(msa, tags_list, aligned_seq_count, 0);
	// time_debug(clock(), "update msa done...");
	// exit(1);

	int64_t p_pp_score, p_pp_score_, pp_ppp_score, global_best_score;
	pp_ppp_score = global_best_score = -10;
	align_tag global_best_p = {0, 0, -1};

	msa_p_d_b_pp_ppp *pp_ppp_m, *pp_ppp_n;

	int factor = READS_TYPE == READS_HIFI ? 4 : 3;

	for (p = 0; p < len; p++){
		for ( d = 0; d < msa[p].max_size; d++ ){
			for (b = 0; b < 5; b++){
				msa_p_d_b *p_base = &msa[p].d_b[d][b];
				p_base->best_score = -10;
				p_base->best_pp.t_pos = -1;
				p_pp_score_ = p_pp_score = INT64_MIN;
				for (m = 0; m < p_base->len; m++){
					pp_ppp_m = &p_base->pp_ppp[m];
					if (pp_ppp_m->pp.t_pos == -1){
						pp_ppp_m->score = 10 * pp_ppp_m->link_count - factor * msa[p].coverage;
					}else{

						msa_p_d_b *pp_base = &msa[pp_ppp_m->pp.t_pos].d_b[pp_ppp_m->pp.delta][pp_ppp_m->pp.q_base];

						for (n = 0; n < pp_base->len; n++){
							pp_ppp_n = &pp_base->pp_ppp[n];

							if (pp_ppp_n->pp.t_pos == pp_ppp_m->ppp.t_pos &&
								pp_ppp_n->pp.delta == pp_ppp_m->ppp.delta &&
								pp_ppp_n->pp.q_base == pp_ppp_m->ppp.q_base){
								pp_ppp_score = pp_ppp_n->score + 10 * pp_ppp_m->link_count - factor * msa[p].coverage;
								if (pp_ppp_score > pp_ppp_m->score) {
									pp_ppp_m->score = pp_ppp_score;
									p_pp_score_ = pp_ppp_n->score;
								}
								if (pp_ppp_n->score > p_pp_score && (pp_ppp_m->pp.q_base == 4 || pp_ppp_m->pp.q_base == b)){//TODO CHECK ONT/PB
								// if (pp_ppp_n->score > p_pp_score && (pp_ppp_m->pp.q_base == 4 || pp_ppp_m->pp.q_base == b
									// || pp_ppp_m->ppp.q_base == b || pp_ppp_m->pp.q_base == pp_ppp_m->ppp.q_base)){
									p_pp_score = pp_ppp_n->score;
									p_base->best_score = pp_ppp_m->score;
									p_base->best_pp = pp_ppp_m->pp;
									p_base->best_link_count = pp_ppp_m->link_count;
								}
							}
						}
					}

					if (pp_ppp_m->score > p_base->best_score || (pp_ppp_m->score == p_base->best_score && pp_ppp_m->pp.q_base != 4)){// find the bese score for a given p
						p_pp_score = p_pp_score_;
						p_base->best_score = pp_ppp_m->score;
						p_base->best_pp = pp_ppp_m->pp;
						p_base->best_link_count = pp_ppp_m->link_count;
					}
				}
				if (p_base->best_score >= global_best_score - 3000){
					global_best_p.t_pos = p;
					global_best_p.delta = d;
					global_best_p.q_base = b;
					if (p_base->best_score > global_best_score) global_best_score = p_base->best_score;
				}
			}
		}
	}
	// free (pp_ppps.pp_ppp);
	// time_debug(clock(), "update score chain done...");

	consensus_trimed *cons_trimed;
	if (fast){
		cons_trimed = generate_cns_from_best_score_fast(msa, tags_list, &global_best_p,
			aligned_seq_count, len, min_cov, min_error_corrected_ratio);
	}else{
		cons_trimed = READS_TYPE == READS_HIFI ? generate_cns_from_best_score_kmer(msa, tags_list, &global_best_p, \
			aligned_seq_count, len, min_cov, min_error_corrected_ratio, split) : generate_cns_from_best_score(msa, \
			tags_list, &global_best_p, aligned_seq_count, len, min_cov, min_error_corrected_ratio, split);
		if (cons_trimed->len > 1000 && cons_trimed->identity > 0.8) trim_terminal_ssr(cons_trimed);
	}
	return cons_trimed;
}

consensus_trimed *nextCorrect(char **seqs, unsigned int *aln_start, unsigned int *aln_end, unsigned int seq_count, \
	unsigned int max_mem_len, unsigned int min_len_aln, unsigned int max_cov_aln, unsigned int min_cov, \
		unsigned int lqseq_max_length_, float min_error_corrected_ratio, unsigned int split, unsigned int fast, int read_type){
	
	//for hifi debug
	// split = 1;
	// max_cov_aln = 130;
	// min_cov = 2;
	// lqseq_max_length_ = 1000;

	// time_debug(clock(), "begin...");
	READS_TYPE = read_type;
	if ((lqseq_max_length = lqseq_max_length_) > DAG_MAX_LENGTH) lqseq_max_length = DAG_MAX_LENGTH;

	alignment aln;
	align_tags_t *tags_list;
	msa_p * msa;
	int aligned_seq_count = 0;
	int total_cov_aln = 0;
	int seed_len = aln_end[0] + 1;
	tags_list = malloc(seq_count * sizeof(align_tags_t));
	msa = calloc(seed_len, sizeof(msa_p));

	//reduce cpu time cost for malloc during reads aligning
	max_mem_len += 2;
	uint64_t max_mem_d = (int) max_mem_len * 0.4 + 1;
	aln.t_aln_str = malloc(max_mem_len * sizeof(char));
	aln.q_aln_str = malloc(max_mem_len * sizeof(char));
	int *V = malloc( max_mem_d * 2 * sizeof(int));
	
	uint8_t **D = NULL;
	for (int i = 0; i < 3 && ! D; i++){
		D = malloc( max_mem_d * sizeof(uint8_t *) + max_mem_d * (max_mem_d + 1)/2 * sizeof(uint8_t));
	}

	if (D == NULL) {
		free(V);
		free (msa);
		free (tags_list);
		free(aln.t_aln_str);
		free(aln.q_aln_str);
		return error_seed(3);
	}

	uint8_t * const _ = (uint8_t *) (D + max_mem_d);
	for (uint64_t d = 0; d < max_mem_d; d ++ ) {
		D[d] = _ + d * (d + 1)/2;
	}

	void (*aln_fun)(char *, int, char *, int, alignment *, int *, uint8_t **);
	aln_fun = READS_TYPE == READS_HIFI ? align_hq : align;

	for (int i = 0; i < seq_count && total_cov_aln/seed_len <= max_cov_aln; i++){
		aln.aln_t_s = aln_start[i];
		aln.aln_t_e = aln_end[i];
		aln.aln_len = 0;
		aln.shift = 0;
		int seq1_len = aln.aln_t_e - aln.aln_t_s + 1;
		int seq2_len = strlen(seqs[i]);

		if (i == 0){
			strcpy(aln.t_aln_str, seqs[0]);
			strcpy(aln.q_aln_str, seqs[i]);
			aln.aln_len = seq2_len;
		}else {
			memset(V, 0, max_mem_d * 2 * sizeof(int));
			aln_fun(seqs[i], seq2_len, seqs[0] + aln.aln_t_s, seq1_len, &aln, V, D);
			get_align_shift(&aln, 8);
			// out_align(&aln);
		}
		if (aln.aln_len >= min_len_aln){
			total_cov_aln += aln.aln_t_e - aln.aln_t_s + 1;
			get_align_tags(&aln, &tags_list[aligned_seq_count++], msa);
		}
	}

	free(aln.t_aln_str);
	free(aln.q_aln_str);
	free(V);
	free(D);
	// time_debug(clock(), "after align");
	consensus_trimed *cons_trimed = get_cns_from_align_tags( tags_list, msa, aligned_seq_count, seed_len, min_cov,\
		min_error_corrected_ratio, split, fast);
	free_msa(msa, seed_len);
	// time_debug(clock(), "done....");
	return cons_trimed;
}

void free_consensus_trimed(consensus_trimed *consensus_trimed){
	free(consensus_trimed->seq);
	free (consensus_trimed);
}

int main(){

	// start = clock();
	// time_i = 0;
	const int t = 500;
	char ids[t][1024];
	unsigned int aln_start[t];
	unsigned int aln_end[t];
	unsigned int min_cov = 4;
	unsigned int count = 0;
	unsigned int max_aln_length = 0;

	char **seqs = calloc(t, sizeof(char *));
	for (int i = 0;i < t; i++){
		seqs[i] = calloc(500000, sizeof(char));
	}

	char *line = NULL;
	size_t len = 0;
	
	while (getline(&line, &len, stdin) != -1){
		sscanf(line, "%s\t%d\t%d\t%s", ids[count], &aln_start[count], &aln_end[count], seqs[count]);
		
		if (strcmp(ids[count], "+") != 0 && count + 1 >= t){
			continue;
		}

		if (strcmp(ids[count ++], "+") == 0){
			// time_i ++;
			consensus_trimed *consensus_trimed = nextCorrect(seqs, aln_start, aln_end, count - 1, \
				max_aln_length, 500, 90, min_cov, DAG_MAX_LENGTH, 0.8, 1, 0, READS_HIFI);
			// time_debug(clock(), " after generate_consenses");
			// printf("%d %f\n",consensus_trimed->len,consensus_trimed->identity );
			if (consensus_trimed->len > 1000 && consensus_trimed->identity > 0.8) {
				printf(">%s %d %f\n", ids[0], consensus_trimed->len, consensus_trimed->identity);
				printf("%s\n", consensus_trimed->seq);
				// char *p = consensus_trimed->seq;
				// while (*p++ > 90) {};
				// p --;
				// while (p < consensus_trimed->seq + consensus_trimed->len){
				// 	putchar(*p++);
				// }
				// putchar('\n');
			}
			free_consensus_trimed(consensus_trimed);
			max_aln_length = count = 0;
		}else{
			unsigned int _ = strlen(seqs[count-1]) + aln_end[count-1] - aln_start[count-1] + 1;
			if (_ > max_aln_length) max_aln_length = _;
		}
	}

	free(line);
	for (int i = 0;i < t; i++){
		free (seqs[i]);
	}
	free (seqs);

	return 0;
}
