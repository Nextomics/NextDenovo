#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"
#include "../lib/ovl.h"
#include "../lib/align.h"

struct mm_tbuf_s {
	void *km;
	int rep_len, frag_gap;
};

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

void *mm_tbuf_get_km(mm_tbuf_t *b)
{
	return b->km;
}

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (int32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && (int32_t)(dreg[u]>>32) < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && (int32_t)(dreg[v]>>32) < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > (int32_t)(dreg[v]>>32)? s : dreg[v]>>32;
				int ee = e < (int32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
	int i, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

typedef struct {
	uint32_t n;
	uint32_t q_pos, q_span;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} mm_match_t;

static mm_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mm_idx_t *mi, const mm128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, uint64_t **mini_pos)
{
	int rep_st = 0, rep_en = 0, n_m;
	size_t i;
	mm_match_t *m;
	*n_mini_pos = 0;
	*mini_pos = (uint64_t*)kmalloc(km, mv->n * sizeof(uint64_t));
	m = (mm_match_t*)kmalloc(km, mv->n * sizeof(mm_match_t));
	for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
		const uint64_t *cr;
		mm128_t *p = &mv->a[i];
		uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
		int t;
		cr = mm_idx_get(mi, p->x>>8, &t);
		if (t >= max_occ) {
			int en = (q_pos >> 1) + 1, st = en - q_span;
			if (st > rep_en) {
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			mm_match_t *q = &m[n_m++];
			q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
			q->is_tandem = 0;
			if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
			if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
			*n_a += q->n;
			(*mini_pos)[(*n_mini_pos)++] = (uint64_t)q_span<<32 | q_pos>>1;
		}
	}
	*rep_len += rep_en - rep_st;
	*_n_m = n_m;
	return m;
}

static inline int skip_seed(int flag, uint64_t r, const mm_match_t *q, const char *qname, int qlen, const mm_idx_t *mi, int *is_self)
{
	*is_self = 0;
	if (qname && (flag & (MM_F_NO_DIAG|MM_F_NO_DUAL))) {
		const mm_idx_seq_t *s = &mi->seq[r>>32];
		int cmp;
		cmp = strcmp(qname, s->name);
		if ((flag&MM_F_NO_DIAG) && cmp == 0 && (int)s->len == qlen) {
			if ((uint32_t)r>>1 == (q->q_pos>>1)) return 1; // avoid the diagnonal anchors
			if ((r&1) == (q->q_pos&1)) *is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag&MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY|MM_F_REV_ONLY)) {
		if ((r&1) == (q->q_pos&1)) { // forward strand
			if (flag & MM_F_REV_ONLY) return 1;
		} else {
			if (flag & MM_F_FOR_ONLY) return 1;
		}
	}
	return 0;
}

static mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m, heap_size = 0;
	int64_t j, n_for = 0, n_rev = 0;
	mm_match_t *m;
	mm128_t *a, *heap;

	m = collect_matches(km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);

	heap = (mm128_t*)kmalloc(km, n_m * sizeof(mm128_t));
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));

	for (i = 0, heap_size = 0; i < n_m; ++i) {
		if (m[i].n > 0) {
			heap[heap_size].x = m[i].cr[0];
			heap[heap_size].y = (uint64_t)i<<32;
			++heap_size;
		}
	}
	ks_heapmake_heap(heap_size, heap);
	while (heap_size > 0) {
		mm_match_t *q = &m[heap->y>>32];
		mm128_t *p;
		uint64_t r = heap->x;
		int32_t is_self, rpos = (uint32_t)r >> 1;
		if (!skip_seed(opt->flag, r, q, qname, qlen, mi, &is_self)) {
			if ((r&1) == (q->q_pos&1)) { // forward strand
				p = &a[n_for++];
				p->x = (r&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else { // reverse strand
				p = &a[(*n_a) - (++n_rev)];
				p->x = 1ULL<<63 | (r&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
		// update the heap
		if ((uint32_t)heap->y < q->n - 1) {
			++heap[0].y;
			heap[0].x = m[heap[0].y>>32].cr[(uint32_t)heap[0].y];
		} else {
			heap[0] = heap[heap_size - 1];
			--heap_size;
		}
		ks_heapdown_heap(0, heap_size, heap);
	}
	kfree(km, m);
	kfree(km, heap);

	// reverse anchors on the reverse strand, as they are in the descending order
	for (j = 0; j < n_rev>>1; ++j) {
		mm128_t t = a[(*n_a) - 1 - j];
		a[(*n_a) - 1 - j] = a[(*n_a) - (n_rev - j)];
		a[(*n_a) - (n_rev - j)] = t;
	}
	if (*n_a > n_for + n_rev) {
		memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mm128_t));
		*n_a = n_for + n_rev;
	}
	return a;
}

static mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m;
	mm_match_t *m;
	mm128_t *a;
	m = collect_matches(km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mm_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		uint32_t k;
		for (k = 0; k < q->n; ++k) {
			int32_t is_self, rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self)) continue;
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) { // forward strand
				p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else { // reverse strand
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
	}
	kfree(km, m);
	radix_sort_128x(a, a + (*n_a));
	return a;
}

static void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL);
		if (n_segs <= 1) mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		else mm_select_sub_multi(km, opt->pri_ratio, 0.2f, 0.7f, max_chain_gap_ref, mi->k*2, opt->best_n, n_segs, qlens, n_regs, regs);
		if (!(opt->flag & (MM_F_SPLICE|MM_F_SR|MM_F_NO_LJOIN))) // long join not working well without primary chains
			mm_join_long(km, opt, qlen, n_regs, regs, a);
	}
}
//////extend/////
//homopolymer-compress SSR to min_hlen
int32_t homo_compres_seq(char *seq, int32_t len, int min_hlen, char *cseq){
	int32_t i, j, h;
	char c, p = '\0', *o = cseq == NULL ? seq : cseq;
	for (i = j = h = 0; i < len; i ++){
		c = seq[i];
		if (c == p) h ++;
		else h = 0;
		if (h < min_hlen) o[j++] = c;
		p = c;
	}
	o[j] = '\0';
	return j;
}

int nd_idx_getbseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, char *seq, int rev)
{	
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1; 
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st; 
	en1 = mi->seq[rid].offset + en;
	if (rev){
		for (i = st1; i < en1; ++i)
			seq[en1 - 1 - i] = "TGCA"[mm_seq4_get(mi->S, i)];
	}else{
		for (i = st1; i < en1; ++i)
			seq[i - st1] = "ACGT"[mm_seq4_get(mi->S, i)];
	}

	seq[i - st1] = '\0';
	return en - st;
}

uint64_t nd_idx_get_homo_compres_bseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, char *seq, int rev, int min_hlen)
{
	uint64_t i, j, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1; 
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st; 
	en1 = mi->seq[rid].offset + en;

	int h;
	char c, p = '\0';
	if (rev){
		for (j = 0, i = en1 - 1; i >= st1 && i < en1; --i){
			c = "TGCA"[mm_seq4_get(mi->S, i)];
			if (c == p) h ++;
			else h = 0;
			if (h < min_hlen) seq[j ++] = c;
			p = c;
		}
	}else{
		for (j = 0, i = st1; i < en1; ++i){
			c = "ACGT"[mm_seq4_get(mi->S, i)];
			if (c == p) h ++;
			else h = 0;
			if (h < min_hlen) seq[j ++] = c;
			p = c;
		}
	}
	seq[j] = '\0';
	return j;
}

static inline void nd_update_coors(mm_reg1_t *r, int32_t qlen, const mm128_t *a) 
{ 
	int32_t k = r->as, q_span = (int32_t)(a[k].y>>32&0xff);
	r->rs = (int32_t)a[k].x + 1 > q_span? (int32_t)a[k].x + 1 - q_span : 0; 
	r->re = (int32_t)a[k + r->cnt - 1].x + 1;
	if (!r->rev) {
		r->qs = (int32_t)a[k].y + 1 - q_span;
		r->qe = (int32_t)a[k + r->cnt - 1].y + 1;
	} else {
		r->qs = qlen - ((int32_t)a[k + r->cnt - 1].y + 1); 
		r->qe = qlen - ((int32_t)a[k].y + 1 - q_span);
	}   
}

static int nd_fix_bad_ends(mm_reg1_t *r, const mm128_t *a, int bw, int min_match)
{
	int32_t i, l, m;
	int32_t as = r->as;
	int32_t cnt = r->cnt;
	if (r->cnt < 3) return 0;
	m = l = a[r->as].y >> 32 & 0xff;//q_span
	for (i = r->as + 1; i < r->as + r->cnt - 1; ++i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i].y >> 32 & 0xff;
		if (a[i].y & MM_SEED_LONG_JOIN) break;
		lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
		lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) as = i;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r->mlen >> 1) break;
	}
	cnt = r->as + r->cnt - as;
	m = l = a[r->as + r->cnt - 1].y >> 32 & 0xff;
	for (i = r->as + r->cnt - 2; i > as; --i) {
		int32_t lq, lr, min, max;
		int32_t q_span = a[i+1].y >> 32 & 0xff;
		if (a[i+1].y & MM_SEED_LONG_JOIN) break;
		lr = (int32_t)a[i+1].x - (int32_t)a[i].x;
		lq = (int32_t)a[i+1].y - (int32_t)a[i].y;
		min = lr < lq? lr : lq;
		max = lr > lq? lr : lq;
		if (max - min > l >> 1) cnt = i + 1 - as;
		l += min;
		m += min < q_span? min : q_span;
		if (l >= bw << 1 || (m >= min_match && m >= bw) || m >= r->mlen >> 1) break;
	}
	int ret = r->as != as || r->cnt != cnt ? 1 : 0;
	r->as = as;
	r->cnt = cnt;
	return ret;
}

static inline int check_realign_nextdenovo(const int rev, const uint32_t qs,
	const uint32_t qe, const uint32_t qlen, const uint32_t ts, const uint32_t te,
	const uint32_t tlen, const int32_t maxhan1, const int32_t maxhan2);

static void nd_extend_ends(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, 
	const char *seq, int *n_regs, mm_reg1_t *regs, int *V, uint8_t **D, int max_mem_d, char *skip_name){
	
	int i;
	int min_clen = 10;
	int max_tseq_len = qlen * 2;
	char *tseq = (char*)kmalloc(km, max_tseq_len);
	int32_t bstx, bsty, max_d, band_size, minlen;

	for (i = 0; i < *n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		if (skip_name && strcmp(skip_name, mi->seq[r->rid].name) == 0) continue;
		if (opt->dvt && !check_realign_nextdenovo(r->rev, r->qs, r->qe, qlen, 
			r->rs, r->re, mi->seq[r->rid].len, opt->maxhan1 * 3, opt->maxhan2 * 3)) continue;//NB:need rewrite

		if (mi->seq[r->rid].len >= max_tseq_len){
			max_tseq_len = mi->seq[r->rid].len;
			tseq = (char*)krealloc(km, tseq, max_tseq_len);
		}
		// NB:should already called nd_fix_bad_ends && nd_update_coors
		// int ret = nd_fix_bad_ends(r, a, opt->bw, opt->min_chain_score * 2);
		// if (ret) nd_update_coors(r, qlen, a);
		if (r->rev == 0){
			int subtlen = r->rs;
			int subqlen = r->qs;
			minlen = subtlen > subqlen ? subqlen : subtlen;
			if (minlen >= min_clen){
				clean_V(V, max_mem_d);
				max_d = minlen/4 > max_mem_d ? max_mem_d : (minlen > 20 ? minlen/4 : minlen);
				band_size = 500;//max_d/3 > 500 ? 500 : max_d/3;//500;TODO: check
				if (subtlen > (minlen << 1)){//only extract partial un-aligned bseqs
					subtlen = minlen << 1;
					nd_idx_getbseq(mi, r->rid, r->rs - (minlen << 1), r->rs, tseq, 0);
				}else{
					nd_idx_getbseq(mi, r->rid, 0, r->rs, tseq, 0);
				}
				extend_rev(seq, subqlen, tseq, subtlen, V, D, max_d, band_size, opt->d_factor, &bstx, &bsty);
				r->qs -= bstx;
				r->rs -= bsty;
			}

			subtlen = mi->seq[r->rid].len - r->re;
			subqlen = qlen - r->qe;
			minlen = subtlen > subqlen ? subqlen : subtlen;
			if (minlen >= min_clen){
				clean_V(V, max_mem_d);
				max_d = minlen/4 > max_mem_d ? max_mem_d : (minlen > 20 ? minlen/4 : minlen);
				band_size = 500;//max_d/3 > 500 ? 500 : max_d/3;//500;
				if (subtlen > (minlen << 1)){//only extract partial un-aligned bseqs
					subtlen = minlen << 1;
					nd_idx_getbseq(mi, r->rid, r->re, r->re + (minlen << 1), tseq, 0);
				}else{
					nd_idx_getbseq(mi, r->rid, r->re, mi->seq[r->rid].len, tseq, 0);
				}
				extend_fwd(seq + r->qe, subqlen, tseq, subtlen, V, D, max_d, band_size, opt->d_factor, &bstx, &bsty);
				r->qe += bstx;
				r->re += bsty;
			}
		}else{
			int subtlen = mi->seq[r->rid].len - r->re;
			int subqlen = r->qs;
			minlen = subtlen > subqlen ? subqlen : subtlen;
			if (minlen >= min_clen){
				clean_V(V, max_mem_d);
				max_d = minlen/4 > max_mem_d ? max_mem_d : (minlen > 20 ? minlen/4 : minlen);
				band_size = 500;//max_d/3 > 500 ? 500 : max_d/3;//500;
				if (subtlen > (minlen << 1)){//only extract partial un-aligned bseqs
					subtlen = minlen << 1;
					nd_idx_getbseq(mi, r->rid, r->re, r->re + (minlen << 1), tseq, 1);
				}else{
					nd_idx_getbseq(mi, r->rid, r->re, mi->seq[r->rid].len, tseq, 1);
				}
				extend_rev(seq, subqlen, tseq, subtlen, V, D, max_d, band_size, opt->d_factor, &bstx, &bsty);
				r->qs -= bstx;
				r->re += bsty;
			}

			subtlen = r->rs;
			subqlen = qlen - r->qe;
			minlen = subtlen > subqlen ? subqlen : subtlen;
			if (minlen >= min_clen){
				clean_V(V, max_mem_d);
				max_d = minlen/4 > max_mem_d ? max_mem_d : (minlen > 20 ? minlen/4 : minlen);
				band_size = 500;//max_d/3 > 500 ? 500 : max_d/3;//500;
				if (subtlen > (minlen << 1)){//only extract partial un-aligned bseqs
					subtlen = minlen << 1;
					nd_idx_getbseq(mi, r->rid, r->rs - (minlen << 1), r->rs, tseq, 1);
				}else{
					nd_idx_getbseq(mi, r->rid, 0, r->rs, tseq, 1);
				}
				extend_fwd(seq + r->qe, subqlen, tseq, subtlen, V, D, max_d, band_size, opt->d_factor, &bstx, &bsty);
				r->qe += bstx;
				r->rs -= bsty;
			}
		}
	}
	kfree (km, tseq);
}
//////end/////
static mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, 
	int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_CIGAR)) {
		if (opt->mode == 3) {//fix bad end for hifi
			int i;
			for (i = 0; i < *n_regs; ++i) {
				mm_reg1_t *r = &regs[i];
				int ret = nd_fix_bad_ends(r, a, opt->bw, opt->min_chain_score * 2);
				if (ret) nd_update_coors(r, qlen, a);
			}
		}
		return regs;
	}
	regs = mm_align_skeleton(km, opt, mi, qlen, seq, n_regs, regs, a); // this calls mm_filter_regs()
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}
	return regs;
}

void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	int i, j, rep_len, qlen_sum, n_regs0, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u, *mini_pos;
	mm128_t *a;
	mm128_v mv = {0,0,0};
	mm_reg1_t *regs0;
	km_stat_t kmst;

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	hash  = qname? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	else a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);

	if (opt->max_occ > opt->mid_occ && rep_len > 0) {
		int rechain = 0;
		if (n_regs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < n_regs0; ++i) { // find the best chain
				if (max < (int)(u[i]>>32)) max = u[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)u[i];
			}
			for (i = 1; i < (int32_t)u[max_i]; ++i) // count the number of segments in the best chain
				if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < n_segs)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(b->km, a);
			kfree(b->km, u);
			kfree(b->km, mini_pos);
			if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}
	b->frag_gap = max_chain_gap_ref;
	b->rep_len = rep_len;

	regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED)
		for (j = 0; j < n_regs0; ++j)
			for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
				fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
	if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);

	if (n_segs == 1) { // uni-segment
		regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], &n_regs0, regs0, a);
		mm_set_mapq(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
		n_regs[0] = n_regs0, regs[0] = regs0;
	} else { // multi-segment
		mm_seg_t *seg;
		seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
		free(regs0);
		for (i = 0; i < n_segs; ++i) {
			mm_set_parent(b->km, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a);
			mm_set_mapq(b->km, n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
		}
		mm_seg_free(b->km, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
	}

	kfree(b->km, mv.a);
	kfree(b->km, a);
	kfree(b->km, u);
	kfree(b->km, mini_pos);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag(mi, 1, &qlen, &seq, n_regs, &regs, b, opt, qname);
	return regs;
}

void mm_map_frag_nextdenovo1(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	int i, j, rep_len, qlen_sum, n_regs0, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	uint32_t hash;
	int64_t n_a;
	uint64_t *u, *mini_pos;
	mm128_t *a;
	mm128_v mv = {0,0,0};
	mm_reg1_t *regs0;
	km_stat_t kmst;

	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	hash  = qname? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	else a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
		fprintf(stdout, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			fprintf(stdout, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
	}

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	a = mm_chain_dp_nextdenovo(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);

	if (opt->max_occ > opt->mid_occ && rep_len > 0) {
		int rechain = 0;
		if (n_regs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < n_regs0; ++i) { // find the best chain
				if (max < (int)(u[i]>>32)) max = u[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)u[i];
			}
			for (i = 1; i < (int32_t)u[max_i]; ++i) // count the number of segments in the best chain
				if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < n_segs)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(b->km, a);
			kfree(b->km, u);
			kfree(b->km, mini_pos);
			if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}
	b->frag_gap = max_chain_gap_ref;
	b->rep_len = rep_len;

	regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED)
		for (j = 0; j < n_regs0; ++j)
			for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
				fprintf(stdout, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
	if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);

	if (n_segs == 1) { // uni-segment
		regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], &n_regs0, regs0, a);
		mm_set_mapq(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
		n_regs[0] = n_regs0, regs[0] = regs0;
	} else { // multi-segment
		mm_seg_t *seg;
		seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
		free(regs0);
		for (i = 0; i < n_segs; ++i) {
			mm_set_parent(b->km, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a);
			mm_set_mapq(b->km, n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
		}
		mm_seg_free(b->km, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
	}

	kfree(b->km, mv.a);
	kfree(b->km, a);
	kfree(b->km, u);
	kfree(b->km, mini_pos);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

mm_reg1_t *mm_map_nextdenovo1(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag_nextdenovo1(mi, 1, &qlen, &seq, n_regs, &regs, b, opt, qname);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int mini_batch_size, n_processed, n_threads, n_fp;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;

	int n_parts;
	uint32_t *rid_shift;
	FILE *fp_split, **fp_parts;
} pipeline_t;

typedef struct {
	uint32_t len;
	uint8_t *seq;
} tseq_nextdenovo;

KHASH_MAP_INIT_INT(rid, int);
typedef struct {
	const pipeline_t *p;
	int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
	tseq_nextdenovo *rseq;
	mm_idx_t **mi_new;
	khash_t(rid) **rids;
	int **V;
	uint8_t ***D;
} step_t;

static int cmpfunc_nextdenovo (const void *a, const void *b)
{
	return (((mm_reg1_t*)a)->rid - ((mm_reg1_t*)b)->rid);
}

static void clean_align_nextdenovo(mm_reg1_t *reg, int n_reg, int k){
	while (--n_reg >= 0){
		if (n_reg != k) reg[n_reg].mlen = 1;
	}
}

static inline int check_realign_nextdenovo(const int rev, const uint32_t qs,
	const uint32_t qe, const uint32_t qlen, const uint32_t ts, const uint32_t te,
	const uint32_t tlen, const int32_t maxhan1, const int32_t maxhan2){

	if (rev){
		if (qs <= maxhan1 && ts <= maxhan1) return 1;
		else if (qlen - qe <= maxhan1 && tlen - te <= maxhan1) return 2;
	}else{
		if (qlen - qe <= maxhan1 && ts <= maxhan1) return 4;
		else if (qs <= maxhan1 && tlen - te <= maxhan1) return 7;
	}

	if (maxhan2 > 0){
		if (qs <= maxhan2 && qe + maxhan2 >= qlen) return 8;
		if (ts <= maxhan2 && te + maxhan2 >= tlen) return 9;
	}
	return 0;
}

static int update_reg_nextdenovo(mm_reg1_t *reg_new, int n_reg_new, mm_reg1_t *reg,
		int n_reg, int s, int e, int t_l, mm_idx_seq_t *mi_new_seq, int32_t maxhan1, int32_t maxhan2){
	int i, c, t, l, pi;
	uint32_t alnlen;
	mm_reg1_t *r, *r_new;
	for (i = c = 0; s < e; s ++){
		r = &reg[s];
		if (r->mlen == 2){
			for ( l = -1, pi = i, alnlen = t = 0; i < n_reg_new && t < 10; i ++){
				r_new = &reg_new[i];
				if (r_new->rid == r->blen){
					if (l == -1){
						l = i;
						alnlen = reg_new[i].blen * 0.8;
						// for (k = i; k < n_reg_new && k < i + 5; k++){
						// 	if (reg_new[k].rid == r->blen && reg_new[k].blen > alnlen){
						// 		l = k;
						// 		alnlen = reg_new[k].blen;
						// 	}
						// }
					}

					if (r_new->blen >= alnlen && check_realign_nextdenovo(r_new->rev, r_new->qs,
						r_new->qe, t_l, r_new->rs, r_new->re, mi_new_seq[r_new->rid].len, maxhan1, maxhan2)){
						l = i;
						break;//TODO SET 0.8 as an option
					}
					// printf("3#%d %d %d %d %d %d %d %d\n", r->blen,l, reg_new[l].rev, reg_new[l].qs, reg_new[l].qe, t_l, reg_new[l].rs,reg_new[l].re );

					t ++;
					if (r_new->qs <= maxhan2 && r_new->qe + maxhan2 >= t_l) {
						l = i;
						c++;
						break;
					}
				}else if (l >= 0) {
					i --;
					break;
				}
			}
			if (l >= 0){
		// printf("4#%d %d %d %d %d %d %d %d\n", r->blen,l, reg_new[l].rev, reg_new[l].qs, reg_new[l].qe, t_l, reg_new[l].rs,reg_new[l].re );
				r_new = &reg_new[l];
				r->rev = r_new->rev;
				r->qs = r_new->qs;
				r->qe = r_new->qe;
				r->rs = r_new->rs;
				r->re = r_new->re;
				r->mlen = r_new->mlen;
				r->blen = r_new->blen;
			}else i = pi;
		}
	}
	return c;
}

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
	step_t *s = (step_t*)_data;
	int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MM_MAX_SEG];
	mm_tbuf_t *b = s->buf[tid];
	assert(s->n_seg[i] <= MM_MAX_SEG);
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
			mm_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
	}
	if (s->p->opt->flag & MM_F_INDEPEND_SEG) {
		for (j = 0; j < s->n_seg[i]; ++j) {
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
			s->rep_len[off + j] = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
		for (j = 0; j < s->n_seg[i]; ++j) {
			s->rep_len[off + j] = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	}
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
			int k, t;
			mm_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mm_reg1_t *r = &s->reg[off + j][k];
				t = r->qs;
				r->qs = qlens[j] - r->qe;
				r->qe = qlens[j] - t;
				r->rev = !r->rev;
			}
		}
	if (s->p->opt->mode == 3 && !(s->p->opt->flag & MM_F_CIGAR)){
		assert (s->n_seg[i] == 1);
		uint32_t k, alnlen, qalen,  talen;
		mm_reg1_t *r;
		mm_bseq1_t *t = &s->seq[off];
		const mm_idx_t *mi = s->p->mi;
		tseq_nextdenovo *rseq = &s->rseq[tid];
		str_toupper(t->seq);
		nd_extend_ends(s->p->opt, mi, b->km, t->l_seq, t->seq, &s->n_reg[off], 
			s->reg[off], s->V[tid], s->D[tid], s->p->opt->ide_ml, t->name);
		if (s->p->opt->step == 1) return;
		char *htseq = (char*)kmalloc(b->km, t->l_seq);
		alignpos apos = {0};
		for (k = 0; k < s->n_reg[off]; ++k) {
			r = &s->reg[off][k];
			if (strcmp(t->name, mi->seq[r->rid].name) == 0) {r->mlen = 0; continue;}
			qalen = r->qe - r->qs;
			talen = r->re - r->rs;
			alnlen = qalen > talen ? qalen : talen;
			if (alnlen >= s->p->opt->minlen && qalen >= t->l_seq/50 && talen >= mi->seq[r->rid].len/50 && 
					check_realign_nextdenovo(r->rev, r->qs, r->qe, t->l_seq, r->rs, r->re, 
					mi->seq[r->rid].len, s->p->opt->maxhan1, s->p->opt->maxhan1)){
				clean_V(s->V[tid], s->p->opt->ide_ml);
				int max_d = alnlen/5 > s->p->opt->ide_ml ? s->p->opt->ide_ml : alnlen/5;
				int band_size = max_d > 1500 ? 500 : max_d/3;
				if (talen > rseq->len){
					rseq->len = talen;
					rseq->seq = (uint8_t *) realloc(rseq->seq, talen + 1);
				}
				talen = nd_idx_get_homo_compres_bseq(mi, r->rid, r->rs, r->re, (char *)rseq->seq, r->rev, 3);
				qalen = homo_compres_seq(t->seq + r->qs, qalen, 3, htseq);
				// ide(htseq, qalen, (char *)rseq->seq, talen,
				// 	s->V[tid], s->D[tid], max_d, band_size, &r->mlen, &r->blen);
				// alignpos apos = {0};
				// clean_V(s->V[tid], s->p->opt->ide_ml);
				apos.aln_len = 0;
				alnpos(htseq, qalen, (char *)rseq->seq, talen,
					s->V[tid], s->D[tid], max_d, band_size, &apos);
				// printf("qalen:%d talen:%d mlen:%d qlen:%d\n", qalen, talen, r->mlen, r->blen);
				// printf("apos-> aln_len:%d aln_mlen:%d aln_q_s:%d aln_q_e:%d aln_t_s:%d aln_t_e:%d %s %s\n",
				// 	apos.aln_len, apos.aln_mlen, apos.aln_q_s, apos.aln_q_e, apos.aln_t_s, apos.aln_t_e,
				// 	mi->seq[r->rid].name, t->name);
				if (apos.aln_len) {
					r->mlen = apos.aln_mlen;
					r->blen = apos.aln_len;
					if (qalen - apos.aln_q_e <= 10 && talen - apos.aln_t_e <= 10){
						r->qs += apos.aln_q_s;
						r->qe -= qalen - apos.aln_q_e;
						if (r->rev){
							r->rs += talen - apos.aln_t_e;
							r->re -= apos.aln_t_s;
						}else{
							r->rs += apos.aln_t_s;
							r->re -= talen - apos.aln_t_e;
						}
					}
				}
				// alignment aln = {0};
				// aln.q_aln_str = malloc(65536);
				// aln.t_aln_str = malloc(65536);
				// clean_V(s->V[tid], s->p->opt->ide_ml);
				// align(htseq, qalen, (char *)rseq->seq, talen,&aln, s->V[tid], s->D[tid]);
				// printf("%s\n", aln.q_aln_str);
				// printf("%s\n", aln.t_aln_str);
				// free(aln.q_aln_str);
				// free(aln.t_aln_str);
			}else r->mlen = 1;
		}
		kfree(b->km, htseq);
	}else if (s->p->opt->step == 2 && !(s->p->opt->flag & MM_F_CIGAR)){
		khint_t rk;
		int absent;
		uint32_t l, tp, alnlen;
		int k, c, n_reg_new, seq_index;
		mm_reg1_t *r, *reg_new;
		mm_bseq1_t *t = &s->seq[off];
		mm_idx_t *mi_new = s->mi_new[tid];
		tseq_nextdenovo *rseq = &s->rseq[tid];
		khash_t(rid) *h = s->rids[tid];
		const mm_idx_t *mi = s->p->mi;
		for (j = 0; j < s->n_seg[i]; ++j){
			c = 0;
			kh_clear(rid, h);
			mm_idx_clean_nextdenovo(mi_new, s->p->opt->cn);
			for (seq_index = k = 0; k < s->n_reg[off + j]; ++k) {
				r = &s->reg[off + j][k];
				if (strcmp(t->name, mi->seq[r->rid].name) == 0) {r->mlen = 0; continue;}
				rk = kh_put(rid, h, r->rid, &absent);
				tp = r->mlen; 
				if (absent) {kh_key(h, rk) = r->rid; kh_val(h, rk) = k;}
				else r->mlen = 1;
				l = kh_val(h, rk);
				alnlen = mi->seq[r->rid].len > t->l_seq ? mi->seq[r->rid].len : t->l_seq;
				if (l != k && (s->reg[off + j][l].mlen == 2 || r->blen < s->reg[off + j][l].blen * 0.8 || r->blen < alnlen/3)) continue;
				else if (r->qe - r->qs >= s->p->opt->minlen && tp >= r->blen * s->p->opt->minide && \
						tp >= s->p->opt->minmatch){
					if (check_realign_nextdenovo(r->rev, r->qs, r->qe, t->l_seq, r->rs, r->re, \
						mi->seq[r->rid].len, s->p->opt->maxhan1, 0)){
						if (s->reg[off + j][l].mlen == 3) c--;
						s->reg[off + j][l].mlen = s->p->opt->mode ? 2 : tp;
						seq_index ++;
					}else if (r->qs <= s->p->opt->maxhan2 && r->qe + s->p->opt->maxhan2 >= t->l_seq) {
						s->reg[off + j][l].mlen = 3;
						if (++c >= MAX_CON) break;
					}else if (s->p->opt->mode == 3 && s->reg[off + j][l].mlen != 3 && (
						(r->qs <= s->p->opt->maxhan1 && r->qe + s->p->opt->maxhan1 >= t->l_seq) ||
						(r->rs <= s->p->opt->maxhan1 && r->re + s->p->opt->maxhan1 >= mi->seq[r->rid].len))) {
						s->reg[off + j][l].mlen = 2;
						seq_index ++;
					}
				}
			}
			if (c >= MAX_CON || !s->p->opt->mode) continue;
			if ((s->p->opt->mode == 2 && seq_index < 200) || seq_index < 20){
				mm_idx_str_nextdenovo3(mi_new, t->seq, t->l_seq);
				for (k = 0; k < s->n_reg[off + j]; ++k) {
					r = &s->reg[off + j][k];
					if (r->mlen != 2) continue;
					if (mi->seq[r->rid].len > rseq->len){
						rseq->len = mi->seq[r->rid].len;
						rseq->seq = (uint8_t *) realloc(rseq->seq, rseq->len + 1);
					}
					mm_idx_getseq(mi, r->rid, 0, mi->seq[r->rid].len, rseq->seq);
					for (l = 0; l < mi->seq[r->rid].len; l ++) rseq->seq[l] = "ACGT"[rseq->seq[l]];
					rseq->seq[l] = '\0';
					if (s->p->opt->mode == 2){
						reg_new = mm_map(mi_new, mi->seq[r->rid].len, (char *)(rseq->seq), &n_reg_new, b, s->p->opt, 0);
					}else{
						reg_new = mm_map_nextdenovo1(mi_new, mi->seq[r->rid].len, (char *)(rseq->seq), &n_reg_new, b, s->p->opt, 0);
					}
					if (reg_new){
						alnlen = reg_new[0].blen * 0.8;
						
						// for (l = 0 ; l < n_reg_new; l++){
						// 	printf("1#%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
						// 		 mi->seq[r->rid].name, reg_new[l].rev, mi->seq[r->rid].len, 
						// 		 reg_new[l].qs, reg_new[l].qe, t->name, t->l_seq, reg_new[l].rs, 
						// 		 reg_new[l].re, reg_new[l].mlen, reg_new[l].blen, reg_new[l].score);
						// }

						for (l = 0 ; l < n_reg_new && l < 10; l++){
						// 	printf("1#%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
						// 		 mi->seq[r->rid].name, reg_new[l].rev, mi->seq[r->rid].len, 
						// 		 reg_new[l].qs, reg_new[l].qe, t->name, t->l_seq, reg_new[l].rs, 
						// 		 reg_new[l].re, reg_new[l].mlen, reg_new[l].blen, reg_new[l].score);
							if (reg_new[l].blen >= alnlen && check_realign_nextdenovo(reg_new[l].rev, reg_new[l].qs,
								reg_new[l].qe, mi->seq[r->rid].len, reg_new[l].rs, reg_new[l].re,t->l_seq, 
								s->p->opt->maxhan1, s->p->opt->maxhan2)) break;
						}
						if (l == 10 || l == n_reg_new) l = 0;
						tp = reg_new[l].qs; reg_new[l].qs = reg_new[l].rs; reg_new[l].rs = tp;
						tp = reg_new[l].qe; reg_new[l].qe = reg_new[l].re; reg_new[l].re = tp;
						reg_new[l].rid = r->rid;
						s->reg[off + j][k] = reg_new[l];
						free (reg_new);
						if (r->qs <= s->p->opt->maxhan2 && r->qe + s->p->opt->maxhan2 >= t->l_seq){
							if (++c >= MAX_CON) break;
						}
					}
				}
			}else{

				int cn = (int)((float) seq_index / ((seq_index + s->p->opt->cn - 1)/s->p->opt->cn) + 0.999);
				for (tp = seq_index = k = 0; k < s->n_reg[off + j]; ++k) {
					r = &s->reg[off + j][k];
					if (r->mlen != 2) continue;
					r->blen = seq_index;
					if (mi->seq[r->rid].len > rseq->len){
						rseq->len = mi->seq[r->rid].len;
						rseq->seq = (uint8_t *) realloc(rseq->seq, rseq->len + 1);
					}
					mm_idx_getseq(mi, r->rid, 0, mi->seq[r->rid].len, rseq->seq);
					mm_idx_str_nextdenovo2(mi_new, mi->seq[r->rid].name, rseq->seq, mi->seq[r->rid].len, seq_index);
					if (++ seq_index >= cn){
						mm_idx_post_nextdenovo(mi_new, 1);
						reg_new = mm_map(mi_new, t->l_seq, t->seq, &n_reg_new, b, s->p->opt, 0);
						qsort(reg_new, n_reg_new, sizeof(mm_reg1_t), cmpfunc_nextdenovo);
						c += update_reg_nextdenovo(reg_new, n_reg_new, s->reg[off + j], s->n_reg[off + j],
							tp, k + 1, t->l_seq, mi_new->seq, s->p->opt->maxhan1, s->p->opt->maxhan2);
						//  for (l = 0 ; l < n_reg_new; l++){
						//  printf("2.1#%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
						//      t->name, t->l_seq, reg_new[l].rev, reg_new[l].qs, reg_new[l].qe, 
						//      mi_new->seq[reg_new[l].rid].name, mi_new->seq[reg_new[l].rid].len, reg_new[l].rs, 
						//      reg_new[l].re, reg_new[l].mlen, reg_new[l].blen, reg_new[l].score);
						// };
						free (reg_new);
						seq_index = 0;
						tp = k + 1;
						mm_idx_clean_nextdenovo(mi_new, s->p->opt->cn);
						if (c >= MAX_CON) break;
					}
				}
				if (c < MAX_CON && seq_index){
					mm_idx_post_nextdenovo(mi_new, 1);
					reg_new = mm_map(mi_new, t->l_seq, t->seq, &n_reg_new, b, s->p->opt, 0);
					qsort(reg_new, n_reg_new, sizeof(mm_reg1_t), cmpfunc_nextdenovo);
					update_reg_nextdenovo(reg_new, n_reg_new, s->reg[off + j], s->n_reg[off + j],
						tp, k, t->l_seq, mi_new->seq, s->p->opt->maxhan1, s->p->opt->maxhan2);
					//  for (l = 0 ; l < n_reg_new; l++){
					//  printf("2.2#%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
					//      t->name, t->l_seq, reg_new[l].rev, reg_new[l].qs, reg_new[l].qe, 
					//      mi_new->seq[reg_new[l].rid].name, mi_new->seq[reg_new[l].rid].len, reg_new[l].rs, 
					//      reg_new[l].re, reg_new[l].mlen, reg_new[l].blen, reg_new[l].score);
					// };
					free (reg_new);
				}
			}
		}
		mm_idx_clean_nextdenovo(mi_new, s->p->opt->cn);
	}
}

static void merge_hits(step_t *s)
{
	int f, i, k0, k, max_seg = 0, *n_reg_part, *rep_len_part, *frag_gap_part, *qlens;
	void *km;
	FILE **fp = s->p->fp_parts;
	const mm_mapopt_t *opt = s->p->opt;

	km = km_init();
	for (f = 0; f < s->n_frag; ++f)
		max_seg = max_seg > s->n_seg[f]? max_seg : s->n_seg[f];
	qlens = CALLOC(int, max_seg + s->p->n_parts * 3);
	n_reg_part = qlens + max_seg;
	rep_len_part = n_reg_part + s->p->n_parts;
	frag_gap_part = rep_len_part + s->p->n_parts;
	for (f = 0, k = k0 = 0; f < s->n_frag; ++f) {
		k0 = k;
		for (i = 0; i < s->n_seg[f]; ++i, ++k) {
			int j, l, t, rep_len = 0;
			qlens[i] = s->seq[k].l_seq;
			for (j = 0, s->n_reg[k] = 0; j < s->p->n_parts; ++j) {
				mm_err_fread(&n_reg_part[j],    sizeof(int), 1, fp[j]);
				mm_err_fread(&rep_len_part[j],  sizeof(int), 1, fp[j]);
				mm_err_fread(&frag_gap_part[j], sizeof(int), 1, fp[j]);
				s->n_reg[k] += n_reg_part[j];
				if (rep_len < rep_len_part[j])
					rep_len = rep_len_part[j];
			}
			s->reg[k] = CALLOC(mm_reg1_t, s->n_reg[k]);
			for (j = 0, l = 0; j < s->p->n_parts; ++j) {
				for (t = 0; t < n_reg_part[j]; ++t, ++l) {
					mm_reg1_t *r = &s->reg[k][l];
					uint32_t capacity;
					mm_err_fread(r, sizeof(mm_reg1_t), 1, fp[j]);
					r->rid += s->p->rid_shift[j];
					if (opt->flag & MM_F_CIGAR) {
						mm_err_fread(&capacity, 4, 1, fp[j]);
						r->p = (mm_extra_t*)calloc(capacity, 4);
						r->p->capacity = capacity;
						mm_err_fread(r->p, r->p->capacity, 4, fp[j]);
					}
				}
			}
			mm_hit_sort(km, &s->n_reg[k], s->reg[k]);
			mm_set_parent(km, opt->mask_level, s->n_reg[k], s->reg[k], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL);
			if (!(opt->flag & MM_F_ALL_CHAINS)) {
				mm_select_sub(km, opt->pri_ratio, s->p->mi->k*2, opt->best_n, &s->n_reg[k], s->reg[k]);
				mm_set_sam_pri(s->n_reg[k], s->reg[k]);
			}
			mm_set_mapq(km, s->n_reg[k], s->reg[k], opt->min_chain_score, opt->a, rep_len, !!(opt->flag & MM_F_SR));
		}
		if (s->n_seg[f] == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(km, frag_gap_part[0], opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, &s->n_reg[k0], &s->reg[k0]);
	}
	free(qlens);
	km_destroy(km);
}

static void *worker_pipeline(void *shared, int step, void *in)
{	
	int i, j, k;
	pipeline_t *p = (pipeline_t*)shared;
	if (step == 0) { // step 0: read sequences
		int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
		int with_comment = !!(p->opt->flag & MM_F_COPY_COMMENT);
		int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
		step_t *s;
		s = (step_t*)calloc(1, sizeof(step_t));
		if (p->n_fp > 1) s->seq = mm_bseq_read_frag2(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mm_bseq_read3(p->fp[0], p->mini_batch_size, with_qual, with_comment, frag_mode, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			if (s->p->opt->step == 2 || s->p->opt->mode == 3){
				s->rseq = (tseq_nextdenovo*) malloc(p->n_threads * sizeof(tseq_nextdenovo));
				s->mi_new = (mm_idx_t**) malloc(p->n_threads * sizeof(mm_idx_t *));
				s->rids = (khash_t(rid) **) malloc(p->n_threads * sizeof(khash_t(rid) *));
				if (s->p->opt->mode == 3){
					s->V = (int **) malloc(p->n_threads * sizeof(int *));
					s->D = (uint8_t ***) malloc(p->n_threads * sizeof(uint8_t **));
				}
			}
			for (i = 0; i < p->n_threads; ++i){
				s->buf[i] = mm_tbuf_init();
				if (s->mi_new){
					s->rseq[i].len = 150000;
					s->rseq[i].seq = (uint8_t *) malloc(s->rseq[i].len + 1);
					s->mi_new[i] = mm_idx_init_nextdenovo(s->p->opt->wn, s->p->opt->kn, s->p->mi->flag & MM_I_HPC, -1, s->p->opt->cn + 1);
					s->rids[i] = kh_init(rid);
					if (s->V) {
						malloc_vd(&s->V[i], &s->D[i], s->p->opt->ide_ml);
						if (!s->V[i] || !s->D[i]){ fprintf(stderr, "memory out!"); exit(EXIT_FAILURE);}
					}
				}
			}
			s->n_reg = (int*)calloc(5 * s->n_seq, sizeof(int));
			s->seg_off = s->n_reg + s->n_seq; // seg_off, n_seg, rep_len and frag_gap are allocated together with n_reg
			s->n_seg = s->seg_off + s->n_seq;
			s->rep_len = s->n_seg + s->n_seq;
			s->frag_gap = s->rep_len + s->n_seq;
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
			for (i = 1, j = 0; i <= s->n_seq; ++i)
				if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					j = i;
				}
			return s;
		} else free(s);
	} else if (step == 1) { // step 1: map
		if (p->n_parts > 0) merge_hits((step_t*)in);
		else kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_frag);
		return in;
	} else if (step == 2) { // step 2: output
		extern bytes_t *encode_tbl;
		extern prev_t pid;
		void *km = 0;
		step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) {
			mm_tbuf_destroy(s->buf[i]); 
			if (s->mi_new) {
				free(s->rseq[i].seq);
				mm_idx_destroy(s->mi_new[i]);
				kh_destroy(rid, s->rids[i]);
				if (s->V){destory_vd(s->V[i], s->D[i]);}
			}
		}
		free(s->buf);
		if (s->mi_new) { 
			free(s->mi_new); 
			free(s->rseq); 
			free (s->rids); 
			if (s->V){free(s->V); free(s->D);}
		}
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		buffer_t *buf = init_buffer(ECBUFSIZE);
		for (k = 0; k < s->n_frag; ++k) {
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				if (p->opt->split_prefix && p->n_parts == 0) { // then write to temporary files
					mm_err_fwrite(&s->n_reg[i],    sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->rep_len[i],  sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->frag_gap[i], sizeof(int), 1, p->fp_split);
					for (j = 0; j < s->n_reg[i]; ++j) {
						mm_reg1_t *r = &s->reg[i][j];
						mm_err_fwrite(r, sizeof(mm_reg1_t), 1, p->fp_split);
						if (p->opt->flag & MM_F_CIGAR) {
							mm_err_fwrite(&r->p->capacity, 4, 1, p->fp_split);
							mm_err_fwrite(r->p, r->p->capacity, 4, p->fp_split);
						}
					}
				} else if (s->n_reg[i] > 0) { // the query has at least one hit
					for (j = 0; j < s->n_reg[i]; ++j) {
						mm_reg1_t *r = &s->reg[i][j];
						assert(!r->sam_pri || r->id == r->parent);
						if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
							continue;
						if (p->opt->flag & MM_F_OUT_SAM){
							mm_write_sam3(&p->str, mi, t, i - seg_st, j, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag, s->rep_len[i]);
							mm_err_puts(p->str.s);
						}else if (p->opt->outraw||!p->opt->step){
							mm_write_paf3(&p->str, mi, t, r, km, p->opt->flag, s->rep_len[i]);
							mm_err_puts(p->str.s);
						}else if (strcmp(t->name, mi->seq[r->rid].name) == 0){//TODO up up
							continue;
						}else if (p->opt->step == 1 && r->qe - r->qs >= p->opt->minlen){
							overlap ovl = { r->rev, strtoul(t->name, NULL, 10), r->qs, r->qe, \
								strtoul(mi->seq[r->rid].name, NULL, 10), r->rs, r->re, r->mlen};
							if ((!p->opt->dvt) || check_realign_nextdenovo(r->rev, r->qs, r->qe, t->l_seq, 
									r->rs, r->re, mi->seq[r->rid].len, p->opt->maxhan1, p->opt->maxhan2)) 
								encode_ovl(stdout, encode_tbl, &pid, &ovl, buf);
						}else if (p->opt->step == 2 && (r->qe - r->qs >= p->opt->minlen || r->mlen == r->blen) && \
								(r->mlen == 3 || (r->mlen >= r->blen * p->opt->minide && r->mlen >= p->opt->minmatch)) && \
								// r->mlen >= r->blen * p->opt->minide && r->mlen >= p->opt->minmatch && 
								r->blen >= t->l_seq/50 && r->blen >= mi->seq[r->rid].len/50){
							overlap_i ovl = {r->rev, strtoul(t->name, NULL, 10), r->qs, r->qe, t->l_seq, \
								strtoul(mi->seq[r->rid].name, NULL, 10), r->rs, r->re, mi->seq[r->rid].len, (uint64_t) r->mlen * 10000/r->blen};
							int lable = filter_ovl(&ovl, (khash_t(ovlh_) *)p->opt->os, p->opt->maxhan1, p->opt->maxhan2);
							// fprintf(stderr, "%u %u %u %u %u %u %u %u %u\n",ovl.qname, ovl.qlen, ovl.qs, ovl.qe, ovl.tname, ovl.tlen, ovl.ts, ovl.te, lable);
							if (lable) encode_ovl_i(stdout, encode_tbl, &pid, &ovl, buf);
							else if (p->opt->outctn){
								khash_t(ovlh_) *os = (khash_t(ovlh_) *)p->opt->os;
								if (ovl.qs <= p->opt->maxhan2 && ovl.qe + p->opt->maxhan2 >= ovl.qlen){
									khint_t k = kh_get(ovlh_, os, ovl.qname);
									ovlinfo_aln *loli = &kh_val(os, k);
									loli->con = 0;
									encode_ovl_i(stdout, encode_tbl, &pid, &ovl, buf);
								}else if(ovl.ts <= p->opt->maxhan2 && ovl.te + p->opt->maxhan2 >= ovl.tlen){
									khint_t k = kh_get(ovlh_, os, ovl.tname);
									ovlinfo_aln *roli = &kh_val(os, k);
									roli->con = 0;
									encode_ovl_i(stdout, encode_tbl, &pid, &ovl, buf);
								}
							}
						}
					}
				} else if ((p->opt->flag & MM_F_PAF_NO_HIT) || ((p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_SAM_HIT_ONLY))) { // output an empty hit, if requested
					if (p->opt->flag & MM_F_OUT_SAM)
						mm_write_sam3(&p->str, mi, t, i - seg_st, -1, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag, s->rep_len[i]);
					else
						mm_write_paf3(&p->str, mi, t, 0, 0, p->opt->flag, s->rep_len[i]);
					mm_err_puts(p->str.s);
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
				free(s->reg[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
				if (s->seq[i].comment) free(s->seq[i].comment);
			}
		}
		flush_buffer(stdout, buf);
		free(s->reg); free(s->n_reg); free(s->seq); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory leak here
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
	return 0;
}

static mm_bseq_file_t **open_bseqs(int n, const char **fn)
{
	mm_bseq_file_t **fp;
	int i, j;
	fp = (mm_bseq_file_t**)calloc(n, sizeof(mm_bseq_file_t*));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = mm_bseq_open(fn[i])) == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s'\n", fn[i]);
			for (j = 0; j < i; ++j)
				mm_bseq_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads)
{
	int i, pl_threads;
	pipeline_t pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	if (opt->step == 3) set_minlen_nextdenovo(pl.fp, pl.n_fp, opt->minlen);
	if (opt->split_prefix)
		pl.fp_split = mm_split_init(opt->split_prefix, idx);
	pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);

	free(pl.str.s);
	if (pl.fp_split) fclose(pl.fp_split);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads)
{
	return mm_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int mm_split_merge(int n_segs, const char **fn, const mm_mapopt_t *opt, int n_split_idx)
{
	int i;
	pipeline_t pl;
	mm_idx_t *mi;
	if (n_segs < 1 || n_split_idx < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt;
	pl.mini_batch_size = opt->mini_batch_size;

	pl.n_parts = n_split_idx;
	pl.fp_parts  = CALLOC(FILE*, pl.n_parts);
	pl.rid_shift = CALLOC(uint32_t, pl.n_parts);
	pl.mi = mi = mm_split_merge_prep(opt->split_prefix, n_split_idx, pl.fp_parts, pl.rid_shift);
	if (pl.mi == 0) {
		free(pl.fp_parts);
		free(pl.rid_shift);
		return -1;
	}
	for (i = n_split_idx - 1; i > 0; --i)
		pl.rid_shift[i] = pl.rid_shift[i - 1];
	for (pl.rid_shift[0] = 0, i = 1; i < n_split_idx; ++i)
		pl.rid_shift[i] += pl.rid_shift[i - 1];
	if (opt->flag & MM_F_OUT_SAM)
		for (i = 0; i < (int32_t)pl.mi->n_seq; ++i)
			printf("@SQ\tSN:%s\tLN:%d\n", pl.mi->seq[i].name, pl.mi->seq[i].len);

	kt_pipeline(2, worker_pipeline, &pl, 3);

	free(pl.str.s);
	mm_idx_destroy(mi);
	free(pl.rid_shift);
	for (i = 0; i < n_split_idx; ++i)
		fclose(pl.fp_parts[i]);
	free(pl.fp_parts);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	mm_split_rm_tmp(opt->split_prefix, n_split_idx);
	return 0;
}
