#include "out.h"
#include "asg.h"
#include "kit.h"
#include "../../lib/bseq.h"

static uint8_t nt_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static uint32_t re_cal_len(ctg_ *s, const idxs *idx){
	uint32_t i, t, len = 0;
	for (i = 0; i < s->i; i++){
		if (s->nd[i].s == UINT32_MAX) s->nd[i].s = idx->idx[s->nd[i].name].length - 1;
		else if (s->nd[i].e == UINT32_MAX) s->nd[i].e = idx->idx[s->nd[i].name].length - 1;

		if (s->nd[i].e < s->nd[i].s){
			t = s->nd[i].e;
			s->nd[i].e = s->nd[i].s;
			s->nd[i].s = t;
		}
		len += s->nd[i].e - s->nd[i].s;
	}
	return len;
}

static int cmpfunc(const void * a_, const void * b_){
	ctg_* a = (ctg_*) a_;
	ctg_* b = (ctg_*) b_;
	// return ( ((ctg_*)b)->len - ((ctg_*)a)->len );
	return (a->len < b->len) - (a->len > b->len);
}

void stat_ctg(ctgs *m, const uint64_t gs){
	qsort(m->ctg, m->i, sizeof(ctg_), cmpfunc);
	uint32_t i, k = 0, j = 1;
	uint64_t l = 0;
	plog(2, "Assembly stat:");
	fprintf(stderr, "%-5s %20s %20s\n", "Type", "Length (bp)", "Count (#)");
	for (i = 0; i < m->i; i++){
		if (!m->ctg[i].len){
			k ++;
			continue;
		}
		l += m->ctg[i].len;
		while (l >= gs * 0.1 * j && j < 10) fprintf(stderr, "N%d0 %20"PRIu32" %20"PRIu32"\n", j++, m->ctg[i].len, i + 1 - k);
	}
	fprintf(stderr, "\n%-5s %18"PRIu32" %20s\n", "Min.", m->ctg[i-1-k].len, "-");
	fprintf(stderr, "%-5s %18"PRIu32" %20s\n", "Max.", m->ctg[0].len, "-");
	fprintf(stderr, "%-5s %18"PRIu32" %20s\n", "Ave.", (uint32_t)((float) gs/(i - k) + 0.5), "-");
	fprintf(stderr, "%-5s %18"PRIu64" %20d\n", "Total", gs, i - k);
}

void out_node(graph *s, const UINTL_T j){
	UINT_T i;
	node *n = &s->nd[j];
	UINT_T id = n->id, od = n->od;

	printf("##########node:%"UINTL_FORMAT" reversed node:%"UINTL_FORMAT"###########\n", j, get_reversed_node(s, j));
	for (i = 0; i < n->idm; i++){
		if (s->eg[n->ie[i]].l & MFLAG_FIT) continue;
		id--;
		printf("node:%"UINTL_FORMAT" inedge:%"UINTL_FORMAT" innode:[id:%"UINTL_FORMAT\
			" name:%"UINTL_FORMAT"]->outnode:[id:%"UINTL_FORMAT" name:%"UINTL_FORMAT"]:lable:%"\
			PRIu16":len:%"PRIu32":sco:%"PRIu32" lable:%"UINT_FORMAT" id:%"UINT_FORMAT" od:%"UINT_FORMAT"\n",\
			j, n->ie[i], s->eg[n->ie[i]].in, s->nd[s->eg[n->ie[i]].in].name, s->eg[n->ie[i]].ou,\
			s->nd[s->eg[n->ie[i]].ou].name, s->eg[n->ie[i]].l, s->eg[n->ie[i]].len, s->eg[n->ie[i]].sco, n->l, n->id, n->od);
	}
	printf("\n");
	for (i = 0; i < n->odm; i++){
		if (s->eg[n->oe[i]].l & MFLAG_FIT) continue;
		od--;
		printf("node:%"UINTL_FORMAT" outedge:%"UINTL_FORMAT" innode:[id:%"UINTL_FORMAT\
			" name:%"UINTL_FORMAT"]->outnode:[id:%"UINTL_FORMAT" name:%"UINTL_FORMAT"]:lable:%"\
			PRIu16":len:%"PRIu32":sco:%"PRIu32" lable:%"UINT_FORMAT" id:%"UINT_FORMAT" od:%"UINT_FORMAT"\n",\
			j, n->oe[i], s->eg[n->oe[i]].in, s->nd[s->eg[n->oe[i]].in].name, s->eg[n->oe[i]].ou,\
			s->nd[s->eg[n->oe[i]].ou].name, s->eg[n->oe[i]].l, s->eg[n->oe[i]].len, s->eg[n->oe[i]].sco, n->l, n->id, n->od);
	}
}

void out_edges(graph *g){
	edge *e;
	uint64_t i;
	for (i = 0; i < g->ei; i++){
		e = &g->eg[i];
		if (e->l & MFLAG_FIT) continue;
		if (e->l & MFLAG_IL) printf("%"UINTL_FORMAT".b", g->nd[e->in].name);
		else printf("%"UINTL_FORMAT".e", g->nd[e->in].name);
		if (e->l & MFLAG_OL) printf("\t%"UINTL_FORMAT".b\n", g->nd[e->ou].name);
		else printf("\t%"UINTL_FORMAT".e\n", g->nd[e->ou].name);
	}
}

uint64_t out_ctg_path(ctgs *m, const idxs *idx, const int min_ctg_len){
	uint32_t i, j, l;
	ctg_ *s = NULL;
	uint64_t gs = 0;
	char *type[5] = {"unknown", "linear", "loop", "breakpoint", "junction"};
	for (i = 0; i < m->i; i++){
		s = &m->ctg[i];
		s->len = re_cal_len(s, idx);
		if (s->len < min_ctg_len) {
			s->len = 0;
			continue;
		}
		gs += s->len;
		l = 0;
		#ifdef RESTRICT
			if (gs > RESTRICT_GS) continue;
		#endif
		printf(">ctg%06"PRIu32" type:s:%s length:i:%"PRIu32" node:i:%"PRIu32"\n", i, type[s->type], s->len, s->i);
		for (j = 0; j < s->i; j++){
			if (s->nd[j].l){
				printf("%"UINTL_FORMAT"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", s->nd[j].name,\
				s->nd[j].l, s->nd[j].s, s->nd[j].e-1, l, s->nd[j].ide, s->nd[j].ort, s->nd[j].irt);
			}else{
				printf("%"UINTL_FORMAT"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\t%"PRIu32"\n", s->nd[j].name,\
				s->nd[j].l, s->nd[j].s+1, s->nd[j].e, l, s->nd[j].ide, s->nd[j].ort, s->nd[j].irt);
			}
			l += s->nd[j].e - s->nd[j].s;
		}
	}
	return gs;
}

uint64_t out_ctg_gfa(ctgs *m, graph *g, const idxs *idx, const int min_ctg_len){//TODO need test
	uint32_t i, j;
	uint64_t l, gs = 0;
	ctg_ *s1 = NULL;
	printf("%c\t%s\n",'H', "VN:Z:2.0");
	for (i = 0; i < m->i; i++){
		s1 = &m->ctg[i];
		s1->len = re_cal_len(s1, idx);

		if (s1->len < min_ctg_len) {
			s1->len = 0;
			continue;
		}
		gs += s1->len;
		#ifdef RESTRICT
			if (gs > RESTRICT_GS) continue;
		#endif

		printf("%c\tctg%06"PRIu32"\t%"PRIu32"\t*\n", 'S', i, s1->len);
		l = 0;
		for (j = 0; j < s1->i; j++){
			if (s1->nd[j].l){
				printf("%c\tctg%06"PRIu32"\t%"UINTL_FORMAT"-\t%"PRIu64"\t%"PRIu64"\t%"PRIu32"\t%"\
					PRIu32"\t*\tNI:i:%"UINTL_FORMAT"\n", 'F', i, s1->nd[j].name, l, \
					l + s1->nd[j].e-s1->nd[j].s-1, s1->nd[j].s, s1->nd[j].e-1, s1->nd[j].nidex);
			}else{
				printf("%c\tctg%06"PRIu32"\t%"UINTL_FORMAT"+\t%"PRIu64"\t%"PRIu64"\t%"PRIu32"\t%"\
					PRIu32"\t*\tNI:i:%"UINTL_FORMAT"\n", 'F', i, s1->nd[j].name, l, \
					l + s1->nd[j].e-s1->nd[j].s-1, s1->nd[j].s+1, s1->nd[j].e, s1->nd[j].nidex);
			}
			l += (s1->nd[j].e - s1->nd[j].s);
			s1->nd[j].s = j ? l + 1 - idx->idx[s1->nd[j].name].length : 0;// recal reads approximate start pos in ctg.
		}
	}
	ctg_ *s2 = NULL;
	int64_t s1s, s1e, s2s, s2e, re, t;
	int k = 0;
	char s1l, s2l;
	node *n1s, *n1e;
	for (i = 0; i < m->i - 1; i++){
		s1 = &m->ctg[i];		
		if (!s1->len) continue;
		n1s = &g->nd[s1->nd[0].nidex];
		n1e = s1->i > 1 ? &g->nd[s1->nd[s1->i - 1].nidex] : NULL;
		for (j = i + 1; j < m->i; j++){
			s2 = &m->ctg[j];
			if (!s2->len) continue;
			for (l = 0; l < n1s->idm; l++){
				if (g->eg[n1s->ie[l]].l & MFLAG_FIT) continue;
				if (g->nd[g->eg[n1s->ie[l]].in].name == s2->nd[0].name){
					re = get_reversed_edge(g, n1s->ie[l]);
					s1s = g->eg[n1s->ie[l]].oe;
					s1e = g->eg[re].ie;
					s1l = (g->eg[n1s->ie[l]].l & MFLAG_OL) == s1->nd[0].l * MFLAG_OL ? '+' : '-';
					s2s = g->eg[n1s->ie[l]].ie;
					s2e = g->eg[re].oe;
					s2l = (g->eg[n1s->ie[l]].l & MFLAG_IL) == s2->nd[0].l * MFLAG_IL ? '+' : '-';
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1->nd[0].l){
						t = s1s;
						s1s = s1->nd[0].s + idx->idx[s1->nd[0].name].length - s1e;
						s1e = s1->nd[0].s + idx->idx[s1->nd[0].name].length - t;
					}else{
						s1s += s1->nd[0].s;
						s1e += s1->nd[0].s;
					}

					if (s2->nd[0].l){
						t = s2s;
						s2s = s2->nd[0].s + idx->idx[s2->nd[0].name].length - s2e;
						s2e = s2->nd[0].s + idx->idx[s2->nd[0].name].length - t;
					}else{
						s2s += s2->nd[0].s;
						s2e += s2->nd[0].s;
					}
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1s < 0) s1s = 0;
					if (s2s < 0) s2s = 0;
					printf("%c\teg%06d\tctg%06"PRIu32"%c\tctg%06"PRIu32"%c\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t*\n",\
					'E', k++, i, s1l, j, s2l, s1s, s1e, s2s, s2e);
				}else if(s2->i > 1 && g->nd[g->eg[n1s->ie[l]].in].name == s2->nd[s2->i - 1].name){

					re = get_reversed_edge(g, n1s->ie[l]);
					s1s = g->eg[n1s->ie[l]].oe;
					s1e = g->eg[re].ie;
					s1l = (g->eg[n1s->ie[l]].l & MFLAG_OL) == s1->nd[0].l * MFLAG_OL ? '+' : '-';
					s2s = g->eg[n1s->ie[l]].ie;
					s2e = g->eg[re].oe;
					s2l = (g->eg[n1s->ie[l]].l & MFLAG_IL) == s2->nd[s2->i - 1].l * MFLAG_IL  ? '+' : '-';
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1->nd[0].l){
						t = s1s;
						s1s = s1->nd[0].s + idx->idx[s1->nd[0].name].length - s1e;
						s1e = s1->nd[0].s + idx->idx[s1->nd[0].name].length - t;
					}else{
						s1s += s1->nd[0].s;
						s1e += s1->nd[0].s;
					}

					if (s2->nd[s2->i - 1].l){
						t = s2s;
						s2s = s2->nd[s2->i - 1].s + idx->idx[s2->nd[s2->i - 1].name].length - s2e;
						s2e = s2->nd[s2->i - 1].s + idx->idx[s2->nd[s2->i - 1].name].length - t;
					}else{
						s2s += s2->nd[s2->i - 1].s;
						s2e += s2->nd[s2->i - 1].s;
					}
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1s < 0) s1s = 0;
					if (s2s < 0) s2s = 0;
					printf("%c\teg%06d\tctg%06"PRIu32"%c\tctg%06"PRIu32"%c\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t*\n",\
					'E', k++, i, s1l, j, s2l, s1s, s1e, s2s, s2e);
				}
			}

			for (l = 0; l < n1s->odm; l++){
				if (g->eg[n1s->oe[l]].l & MFLAG_FIT) continue;
				if (g->nd[g->eg[n1s->oe[l]].ou].name == s2->nd[0].name){
					re = get_reversed_edge(g, n1s->oe[l]);
					s1s = g->eg[n1s->oe[l]].ie;
					s1e = g->eg[re].oe;
					s1l = (g->eg[n1s->oe[l]].l & MFLAG_IL) == s1->nd[0].l * MFLAG_IL? '+' : '-';
					s2s = g->eg[n1s->oe[l]].oe;
					s2e = g->eg[re].ie;
					s2l = (g->eg[n1s->oe[l]].l & MFLAG_OL) == s2->nd[0].l * MFLAG_OL ? '+' : '-';
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1->nd[0].l){
						t = s1s;
						s1s = s1->nd[0].s + idx->idx[s1->nd[0].name].length - s1e;
						s1e = s1->nd[0].s + idx->idx[s1->nd[0].name].length - t;
					}else{
						s1s += s1->nd[0].s;
						s1e += s1->nd[0].s;
					}

					if (s2->nd[0].l){
						t = s2s;
						s2s = s2->nd[0].s + idx->idx[s2->nd[0].name].length - s2e;
						s2e = s2->nd[0].s + idx->idx[s2->nd[0].name].length - t;
					}else{
						s2s += s2->nd[0].s;
						s2e += s2->nd[0].s;
					}
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1s < 0) s1s = 0;
					if (s2s < 0) s2s = 0;
					printf("%c\teg%06d\tctg%06"PRIu32"%c\tctg%06"PRIu32"%c\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t*\n",\
					'E', k++, i, s1l, j, s2l, s1s, s1e, s2s, s2e);
				}else if(s2->i > 1 && g->nd[g->eg[n1s->oe[l]].ou].name == s2->nd[s2->i - 1].name){
					re = get_reversed_edge(g, n1s->oe[l]);
					s1s = g->eg[n1s->oe[l]].ie;
					s1e = g->eg[re].oe;
					s1l = (g->eg[n1s->oe[l]].l & MFLAG_IL) == s1->nd[0].l * MFLAG_IL ? '+' : '-';
					s2s = g->eg[n1s->oe[l]].oe;
					s2e = g->eg[re].ie;
					s2l = (g->eg[n1s->oe[l]].l & MFLAG_OL) == s2->nd[s2->i - 1].l * MFLAG_OL ? '+' : '-';
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					
					if (s1->nd[0].l){
						t = s1s;
						s1s = s1->nd[0].s + idx->idx[s1->nd[0].name].length - s1e;
						s1e = s1->nd[0].s + idx->idx[s1->nd[0].name].length - t;
					}else{
						s1s += s1->nd[0].s;
						s1e += s1->nd[0].s;
					}

					if (s2->nd[s2->i - 1].l){
						t = s2s;
						s2s = s2->nd[s2->i - 1].s + idx->idx[s2->nd[s2->i - 1].name].length - s2e;
						s2e = s2->nd[s2->i - 1].s + idx->idx[s2->nd[s2->i - 1].name].length - t;
					}else{
						s2s += s2->nd[s2->i - 1].s;
						s2e += s2->nd[s2->i - 1].s;
					}
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1s < 0) s1s = 0;
					if (s2s < 0) s2s = 0;
					printf("%c\teg%06d\tctg%06"PRIu32"%c\tctg%06"PRIu32"%c\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t*\n",\
					'E', k++, i, s1l, j, s2l, s1s, s1e, s2s, s2e);
				}
			}

			for (l = 0; n1e && l < n1e->idm; l++){
				if (g->eg[n1e->ie[l]].l & MFLAG_FIT) continue;
				if (g->nd[g->eg[n1e->ie[l]].in].name == s2->nd[0].name){
					re = get_reversed_edge(g, n1e->ie[l]);
					s1s = g->eg[n1e->ie[l]].oe;
					s1e = g->eg[re].ie;
					s1l = (g->eg[n1e->ie[l]].l & MFLAG_OL) == s1->nd[s1->i - 1].l * MFLAG_OL ? '+' : '-';
					s2s = g->eg[n1e->ie[l]].ie;
					s2e = g->eg[re].oe;
					s2l = (g->eg[n1e->ie[l]].l & MFLAG_IL) == s2->nd[0].l * MFLAG_IL ? '+' : '-';
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1->nd[s1->i - 1].l){
						t = s1s;
						s1s = s1->nd[s1->i - 1].s + idx->idx[s1->nd[s1->i - 1].name].length - s1e;
						s1e = s1->nd[s1->i - 1].s + idx->idx[s1->nd[s1->i - 1].name].length - t;
					}else{
						s1s += s1->nd[s1->i - 1].s;
						s1e += s1->nd[s1->i - 1].s;
					}

					if (s2->nd[0].l){
						t = s2s;
						s2s = s2->nd[0].s + idx->idx[s2->nd[0].name].length - s2e;
						s2e = s2->nd[0].s + idx->idx[s2->nd[0].name].length - t;
					}else{
						s2s += s2->nd[0].s;
						s2e += s2->nd[0].s;
					}
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1s < 0) s1s = 0;
					if (s2s < 0) s2s = 0;
					printf("%c\teg%06d\tctg%06"PRIu32"%c\tctg%06"PRIu32"%c\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t*\n",\
					'E', k++, i, s1l, j, s2l, s1s, s1e, s2s, s2e);
				}else if(s2->i > 1 && g->nd[g->eg[n1e->ie[l]].in].name == s2->nd[s2->i - 1].name){
					re = get_reversed_edge(g, n1e->ie[l]);
					s1s = g->eg[n1e->ie[l]].oe;
					s1e = g->eg[re].ie;
					s1l = (g->eg[n1e->ie[l]].l & MFLAG_OL) == s1->nd[s1->i - 1].l * MFLAG_OL ? '+' : '-';
					s2s = g->eg[n1e->ie[l]].ie;
					s2e = g->eg[re].oe;
					s2l = (g->eg[n1e->ie[l]].l & MFLAG_IL) == s2->nd[s2->i - 1].l * MFLAG_IL ? '+' : '-';
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1->nd[s1->i - 1].l){
						t = s1s;
						s1s = s1->nd[s1->i - 1].s + idx->idx[s1->nd[s1->i - 1].name].length - s1e;
						s1e = s1->nd[s1->i - 1].s + idx->idx[s1->nd[s1->i - 1].name].length - t;
					}else{
						s1s += s1->nd[s1->i - 1].s;
						s1e += s1->nd[s1->i - 1].s;
					}

					if (s2->nd[s2->i - 1].l){
						t = s2s;
						s2s = s2->nd[s2->i - 1].s + idx->idx[s2->nd[s2->i - 1].name].length - s2e;
						s2e = s2->nd[s2->i - 1].s + idx->idx[s2->nd[s2->i - 1].name].length - t;
					}else{
						s2s += s2->nd[s2->i - 1].s;
						s2e += s2->nd[s2->i - 1].s;
					}
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1s < 0) s1s = 0;
					if (s2s < 0) s2s = 0;
					printf("%c\teg%06d\tctg%06"PRIu32"%c\tctg%06"PRIu32"%c\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t*\n",\
					'E', k++, i, s1l, j, s2l, s1s, s1e, s2s, s2e);
				}
			}
			for (l = 0; n1e && l < n1e->odm; l++){
				if (g->eg[n1e->oe[l]].l & MFLAG_FIT) continue;
				if (g->nd[g->eg[n1e->oe[l]].ou].name == s2->nd[0].name){
					re = get_reversed_edge(g, n1e->oe[l]);
					s1s = g->eg[n1e->oe[l]].ie;
					s1e = g->eg[re].oe;
					s1l = (g->eg[n1e->oe[l]].l & MFLAG_IL) == s1->nd[s1->i - 1].l * MFLAG_IL ? '+' : '-';
					s2s = g->eg[n1e->oe[l]].oe;
					s2e = g->eg[re].ie;
					s2l = (g->eg[n1e->oe[l]].l & MFLAG_OL) == s2->nd[0].l * MFLAG_OL ? '+' : '-';
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1->nd[s1->i - 1].l){
						t = s1s;
						s1s = s1->nd[s1->i - 1].s + idx->idx[s1->nd[s1->i - 1].name].length - s1e;
						s1e = s1->nd[s1->i - 1].s + idx->idx[s1->nd[s1->i - 1].name].length - t;
					}else{
						s1s += s1->nd[s1->i - 1].s;
						s1e += s1->nd[s1->i - 1].s;
					}

					if (s2->nd[0].l){
						t = s2s;
						s2s = s2->nd[0].s + idx->idx[s2->nd[0].name].length - s2e;
						s2e = s2->nd[0].s + idx->idx[s2->nd[0].name].length - t;
					}else{
						s2s += s2->nd[0].s;
						s2e += s2->nd[0].s;
					}
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);					
					if (s1s < 0) s1s = 0;
					if (s2s < 0) s2s = 0;
					printf("%c\teg%06d\tctg%06"PRIu32"%c\tctg%06"PRIu32"%c\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t*\n",\
					'E', k++, i, s1l, j, s2l, s1s, s1e, s2s, s2e);
				}else if(s2->i > 1 && g->nd[g->eg[n1e->oe[l]].ou].name == s2->nd[s2->i - 1].name){
					re = get_reversed_edge(g, n1e->oe[l]);
					s1s = g->eg[n1e->oe[l]].ie;
					s1e = g->eg[re].oe;
					s1l = (g->eg[n1e->oe[l]].l & MFLAG_IL) == s1->nd[s1->i - 1].l * MFLAG_IL ? '+' : '-';
					s2s = g->eg[n1e->oe[l]].oe;
					s2e = g->eg[re].ie;
					s2l = (g->eg[n1e->oe[l]].l & MFLAG_OL) == s2->nd[s2->i - 1].l * MFLAG_OL ? '+' : '-';
					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);
					if (s1->nd[s1->i - 1].l){
						t = s1s;
						s1s = s1->nd[s1->i - 1].s + idx->idx[s1->nd[s1->i - 1].name].length - s1e;
						s1e = s1->nd[s1->i - 1].s + idx->idx[s1->nd[s1->i - 1].name].length - t;
					}else{
						s1s += s1->nd[s1->i - 1].s;
						s1e += s1->nd[s1->i - 1].s;
					}

					if (s2->nd[s2->i - 1].l){
						t = s2s;
						s2s = s2->nd[s2->i - 1].s + idx->idx[s2->nd[s2->i - 1].name].length - s2e;
						s2e = s2->nd[s2->i - 1].s + idx->idx[s2->nd[s2->i - 1].name].length - t;
					}else{
						s2s += s2->nd[s2->i - 1].s;
						s2e += s2->nd[s2->i - 1].s;
					}

					if (s1s > s1e) SWAP(s1s, s1e, int64_t);
					if (s2s > s2e) SWAP(s2s, s2e, int64_t);					
					if (s1s < 0) s1s = 0;
					if (s2s < 0) s2s = 0;
					printf("%c\teg%06d\tctg%06"PRIu32"%c\tctg%06"PRIu32"%c\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t*\n",\
					'E', k++, i, s1l, j, s2l, s1s, s1e, s2s, s2e);
				}
			}
		}
	}
	return gs;
}

uint64_t out_ctg_fasta(ctgs *m, const idxs *idx, const int min_ctg_len, const int out_seq){
	uint32_t i, j, k, l;
	uint64_t qv, gs = 0;
	ctg_ *s = NULL;
	char *type[5] = {"unknown", "linear", "loop", "breakpoint", "junction"};
	sbbuf *buf = init_sbbuf(MAX_INIT_N);
	for (i = 0; i < m->i; i++){
		s = &m->ctg[i];
		s->len = re_cal_len(s, idx);
		if (s->len < min_ctg_len) {
			s->len = 0;
			continue;
		}
		gs += s->len;
		#ifdef RESTRICT
			if (gs > RESTRICT_GS) continue;
		#endif
		if (out_seq){
			printf(">ctg%06"PRIu32" type:s:%s length:i:%"PRIu32" node:i:%"PRIu32" qv:i", i, type[s->type], s->len, s->i);
			
			for (l = j = 0; j < s->i; j++){
				qv = 0;
				qv |= l;
				qv <<= 12;
				qv |= s->nd[j].ide;
				qv <<= 10;
				qv |= s->nd[j].ort;
				qv <<= 10;
				qv |= s->nd[j].irt;
				printf(":%"PRIx64, qv);//PRIu64
				l += s->nd[j].e - s->nd[j].s;
			}
			printf("\n");

			for (j = 0; j < s->i; j++){
				if (s->nd[j].l){
					subfa(buf, idx->idx[s->nd[j].name].f2bit, s->nd[j].l, idx->idx[s->nd[j].name].offset, s->nd[j].s, s->nd[j].e-1);
				}else{
					subfa(buf, idx->idx[s->nd[j].name].f2bit, s->nd[j].l, idx->idx[s->nd[j].name].offset, s->nd[j].s+1, s->nd[j].e);
				}
				if (s->nd[j].lq){
					for (k = 0; buf->seq[k] != '\0'; k++) putchar("acgt"[nt_table[(uint8_t)buf->seq[k]]]);
				}else printf("%s", buf->seq);
				// printf("%s", buf->seq);
			}
			printf("\n");
		}
	}
	destroy_sbbbuf(buf);
	return gs;
}

void out_graph_raw(graph *s){
	UINTL_T i;
	edge *e = NULL;
	printf("##########graph###########\n");
	for (i = 0; i < s->ei; i++){
		e = &s->eg[i];
		if (e->l & MFLAG_FIT) continue;
		printf("edge:%"UINTL_FORMAT" innode:[id:%"UINTL_FORMAT" name:%"UINTL_FORMAT"]->outnode:[id:%"\
			UINTL_FORMAT" name:%"UINTL_FORMAT"] lable:%"PRIu16" ide:%"PRIu16" ie:%"PRIu32" oe:%"PRIu32\
			" len:%"PRIu32" sco:%"PRIu32"\n", i,e->in,s->nd[e->in].name, e->ou, s->nd[e->ou].name, e->l,\
			e->ide, e->ie, e->oe, e->len, e->sco);
	}
}

void out_graph_graphml(graph *s){
	UINTL_T i;
	edge *e = NULL;
	printf("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	printf("<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">\n");
	printf("<key attr.name=\"l\" attr.type=\"int\" for=\"edge\" id=\"e0\"/>\n");
	printf("<key attr.name=\"ide\" attr.type=\"int\" for=\"edge\" id=\"e1\"/>\n");
	#ifndef RESTRICT
	printf("<key attr.name=\"ie\" attr.type=\"int\" for=\"edge\" id=\"e2\"/>\n");
	printf("<key attr.name=\"oe\" attr.type=\"int\" for=\"edge\" id=\"e3\"/>\n");
	#endif
	printf("<key attr.name=\"len\" attr.type=\"int\" for=\"edge\" id=\"e4\"/>\n");
	printf("<key attr.name=\"sco\" attr.type=\"int\" for=\"edge\" id=\"e5\"/>\n");

	printf("<key attr.name=\"nodeid\" attr.type=\"int\" for=\"node\" id=\"n0\"/>\n");
	printf("<key attr.name=\"revnodeid\" attr.type=\"int\" for=\"node\" id=\"n1\"/>\n");

	printf("<graph edgedefault=\"directed\">\n");
	for (i = 0; i < s->ni; i++){
		if (s->nd[i].id + s->nd[i].od == 0) continue;
		if (i > get_reversed_node(s, i)) continue;
		printf("<node id=\"%"UINTL_FORMAT".b\">\n", s->nd[i].name);
		printf("<data key=\"n0\">%"UINTL_FORMAT"</data>\n", i);
		printf("<data key=\"n1\">%"UINTL_FORMAT"</data>\n", get_reversed_node(s, i));
		printf("</node>\n");

		printf("<node id=\"%"UINTL_FORMAT".e\">\n", s->nd[i].name);
		printf("<data key=\"n0\">%"UINTL_FORMAT"</data>\n", i);
		printf("<data key=\"n1\">%"UINTL_FORMAT"</data>\n", get_reversed_node(s, i));
		printf("</node>\n");
	}

	for (i = 0; i < s->ei; i++){
		e = &s->eg[i];
		if (e->l & MFLAG_FIT) continue;
		if (e->l & MFLAG_IL) printf("<edge source=\"%"UINTL_FORMAT".b\"", s->nd[e->in].name);
		else printf("<edge source=\"%"UINTL_FORMAT".e\"", s->nd[e->in].name);

		if (e->l & MFLAG_OL) printf(" target=\"%"UINTL_FORMAT".b\">\n", s->nd[e->ou].name);
		else printf(" target=\"%"UINTL_FORMAT".e\">\n", s->nd[e->ou].name);
		printf("<data key=\"e0\">%"PRIu16"</data>\n", e->l);
		printf("<data key=\"e1\">%"PRIu16"</data>\n", e->ide);
		#ifndef RESTRICT
		printf("<data key=\"e2\">%"PRIu32"</data>\n", e->ie);
		printf("<data key=\"e3\">%"PRIu32"</data>\n", e->oe);
		#endif
		printf("<data key=\"e4\">%"PRIu32"</data>\n", e->len);
		printf("<data key=\"e5\">%"PRIu32"</data>\n", e->sco);
		printf("</edge>\n");

	}
	printf("</graph></graphml>\n");
}
