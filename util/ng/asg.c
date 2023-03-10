#include <assert.h>
#include "asg.h"
#include "opt.h"

#define BFLAG_N   0x1
#define BFLAG_P1  0x2
#define BFLAG_P2  0x4
#define BFLAG_LP  0x8
#define BFLAG_U   0x10

typedef struct {
	UINTL_T *n;
	int64_t sco;
	// double sco;
	int i;
	int im;
} ph_;

typedef struct {
	ph_ *ph;
	// int i;
	// int im;
	int64_t i;
	int64_t im;
} phs;

typedef struct {
	UINTL_T *n;
	int i;
	int im;
	uint16_t ide;
	uint16_t perc;// should x 10000
	uint32_t sco;
} pl_;

typedef struct {
	pl_ *pl;
	int i;
	int im;
} pls;

typedef struct {
    UINT_T id;
    UINT_T l;
    UINTL_T len; 
    UINTL_T pnode;
} cpath_info;

typedef struct bfs_node_info{ 
	UINTL_T start;
	UINTL_T predecessor;
	uint32_t edge_num;
	int64_t prop_diff_sum;
} bfs_node_info;

typedef struct bfs{
	uint16_t l;
	UINT_T unvisited;
	int bfs_depth;
	int info_len;
	union {
		bfs_node_info *info;
		uint64_t pnode;
	};
} bfs;

#if GENOME_SIZE == 2
	KHASH_MAP_INIT_INT64(bfs, bfs);
#else
	KHASH_MAP_INIT_INT(bfs, bfs);
#endif

typedef struct bfs_r{
	void *data1;
	void *data2;
	khash_t(bfs) *bfs_info;
} bfs_r;

void init_graph(graph *s, const uint32_t len){
	s->ni = 1;
	s->ei = 0;
	s->nm = len;
	s->em = 2 * len;
	s->nd = calloc (s->nm, sizeof(node));
	s->eg = calloc (s->em, sizeof(edge));
}

void destroy_graph(graph *s){
	UINTL_T i;
	for (i = 0; i < s->ni; i++){
		if (s->nd[i].ie) free (s->nd[i].ie);
		if (s->nd[i].oe) free (s->nd[i].oe);
	}
	free (s->nd);
	free (s->eg);
}

static void init_ph(ph_ *s, const uint32_t len){
	s->im = len;
	s->sco = s->i = 0;
	s->n = malloc(s->im * sizeof(UINTL_T));
}

static void reallocate_ph(ph_ *s, const uint32_t len){
	s->im += len;
	s->n = realloc(s->n, s->im * sizeof(UINTL_T));
}

static void init_phs(phs *s, const uint32_t len, const uint32_t ph_len){
	int i;
	s->im = len;
	s->i = 0;
	s->ph = malloc(s->im * sizeof(ph_));
	for (i = 0; i < s->im; i++) init_ph(&s->ph[i], ph_len);
}

static void reallocate_phs(phs *s, const uint32_t len, const uint32_t ph_len){
	int i = s->im;
	s->im += len;
	s->ph = realloc(s->ph, sizeof(ph_) * s->im);
	for (;i < s->im; i++) init_ph(&s->ph[i], ph_len);
}

static void destroy_phs(phs *s){
	int i;
	for (i = 0; i < s->im; i++) free(s->ph[i].n);
	//{printf("%p %d %d %p\n",s, i,s->im,s->ph[i].n ); 
	free (s->ph);
}

static void clear_phs(phs *s){
	int i;
	for (s->i = i = 0; i < s->im; i++) s->ph[i].sco = s->ph[i].i = 0;
}

static void sort_phs_g(phs *s, graph *g){
	int i, j;
	ph_ p;
	for (i = 1; i < s->i; i ++){// insert sort
		p = s->ph[i];
		for (j = i; j > 0 && (s->ph[j - 1].sco > p.sco || (s->ph[j - 1].sco == p.sco 
			&& g->nd[g->eg[s->ph[j - 1].n[0]].in].name > g->nd[g->eg[p.n[0]].in].name )); j--) s->ph[j] = s->ph[j - 1];
		s->ph[j] = p;
	}
}

static int cmp_sco(const void * a_, const void * b_){
	ph_ * a = (ph_ *)a_;
	ph_ * b = (ph_ *)b_;
	return (a->sco > b->sco) - (a->sco < b->sco);
	// return ((ph_ *)a)->sco - ((ph_ *)b)->sco;
}

static void qsort_phs(phs *s){// for bubble z-path
	qsort(s->ph, s->i, sizeof(ph_), cmp_sco);
}

static void msort_phs(phs *s){// for bubble z-path
	int i, j;
	ph_ p;
	for (i = 1; i < s->i; i ++){// insert sort
		p = s->ph[i];
		for (j = i; j > 0 && s->ph[j - 1].sco > p.sco; j--) s->ph[j] = s->ph[j - 1];
		s->ph[j] = p;
	}
}

static void init_pl(pl_ *s, const uint32_t len){
	s->im = len;
	s->i = 0;
	s->n = malloc(s->im * sizeof(UINTL_T));
}

static void reallocate_pl(pl_ *s, const uint32_t len){
	s->im += len;
	s->n = realloc(s->n, s->im * sizeof(UINTL_T));
}

static void init_pls(pls *s, const uint32_t len, const uint32_t pl_len){
	int i;	
	s->im = len;
	s->i = 0;
	s->pl = malloc(s->im * sizeof(pl_));
	for (i = 0; i < s->im; i++) init_pl(&s->pl[i], pl_len);
}

static void reallocate_pls(pls *s, const uint32_t len, const uint32_t pl_len){
	int i = s->im;
	s->im += len;
	s->pl = realloc(s->pl, sizeof(pl_) * s->im);
	for (;i < s->im; i++) init_pl(&s->pl[i], pl_len);
}

static void destroy_pls(pls *s){
	int i;
	for (i = 0; i < s->im; i++) free(s->pl[i].n);
	free (s->pl);
}

static int sort_pls(const void *p1, const void *p2){// for bubble z-path
	pl_ *a = (pl_ *) p1;
	pl_ *b = (pl_ *) p2;

	if (a->perc != b->perc){
		return b->perc - a->perc;
	}else if (a->ide > b->ide * 5 / 4 || a->ide < b->ide * 4 / 5){
		return a->ide - b->ide;
	}else{
		// return a->sco - b->sco;
		return (a->sco > b->sco) - (a->sco < b->sco);
	}
}

static inline UINTL_T add_edge(graph *g, const UINTL_T in, const UINTL_T ou, const uint32_t ie, \
	const uint32_t oe, const uint32_t sco, const uint16_t ide, const uint32_t len, const uint32_t l) {
	if (g->ei >= g->em) reallocate_edge(g, MAX_INIT_N);
	UINTL_T e;
	if (l){
		e = rp_exited_edge(g, in, ou);
		g->eg[e].l = 0;
	}else{
		e = g->ei++;
	}
	g->eg[e].ide = ide;
	g->eg[e].in = in;
	g->eg[e].ou = ou;
	g->eg[e].ie = ie;
	g->eg[e].oe = oe;
	g->eg[e].len = len;
	g->eg[e].sco = sco;
	add_oe(g, in, e);
	add_ie(g, ou, e);
	return e;
}

static inline void rm_edge(graph *g, const UINTL_T n){ 
	if (!(g->eg[n].l & MFLAG_FIT)){
		g->eg[n].l |= MFLAG_FIT;
		g->nd[g->eg[n].in].od --; 
		g->nd[g->eg[n].ou].id --; 
	}   
}

void rm_node_con(graph *g, khash_t(ovlinfo_) *os, opt *p, list *li){
	UINTL_T i;
	for (i = 1; i < g->ni; i++){//should re-calculate odm and idm firstly
		if (p->out_alt_ctg){
			if (g->nd[i].od >= g->nd[i].odm){// used to save an extra attribute
				g->nd[i].oe = realloc(g->nd[i].oe, (g->nd[i].od + 1) * sizeof(UINTL_T));
			}
			if (g->nd[i].id >= g->nd[i].idm){// used to save an extra attribute
				g->nd[i].ie = realloc(g->nd[i].ie, (g->nd[i].id + 1) * sizeof(UINTL_T));
			}
		}
		g->nd[i].odm = g->nd[i].od;
		g->nd[i].idm = g->nd[i].id;
	}
	for (i = 1; i < g->ni; i++){
		if (kh_val(os, kh_get(ovlinfo_, os, li->data[g->nd[i].name])).con >= p->min_con_count) {
			rm_node(g, i);
			rm_node(g, get_reversed_node(g, i));
		}
	}
}

void sort_stat_oe(graph *g, det *d, opt *p){
	UINTL_T i, j, k, e;
	node *n = NULL;
	int depth_i = 0, *depth = malloc(RANDOM_COUNT * sizeof(int));
	for (i = 1; i < g->ni; i++){
		n = &g->nd[i];

		if (n->od < p->min_aln_depth) d->l[1] ++;
		else if (n->od < p->max_aln_depth) d->n[1] ++;
		else d->h[1] ++;

		if (depth_i < RANDOM_COUNT){
			depth[depth_i++] = n->od;
		}

		if (p->sort){// insert sort
			for(j = 1; j < n->odm; j++){
				e = n->oe[j];
				for(k = j; k > 0 && (g->eg[n->oe[k-1]].len > g->eg[e].len || (g->eg[n->oe[k-1]].len == \
					g->eg[e].len && g->eg[n->oe[k-1]].sco < g->eg[e].sco)); k--) n->oe[k] = n->oe[k-1];
				n->oe[k] = e;
			}
		}

		// printf("\n%u ", i);
		// for (j = 0; j < n->odm; j ++){
		// 	printf("(%d %u %u) ",j, g->eg[n->oe[j]].len, g->eg[n->oe[j]].sco);
		// }
	}
	p->median_outdegree = quick_select(depth, 0, depth_i - 1, depth_i/2);
	free (depth);
}

void rm_edge_lq(graph *g, khash_t(ovlinfo_) *os, opt *p, list *li){
	UINTL_T i;
	ovlinfo *loli = NULL, *roli = NULL;
	for (i = 0; i < g->ei; i++){
		if (g->eg[i].l & MFLAG_FIT) continue;
		loli = &kh_val(os, kh_get(ovlinfo_, os, li->data[g->nd[g->eg[i].in].name]));
		roli = &kh_val(os, kh_get(ovlinfo_, os, li->data[g->nd[g->eg[i].ou].name]));
		// printf("\n%u %u %u",i,g->nd[g->eg[i].in].name, g->nd[g->eg[i].ou].name);
		if (check_valid_edge(&g->eg[i], p, loli, roli) < p->min_node_count){
			rm_edge(g, i);
			rm_edge(g, get_reversed_edge(g, i));
		}
	}
}


static void sort_sco(UINTL_T *st, const graph *g, const UINT_T odm){//from large to small
	UINT_T i, j;
	UINTL_T e;
	for (i = 1; i < odm; i ++){
		e = st[i];
		for (j = i; j > 0 && (g->eg[st[j-1]].sco < g->eg[e].sco || (g->eg[st[j-1]].sco == g->eg[e].sco && \
			g->eg[st[j-1]].ide < g->eg[e].ide)); j--) st[j] = st[j-1];
		st[j] = e;
	}
	// printf("\n%u ", i);
	// for (j = 0; j < odm; j ++){
	// 	printf("(%d %u %d) ",j, g->eg[st[j]].sco, g->eg[st[j]].ide);
	// }
}

static void sort_ide(UINTL_T *st, const graph *g, const UINT_T odm){//from large to small
	UINT_T i, j;
	UINTL_T e;
	for (i = 1; i < odm; i ++){
		e = st[i];
		for (j = i; j > 0 && (g->eg[st[j-1]].ide < g->eg[e].ide || (g->eg[st[j-1]].ide == g->eg[e].ide && \
			g->eg[st[j-1]].sco < g->eg[e].sco)); j--) st[j] = st[j-1];
		st[j] = e;
	}
}

void mark_edge_rep(graph *g, khash_t(ovlinfo_) *os, opt *p, const det *d, list *li){
	float d0 = (float) p->median_aln_depth;
	float d1 = (float) p->median_outdegree;
	// float maxd11 = d0 * p->min_depth_multi * (1 + (float)d->h[0]/d->n[0]);
	// float maxd21 = d1 * p->min_depth_multi * (1 + (float)d->h[1]/d->n[1]);
	float maxd11 = d0 * p->min_depth_multi;
	float maxd21 = d1 * p->min_depth_multi;

	float maxd12 = d0 * p->max_depth_multi;
	float maxd22 = d1 * p->max_depth_multi; //TODO: should node? not edge
	int maxd23 = (int) d1 * 10;

	double total1 = (double) d->l[0] + d->n[0] + d->h[0];
	double total2 = (double) d->l[1] + d->n[1] + d->h[1];
	plog(2, "Depth stat, Mid: %.3f Max: %.3f Repeat: %.3f L:N:H: %.3f:%.3f:%.3f", d0, maxd12, maxd11, 
		d->l[0]/total1, d->n[0]/total1, d->h[0]/total1);
	plog(2, "Outdegree stat, Mid: %.3f Max: %.3f Repeat: %.3f L:N:H: %.3f:%.3f:%.3f", d1, maxd22, maxd21, 
		d->l[1]/total2, d->n[1]/total2, d->h[1]/total2);
	UINTL_T i, j;
	node *n = NULL;
	UINT_T dp;

	UINTL_T st_t, st_i = 0;
	UINTL_T st_m = d1 * 1000;
	UINTL_T *st = malloc(sizeof(UINTL_T) * st_m);

	for (i = 1; i < g->ni; i++){
		n = &g->nd[i];
		j = 0;
		if (n->od >= maxd22) {
			rm_node(g, i);
			rm_node(g, get_reversed_node(g, i));
			continue;
		}else if (n->od >= maxd21) {
			// printf("1 %u %u %u %f\n",i,n->name, n->od, maxd21);
			j = 1;
		}else if (n->od){
			dp = g->eg[n->oe[0]].l & MFLAG_IL ? kh_val(os, kh_get(ovlinfo_, os, li->data[n->name])).lc : \
				kh_val(os, kh_get(ovlinfo_, os, li->data[n->name])).rc;
			if (dp >= maxd12){
				rm_node(g, i);
				rm_node(g, get_reversed_node(g, i));
				continue;
			}else if (dp >= maxd11) {
				j = 1;
				// printf("2 %u %u %u %f\n",i,n->name, dp, maxd11);
			}
		}

		if (n->od >= maxd23){
			st_i = 0;
			if (n->odm >= st_m) {
				st_m = n->odm + 1;
				st = realloc(st, sizeof(UINTL_T) * st_m);
			}
			memcpy(st, n->oe, sizeof(UINTL_T) * n->odm);
			sort_sco(st, g, n->odm);
			for (st_t = st_i = 0; st_i < n->odm; st_i++){
				if (g->eg[st[st_i]].l & MFLAG_FIT) continue;
				if (st_t ++ >= maxd23){
					rm_edge(g, st[st_i]);
					rm_edge(g, get_reversed_edge(g, st[st_i]));
				}
			}
		}

		if (j){
			for (j = 0; j < n->odm; j++) g->eg[n->oe[j]].l |= MFLAG_REP1;
			j = get_reversed_node(g, i);
			n = &g->nd[j];
			for (j =0; j < n->idm; j++) g->eg[n->ie[j]].l |= MFLAG_REP1;
		}
	}
	free (st);
}

void mark_edge_tr(graph *g, const int fuzz){
	int l;
	uint32_t mlen = 0;
	UINTL_T i, j, k;
	node *v = NULL, *w = NULL;	
	for (i = 1; i < g->ni; i++){
		v = &g->nd[i];
		if (!v->od) continue;
		mlen = 0;
		for (j = 0; j < v->odm; j++){
			if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
			g->nd[g->eg[v->oe[j]].ou].l = 1;
			if (g->eg[v->oe[j]].len > mlen) {
				mlen = g->eg[v->oe[j]].len;
			}
		}
		mlen += fuzz;

		for (j = 0; j < v->odm; j++){
			if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
			w = &g->nd[g->eg[v->oe[j]].ou];
			if (w->l){
				for (k = 0; k < w->odm; k++){
					if (g->eg[w->oe[k]].l & MFLAG_FIT) continue;
					if (g->eg[v->oe[j]].len + g->eg[w->oe[k]].len <= mlen) \
						// && g->nd[g->eg[w->oe[k]].ou].l) 
						g->nd[g->eg[w->oe[k]].ou].l = 0;
				}
			}
		}

		for (j = 0; j < v->odm; j++){
			if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
			w = &g->nd[g->eg[v->oe[j]].ou];
			for (l = 1, k = 0; k < w->odm; k++){
				if (g->eg[w->oe[k]].l & MFLAG_FIT) continue;
				if (l || g->eg[w->oe[k]].len < fuzz) g->nd[g->eg[w->oe[k]].ou].l = 0;
				if (l) l = 0;
			}
		}
		
		for (j = 0; j < v->odm; j++){
			if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
			if (g->nd[g->eg[v->oe[j]].ou].l == 0){
				rm_edge(g, v->oe[j]);
				rm_edge(g, get_reversed_edge(g, v->oe[j]));
			}else g->nd[g->eg[v->oe[j]].ou].l = 0;
		}
	}
}

void rm_edge_spur(graph *g){
	UINTL_T i;
	UINT_T j;
	node *v = NULL;
	for (i = 1; i < g->ni; i++){
		v = &g->nd[i];
		if (v->od > 1){
			for (j = 0; j < v->odm; j++){
				if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
				// printf("%u %u %d %u %d len:%u\n",i, v->name, v->od, g->nd[g->eg[v->oe[j]].ou].name, g->nd[g->eg[v->oe[j]].ou].od, g->eg[v->oe[j]].len);
				if (v->od > 1 && !g->nd[g->eg[v->oe[j]].ou].od){
					rm_edge(g, v->oe[j]);
					rm_edge(g, get_reversed_edge(g, v->oe[j]));
				}
			}
		}
	}
}

void mark_edge_hli(graph *g, const float max_ide_ratio){
	UINTL_T i;
	UINT_T j;
	double ms = 0;
	node *v = NULL;
	for (i = 1; i < g->ni; i++){
		v = &g->nd[i];
		if (v->od){
			ms = 0;
			for (j = 0; j < v->odm; j++){
				if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->oe[j]].ide > ms) ms = g->eg[v->oe[j]].ide;
			}
			ms *= max_ide_ratio;
			for (j = 0; j < v->odm; j++){
				if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->oe[j]].ide >= ms){
					g->eg[v->oe[j]].l |= MFLAG_HI;
				}else{
					g->eg[v->oe[j]].l |= MFLAG_LI;
				}
			}
		}

		if (v->id){
			ms = 0;
			for (j = 0; j < v->idm; j++){
				if (g->eg[v->ie[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->ie[j]].ide > ms) ms = g->eg[v->ie[j]].ide;
			}
			ms *= max_ide_ratio;
			for (j = 0; j < v->idm; j++){
				if (g->eg[v->ie[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->ie[j]].ide >= ms){
					g->eg[v->ie[j]].l |= MFLAG_HI;
				}else{
					g->eg[v->ie[j]].l |= MFLAG_LI;
				}
			}
		}
	}

	for (i = 0; i < g->ei; i++){
		// if (g->eg[i].l & MFLAG_FIT) continue;
		if (!(g->eg[i].l & MFLAG_HI)){
			rm_edge(g, i);
			rm_edge(g, get_reversed_edge(g, i));
		}
	}
	// out_graph(g);
}

void rm_edge_li(graph *g){
	UINTL_T i;
	node *v = NULL;

	UINTL_T st_i = 0, st_m = 5000;
	UINTL_T *st = malloc(sizeof(UINTL_T) * st_m);

	for (i = 1; i < g->ni; i++){
		v = &g->nd[i];
		st_i = 0;
		if (v->od <= 1) continue;
		if (v->odm >= st_m) {
			st_m = v->odm + 1;
			st = realloc(st, sizeof(UINTL_T) * st_m);
		}
		memcpy(st, v->oe, sizeof(UINTL_T) * v->odm);
		sort_ide(st, g, v->odm);
		for (st_i = v->odm; st_i && v->od > 1; st_i--){
			if (g->eg[st[st_i-1]].l & MFLAG_FIT) continue;
			if(g->eg[st[st_i-1]].l & MFLAG_LI && g->nd[g->eg[st[st_i-1]].ou].id > 1){
				rm_edge(g, st[st_i-1]);
				rm_edge(g, get_reversed_edge(g, st[st_i-1]));
			}
		}
	}
	free (st);
}

void mark_edge_hls(graph *g, const float max_sco_ratio){
	UINTL_T i;
	UINT_T j;
	double ms = 0;
	node *v = NULL;
	for (i = 1; i < g->ni; i++){
		v = &g->nd[i];
		if (v->od){
			ms = 0;
			for (j = 0; j < v->odm; j++){
				if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->oe[j]].sco > ms) ms = g->eg[v->oe[j]].sco;
			}
			ms *= max_sco_ratio;
			for (j = 0; j < v->odm; j++){
				if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->oe[j]].sco >= ms){
					g->eg[v->oe[j]].l |= MFLAG_HS;
				}else{
					g->eg[v->oe[j]].l |= MFLAG_LS;
				}
			}
		}

		if (v->id){
			ms = 0;
			for (j = 0; j < v->idm; j++){
				if (g->eg[v->ie[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->ie[j]].sco > ms) ms = g->eg[v->ie[j]].sco;
			}
			ms *= max_sco_ratio;
			for (j = 0; j < v->idm; j++){
				if (g->eg[v->ie[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->ie[j]].sco >= ms){
					g->eg[v->ie[j]].l |= MFLAG_HS;
				}else{
					g->eg[v->ie[j]].l |= MFLAG_LS;
				}
			}
		}
	}

	for (i = 0; i < g->ei; i++){
		// if (g->eg[i].l & MFLAG_FIT) continue;
		if (!(g->eg[i].l & MFLAG_HS)){
			rm_edge(g, i);
			rm_edge(g, get_reversed_edge(g, i));
		}
	}
	// out_graph(g);
}

void rm_edge_ls(graph *g){
	UINTL_T i;
	node *v = NULL;

	UINTL_T st_i = 0, st_m = 5000;
	UINTL_T *st = malloc(sizeof(UINTL_T) * st_m);

	for (i = 1; i < g->ni; i++){
		v = &g->nd[i];
		st_i = 0;
		if (v->od <= 1) continue;
		if (v->odm >= st_m) {
			st_m = v->odm + 1;
			st = realloc(st, sizeof(UINTL_T) * st_m);
		}
		memcpy(st, v->oe, sizeof(UINTL_T) * v->odm);
		sort_sco(st, g, v->odm);
		for (st_i = v->odm; st_i && v->od > 1; st_i--){
			if (g->eg[st[st_i-1]].l & MFLAG_FIT) continue;
			// printf("%u %u %u\n",i, st_i, g->eg[st[st_i-1]].sco);
			if(g->eg[st[st_i-1]].l & MFLAG_LS && g->nd[g->eg[st[st_i-1]].ou].id > 1){
				rm_edge(g, st[st_i-1]);
				rm_edge(g, get_reversed_edge(g, st[st_i-1]));
			}
		}
		// for (st_i = 0; st_i < v->odm && v->od > 1; st_i++){
		// 	if (g->eg[v->oe[st_i]].l & MFLAG_FIT) continue;
		// 	if (g->eg[v->oe[st_i]].l & MFLAG_LS && g->nd[g->eg[v->oe[st_i]].ou].id > 1){
		// 		rm_edge(g, v->oe[st_i]);
		// 		rm_edge(g, get_reversed_edge(g, v->oe[st_i]));
		// 	}
		// }
	}
	free (st);
}

void mark_edge_bs(graph *g){
	UINT_T j;
	UINTL_T i, k = 0;
	uint32_t bs, bide;
	node *v = NULL;
	for (i = 1; i < g->ni; i++){
		v = &g->nd[i];
		if (v->od){
			for (bs = bide = j = 0; j < v->odm; j++){
				if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->oe[j]].l & MFLAG_REP1) g->eg[v->oe[j]].l |= MFLAG_BS;
				if (g->eg[v->oe[j]].sco > bs || (g->eg[v->oe[j]].sco == bs && g->eg[v->oe[j]].ide > bide)) {
					k = v->oe[j];
					bs = g->eg[k].sco;
					bide = g->eg[k].ide;
				}
			}
			if (bs) g->eg[k].l |= MFLAG_BS;
		}
		if (v->id){
			for (bs = bide = j = 0; j < v->idm; j++){
				if (g->eg[v->ie[j]].l & MFLAG_FIT) continue;
				if (g->eg[v->ie[j]].l & MFLAG_REP1) g->eg[v->ie[j]].l |= MFLAG_BS;
				if (g->eg[v->ie[j]].sco > bs || (g->eg[v->ie[j]].sco == bs && g->eg[v->ie[j]].ide > bide)) {
					k = v->ie[j];
					bs = g->eg[k].sco;
					bide = g->eg[k].ide;
				}
			}
			if (bs) g->eg[k].l |= MFLAG_BS;
		}
	}
	
	for (i = 0; i < g->ei; i++){
		if (g->eg[i].l & MFLAG_FIT) continue;
		if (!(g->eg[i].l & MFLAG_BS)){
			rm_edge(g, i);
			rm_edge(g, get_reversed_edge(g, i));
		}
	}
}

void rm_sht_brh(graph *g, const int s, const int l){
	UINTL_T i, n;
	UINT_T t;
	ph p;
	init_sn(&p, s + 1);
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].od == 1 && !g->nd[i].id){
			p.i = 0;
			t = get_validly_oe(g, i, 0);
			p.n[p.i++] = g->nd[i].oe[t];
			n = g->eg[g->nd[i].oe[t]].ou;
			while (g->nd[n].id == 1 && g->nd[n].od == 1 && p.i < s){
				t = get_validly_oe(g, n, 0);
				p.n[p.i++] = g->nd[n].oe[t];
				n = g->eg[g->nd[n].oe[t]].ou;
			}

			if (g->nd[n].id > 1 || !g->nd[n].od){
				if (l && p.i >= l){
					if (g->nd[n].od){
						rm_edge(g, p.n[p.i-1]);
						rm_edge(g, get_reversed_edge(g, p.n[p.i-1]));
					}
				}else{
					while (--p.i >= 0){
						rm_edge(g, p.n[p.i]);
						rm_edge(g, get_reversed_edge(g, p.n[p.i]));
					}
				}
			}
		}
	}
	destroy_sn(&p);
}

// void rm_z_clip_lable(graph *g, const int s, const int m, const uint16_t l){
// 	UINTL_T i, n;
// 	UINT_T j, t;
// 	int rm;
// 	ph p;
// 	init_sn(&p, s + 1);
// 	for (i = 1; i < g->ni; i++){
// 		for (j = 0; g->nd[i].od > 1 && j < g->nd[i].odm; j++){
// 			if (g->eg[g->nd[i].oe[j]].l & MFLAG_FIT) continue;
// 			rm = p.i = 0;		
// 			p.n[p.i++] = g->nd[i].oe[j];
// 			if (g->eg[g->nd[i].oe[j]].l & l) rm = 1;
// 			n = g->eg[g->nd[i].oe[j]].ou;
// 			while (g->nd[n].id == 1 && g->nd[n].od == 1 && p.i <= s){
// 				t = get_validly_oe(g, n, 0);
// 				p.n[p.i++] = g->nd[n].oe[t];
// 				if (g->eg[g->nd[n].oe[t]].l & l) rm = 1;
// 				n = g->eg[g->nd[n].oe[t]].ou;
// 			}

// 			if (g->nd[n].id > 1){
// 				if (l && !rm) continue;
// 				if (m && p.i >= m){
// 					rm_edge(g, p.n[p.i-1]);
// 					rm_edge(g, get_reversed_edge(g, p.n[p.i-1]));
// 					rm_edge(g, p.n[0]);
// 					rm_edge(g, get_reversed_edge(g, p.n[0]));
// 				}else{
// 					while (--p.i >= 0){
// 						rm_edge(g, p.n[p.i]);
// 						rm_edge(g, get_reversed_edge(g, p.n[p.i]));
// 					}
// 				}
// 			}
// 		}
// 	}
// 	destroy_sn(&p);
// }

void rm_z_clip_lable(graph *g, const int s, const int m, const uint16_t l){	
	UINT_T t;
	UINTL_T i, j, n;
	pls ps;
	pl_ *p = NULL;
	init_pls(&ps, 100, s + 1);
	for (i = 1; i < g->ni; i++){
		ps.i = 0;
		for (j = 0; g->nd[i].od > 1 && j < g->nd[i].odm; j++){
			if (g->eg[g->nd[i].oe[j]].l & MFLAG_FIT) continue;
			p = &ps.pl[ps.i];
			p->i = p->perc = 0;
			p->ide = g->eg[g->nd[i].oe[j]].ide;
			p->sco = g->eg[g->nd[i].oe[j]].sco;
			if (g->eg[g->nd[i].oe[j]].l & l) p->perc++;
			
			p->n[p->i++] = g->nd[i].oe[j];
			n = g->eg[g->nd[i].oe[j]].ou;
			
			while (g->nd[n].id == 1 && g->nd[n].od == 1 && p->i <= s){
				t = get_validly_oe(g, n, 0);
				p->n[p->i++] = g->nd[n].oe[t];
				if (g->eg[g->nd[n].oe[t]].l & l) p->perc++;
				n = g->eg[g->nd[n].oe[t]].ou;
			}
			if (g->nd[n].id > 1 && p->perc) {
				// printf("#%u %d %u\n", p->perc, p->i,  p->perc * 10000 / p->i);
				p->perc = (uint32_t) p->perc * 10000 / p->i;
				ps.i++;
				if (ps.i >= ps.im) reallocate_pls(&ps, 100, s + 1);
			}
		}

		if (!ps.i) continue;
		qsort(ps.pl, ps.i, sizeof(pl_), sort_pls);
		// for (j = 0; j < ps.i; j++){
		// 	printf("%u %u %u %u\n",j, ps.pl[j].perc,ps.pl[j].ide,ps.pl[j].sco);
		// }
		for (j = 0; j < ps.i; j ++){
			p = &ps.pl[j];
			if (g->nd[g->eg[p->n[0]].in].od < 2 || g->nd[g->eg[p->n[p->i-1]].ou].id < 2) continue;

			if (m && p->i >= m){
				rm_edge(g, p->n[p->i-1]);
				rm_edge(g, get_reversed_edge(g, p->n[p->i-1]));
				rm_edge(g, p->n[0]);
				rm_edge(g, get_reversed_edge(g, p->n[0]));
			}else{
				while (--p->i >= 0){
					rm_edge(g, p->n[p->i]);
					rm_edge(g, get_reversed_edge(g, p->n[p->i]));
				}
			}
			// printf("del:%u\n",j );
		}
	}
	destroy_pls(&ps);
}

void cal_node_io_bstsc(graph *g, const int s){
	UINTL_T i;
	UINT_T j;
	uint64_t score, temp;
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].od >= s){
			for (score = j = 0; j < g->nd[i].odm; j++){
				if (g->eg[g->nd[i].oe[j]].l & MFLAG_FIT) continue;
				temp = (uint64_t) g->eg[g->nd[i].oe[j]].sco * g->eg[g->nd[i].oe[j]].ide;
				if (temp > score) score = temp;
			}
			g->nd[i].oe[g->nd[i].odm] = score / 10000;
			// printf("oe:%u %u %u %lu\n",i,g->nd[i].odm, g->nd[i].oe[g->nd[i].odm], score );
		}
		
		if (g->nd[i].id >= s){
			for (score = j = 0; j < g->nd[i].idm; j++){
				if (g->eg[g->nd[i].ie[j]].l & MFLAG_FIT) continue;
				temp = (uint64_t) g->eg[g->nd[i].ie[j]].sco * g->eg[g->nd[i].ie[j]].ide;
				if (temp > score) score = temp;
			}
			g->nd[i].ie[g->nd[i].idm] = score / 10000;
			// printf("ie:%u %u %u %lu\n",i, g->nd[i].idm, g->nd[i].ie[g->nd[i].idm], score );
		}
	}
}

void rm_z_clip_score(graph *g, const int s, const int m){
	
	UINT_T t;
	UINTL_T i, j, n;
	phs ps;
	ph_ *p = NULL;
	init_phs(&ps, 100, s + 1);
	for (i = 1; i < g->ni; i++){
		// if (i > get_reversed_node(g, i)) continue;
		// int k =0;
		// while (k++ < 2){
		// i=get_reversed_node(g, i);
		// printf("\nnode:%u %u %d\n",i,g->nd[i].name,s );
		ps.i = 0;
		for (j = 0; g->nd[i].od > 1 && j < g->nd[i].odm; j++){
			if (g->eg[g->nd[i].oe[j]].l & MFLAG_FIT) continue;
			p = &ps.ph[ps.i];
			p->i = p->sco = 0;
			p->n[p->i++] = g->nd[i].oe[j];
			n = g->eg[g->nd[i].oe[j]].ou;
			
			while (g->nd[n].id == 1 && g->nd[n].od == 1 && p->i < s){
				t = get_validly_oe(g, n, 0);
				p->n[p->i++] = g->nd[n].oe[t];
				n = g->eg[g->nd[n].oe[t]].ou;
				// printf("#%u(%d)#\n", g->nd[n].name,p->i);
			}
			if (g->nd[n].id > 1) {
				ps.i++;
				if (ps.i >= ps.im) reallocate_phs(&ps, 100, s + 1);
			}
		}

		if (ps.i > 1){
			for (j = 0; j < ps.i; j ++){
				p = &ps.ph[j];
				n = g->eg[p->n[0]].in;
				p->sco = (double) g->eg[p->n[0]].sco * g->eg[p->n[0]].ide / \
					g->nd[n].oe[g->nd[n].odm] / 100;// /10000 * 100
				n = g->eg[p->n[p->i-1]].ou;
				p->sco += (double) g->eg[p->n[p->i-1]].sco * g->eg[p->n[p->i-1]].ide / \
					g->nd[n].ie[g->nd[n].idm] / 100; // /10000 * 100
				p->sco -= g->nd[n].id * 5;
			}
			msort_phs(&ps);
		}

		// for (j = 0; j < ps.i; j ++){
		// 	p = &ps.ph[j];
		// 	printf("\n%u %f",j, p->sco);
		// 	for (n =0; n < p->i; n++){
		// 		printf(" [%u %u]", g->nd[g->eg[p->n[n]].in].name, g->nd[g->eg[p->n[n]].ou].name);
		// 	}
		// 	if (j == ps.i - 1) printf("\n");
		// }

		for (j = 0; j < ps.i; j ++){
			p = &ps.ph[j];
			if (g->nd[g->eg[p->n[0]].in].od < 2 || g->nd[g->eg[p->n[p->i-1]].ou].id < 2) continue;

			// printf("\n%u %ld",j, p->sco);
			// for (n =0; n < p->i; n++){
			// 	printf(" [%u %u]", g->nd[g->eg[p->n[n]].in].name, g->nd[g->eg[p->n[n]].ou].name);
			// }
			// if (j == ps.i - 1) printf("\n");

			if (m && p->i >= m){
				rm_edge(g, p->n[p->i-1]);
				rm_edge(g, get_reversed_edge(g, p->n[p->i-1]));
				rm_edge(g, p->n[0]);
				rm_edge(g, get_reversed_edge(g, p->n[0]));
			}else{
				while (--p->i >= 0){
					rm_edge(g, p->n[p->i]);
					rm_edge(g, get_reversed_edge(g, p->n[p->i]));
				}
			}
		}
	}
	// }
	destroy_phs(&ps);
	// exit(1);
}

void rm_sht_loop(graph *g, const int s){
	UINTL_T i, n;
	UINT_T t;
	ph p;
	init_sn(&p, s + 1);
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].od != 1 || g->nd[i].id != 1) continue;
		p.i = 0;
		t = get_validly_oe(g, i, 0);
		p.n[p.i++] = g->nd[i].oe[t];
		n = g->eg[g->nd[i].oe[t]].ou;
		while (g->nd[n].id == 1 && g->nd[n].od == 1 && p.i < s && n != i){
				t = get_validly_oe(g, n, 0);
				p.n[p.i++] = g->nd[n].oe[t];
				n = g->eg[g->nd[n].oe[t]].ou;
			}

		if (n == i){
			while (--p.i >= 0){
				rm_edge(g, p.n[p.i]);
				rm_edge(g, get_reversed_edge(g, p.n[p.i]));
			}
		}
	}
	destroy_sn(&p);
}

void rm_sht_bubble(graph *g, const int s){
	UINTL_T i, n, m;
	UINT_T t;
	ph p1, p2;
	init_sn(&p1, s + 1);
	init_sn(&p2, s + 1);
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].od != 2) continue;//TODO: need change for polyploid
		p1.i = 0;
		t = get_validly_oe(g, i, 0);
		p1.n[p1.i++] = g->nd[i].oe[t];
		n = g->eg[g->nd[i].oe[t]].ou;
		while (g->nd[n].id == 1 && g->nd[n].od == 1 && p1.i < s){
			t = get_validly_oe(g, n, 0);
			p1.n[p1.i++] = g->nd[n].oe[t];
			n = g->eg[g->nd[n].oe[t]].ou;
		}

		p2.i = 0;
		t = get_validly_oe(g, i, 1);
		p2.n[p2.i++] = g->nd[i].oe[t];
		m = g->eg[g->nd[i].oe[t]].ou;
		while (g->nd[m].id == 1 && g->nd[m].od == 1 && p2.i < s){
			t = get_validly_oe(g, m, 0);
			p2.n[p2.i++] = g->nd[m].oe[t];
			m = g->eg[g->nd[m].oe[t]].ou;
		}

		if (n == m){
			ph *p = p1.i > p2.i ? &p2 : &p1;
			while (--p->i >= 0){
				rm_edge(g, p->n[p->i]);
				rm_edge(g, get_reversed_edge(g, p->n[p->i]));
			}
		}
	}
	destroy_sn(&p1);
	destroy_sn(&p2);
}

void rm_end_loop(graph *g, const int s){
	UINTL_T i, n, m;
	UINT_T k, t;
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].id == 2 && g->nd[i].od == 1){
			k = 1;
			t = get_validly_oe(g, i, 0);
			m = n = g->eg[g->nd[i].oe[t]].ou;
			while (g->nd[n].id == 1 && g->nd[n].od == 1 && k++ <= s && n != i){
				m = n;
				t = get_validly_oe(g, n, 0);
				n = g->eg[g->nd[m].oe[t]].ou;
			}
			if (n == i){
				rm_edge(g, g->nd[m].oe[t]);
				rm_edge(g, get_reversed_edge(g, g->nd[m].oe[t]));
			}
		}
	}
}

static void mark_node_lable(graph *g, const uint32_t s, const uint16_t l){
	UINT_T i;
	node *v = &g->nd[s];
	for (i = 0; i < v->idm; i++){
		g->eg[v->ie[i]].l |= l;
	}
	for (i = 0; i < v->odm; i++){
		g->eg[v->oe[i]].l |= l;
	}
}

static void find_candnode_chim(graph *g, ph *p){
	UINTL_T i, ri;
	UINT_T j, k, l, n;
	node *v = NULL, *w = NULL;
	for (i = 1; i < g->ni; i++){

		if (i > get_reversed_node(g, i)) continue; // skip checked reversed node
		v = &g->nd[i];
		if (!v->od || !v->id) continue;
		n = 0;
		for (j = 0; j < v->odm; j++){
			if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
			if (g->nd[g->eg[v->oe[j]].ou].id >= 2) n = 1;
			g->nd[g->eg[v->oe[j]].ou].l = 1;
		}

		for (l = j = 0; j < v->idm && !l && n; j++){
			if (g->eg[v->ie[j]].l & MFLAG_FIT) continue;
			if (g->nd[g->eg[v->ie[j]].in].od >= 2) n = 2;
			w = &g->nd[g->eg[v->ie[j]].in];
			for (k = 0; k < w->odm && !l; k++){
			// for (k = 0; k < w->odm && l < 2; k++){
				if (g->eg[w->oe[k]].l & MFLAG_FIT) continue;
				if (g->nd[g->eg[w->oe[k]].ou].l){
					l ++;
				}
			}
		}

		// if (l < 2) {
		#ifdef P_VERSION
			if (n == 2 && !l){
		#else
			if (!l){
		#endif
			mark_node_lable(g, i, MFLAG_LQ);
			ri = get_reversed_node(g, i);
			mark_node_lable(g, ri, MFLAG_LQ);
			// printf("%u\n", g->nd[i].name);
		}
		
		#ifdef P_VERSION
			if (n == 2) {
		#else
			if (n == 2 && !l) {
		#endif
			p->n[p->i++] = i;
			if (p->i >= p->im) reallocate_sn(p, 1000);
			mark_node_lable(g, i, MFLAG_CC);
			ri = get_reversed_node(g, i);
			mark_node_lable(g, ri, MFLAG_CC);
			// printf("%u\n", g->nd[i].name);
		}

		for (j = 0; j < v->odm; j++){
			if (g->eg[v->oe[j]].l & MFLAG_FIT) continue;
			g->nd[g->eg[v->oe[j]].ou].l = 0;
		}
	}
}

#if GENOME_SIZE == 2
	KHASH_SET_INIT_INT64(set_)
#else
	KHASH_SET_INIT_INT(set_)
#endif
void mark_node_chim(graph *g, const int s, const int m, const int md, int *kc){
	UINTL_T i, t = 0;
	UINT_T j, n;
	// float maxd = (float) d->d[1]/d->c[1] * m * (1 + (float)d->u[1]/d->c[1]);
	int maxd = md * m;
	node *v = NULL;
	edge *e = NULL;
	ph p, p1, p2, *pp1, *pp2;
	init_sn(&p, 5000);
	find_candnode_chim(g, &p);
	int64_t pi = p.i;

	khint_t k;
	int absent;
	khash_t(set_) *set1 = kh_init(set_);
	khash_t(set_) *set2 = kh_init(set_);
	init_sn(&p1, 500);
	init_sn(&p2, 500);
	pp1 = &p1;
	pp2 = &p2;
	while (--p.i >= 0){
		pp1->i = 0;
		v = &g->nd[p.n[p.i]];
		for (i = 0; i < v->idm; i++){
			if (g->eg[v->ie[i]].l & MFLAG_FIT) continue;
			pp1->n[pp1->i++] = g->eg[v->ie[i]].in;
			if (pp1->i >= pp1->im) reallocate_sn(pp1, 100);
			k = kh_put(set_, set1, g->eg[v->ie[i]].in, &absent);
			kh_key(set1, k) = g->eg[v->ie[i]].in;
		}

		for (i = 0; i < s; i++){
			pp2->i = 0;
			while (--pp1->i >= 0){
				v = &g->nd[pp1->n[pp1->i]];
				for (j = 0; j < v->odm; j++){
					e = &g->eg[v->oe[j]];
					if (e->l & MFLAG_FIT || e->ou == p.n[p.i]) continue;
					if (g->nd[e->ou].od >= maxd) continue;
					k = kh_put(set_, set1, e->ou, &absent);
					if (absent) {
						kh_key(set1, k) = pp2->n[pp2->i++] = e->ou;
						if (pp2->i >= pp2->im) reallocate_sn(pp2, 100);
					}
				}
			}
			SWAP(pp1, pp2, ph*);
		}
		n = pp1->i = 0;
		v = &g->nd[p.n[p.i]];
		for (i = 0; i < v->odm; i++){
			if (g->eg[v->oe[i]].l & MFLAG_FIT) continue;
			pp1->n[pp1->i++] = g->eg[v->oe[i]].ou;
			if (pp1->i >= pp1->im) reallocate_sn(pp1, 100);
			if (kh_get(set_, set1, g->eg[v->oe[i]].ou) != kh_end(set1)) {
				n = 1;
				break;
			}
			kh_put(set_, set2, g->eg[v->oe[i]].ou, &absent);
		}

		for (i = 0; i < s && !n; i++){
			pp2->i = 0;
			while (--pp1->i >= 0 && !n){
				v = &g->nd[pp1->n[pp1->i]];
				for (j = 0; j < v->idm; j++){
					e = &g->eg[v->ie[j]];
					if (e->l & MFLAG_FIT || e->in == p.n[p.i]) continue;
					if (kh_get(set_, set1, e->in) != kh_end(set1)) {
						n = 1;
						break;
					}
					if (g->nd[e->in].id >= maxd) continue;
					kh_put(set_, set2, e->in, &absent);
					if (absent) {
						pp2->n[pp2->i++] = e->in;
						if (pp2->i >= pp2->im) reallocate_sn(pp2, 100);
					}
				}
			}
			SWAP(pp1, pp2, ph*);
		}
		kh_clear(set_, set1);
		kh_clear(set_, set2);
		if (!n) {
			// printf("%u\n", g->nd[p.n[p.i]].name);
			t ++;
			mark_node_lable(g, p.n[p.i], MFLAG_CN);
			i = get_reversed_node(g, p.n[p.i]);
			mark_node_lable(g, i, MFLAG_CN);
		}
	}
	plog(2, "Chimeric node ratio: %.3f%% (candidate: %.3f%%)", (float) t/g->ni/2 * 100, (float) pi/g->ni/2 * 100);
	// if ((float) pi/g->ni/2 > 0.05) *kc = 1;//TODO CHECK
	kh_destroy(set_, set1);
	kh_destroy(set_, set2);
	destroy_sn(&p1);
	destroy_sn(&p2);
	destroy_sn(&p);
}

void calc_edge_tc(graph *g, const int s){
	UINTL_T i, j;
	edge *e;
	node *v;
	for (i = 0; i < g->ei; i++){
		e = &g->eg[i];
		if (e->tc || e->l & MFLAG_FIT) continue;
		v = &g->nd[e->in];
		for (j = 0; j < v->idm; j ++) g->nd[g->eg[v->ie[j]].in].l = 1;
		for (j = 0; j < v->odm; j ++) g->nd[g->eg[v->oe[j]].ou].l = 2;

		v = &g->nd[e->ou];
		for (j = 0; j < v->idm && e->tc <= s; j ++){
			if (g->nd[g->eg[v->ie[j]].in].l) e->tc ++;
		}
		for (j = 0; j < v->odm && e->tc <= s; j ++){
			if (g->nd[g->eg[v->oe[j]].ou].l == 2) e->tc ++;
		}

		v = &g->nd[e->in];
		for (j = 0; j < v->idm; j ++) g->nd[g->eg[v->ie[j]].in].l = 0;
		for (j = 0; j < v->odm; j ++) g->nd[g->eg[v->oe[j]].ou].l = 0;
		// printf("%u %u %d\n",g->nd[e->in].name, g->nd[e->ou].name, e->tc);
	}
}

void rm_edge_chim(graph *g, const int s, const float sco, const uint32_t ide, const uint16_t l){
	UINTL_T i, n;
	UINT_T j, t;
	uint32_t rlen;
	node *v = NULL;
	ph p;
	init_sn(&p, 5000);
	for (i = 1; i < g->ni; i++){
		v = &g->nd[i];
		if (!v->od) continue;
		// t = get_validly_oe(g, i, 0);
		if(!check_node_lable(g, i, l)) continue;

		j = 0; n = i;
		while (g->nd[n].id == 1 && g->nd[n].od == 1 && ++j < s){
			t = get_validly_oe(g, n, 0);
			n = g->eg[g->nd[n].oe[t]].ou;
		}
		// printf("%u %d %d\n", v->name, j, s);
		if (j >= s) continue;
		
		for (n = 0; n < v->odm; n++){
			if (g->eg[v->oe[n]].l & MFLAG_FIT) continue;
			rlen = min(get_ireadlen(g, v->oe[n]), get_oreadlen(g, v->oe[n]));
	// printf("%u %u %u %f %d\n",g->nd[g->eg[v->oe[n]].in].name,g->nd[g->eg[v->oe[n]].ou].name, g->eg[v->oe[n]].ide,(double)g->eg[v->oe[n]].sco/rlen,  g->eg[v->oe[n]].tc);
			if (((g->eg[v->oe[n]].ide < ide || g->eg[v->oe[n]].l & MFLAG_LS) && g->eg[v->oe[n]].sco < sco * rlen) \
					|| !g->eg[v->oe[n]].tc){
				p.n[p.i++] = v->oe[n];
				if (p.i >= p.im) reallocate_sn(&p, 1000);
				// rm_edge(g, v->oe[n]);
				// rm_edge(g, get_reversed_edge(g, v->oe[n]));
			}
		}

		for (n = 0; n < v->idm; n++){
			if (g->eg[v->ie[n]].l & MFLAG_FIT) continue;
			rlen = min(get_ireadlen(g, v->ie[n]), get_oreadlen(g, v->ie[n]));
		// printf("%u %u %u %f %d\n",g->nd[g->eg[v->ie[n]].in].name,g->nd[g->eg[v->ie[n]].ou].name, g->eg[v->ie[n]].ide,(double)g->eg[v->ie[n]].sco/rlen,  g->eg[v->ie[n]].tc);
			if (((g->eg[v->ie[n]].ide < ide || g->eg[v->ie[n]].l & MFLAG_LS) && g->eg[v->ie[n]].sco < sco * rlen) \
					|| !g->eg[v->ie[n]].tc){
				p.n[p.i++] = v->ie[n];
				if (p.i >= p.im) reallocate_sn(&p, 1000);
				// rm_edge(g, v->ie[n]);
				// rm_edge(g, get_reversed_edge(g, v->ie[n]));
			}
		}

		// if (j == 0 && (v->od != 1 || v->id != 1)){
		// 	rm_node(g, i);
		// 	rm_node(g, get_reversed_node(g, i));
		// }
	}
	for (i = 0; i < p.i; i ++){
		rm_edge(g, p.n[i]);
		rm_edge(g, get_reversed_edge(g, p.n[i]));
	}
	destroy_sn(&p);
	// exit(1);
}

void rm_edge_ltc(graph *g, const int s, const float sco, const float idt, const int l){
	UINTL_T i, n;
	UINT_T j, t;
	uint32_t rlen, mide, msco, mide_, msco_;
	edge *e;
	node *v;
	ph p;
	init_sn(&p, 5000);
	for (i = 0; i < g->ei; i++){
		e = &g->eg[i];
		if (e->tc > l || e->l & MFLAG_FIT) continue;

		j = 0; n = e->in;
		while (g->nd[n].id == 1 && g->nd[n].od == 1 && ++j < s){
			t = get_validly_oe(g, n, 0);
			n = g->eg[g->nd[n].oe[t]].ou;
		}
		if (j >= s) continue;

		rlen = min(get_ireadlen(g, i), get_oreadlen(g, i));
		if (e->sco < sco * rlen){
			p.n[p.i++] = i;
			if (p.i >= p.im) reallocate_sn(&p, 1000);
			continue;
		}

		v = &g->nd[e->in];
		for (mide = msco = j = 0; j < v->odm; j ++) {
			if (g->eg[v->oe[j]].ide > mide) mide = g->eg[v->oe[j]].ide;
			if (g->eg[v->oe[j]].sco > msco) msco = g->eg[v->oe[j]].sco;
		}
		v = &g->nd[e->ou];
		for (mide_ = msco_ = j = 0; j < v->idm; j ++){
			if (g->eg[v->ie[j]].ide > mide_) mide_ = g->eg[v->ie[j]].ide;
			if (g->eg[v->ie[j]].sco > msco_) msco_ = g->eg[v->ie[j]].sco;
		}
		if (mide_ < mide) mide = mide_;
		if (msco_ < msco) msco = msco_;

		if (e->ide < mide * idt && e->sco < msco * idt){
			p.n[p.i++] = i;
			if (p.i >= p.im) reallocate_sn(&p, 1000);
		}
	}
	for (i = 0; i < p.i; i ++){
		rm_edge(g, p.n[i]);
		rm_edge(g, get_reversed_edge(g, p.n[i]));
	}
	destroy_sn(&p);
}

#if GENOME_SIZE == 2
	KHASH_MAP_INIT_INT64(cpath_info, cpath_info);
#else
	KHASH_MAP_INIT_INT(cpath_info, cpath_info);
#endif
void clean_complex_path(graph *g, const int s, const int r, const uint16_t l){
	khint_t k;
	int absent;
	cpath_info *bsv = NULL;

	ph p; //visited edge
	init_sn(&p, 100);

	UINTL_T i, d, v, w, pending;
	UINT_T j, m, loop, exclen;

	UINTL_T qi, ql, qm = 200;
	UINTL_T *q = malloc(qm * sizeof(UINTL_T));
	khash_t(cpath_info) *bs = kh_init(cpath_info);
	size_t bss = sizeof(cpath_info);

	for (i = 1; i < g->ni; i++){
		if (g->nd[i].od > 1){
			kh_clear(cpath_info, bs);
			loop = exclen = pending = p.i = qi = ql = 0;			
			q[ql++] = i;
			k = kh_put(cpath_info, bs, i, &absent);
			memset(&kh_val(bs, k), 0, bss);
			kh_key(bs, k) = i;
			while (1){
				v = q[qi++];
				bsv = &kh_val(bs, kh_get(cpath_info, bs, v));
				d = bsv->len;
				for (j = 0; j < g->nd[v].odm; j ++){
					if (g->eg[g->nd[v].oe[j]].l & MFLAG_FIT) continue;
					if (g->eg[g->nd[v].oe[j]].l & l){
						p.n[p.i++] = g->nd[v].oe[j];
						if (p.i >= p.im) reallocate_sn(&p, 100);
						continue;
					}
					w = g->eg[g->nd[v].oe[j]].ou;
					if (w == i) {loop = 1; break;}
					// if (g->eg[g->nd[v].oe[j]].len + d >= s){exclen = 1; break;}
					if (d + 1 >= s){exclen = 1; break;}

					p.n[p.i++] = g->nd[v].oe[j];
					if (p.i >= p.im) reallocate_sn(&p, 100);

					k = kh_get(cpath_info, bs, w);
					bsv = NULL;
					if (k == kh_end(bs)){
						k = kh_put(cpath_info, bs, w, &absent);
						bsv = &kh_val(bs, k);
						memset(bsv, 0, bss);
						kh_key(bs, k) = w;
						pending ++;
						for (m = 0; m < g->nd[w].idm; m++){
							if (g->eg[g->nd[w].ie[m]].l & MFLAG_FIT) continue;
							if (g->eg[g->nd[w].ie[m]].l & l) continue;
							bsv->id ++;
						}
					}

					if (!bsv) bsv = &kh_val(bs, k);
					bsv->pnode = v;
					bsv->id --;
					bsv->len = d + 1;
					if (!bsv->id){
						if (g->nd[w].od) {
							q[ql++] = w;
							if (ql >= qm) {
								qm += 200;
								q = realloc(q, qm * sizeof(UINTL_T));
							}
						}
						pending --;
					}
				}
				if (loop || exclen || ql == qi) break;
				if (ql == qi + 1 && !pending) break;
			}
			if (loop || exclen || ql == qi) continue;

			v = q[qi++];
			while (1){
				bsv = &kh_val(bs, kh_get(cpath_info, bs, v));
				bsv->l = 1;
				w = bsv->pnode;
				if (r){
					for (m = 0; m < g->nd[v].idm; m++){
						if (g->eg[g->nd[v].ie[m]].l & MFLAG_FIT) continue;
						d = g->nd[v].ie[m];
						if (g->eg[d].in != w) {
							rm_edge(g, d);
							rm_edge(g, get_reversed_edge(g, d));
						}
					}

					for (m = 0; m < g->nd[w].odm; m++){
						if (g->eg[g->nd[w].oe[m]].l & MFLAG_FIT) continue;
						d = g->nd[w].oe[m];
						if (g->eg[d].ou != v){
							rm_edge(g, d);
							rm_edge(g, get_reversed_edge(g, d));
						}
					}
				}
				v = w;
				if (v == i) break;
			}
			if (!r){
				while (--p.i >= 0){
					if (kh_val(bs, kh_get(cpath_info, bs, g->eg[p.n[p.i]].in)).l) continue;
					rm_edge(g, p.n[p.i]);
					rm_edge(g, get_reversed_edge(g, p.n[p.i]));
				}
			}
		}
	}
	kh_destroy(cpath_info, bs);
	destroy_sn(&p);
	free (q);
}

void update_graph(char *file, graph *g, khash_t(ovlinfo_) *os, opt *p, uint32_t *arr, list *li){
	khiter_t k;
	FILE *fp = fopen(file, "r");
	ovlinfo *loli = NULL, *roli = NULL;
	buffer_t *buf = init_buffer(DCBUFSIZE);
	prev_t id = {0, 0};
	int v, mode = find_ovlb_mode(fp);
	uint32_t l, qlen, tlen, ovl[10], alnlen;
	UINTL_T edge;
	qlen = tlen = 0;
	assert (mode == 10);

	while(decode_ovl(fp, arr, &id, ovl, buf, mode) >= 0) {
		// printf("%u %u %u %u %u %u %u %u\n", ovl[0], ovl[1],ovl[2],ovl[3],ovl[4],ovl[5],ovl[6], g->ni);
		if (ovl[7]) qlen = ovl[7];
		else ovl[7] = qlen;
		if (ovl[8]) tlen = ovl[8];
		else ovl[8] = tlen;
		if (ovl[9] == 0) ovl[9] ++;

		k = kh_get(ovlinfo_, os, li->data[ovl[0]]);
		if (k == kh_end(os)) continue;
		else loli = &kh_val(os, k);

		k = kh_get(ovlinfo_, os, li->data[ovl[4]]);
		if (k == kh_end(os)) continue;
		else roli = &kh_val(os, k);

		if (loli->con < MAX_CON && ovl[2] <= p->max_hang_len + loli->le && \
			ovl[3] >= ovl[7] - loli->re - p->max_hang_len) {//fixed
			loli->con ++;
			continue;
		}else if (roli->con < MAX_CON && ovl[5] <= p->max_hang_len + roli->le && \
			ovl[6] >= ovl[8] - roli->re - p->max_hang_len){//fixed
			roli->con ++;
			continue;
		}
		if (loli->con >= p->min_con_count || roli->con >= p->min_con_count) continue;

		l = v = 0;
		alnlen = max(ovl[3] - ovl[2], ovl[6] - ovl[5]);

		if (ovl[1]){
			if (loli->rid && roli->lid && (l = check_exited_edge(g, loli->rid, roli->lid)) >= alnlen) continue;// skip supp. alignment
			// printf("%d\n", l);
			if (ovl[2] <= p->max_hang_len + loli->le && ovl[5] <= p->max_hang_len + roli->le){
				if (alnlen >= loli->llm) v = 2;
				else{ 
					if (loli->lim >= p->min_ide){
						if (ovl[9] >= loli->lim * p->min_ide_ratio) v++;
					}else{
						if (alnlen >= loli->llm * p->min_sco_ratio) v++;
					}
				}
				if (alnlen >= roli->llm) v = 2;
				else{
					if (roli->lim >= p->min_ide){
						if (ovl[9] >= roli->lim * p->min_ide_ratio) v++;
					}else{
						if (alnlen > roli->llm * p->min_sco_ratio) v++;
					}
				}


				// if (ovl[0] == 676631 && ovl[4] == 2427435){
				// 	printf("%u %u %u %u %u %u %f %f %u %f %f %d\n",ovl[0],ovl[4],ovl[9],alnlen,
				// 		p->min_ide, loli->lim, loli->lim * p->min_ide_ratio, 
				// 		loli->llm * p->min_sco_ratio,roli->lim, 
				// 		roli->lim * p->min_ide_ratio,roli->llm * p->min_sco_ratio,v  );exit(1);
				// }

				if (v < p->min_node_count) continue;
				// should update lm & im after updating align edges
				if (alnlen > loli->llm) loli->llm = alnlen;
				if (alnlen > roli->llm) roli->llm = alnlen;
				if (ovl[9] > loli->lim) loli->lim = ovl[9];
				if (ovl[9] > roli->lim) roli->lim = ovl[9];

				if (!loli->rid) {
					loli->rid = add_node(g, ovl[0], 8);
					loli->lid = add_node(g, ovl[0], 8);
				}
				if (!roli->rid) {
					roli->rid = add_node(g, ovl[4], 8);
					roli->lid = add_node(g, ovl[4], 8);
				}
				
				edge = add_edge(g, loli->lid, roli->rid, ovl[2], ovl[6]-1, alnlen, ovl[9], ovl[8] - roli->re - ovl[6], l);//fixed
				g->eg[edge].l |= MFLAG_IL;	
				edge = add_edge(g, roli->lid, loli->rid, ovl[5], ovl[3]-1, alnlen, ovl[9], ovl[7] - loli->re - ovl[3], l);//fixed
				g->eg[edge].l |= MFLAG_IL;
			}else if (ovl[3] >= ovl[7] - loli->re - p->max_hang_len && ovl[6] >= ovl[8] - roli->re - p->max_hang_len){//fixed

				if (alnlen >= loli->rlm) v = 2;
				else{ 
					if (loli->rim >= p->min_ide){
						if (ovl[9] >= loli->rim * p->min_ide_ratio) v++;
					}else{
						if (alnlen >= loli->rlm * p->min_sco_ratio) v++;
					}
				}

				if (alnlen >= roli->rlm) v = 2;
				else{ 
					if (roli->rim >= p->min_ide){
						if (ovl[9] >= roli->rim * p->min_ide_ratio) v++;
					}else{
						if (alnlen >= roli->rlm * p->min_sco_ratio) v++;
					}
				}
				if (v < p->min_node_count) continue;

				if (alnlen > loli->rlm) loli->rlm = alnlen;
				if (alnlen > roli->rlm) roli->rlm = alnlen;
				if (ovl[9] > loli->rim) loli->rim = ovl[9];
				if (ovl[9] > roli->rim) roli->rim = ovl[9];

				if (!loli->rid) {
					loli->rid = add_node(g, ovl[0], 8);
					loli->lid = add_node(g, ovl[0], 8);
				}
				if (!roli->rid) {
					roli->rid = add_node(g, ovl[4], 8);
					roli->lid = add_node(g, ovl[4], 8);
				}
				edge = add_edge(g, loli->rid, roli->lid, ovl[3]-1, ovl[5], alnlen, ovl[9], ovl[5] - roli->le, l);
				g->eg[edge].l |= MFLAG_OL;
				edge = add_edge(g, roli->rid, loli->lid, ovl[6]-1, ovl[2], alnlen, ovl[9], ovl[2] - loli->le, l);
				g->eg[edge].l |= MFLAG_OL;
			}
		}else{
			if (loli->rid && roli->rid && (l = check_exited_edge(g, loli->rid, roli->rid)) >= alnlen) continue;// skip supp. alignment

			if (ovl[3] >= ovl[7] - loli->re - p->max_hang_len && ovl[5] <= p->max_hang_len + roli->le){//fixed
				if (alnlen >= loli->rlm) v = 2;
				else {
					if (loli->rim >= p->min_ide){
						if (ovl[9] >= loli->rim * p->min_ide_ratio) v++;
					}else{
						if (alnlen >= loli->rlm * p->min_sco_ratio) v++;
					}
				}
				if (alnlen >= roli->llm) v = 2;
				else{ 
					if (roli->lim >= p->min_ide){
						if (ovl[9] >= roli->lim * p->min_ide_ratio) v++;
					}else{
						if (alnlen >= roli->llm * p->min_sco_ratio) v++;
					}
				}
				if (v < p->min_node_count) continue;

				if (alnlen > loli->rlm) loli->rlm = alnlen;
				if (alnlen > roli->llm) roli->llm = alnlen;
				if (ovl[9] > loli->rim) loli->rim = ovl[9];
				if (ovl[9] > roli->lim) roli->lim = ovl[9];

				if (!loli->rid) {
					loli->rid = add_node(g, ovl[0], 8);
					loli->lid = add_node(g, ovl[0], 8);
				}
				if (!roli->rid) {
					roli->rid = add_node(g, ovl[4], 8);
					roli->lid = add_node(g, ovl[4], 8);
				}
				add_edge(g, loli->rid, roli->rid, ovl[3]-1, ovl[6]-1, alnlen, ovl[9], ovl[8] - roli->re - ovl[6], l);//fixed
				edge = add_edge(g, roli->lid, loli->lid, ovl[5], ovl[2], alnlen, ovl[9], ovl[2] - loli->le, l);
				g->eg[edge].l |= MFLAG_IL;
				g->eg[edge].l |= MFLAG_OL;
			}else if (ovl[2] <= p->max_hang_len + loli->le && ovl[6] >= ovl[8] - roli->re - p->max_hang_len){//fixed
				if (alnlen >= loli->llm) v = 2;
				else{ 
					if (loli->lim >= p->min_ide){
						if (ovl[9] >= loli->lim * p->min_ide_ratio) v++;
					}else{
						if (alnlen >= loli->llm * p->min_sco_ratio) v++;
					}
				}

				if (alnlen >= roli->rlm) v = 2;
				else{ 
					if (roli->rim >= p->min_ide){
						if (ovl[9] >= roli->rim * p->min_ide_ratio) v++;
					}else{
						if (alnlen >= roli->rlm * p->min_sco_ratio) v++;
					}
				}
				if (v < p->min_node_count) continue;

				if (alnlen > loli->llm) loli->llm = alnlen;
				if (alnlen > roli->rlm) roli->rlm = alnlen;
				if (ovl[9] > loli->lim) loli->lim = ovl[9];
				if (ovl[9] > roli->rim) roli->rim = ovl[9];

				if (!loli->rid) {
					loli->rid = add_node(g, ovl[0], 8);
					loli->lid = add_node(g, ovl[0], 8);
				}
				if (!roli->rid) {
					roli->rid = add_node(g, ovl[4], 8);
					roli->lid = add_node(g, ovl[4], 8);
				}
				add_edge(g, roli->rid, loli->rid, ovl[6]-1, ovl[3]-1, alnlen, ovl[9], ovl[7] - loli->re - ovl[3], l);//fixed
				edge = add_edge(g, loli->lid, roli->lid, ovl[2], ovl[5], alnlen, ovl[9], ovl[5] - roli->le, l);
				g->eg[edge].l |= MFLAG_IL;
				g->eg[edge].l |= MFLAG_OL;
			}
		}
	}
	fclose (fp);
	flush_buffer(NULL, buf);
}


////////////////////////


static UINTL_T stat_bfs_l(khash_t(bfs) *h, uint16_t l, UINTL_T *q, UINTL_T ql){
	UINTL_T c = 0;
	khiter_t k;
	if (!q){
		for (k = kh_begin(h); k != kh_end(h); k++){
			if (kh_exist(h, k) && kh_val(h, k).l & l) c++;
		}
	}else{
		UINTL_T i;
		for (i = 0; i < ql; i ++){
			k = kh_get(bfs, h, q[i]);
			if (k != kh_end(h) && kh_val(h, k).l & l) c++;
		}
	}
	return c;
}


static int cmp_int(const void * a, const void * b){
    UINTL_T *_a = (UINTL_T *) a;
    UINTL_T *_b = (UINTL_T *) b;
    return (*_a > *_b) - (*_a < *_b);
    // return *_a - *_b;
}

static int binary_search(UINTL_T *n, int64_t n_l, int64_t d) { 
	int64_t low = 0, high = n_l, middle = 0;
	while (low < high) {
		middle = (low + high)/2;
		// printf("\n%d %d %d\n", d,n[middle],middle);
		if (d == n[middle]){
			return 1;
		}else if (d < n[middle]) { 
			high = middle;
		}else if (d > n[middle]) {
			low = middle + 1;
		}
	}
	return 0;
}

static int check_in_exclude(ph *h, UINTL_T n){
	if(!h) return 0;
	else if (h->im){//TODO CHECK 
		qsort(h->n, h->i, sizeof(UINTL_T), cmp_int);
		h->im = 0;
	}
	int ret = binary_search(h->n, h->i, n);
	return ret;
}

static int check_in_exclude_not_sort(ph *h, UINTL_T n){
	if(!h) return 0;
	UINTL_T i;
	for (i = 0; i < h->i; i ++) {
		if (h->n[i] == n) return 1;
	}
	return 0;
}

void *bfs_nodes_compound_orig(graph *g, ph *n, ph *exclude, const int depth,  
		const int max_child_num, const int lastest, const int stop_at_merged, 
		void *(*callback)(UINTL_T, UINTL_T, khash_t(bfs) *, void *), void *callback_opt){
	
	UINTL_T i, j, v, w;
	UINTL_T qi = 0, ql = 0, qm = n->i + 200;
	UINTL_T *q = malloc(qm * sizeof(UINTL_T));
	
	khint_t k;
	bfs *bfs_v, *bfs_w;
	bfs_r *bfs_ret = malloc(sizeof(bfs_r));
	size_t bfs_size = sizeof(bfs);
	khash_t(bfs) *bfs_info = bfs_ret->bfs_info = kh_init(bfs);
	ph *visited_edges = bfs_ret->data1 = malloc(sizeof(ph));
	init_sn(visited_edges, 100);
	
	int absent;
	for (i = 0; i < n->i; i++){
		q[ql++] = n->n[i];
		// printf("###########s:%d##############\n", g->nd[n->n[i]].name);
		k = kh_put(bfs, bfs_info, n->n[i], &absent);
		kh_key(bfs_info, k) = n->n[i];
		bfs_v = &kh_val(bfs_info, k);
		memset(bfs_v, 0, bfs_size);
		bfs_v->l |= BFLAG_N;
	}

	while (ql > qi && (!stop_at_merged || stat_bfs_l(bfs_info, BFLAG_N, q + qi, ql - qi) || 
			(ql - qi + stat_bfs_l(bfs_info, BFLAG_P2, NULL, 0)) > 1)){
		// printf("%d %d %d %d\n",ql > qi,!stop_at_merged,stat_bfs_l(bfs_info, BFLAG_N, q + qi, ql-qi+1),(ql - qi + stat_bfs_l(bfs_info, BFLAG_P2, NULL, 0)) > 1  );
		v = q[qi++];
		bfs_v = &kh_val(bfs_info, kh_get(bfs, bfs_info, v));
		for (i = 0; i < g->nd[v].odm; i ++){
			if (g->eg[g->nd[v].oe[i]].l & MFLAG_FIT) continue;
			w = g->eg[g->nd[v].oe[i]].ou;

// int tt = 0;
// printf("exclude: ");
// for (tt=0;exclude && tt<exclude->i;tt++){printf(" %u(%u,%d)",g->nd[exclude->n[tt]].name,exclude->n[tt],exclude->im );}
// printf("\n%u(%u) %u(%u) %d %d\n",g->nd[v].name,v,  g->nd[w].name,w, check_in_exclude(exclude, w));

			if (check_in_exclude(exclude, w)) continue;

			visited_edges->n[visited_edges->i++] = g->nd[v].oe[i];
			if (visited_edges->i >= visited_edges->im) reallocate_sn(visited_edges, 100);
			
			k = kh_get(bfs, bfs_info, w);
			if (k == kh_end(bfs_info) || (kh_val(bfs_info, k).l & BFLAG_N && !(kh_val(bfs_info, k).l & BFLAG_LP))){
				UINT_T in_degree_w = 0;
				for (j = 0; j < g->nd[w].idm; j++){
					if (g->eg[g->nd[w].ie[j]].l & MFLAG_FIT) continue;
					if (!check_in_exclude(exclude, g->eg[g->nd[w].ie[j]].in)) in_degree_w++;
				}

				if (k == kh_end(bfs_info)){
					k = kh_put(bfs, bfs_info, w, &absent);
					kh_key(bfs_info, k) = w;
					bfs_w = &kh_val(bfs_info, k);
					memset(bfs_w, 0, bfs_size);
					bfs_v = &kh_val(bfs_info, kh_get(bfs, bfs_info, v));
				}else bfs_w = &kh_val(bfs_info, k);

				bfs_w->bfs_depth = bfs_v->bfs_depth + 1;
				bfs_w->unvisited = in_degree_w;
				if (bfs_w->l & BFLAG_N) bfs_w->l |= BFLAG_LP;
			}else{
				bfs_w = &kh_val(bfs_info, k);
				if (bfs_w->l & BFLAG_N) bfs_w->l |= BFLAG_LP;
				if (lastest) bfs_w->bfs_depth = bfs_v->bfs_depth + 1;
			}

			if (callback){
				void *callback_ret = callback(v, w, bfs_info, callback_opt);
				if (callback_ret) {
					bfs_ret->data1 = callback_ret;
					destroy_sn(visited_edges);
					free (visited_edges);
					goto END;
				}
			}

			bfs_w->unvisited --;
			UINT_T out_degree_w = 0;
			for (j = 0; j < g->nd[w].odm; j++){
				if (g->eg[g->nd[w].oe[j]].l & MFLAG_FIT) continue;
				if (!check_in_exclude(exclude, g->eg[g->nd[w].oe[j]].ou)) out_degree_w++;
			}
			// printf("1:%d %d %d %d %d\n",g->nd[w].name, bfs_w->unvisited,out_degree_w,!(bfs_w->l & BFLAG_N),bfs_w->bfs_depth);	
			if (bfs_w->unvisited == 0 && out_degree_w && (out_degree_w < max_child_num || max_child_num < 0) &&
				(!(bfs_w->l & BFLAG_N)) && (bfs_w->bfs_depth < depth || depth < 0)){
				q[ql++] = w;
				if (ql >= qm) {
					qm += 200;
					q = realloc(q, qm * sizeof(UINTL_T));
				}
				
				if (bfs_w->l & BFLAG_P2) bfs_w->l ^= BFLAG_P2;
			}else if (bfs_w->unvisited == 0 && out_degree_w == 0) bfs_w->l |= BFLAG_P1;
			else bfs_w->l |= BFLAG_P2;
		}
		// printf("2:%d %d %d %d\n", ql - qi,!stop_at_merged,stat_bfs_l(bfs_info, BFLAG_N, q + qi, ql - qi+1),ql - qi + stat_bfs_l(bfs_info, BFLAG_P2, NULL, 0));
	}//1 False False 1

	while (ql > qi){
		v = q[qi++];
		kh_val(bfs_info, kh_get(bfs, bfs_info, v)).l |= BFLAG_U;
	}
END:
	// kh_destroy(bfs, bfs_info);
	free (q);
	return bfs_ret;
}

// void reverse_ph(ph *p){
// 	UINTL_T tmp;
// 	UINTL_T *p1 = p->n;
// 	UINTL_T *p2 = p->n + p->i - 1;
// 	while (p1 < p2){
// 		tmp = *p1;
// 		*p1++ = *p2;
// 		*p2-- = tmp;
// 	}
// }

void reverse_ph2(ph_ *p){
    UINTL_T tmp; 
    UINTL_T *p1 = p->n;
    UINTL_T *p2 = p->n + p->i - 1; 
    while (p1 < p2){ 
        tmp = *p1;
        *p1++ = *p2;
        *p2-- = tmp;
    }
}

phs *bfs_nodes_compound_get_path(khash_t(bfs) *h, int s){
	ph_ *p;
	phs *ps = malloc(sizeof(phs));
	
	if (s <= 0) s = 100;
	init_phs(ps, 20, s + 1);

// if (kh_get(bfs, h, 73357) != kh_end(h) && kh_val(h, kh_get(bfs, h, 73357)).info){ 
// 	printf ("461100(73357)_p:%u",g->nd[kh_val(h, kh_get(bfs, h, 73357)).info[2].predecessor].name);}
// else printf("######\n");

	uint16_t l = 0;
	l |= (BFLAG_P1 | BFLAG_P2 | BFLAG_U);
	bfs *h_v, *h_w;
	UINTL_T i, v, w, start;
	khiter_t k;
	for (k = kh_begin(h); k != kh_end(h); k++){
		if (kh_exist(h, k)){
			v = kh_key(h, k);
			h_v = &kh_val(h, k);
			if (h_v->l & l) {
				for (i = 0; i < h_v->info_len; i ++){
					if (!h_v->info[i].edge_num) continue;
				
				// printf("i:%d v:%u(%u) v_p:%u start:%u\n", i, g->nd[v].name,v, g->nd[h_v->info[i].predecessor].name, g->nd[h_v->info[i].start].name);

					start = h_v->info[i].start;
					if (ps->i >= ps->im) reallocate_phs(ps, 20, s + 1);
					p = &ps->ph[ps->i++];
					p->i = 0;
					p->n[p->i++] = v;
					p->n[p->i++] = w = h_v->info[i].predecessor;
					while (w != start) {
						h_w = &kh_val(h, kh_get(bfs, h, w));
						w = h_w->info[i].predecessor;
						p->n[p->i++] = w;
						if (p->i >= p->im) reallocate_ph(p, 100);
					}
					reverse_ph2(p);
				}
			}
		}
	}
	return ps;
}
static int bfs_node_info_cmp(const void * a_, const void * b_){
	bfs_node_info * a = (bfs_node_info *) a_;
	bfs_node_info * b = (bfs_node_info *) b_;
	return (a->start > b->start) - (a->start < b->start);
	// return ((bfs_node_info *) a)->start - ((bfs_node_info *) b)->start;
}

bfs_node_info *init_bfs_node_info(uint32_t n, khash_t(bfs) *h){
	khiter_t k;
	bfs_node_info *info = calloc(n, sizeof(bfs_node_info));
	// printf("\n#");
	for (n = 0, k = kh_begin(h); k != kh_end(h); k++){
		if (kh_exist(h, k) && kh_val(h, k).l & BFLAG_N){
			info[n++].start = kh_key(h, k);
			// printf("kh_key:%u\n", kh_key(h, k));
		}
	}

	qsort(info, n, sizeof(bfs_node_info), bfs_node_info_cmp);
	return info;
}

uint32_t find_v_index_without_sort(bfs *bfs_v, UINTL_T v){
	uint32_t i, j = 0;
	for (i = 0; i < bfs_v->info_len && bfs_v->info[i].start; i ++){
		if (bfs_v->info[i].start == v) return i;
	}
	if (i == bfs_v->info_len) {
		j = bfs_v->info_len;
		bfs_v->info_len += 30;
		bfs_v->info = realloc(bfs_v->info, sizeof(bfs_node_info) * bfs_v->info_len);
		memset(bfs_v->info + j, 0, sizeof(bfs_node_info) * (bfs_v->info_len - j));
	}else j = i;

	bfs_v->info[j].start = v;
	return j;
}

typedef struct bfs_c_opt{
	float matches_ratio;
	uint32_t total_count;
	graph *g;
}bfs_c_opt;

void *bfs_nodes_compound_callback1(UINTL_T v, UINTL_T w, khash_t(bfs) *h, void *opt){
	bfs_c_opt *c_opt = (bfs_c_opt *)opt;
	bfs *bfs_v = &kh_val(h, kh_get(bfs, h, v));
	bfs *bfs_w = &kh_val(h, kh_get(bfs, h, w));
	
	uint32_t i;
	int64_t curr_matches, max_matches, temp;
	graph *g = c_opt->g;
	for (curr_matches = max_matches = i = 0; i < g->nd[v].odm; i++){
		if (g->eg[g->nd[v].oe[i]].l & MFLAG_FIT) continue;
		temp = (int64_t) g->eg[g->nd[v].oe[i]].sco * g->eg[g->nd[v].oe[i]].ide;
		if (temp > max_matches) max_matches = temp;
		if (g->eg[g->nd[v].oe[i]].ou == w) curr_matches = temp;
	}
	// printf("%u %u\n",g->nd[v].name, g->nd[w].name);
	// printf("curr_matches:%ld %ld %f\n", curr_matches, max_matches, c_opt->matches_ratio);
	curr_matches = curr_matches - max_matches * c_opt->matches_ratio;
	if (!bfs_w->info) {
		bfs_w->info_len = c_opt->total_count;
		bfs_w->info = init_bfs_node_info(bfs_w->info_len, h);
	}

	if (bfs_v->l & BFLAG_N){
		i = find_v_index_without_sort(bfs_w, v);
		bfs_w->info[i].predecessor = v;
		bfs_w->info[i].prop_diff_sum = curr_matches;
		bfs_w->info[i].edge_num = 1;
	}

	if (bfs_v->info){
		for (i = 0; i < bfs_v->info_len; i ++){
			if (bfs_v->info[i].start == v || !bfs_v->info[i].edge_num) continue;
			// printf("pmatch: %zd %zd %zd ",curr_matches,bfs_v->info[i].prop_diff_sum,curr_matches + bfs_v->info[i].prop_diff_sum );
			temp = curr_matches + bfs_v->info[i].prop_diff_sum;
		// printf("matches: v:%d w:%d start:%d %lld %lld\n", g->nd[v].name, g->nd[w].name,g->nd[bfs_w->info[i].start].name, temp,bfs_w->info[i].prop_diff_sum);
			if (temp > bfs_w->info[i].prop_diff_sum || !bfs_w->info[i].edge_num){
				// bfs_w->info[i].start = bfs_v->info[i].start;
				assert (bfs_w->info[i].start == bfs_v->info[i].start);
				bfs_w->info[i].predecessor = v;
				bfs_w->info[i].prop_diff_sum = temp;
				bfs_w->info[i].edge_num = bfs_v->info[i].edge_num + 1;
			}
		}
	}
	return NULL;
}



// ph *get_pending_nodes(khash_t(bfs) *h){
// 	ph *pd = malloc(sizeof(ph));
// 	init_sn(pd, 100);
	
// 	uint16_t l = 0;
// 	l |= (BFLAG_P1 | BFLAG_P2 | BFLAG_U);
	
// 	khiter_t k;
// 	for (k = kh_begin(h); k != kh_end(h); k++){
// 		if (kh_exist(h, k)){
// 			if (kh_val(h, k).l & l) {
// 				pd->n[pd->i++] = kh_key(h, k);
// 				if (pd->i >= pd->im) reallocate_sn(pd, 100);
// 			}
// 		}
// 	}
// 	return pd;
// }

bfs_r *bfs_nodes_compound(graph *g, ph *n, ph *exclude, const int depth,
		const int max_child_num, const int stop_at_merged, const float perc){
	
	bfs_c_opt c_opt = {
		.matches_ratio = perc,
		.total_count = n->i,
		.g = g
	};
	bfs_r *bfs_ret = bfs_nodes_compound_orig(g, n, exclude, depth, max_child_num, 1,
		stop_at_merged, bfs_nodes_compound_callback1, &c_opt);





	// ph *pending_nodes = get_pending_nodes(bfs_ret->bfs_info);//TODO test duplication
	// int i;
	// printf("\npending_nodes: ");
	// for (i=0;i<pending_nodes->i;i++){printf("%u ", g->nd[pending_nodes->n[i]].name);}
	// printf("\n");




	bfs_ret->data2 = bfs_nodes_compound_get_path(bfs_ret->bfs_info, depth);
	return bfs_ret;
}

void destroy_bfs_r(bfs_r *b, const int lable){
	khiter_t k;
	khash_t(bfs) *h = b->bfs_info;
	if (lable == 1){
		ph *visited_edges = (ph *) b->data1;
		phs *ps = (phs *) b->data2;
		destroy_sn(visited_edges); free (visited_edges);
		destroy_phs(ps); free(ps);
     	for (k = kh_begin(h); k != kh_end(h); k++){
        	if (kh_exist(h, k) && kh_val(h, k).info) free(kh_val(h, k).info);
     	}
		kh_destroy(bfs, b->bfs_info);
	}
	free (b);
}


void mark_path_with_lable(graph *g, ph_ *path, const uint16_t l){//TODO save edge in path
	UINTL_T i, j, d, v, w;
	// printf("lable edge:");
	for (i = 0; i < path->i - 1; i ++){
		v = path->n[i];
		w = path->n[i + 1];
		// printf("%d ",g->nd[v].name );
		for (j = 0; j < g->nd[v].odm; j++){
			d = g->nd[v].oe[j];
			if (g->eg[d].ou == w) {
				// printf(" (%u %u %u)", g->nd[g->eg[d].in].name, g->nd[g->eg[d].ou].name, d);
				g->eg[d].l |= MFLAG_TT;
				break;
			}
		}
	}
	// printf("\n");
	// printf("%d\n", g->nd[w].name);
}

static int check_paths_exist(graph *g, ph_ *path){
	int i;
	for (i = 0; i < path->i; i ++){
		if (g->eg[path->n[i]].l & MFLAG_FIT) return 0;
	}
	return i;
}

void rm_visited_edges(graph *g, phs *ps, ph *visited_edges, const int m, const int skip_rev){
	UINTL_T j, k, w, v, d, rd;
	ph_ *path;
	for (j = 0; j < ps->i; j ++) {
		path = &ps->ph[j];
		
		// printf("\npath%d %u %u:", j,g->nd[path->n[0]].name,g->nd[path->n[path->i - 1]].name);
		// for (k=0;k<path->i;k++){printf(" %u", g->nd[path->n[k]].name);}
		
		if (skip_rev && path->n[0] == get_reversed_node(g, path->n[path->i - 1])) continue;
		mark_path_with_lable(g, path, MFLAG_TT);
	}

	if (!m){
		// printf("\ndel:");
		for (j = 0; j < visited_edges->i; j ++){
			d = visited_edges->n[j];
			rd = get_reversed_edge(g, d);
			if (g->eg[d].l & MFLAG_TT || g->eg[rd].l & MFLAG_TT) {
				continue;
			}else{
				// printf("(%u %u) ", g->nd[g->eg[d].in].name, g->nd[g->eg[d].ou].name);
				rm_edge(g, d);
				rm_edge(g, rd);
			}
		}
		// printf("\n");
	}else{
		UINTL_T n;
		for (j = 0; j < ps->i; j ++){
			path = &ps->ph[j];
			for (k = 0; k < path->i - 1; k ++){
				v = path->n[k];
				w = path->n[k + 1];
				for (n = 0; n < g->nd[w].idm; n++){
					if (g->eg[g->nd[w].ie[n]].l & MFLAG_FIT) continue;
					d = g->nd[w].ie[n];
					rd = get_reversed_edge(g, d);
					if (g->eg[d].l & MFLAG_TT || g->eg[rd].l & MFLAG_TT) {
						continue;
					}else if (check_in_exclude(visited_edges, d)){
						rm_edge(g, d);
						rm_edge(g, rd);
					}
				}

				for (n = 0; n < g->nd[v].odm; n++){
					if (g->eg[g->nd[v].oe[n]].l & MFLAG_FIT) continue;
					d = g->nd[v].oe[n];
					rd = get_reversed_edge(g, d);
					if (g->eg[d].l & MFLAG_TT || g->eg[rd].l & MFLAG_TT) {
						continue;
					}else if(check_in_exclude(visited_edges, d)){
						rm_edge(g, d);
						rm_edge(g, rd);
					}
				}
			}
		}
	}

	for (j = 0; j < visited_edges->i; j ++){
		d = visited_edges->n[j];
		if (g->eg[d].l & MFLAG_TT) g->eg[d].l ^= MFLAG_TT;
	}
}

void clean_complex_single_path(graph *g, const int s, const int m, const float perc){
									//s:max_node_num, m:remove_method, perc:matches_ratio
	UINTL_T i;

	phs *ps;
	ph p;
	init_sn(&p, 2);

	bfs_r *bfs_ret;
	ph *visited_edges;
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].od < 2) continue;
		p.i = 0;
		p.n[p.i++] = i;
		bfs_ret = bfs_nodes_compound(g, &p, NULL, s, -1, 1, perc);
		ps = (phs *) bfs_ret->data2;
		visited_edges = (ph *) bfs_ret->data1;
		
		// printf("%d %d\n",g->nd[i].name, ps->i);
		// for (j = 0; j < ps->i; j ++){
		// 	path = &ps->ph[j];
		// 	for (k = 0; k < path->i - 1; k ++){
		// 		v = path->n[k];
		// 		w = path->n[k + 1];
		// 		printf("(%d %d) ",g->nd[v].name ,g->nd[w].name);
		// 	}
		// 	printf("\n");
		// }
		
		if (ps->i != 1) {
			destroy_bfs_r(bfs_ret, 1);
			continue;
		}
		// printf("%d\n",g->nd[i].name);
		rm_visited_edges(g, ps, visited_edges, m, 0);
		destroy_bfs_r(bfs_ret, 1);
	}
	destroy_sn(&p);
	// exit(1);
}

void *bfs_nodes_compound_common_callback3(UINTL_T v, UINTL_T w, khash_t(bfs) *h, void *opt){
	bfs_c_opt *c_opt = (bfs_c_opt *)opt;
	bfs *bfs_v = &kh_val(h, kh_get(bfs, h, v));
	bfs *bfs_w = &kh_val(h, kh_get(bfs, h, w));
	uint32_t i, j;
	int64_t curr_matches, max_matches, temp;

	graph *g = c_opt->g;
	for (curr_matches = max_matches = i = 0; i < g->nd[v].odm; i++){
		if (g->eg[g->nd[v].oe[i]].l & MFLAG_FIT) continue;
		temp = (int64_t) g->eg[g->nd[v].oe[i]].sco * g->eg[g->nd[v].oe[i]].ide;
		if (temp > max_matches) max_matches = temp;
		if (g->eg[g->nd[v].oe[i]].ou == w) curr_matches = temp;
	}
	// printf("%ld %ld %ld\n",curr_matches,max_matches,curr_matches - max_matches * c_opt->matches_ratio );
	curr_matches = curr_matches - max_matches * c_opt->matches_ratio;
	// printf("curr_matches:%ld\n", curr_matches);
	if (!bfs_w->info) {
		bfs_w->info_len = c_opt->total_count;
		bfs_w->info = init_bfs_node_info(bfs_w->info_len, h);
	}

	if (g->nd[v].od > 1){
		i = find_v_index_without_sort(bfs_w, v);
		bfs_w->info[i].predecessor = v;
		bfs_w->info[i].prop_diff_sum = curr_matches;
		bfs_w->info[i].edge_num = 1;
	}

	if (bfs_v->info){
		for (i = 0; i < bfs_v->info_len; i ++){
			if (bfs_v->info[i].start == v || !bfs_v->info[i].edge_num) continue;
			temp = curr_matches + bfs_v->info[i].prop_diff_sum;
			j = i < bfs_w->info_len && bfs_w->info[i].start == bfs_v->info[i].start ? \
				i : find_v_index_without_sort(bfs_w, bfs_v->info[i].start);
			if (temp > bfs_w->info[j].prop_diff_sum || !bfs_w->info[j].edge_num){
				bfs_w->info[j].predecessor = v;
				bfs_w->info[j].prop_diff_sum = temp;
				bfs_w->info[j].edge_num = bfs_v->info[i].edge_num + 1;
			}
	// printf("#(%u %u s:%u curr_matches:%ld, v_match, %ld)\n",g->nd[v].name, 
		// g->nd[w].name,g->nd[bfs_w->info[j].start].name, temp,bfs_v->info[i].prop_diff_sum  );
		}
	}
	return NULL;
}

// void bfs_nodes_compound_common_new_exclude(graph *g, UINTL_T n, ph *pending, ph *vis_nodes){
// 	ph *exclude = malloc(sizeof(ph));
// 	init_sn(exclude, 100);

// 	UINTL_T i, j, w;
// 	node *v = &g->nd[n];
// 	for (i = 0; i < v->idm; i++){
// 		if (g->eg[v->ie[i]].l & MFLAG_FIT) continue;
// 		w = g->eg[v->ie[i]].in;
// 		if (!check_in_exclude(vis_nodes, w)){
// 			exclude->n[exclude->i++] = w;
// 			if (exclude->i >= exclude->im) reallocate_sn(exclude, 100);
// 		}
// 	}

// 	for (i = 0; i < pending->i; i++){
// 		v = &g->nd[pending->n[i]];
// 		for (i = 0; i < v->idm; i++){
// 			if (g->eg[v->ie[i]].l & MFLAG_FIT) continue;
// 			w = g->eg[v->ie[i]].in;
// 			if (!check_in_exclude(vis_nodes, w)){
// 				exclude->n[exclude->i++] = w;
// 				if (exclude->i >= exclude->im) reallocate_sn(exclude, 100);
// 			}
// 		}
// 		for (i = 0; i < v->odm; i++){
// 			if (g->eg[v->oe[i]].l & MFLAG_FIT) continue;
// 			w = g->eg[v->oe[i]].ou;
// 			if (!check_in_exclude(vis_nodes, w)){
// 				exclude->n[exclude->i++] = w;
// 				if (exclude->i >= exclude->im) reallocate_sn(exclude, 100);
// 			}
// 		}
// 	}
// 	return exclude;
// }

// int bfs_nodes_compound_common_pred_nodes_callback(UINTL_T v, UINTL_T w, khash_t(map_) *r, void *opt){
// 	ph *p = (ph *) opt;
// 	if (!check_in_exclude_not_sort(p, w)){
// 		p->n[p->i++] = w;
// 		if (p->i >= p->im) reallocate_sn(p, 100);
// 	}
// 	return 0;
// }

// ph *bfs_nodes_compound_common_pred_nodes(graph *g, UINTL_T n, UINTL_T v,  ph *exclude){
// 	ph p, *prev_nodes_v;
// 	init_sn(&p, 1);
// 	init_sn(prev_nodes_v, 100);
// 	p.n[p.i++] = v;
// 	exclude->n[exclude->i++] = n;
// 	bfs_nodes(g, &p, exclude, -1, -1, 2, bfs_nodes_compound_common_pred_nodes_callback, prev_nodes_v);
// 	return prev_nodes_v;
// }

void bfs_nodes_compound_common_get_path(UINTL_T n, UINTL_T w, khash_t(bfs) *h, ph_ *p){
	UINTL_T v;
	bfs *bfs_v;
	p->n[p->i++] = w;
	bfs_v = &kh_val(h, kh_get(bfs, h, w));
	v = find_v_index_without_sort(bfs_v, n);
	v = bfs_v->info[v].predecessor;

	p->n[p->i++] = v;
	while (v != n) {
		bfs_v = &kh_val(h, kh_get(bfs, h, v));
		v = find_v_index_without_sort(bfs_v, n);
		v = bfs_v->info[v].predecessor;
		p->n[p->i++] = v;
		if (p->i >= p->im) reallocate_ph(p, 100);
	}
	reverse_ph2(p);
}

ph *get_pending_nodes(khash_t(bfs) *h){
	ph *pd = malloc(sizeof(ph));
	init_sn(pd, 100);
	
	uint16_t l = 0;
	l |= (BFLAG_P1 | BFLAG_P2 | BFLAG_U);
	
	khiter_t k;
	for (k = kh_begin(h); k != kh_end(h); k++){
		if (kh_exist(h, k)){
			// printf("key: %u\n", kh_key(h, k));
			if (kh_val(h, k).l & l) {
				pd->n[pd->i++] = kh_key(h, k);
				if (pd->i >= pd->im) reallocate_sn(pd, 100);
			}
		}
	}
	return pd;
}

ph *get_visited_nodes(khash_t(bfs) *h){
	ph *n = malloc(sizeof(ph));
	init_sn(n, 100);
	khiter_t k;
	for (k = kh_begin(h); k != kh_end(h); k++){
		if (kh_exist(h, k)){
			n->n[n->i++] = kh_key(h, k);
			if (n->i >= n->im) reallocate_sn(n, 100);
		}
	}
	return n;
}

bfs_r *bfs_nodes_compound_common(graph *g, ph *n, const int depth, const int max_child_num, const float perc){

	bfs_c_opt c_opt = {
		.matches_ratio = perc,
		.total_count = n->i,
		.g = g
	};

	bfs_r *bfs_ret = bfs_nodes_compound_orig(g, n, NULL, depth, max_child_num, 1, 
		1, bfs_nodes_compound_common_callback3, &c_opt);
	ph *pending_nodes = get_pending_nodes(bfs_ret->bfs_info);

	UINTL_T i, j, v, w;
	ph tmp_exclude1, tmp_exclude2;
	init_sn(&tmp_exclude1, 100);
	init_sn(&tmp_exclude2, 100);

	// for (i =0;i<pending_nodes->i; i ++){
	// 	printf("%u %u\n",n->n[0], pending_nodes->n[i]);
	// }

	ph *common_nodes = &tmp_exclude1;
	ph *old_common_nodes = &tmp_exclude2;
	common_nodes->i = old_common_nodes->i = 0;
	
	bfs *bfs_v = &kh_val(bfs_ret->bfs_info, kh_get(bfs, bfs_ret->bfs_info, pending_nodes->n[0]));
	for (i = 0; i < bfs_v->info_len && bfs_v->info[i].start; i++){
		if (!bfs_v->info[i].edge_num) continue;
        common_nodes->n[common_nodes->i++] = bfs_v->info[i].start;
        if (common_nodes->i >= common_nodes->im) reallocate_sn(common_nodes, 100);
    }
	
	// printf(" common_nodes1: ");
	// for (i = 0; i < common_nodes->i; i ++){
	// 	printf(" %d", g->nd[common_nodes->n[i]].name);
	// }

	for (i = 1; i < pending_nodes->i; i ++){
		bfs_v = &kh_val(bfs_ret->bfs_info, kh_get(bfs, bfs_ret->bfs_info, pending_nodes->n[i]));
		SWAP(common_nodes, old_common_nodes, ph*);
		for (common_nodes->i = j = 0; j < bfs_v->info_len && bfs_v->info[j].start; j ++){
			if (!bfs_v->info[j].edge_num) continue;
			w = bfs_v->info[j].start;

		// printf(" \n%d \n", g->nd[w].name);
		// for (v = 0; v < old_common_nodes->i; v++){
		// 	printf(" %d", g->nd[old_common_nodes->n[v]].name);
		// }
		// printf(" %d(%d)\n",check_in_exclude(old_common_nodes, w),check_in_exclude_not_sort(old_common_nodes, w) );
			if (check_in_exclude_not_sort(old_common_nodes, w)){
				common_nodes->n[common_nodes->i++] = w;
				if (common_nodes->i >= common_nodes->im) reallocate_sn(common_nodes, 100);
			}
		}
	}

	// printf(" common_nodes2: ");
	// for (i = 0; i < common_nodes->i; i ++){
	// 	printf(" %d", g->nd[common_nodes->n[i]].name);
	// }

	UINTL_T common_node = n->n[0];
	if (common_nodes->i > 1){
		int v_depth, max_depth = INT_MIN;
		for (i = j = 0; i < common_nodes->i; i ++){
			v = common_nodes->n[i];
			if (v == n->n[0]) continue;
			v_depth = kh_val(bfs_ret->bfs_info, kh_get(bfs, bfs_ret->bfs_info, v)).bfs_depth;
			if (v_depth > max_depth){
				max_depth = v_depth;
				common_node = v;
				j = 1;
			}else if (v_depth == max_depth) j ++;
			// printf("\nlog1:%u %d %d %u %d\n", g->nd[v].name, v_depth, max_depth, g->nd[common_node].name, j);
		}

		if (j > 1){
			int64_t prop_val, max_prop_val = INT32_MIN;
			for (i = 0; i < common_nodes->i; i ++){
				v = common_nodes->n[i];
				bfs_v = &kh_val(bfs_ret->bfs_info, kh_get(bfs, bfs_ret->bfs_info, v));
				if (bfs_v->bfs_depth == max_depth){
					j = find_v_index_without_sort(bfs_v, n->n[0]);
					prop_val = bfs_v->info[j].prop_diff_sum;
					// printf("\n%u s:%u prop_val: %ld\n",g->nd[v].name, g->nd[n->n[0]].name,  prop_val);
					assert(bfs_v->info[j].edge_num);
					for (j = 0; j < pending_nodes->i; j ++){
						bfs_v = &kh_val(bfs_ret->bfs_info, kh_get(bfs, bfs_ret->bfs_info, pending_nodes->n[j]));
						w = find_v_index_without_sort(bfs_v, v);//TODO optimize with sort
						prop_val += bfs_v->info[w].prop_diff_sum;
					}
					if (prop_val > max_prop_val){
						max_prop_val = prop_val;
						common_node = v;
					}
					// printf("log2:%u %d %ld %u\n", g->nd[v].name, bfs_v->bfs_depth, prop_val, g->nd[common_node].name);
				}
			}
		}
	}else assert(common_nodes->n[0] == n->n[0]);
	

	phs *ps = malloc(sizeof(phs));
	init_phs(ps, pending_nodes->i, 31);
	ph_ *p;
	// printf(" final common_node: %d\n", g->nd[common_node].name);
	if (common_node == n->n[0]){
		for (i = 0; i < pending_nodes->i; i ++){
			p = &ps->ph[ps->i++];
			bfs_nodes_compound_common_get_path(n->n[0], pending_nodes->n[i], bfs_ret->bfs_info, p);
		}
	}else{
		ph_ path_nc;
		init_ph(&path_nc, 31);
		bfs_nodes_compound_common_get_path(n->n[0], common_node, bfs_ret->bfs_info, &path_nc);

// printf("path_nc:\n"); for (j=0;j<path_nc.i;j++){printf(" %u",g->nd[path_nc.n[j]].name);}
		for (i = 0; i < pending_nodes->i; i ++){
			p = &ps->ph[ps->i++];
			bfs_nodes_compound_common_get_path(common_node, pending_nodes->n[i], bfs_ret->bfs_info, p);

// printf("\n%u: ", pending_nodes->n[i]); for (j=0;j<p->i;j++){printf(" %u",g->nd[p->n[j]].name);}

			if (p->i + path_nc.i - 1 >= p->im) reallocate_ph(p, p->i + path_nc.i - 1);
			memmove(p->n + path_nc.i - 1, p->n, sizeof(UINTL_T) * p->i);
			memcpy(p->n, path_nc.n, sizeof(UINTL_T) * (path_nc.i - 1));
			p->i += path_nc.i - 1;
		}
		free (path_nc.n);
	}

	bfs_ret->data2 = ps;
	destroy_sn(&tmp_exclude1);
	destroy_sn(&tmp_exclude2);
	destroy_sn(pending_nodes);
	free (pending_nodes);
	return bfs_ret;
}

// bfs_ret *bfs_nodes_compound_common(graph *g, ph *n, const int depth, const int max_child_num, const float perc){

// 	bfs_c_opt c_opt = {
// 		.matches_ratio = perc,
// 		.total_count = n->i,
// 		.g = g
// 	};

// 	bfs_r *bfs_ret = bfs_nodes_compound_orig(g, n, depth, NULL, max_child_num, 1, 
// 		1, bfs_nodes_compound_common_callback3, &c_opt);
// 	ph *vis_nodes = get_visited_nodes(bfs_ret->bfs_info);
// 	ph *pending_nodes = get_pending_nodes(bfs_ret->bfs_info);//TODO test duplication
// 	ph *new_exclude = bfs_nodes_compound_common_new_exclude(n->n[0], pending_nodes, vis_nodes);

// 	UINTL_T i, j, v, w;
// 	ph tmp_exclude1;
// 	init_sn(&tmp_exclude1, new_exclude->i + 100);

// 	ph_s *pred_nodes = malloc(sizeof(ph_s));
// 	pred_nodes->i = pred_nodes->im = pending_nodes->i;
// 	pred_nodes->ph = malloc(pred_nodes->im * sizeof(ph *));
// 	for (i = 0; i < pending_nodes->i; i ++){
// 		v = pending_nodes->n[i];
// 		memcpy(&tmp_exclude1, new_exclude, sizeof(UINTL_T) * new_exclude->i);
// 		tmp_exclude1.i = new_exclude->i;
// 		for (j = 0; j < pending_nodes->i; j ++){
// 			w = pending_nodes->n[j];
// 			if (v != w){
// 				tmp_exclude1.n[tmp_exclude1.i++] = w;
// 				if (tmp_exclude1.i >= tmp_exclude1.im) reallocate_sn(&tmp_exclude1, 100);
// 			}
// 		}
// 		pred_nodes->ph[i] = bfs_nodes_compound_common_pred_nodes(g, n->n[0], v, tmp_exclude);
// 	}

// 	ph *common_nodes = &tmp_exclude1;
// 	ph *old_common_nodes = vis_nodes;
// 	common_nodes->i = old_common_nodes->i = 0;
	
// 	for (j = 0; j < pred_nodes->ph[0].i; j ++){
// 		common_nodes->n[common_nodes->i++] = pred_nodes->ph[0].n[j];
// 		if (common_nodes->i >= common_nodes->im) reallocate_sn(common_nodes, 100);
// 	}
// 	for (i = 1; i < pending_nodes->i && common_nodes->i; i ++){
// 		v = pending_nodes->n[i];
// 		SWAP(common_nodes, old_common_nodes, ph*);
// 		for (common_nodes->i = j = 0; j < pred_nodes->ph[i].i; j ++){
// 			w = pred_nodes->ph[i].n[j];
// 			if (!check_in_exclude(old_common_nodes, w)){
// 				common_nodes->n[common_nodes->i++] = w;
// 				if (common_nodes->i >= common_nodes->im) reallocate_sn(common_nodes, 100);
// 			}
// 		}
// 	}

// 	UINTL_T common_node = n->n[0];
// 	if (common_nodes->i){
// 		int v_depth, max_depth = -1;
// 		for (i = 0; i < common_nodes->i; i ++){
// 			v = common_nodes->n[i];
// 			v_depth != common_node ? kh_val(bfs_ret->bfs_info, kh_get(bfs, bfs_ret->bfs_info, v)).bfs_depth : 0;
// 			if (v_depth > max_depth){
// 				v_depth = max_depth;
// 				common_node = v;
// 			}
// 		}
// 	}
	
// 	ph *pred_node;
// 	if (common_node == n->n[0]){
// 		for (i = 0; i < pending_nodes->i; i ++){
// 			pred_node = &pred_nodes->ph[i];
// 			pred_node->i = 0;
// 			bfs_nodes_compound_common_get_path(n->n[0], pending_nodes->n[i], bfs_ret->bfs_info, pred_node);
// 		}
// 	}else{
// 		ph *path_nc = &tmp_exclude1;
// 		path_nc->i = 0;
// 		bfs_nodes_compound_common_find_path(g, n->n[0], common_node, new_exclude, path_nc);
// 		for (i = 0; i < pending_nodes->i; i ++){
// 			pred_node = &pred_nodes->ph[i];
// 			pred_node->i = path_nc->i - 1;
// 			if (pred_node->i < pred_node->im) reallocate_sn(pred_node, pred_node->i + 100);
// 			memcpy(path_nc->n, pred_node->n, sizeof(UINTL_T) * pred_node->i);
// 			bfs_nodes_compound_common_find_path(n->n[0], pending_nodes->n[i], bfs_ret->bfs_info, pred_node);
// 		}
// 	}

// 	bfs_ret->data2 = pred_nodes;
// 	destroy_sn(&tmp_exclude1);
// 	destroy_sn(vis_nodes);
// 	free (vis_nodes);
// 	return bfs_ret;
// }

void clean_complex_multi_path(graph *g, const int s, const int m, const float perc){
	UINTL_T i;

	ph p;
	init_sn(&p, 2);

	phs *ps;
	bfs_r *bfs_ret;
	ph *visited_edges;
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].od < 2) continue;
		p.i = 0;
		p.n[p.i++] = i;
		bfs_ret = bfs_nodes_compound_common(g, &p, s, 30, perc);
		ps = (phs *) bfs_ret->data2;
		visited_edges = bfs_ret->data1;
		
		// printf("save path:%d %d\n",g->nd[i].name, ps->i);
		// UINTL_T j, v, k;
		// for (j = 0; j < ps->i; j ++){
		// 	ph_ *path = &ps->ph[j];
		// 	for (k = 0; k < path->i; k ++){
		// 		v = path->n[k];
		// 		printf("%d ",g->nd[v].name);
		// 	}
		// 	printf("\n");
		// }

		rm_visited_edges(g, ps, visited_edges, m, 0);
		destroy_bfs_r(bfs_ret, 1);
	}
	destroy_sn(&p);
}


uint64_t get_max_score(graph *g, UINTL_T n, int ou){

	UINTL_T i;
	uint64_t score, temp;
	if (ou){
		for (score = i = 0; i < g->nd[n].odm; i++){
			if (g->eg[g->nd[n].oe[i]].l & MFLAG_FIT) continue;
			temp = (uint64_t) g->eg[g->nd[n].oe[i]].sco * g->eg[g->nd[n].oe[i]].ide;
			if (temp > score) score = temp;
		}
	}else{
		for (score = i = 0; i < g->nd[n].idm; i++){
			if (g->eg[g->nd[n].ie[i]].l & MFLAG_FIT) continue;
			temp = (uint64_t) g->eg[g->nd[n].ie[i]].sco * g->eg[g->nd[n].ie[i]].ide;
			if (temp > score) score = temp;
		}
	}
	return score;//should score/10000
}

uint32_t get_max_tc(graph *g, UINTL_T n, int ou){

	UINTL_T i;
	uint32_t tc;
	if (ou){
		for (tc = i = 0; i < g->nd[n].odm; i++){
			if (g->eg[g->nd[n].oe[i]].l & MFLAG_FIT) continue;
			if (g->eg[g->nd[n].oe[i]].tc > tc) tc = g->eg[g->nd[n].oe[i]].tc ;
		}
	}else{
		for (tc = i = 0; i < g->nd[n].idm; i++){
			if (g->eg[g->nd[n].ie[i]].l & MFLAG_FIT) continue;
			if (g->eg[g->nd[n].ie[i]].tc > tc) tc = g->eg[g->nd[n].ie[i]].tc;
		}
	}
	return tc;
}


void rm_z_clip_score3(graph *g, const int s, const int m, const int perc){

	UINT_T t;
	UINTL_T i, j, n, tc, max_tc_ou, max_tc_in;
	int64_t max_sco_ou, max_sco_in;
	phs ps;
	ph_ *p = NULL;
	init_phs(&ps, 100, s + 1);
	for (i = 1; i < g->ni; i++){
		ps.i = 0;
		for (j = 0; g->nd[i].od > 1 && j < g->nd[i].odm; j++){
			if (g->eg[g->nd[i].oe[j]].l & MFLAG_FIT) continue;
			p = &ps.ph[ps.i];
			p->i = p->sco = 0;
			p->n[p->i++] = g->nd[i].oe[j];
			n = g->eg[g->nd[i].oe[j]].ou;
			
			while (g->nd[n].id == 1 && g->nd[n].od == 1 && p->i < s){
				t = get_validly_oe(g, n, 0);
				p->n[p->i++] = g->nd[n].oe[t];
				n = g->eg[g->nd[n].oe[t]].ou;
				// printf("#%u(%d)#\n", g->nd[n].name,p->i);
			}
			if (g->nd[n].id > 1) {
				ps.i++;
				if (ps.i >= ps.im) reallocate_phs(&ps, 100, s + 1);
			}
		}

		if (ps.i){
			
			// printf("%u\n", g->nd[i].name);
			// for (j = 0; j < ps.i; j ++){
			// 	p = &ps.ph[j];
			// 	printf("\n%u score:%ld tc:%ld",j, p->sco>>8, p->sco & 255);
			// 	// for (n =0; n < p->i; n++){
			// 	// 	printf(" [%u %u]", g->nd[g->eg[p->n[n]].in].name, g->nd[g->eg[p->n[n]].ou].name);
			// 	// }
			// 	// if (j == ps.i - 1) printf("\n");
			// }

			max_sco_ou = get_max_score(g, i, 1);
			max_tc_ou = get_max_tc(g, i, 1);
			for (j = 0; j < ps.i; j ++){
				p = &ps.ph[j];
				n = g->eg[p->n[p->i-1]].ou;
				max_sco_in = get_max_score(g, n, 0);
				p->sco = max_sco_ou ? (int64_t) g->eg[p->n[0]].sco * g->eg[p->n[0]].ide * 50 / max_sco_ou : 0;
				p->sco += max_sco_in ? (int64_t) g->eg[p->n[p->i-1]].sco * g->eg[p->n[p->i-1]].ide * 50 / max_sco_in : 0;

			// for (n =0; n < p->i; n++){printf(" [%u %u]", g->nd[g->eg[p->n[n]].in].name, g->nd[g->eg[p->n[n]].ou].name);}
			// printf("score:%ld ",p->sco);
				if (p->sco >= perc) p->sco = perc;

				max_tc_in = get_max_tc(g, n, 0);
				tc = max_tc_ou ? (int64_t) g->eg[p->n[0]].tc * 50 / max_tc_ou : 0;
				tc += max_tc_in ? (int64_t) g->eg[p->n[p->i-1]].tc * 50 / max_tc_in : 0;
			// printf("(%d %d %d %d %d)",max_tc_ou ? (int64_t) g->eg[p->n[0]].tc * 50 / max_tc_ou : 0,
				// max_tc_in ? (int64_t) g->eg[p->n[0]].tc * 50 / max_tc_in : 0, max_tc_ou, max_tc_in, g->eg[p->n[0]].tc );
			// printf("tc:%ld\n",tc);

				if (tc >= perc) tc = perc;
				p->sco = p->sco << 8 | tc;
			}
			msort_phs(&ps);
		}

		for (j = 0; j < ps.i; j ++){
			p = &ps.ph[j];
			if (g->nd[g->eg[p->n[0]].in].od < 2 || g->nd[g->eg[p->n[p->i-1]].ou].id < 2) continue;
			if ((p->sco >> 8) >= perc && (p->sco & 255) >= perc) break;

			// printf("\n%u %ld",j, p->sco);
			// for (n =0; n < p->i; n++){
			// 	printf(" [%u %u]", g->nd[g->eg[p->n[n]].in].name, g->nd[g->eg[p->n[n]].ou].name);
			// }
			// if (j == ps.i - 1) printf("\n");

			if (m && p->i >= m){
				rm_edge(g, p->n[p->i-1]);
				rm_edge(g, get_reversed_edge(g, p->n[p->i-1]));
				rm_edge(g, p->n[0]);
				rm_edge(g, get_reversed_edge(g, p->n[0]));
			}else{
				while (--p->i >= 0){
					rm_edge(g, p->n[p->i]);
					rm_edge(g, get_reversed_edge(g, p->n[p->i]));
				}
			}
		}
	}
	// }
	destroy_phs(&ps);
	// exit(1);
}

void cp_ph_(ph_ *dest, ph_ *src, int deep){
	if (!deep){
		SWAP(*dest, *src, ph_);
		if (!src->im){
			src->im = 1;
			src->n = malloc(sizeof(UINTL_T));
		}
	}else{
		dest->sco = src->sco;
		dest->i = src->i;
		if (dest->im){
			if (dest->i >= dest->im){
				dest->im = dest->i + 1;
				dest->n = realloc(dest->n, dest->im * sizeof(UINTL_T));
			}
		}else{
			dest->im = dest->i + 1;
			dest->n = malloc(dest->im * sizeof(UINTL_T));
		}
		memcpy(dest->n, src->n, dest->i * sizeof(UINTL_T));
	}
}

static void cal_z_path_score(graph *g, ph_ *p){
	int64_t max_sco_ou, max_sco_in;
	max_sco_ou = get_max_score(g, g->eg[p->n[0]].in, 1);
	p->sco = max_sco_ou ? (int64_t) g->eg[p->n[0]].sco * g->eg[p->n[0]].ide * 10000 / max_sco_ou : 0;
	max_sco_in = get_max_score(g, g->eg[p->n[p->i-1]].ou, 0);
	p->sco += max_sco_in ? (int64_t) g->eg[p->n[p->i-1]].sco * g->eg[p->n[p->i-1]].ide * 10000 / max_sco_in : 0;
}

void add_path(graph *g, phs *new_paths, phs *left_z_paths, ph_ *p, const int s, const int prop_thres){
	cal_z_path_score(g, p);
	// printf("p:new_paths: %p %p %p\n", new_paths, new_paths->ph[0].n,new_paths->ph[1].n);
	// printf("%ld\n", p->sco);
	if (p->sco >= prop_thres && prop_thres >= 0){
		cp_ph_(&left_z_paths->ph[left_z_paths->i ++], p, 0);
		if (left_z_paths->i >= left_z_paths->im) reallocate_phs(left_z_paths, 100, s + 1);
	}else if (!new_paths->i){
		// printf("new_paths: %d %p %p\n", new_paths->i, new_paths->ph[new_paths->i].n,new_paths);
		cp_ph_(&new_paths->ph[new_paths->i ++], p, 0);
		// printf("%d %p\n", new_paths->i, new_paths->ph[new_paths->i-1].n);
	}else{
		int64_t i = 0;
		for (i = 0; i < new_paths->i; i ++){
			if (p->sco < new_paths->ph[i].sco) SWAP(new_paths->ph[i], *p, ph_);
		}
		cp_ph_(&new_paths->ph[new_paths->i ++], p, 0);

		if (new_paths->i >= new_paths->im) reallocate_phs(new_paths, 100, s + 1);
	}
	// printf("f:new_paths: %p %p %p\n", new_paths, new_paths->ph[0].n,new_paths->ph[1].n);
}


void find_z_path_from(graph *g, phs *new_z_paths, phs *left_z_paths, UINTL_T node, const int s, 
		const int prop_thres){//TODO optimize

	UINTL_T n = node, t;
	ph_ cur_path1;
	init_ph(&cur_path1, s + 1);
	while (g->nd[n].id == 1 && g->nd[n].od == 1 && cur_path1.i < s){
		t = get_validly_oe(g, n, 0);
		cur_path1.n[cur_path1.i++] = g->nd[n].oe[t];
		// printf("(1:%u %u)\n",g->nd[g->eg[g->nd[n].oe[t]].in].name,g->nd[g->eg[g->nd[n].oe[t]].ou].name);
		n = g->eg[g->nd[n].oe[t]].ou;
	}
	if (cur_path1.i && g->nd[n].id > 1) {
		n = node;
		ph_ cur_path2;
		init_ph(&cur_path2, s + 2);
		while (g->nd[n].id == 1 && g->nd[n].od == 1 && cur_path1.i + cur_path2.i <= s + 1){
			t = get_validly_ie(g, n, 0);
			cur_path2.n[cur_path2.i++] = g->nd[n].ie[t];
			// printf("(2:%u %u)\n",g->nd[g->eg[g->nd[n].ie[t]].in].name,g->nd[g->eg[g->nd[n].ie[t]].ou].name);
			n = g->eg[g->nd[n].ie[t]].in;
		}
		if (g->nd[n].od > 1 && cur_path2.i + cur_path1.i <= s){
			// for (t=0;t<cur_path2.i;t++){printf("(3.0:%u %u)\n",g->eg[cur_path2.n[t]].in,g->eg[cur_path2.n[t]].ou);}
			reverse_ph2(&cur_path2);
			// for (t=0;t<cur_path2.i;t++){printf("(3.1:%u %u)\n",g->eg[cur_path2.n[t]].in,g->eg[cur_path2.n[t]].ou);}
			memcpy(cur_path2.n + cur_path2.i, cur_path1.n, cur_path1.i * sizeof(UINTL_T));
			cur_path2.i += cur_path1.i;
			// assert (cur_path2.i);
			// printf("%ld %u %u\n", cur_path2.sco, cur_path1.i,cur_path2.i);
			// printf("left_z_paths4:");
			// for (t=0;t<cur_path2.i;t++){printf("(%u %u)",g->nd[g->eg[cur_path2.n[t]].in].name,g->nd[g->eg[cur_path2.n[t]].ou].name);}
			// printf("\n");	
			add_path(g, new_z_paths, left_z_paths, &cur_path2, s, prop_thres);
			// printf("done\n");
		}
		free (cur_path2.n);
	}
	free (cur_path1.n);
	
	// int i;for (i = 0; i < left_z_paths->i; i ++){
	// 	ph_ *p = &left_z_paths->ph[i];
	// 	printf("left_z_paths1:%ld %u(%d) %u(%d) %ld %d\n",i,  g->nd[g->eg[p->n[0]].in].name, g->nd[g->eg[p->n[0]].in].od, g->nd[g->eg[p->n[p->i-1]].ou].name, g->nd[g->eg[p->n[p->i-1]].ou].id, p->sco, p->i);
	// }
	// printf("%p %p\n",new_z_paths,new_z_paths->ph );
}


void merge_z_paths(phs *p1, int64_t p1_s, phs *p2, int64_t p2_s, const int s){//TODO optimize
	int64_t i;
	if (p1_s){
		for (i = 0; i < p1->i - p1_s; i ++) SWAP(p1->ph[i], p1->ph[p1_s + i], ph_);	
		p1->i -= p1_s;
	}

	if (p2->i - p2_s > 0){
		int64_t j, im = p1->i + p2->i - p2_s;
		if (im >= p1->im) reallocate_phs(p1, im + 1 - p1->im, s + 1);
		// for (i = 0; i < p1->i - p1_s; i ++) SWAP(p1->ph[i], p1->ph[p1_s + i], ph_);	
		// printf("0: %d %d %d %d %d\n",p1->i,p1_s,p2_s,p2->i,im );
		// p1->i -= p1_s;
		i = p1->i - 1;
		j = p2->i - 1;
		im --;
		while (i >= 0 || j >= p2_s){
			// printf("%d %d\n",i,j );
			if (i >= 0 && (j < p2_s || p1->ph[i].sco >= p2->ph[j].sco)) {
				if (im != i) SWAP(p1->ph[im], p1->ph[i], ph_);
				i --;
				// printf("1: %ld %ld %ld  %ld\n", i+1,p1->ph[i+1].sco,im,p1->ph[0].n[0]);
			}else cp_ph_(&p1->ph[im], &p2->ph[j--], 0);
				// printf("2: %ld %ld %ld %ld\n", j+1,p2->ph[j+1].sco,im,p1->ph[0].n[0]);}
				// printf("2 %d %ld %d %ld\n",j+1,p2->ph[j].sco,im, p1->ph[0].n[0]);}
			im --;
		}
		// printf("i:%d j:%d im:%d\n",i,j,im );
		p1->i += p2->i - p2_s;
	// }else{
	// 	// assert (p1_s == 0);
	// 	printf("errorerrorerrorerror %ld %ld\n", p1->i, p1_s);
	}
}


void z_clipping_by_updating(graph *g, phs *ps, phs *ps_new, phs *left_z_paths, const int s, const int m, const int prop_thres){
	if (ps->i){
		clear_phs(ps_new);
		clear_phs(left_z_paths);

		int64_t curr_i = 0;
		int i_j_flag;
		ph_ path = {0};
		UINTL_T start, end;
		while (curr_i < ps->i || ps_new->i){
			i_j_flag = 1;
			if (curr_i >= ps->i){
				i_j_flag = 0;
				cp_ph_(&path, &ps_new->ph[0], 0);
			}else if (!ps_new->i){
				cp_ph_(&path, &ps->ph[curr_i], 0);
			}else{
				if (ps->ph[curr_i].sco <= ps_new->ph[0].sco){
					cp_ph_(&path, &ps->ph[curr_i], 0);
				}else{
					i_j_flag = 0;
					cp_ph_(&path, &ps_new->ph[0], 0);
				}
			}
			if (path.sco >= prop_thres && prop_thres >= 0) break;
			// if (path.sco > prop_thres && prop_thres >= 0) break;
			if (i_j_flag) curr_i ++;
			else{
				// printf("p1:new_paths: %p %p %p %p %d\n", ps_new,ps_new->ph, ps_new->ph[0].n,ps_new->ph[1].n, ps_new->im);
				ph_ _ = {.n = ps_new->ph[0].n, .sco = 0, .i = 0, .im = ps_new->ph[0].im };
				ps_new->i --;
				memmove(ps_new->ph, ps_new->ph + 1, ps_new->i * sizeof(ph_)); //TODO optimize
				ps_new->ph[ps_new->i] = _;	
				// printf("p2:new_paths: %p %p %p %p\n", ps_new, ps_new->ph, ps_new->ph[0].n,ps_new->ph[1].n);
			}

			start = g->eg[path.n[0]].in;
			end = g->eg[path.n[path.i-1]].ou;
			if (g->nd[start].od < 2 || g->nd[end].id < 2 || !check_paths_exist(g, &path)) continue;

			if (m && path.i >= m){
				rm_edge(g, path.n[path.i-1]);
				rm_edge(g, get_reversed_edge(g, path.n[path.i-1]));
				rm_edge(g, path.n[0]);
				rm_edge(g, get_reversed_edge(g, path.n[0]));
			}else{
				// printf("del:");
				while (--path.i >= 0){
					// printf(" %ld", g->nd[g->eg[path.n[path.i]].in].name);
					rm_edge(g, path.n[path.i]);
					rm_edge(g, get_reversed_edge(g, path.n[path.i]));
				}
			}
			// printf("%ld\n",ps_new.i);
			// printf("\nstart:%u end:%u sco:%ld index:%d prop_thres %d\n",g->nd[start].name, g->nd[end].name,path.sco,curr_i,prop_thres  );
			find_z_path_from(g, ps_new, left_z_paths, start, s, prop_thres);
			find_z_path_from(g, ps_new, left_z_paths, end, s, prop_thres);
			find_z_path_from(g, ps_new, left_z_paths, get_reversed_node(g, start), s, prop_thres);
			find_z_path_from(g, ps_new, left_z_paths, get_reversed_node(g, end), s, prop_thres);

		// int i;for (i = 0; i < left_z_paths.i; i ++){
		// ph_ *p = &left_z_paths.ph[i];
		// printf("left_z_paths3:%ld %u(%d) %u(%d) %ld %d\n",i,  g->nd[g->eg[p->n[0]].in].name, g->nd[g->eg[p->n[0]].in].od, g->nd[g->eg[p->n[p->i-1]].ou].name, g->nd[g->eg[p->n[p->i-1]].ou].id, p->sco, p->i);
		// }
		// for (i = 0; i < ps_new.i; i ++){
		// ph_ *p = &ps_new.ph[i];
		// printf("ps_new:%ld %u(%d) %u(%d) %ld %d\n",i,  g->nd[g->eg[p->n[0]].in].name, g->nd[g->eg[p->n[0]].in].od, g->nd[g->eg[p->n[p->i-1]].ou].name, g->nd[g->eg[p->n[p->i-1]].ou].id, p->sco, p->i);
		// }
		}
		
		free (path.n);
		assert (ps_new->i == 0);
		msort_phs(left_z_paths);
		// sort_phs_g(&left_z_paths, g);//TODO REPLACE sort_phs

	// int i;for (i = 0; i < left_z_paths->i; i ++){
	// 	ph_ *p = &left_z_paths->ph[i];
	// 	printf("left_z_paths2:%ld %u(%d) %u(%d) %ld %d\n",i,  g->nd[g->eg[p->n[0]].in].name, g->nd[g->eg[p->n[0]].in].od, g->nd[g->eg[p->n[p->i-1]].ou].name, g->nd[g->eg[p->n[p->i-1]].ou].id, p->sco, p->i);
	// }
	// printf("#########\n");
	// printf("%ld %ld\n",ps_new->i,  curr_i);
	// for (i = 0; i < ps->i; i ++){
	// 	ph_ *p = &ps->ph[i];
	// 	printf("merge1:%ld %u-%u %u %ld %d\n",i,  g->nd[g->eg[p->n[0]].in].name, g->nd[g->eg[p->n[0]].ou].name, g->nd[g->eg[p->n[p->i-1]].ou].name, p->sco, p->i);
	// }
		// printf("ps:%p ps->i:%d curr_i:%d  left_z_paths->i:%d\n",ps,ps->i,curr_i,left_z_paths->i);
		merge_z_paths(ps, curr_i, left_z_paths, 0, s);

	// sort_phs_g(ps, g);//TODO del
	// for (i = 0; i < ps->i; i ++){
	// 	ph_ *p = &ps->ph[i];
	// 	printf("merge2:%ld %u-%u %u %ld %d\n",i,  g->nd[g->eg[p->n[0]].in].name, g->nd[g->eg[p->n[0]].ou].name, g->nd[g->eg[p->n[p->i-1]].ou].name, p->sco, p->i);
	// }
	// exit(1);

		// printf("f:%p %p\n",ps_new,ps_new->ph );
		// return ps_new;
	}
}

void rm_z_clip_score2(graph *g, const int s, const int m, const int perc){

	UINT_T t;
	UINTL_T n;
	int64_t i, j, sco;
	phs ps;
	ps.i = 0;

	ph_ *p = NULL;
	init_phs(&ps, 1000, s + 1);
	for (i = 1; i < g->ni; i++){
		for (j = 0; g->nd[i].od > 1 && j < g->nd[i].odm; j++){
			if (g->eg[g->nd[i].oe[j]].l & MFLAG_FIT) continue;
			p = &ps.ph[ps.i];
			p->i = p->sco = 0;
			p->n[p->i++] = g->nd[i].oe[j];
			n = g->eg[g->nd[i].oe[j]].ou;
			
			while (g->nd[n].id == 1 && g->nd[n].od == 1 && p->i < s - 1){
				t = get_validly_oe(g, n, 0);
				p->n[p->i++] = g->nd[n].oe[t];
				n = g->eg[g->nd[n].oe[t]].ou;
				// printf("#%u(%d)#\n", g->nd[n].name,p->i);
			}
			if (g->nd[n].id > 1) {
				cal_z_path_score(g, p);
				ps.i++;
				if (ps.i >= ps.im) reallocate_phs(&ps, 1000, s + 1);
			}
		}
	}
	// msort_phs(&ps);
	qsort_phs(&ps);
	// sort_phs_g(&ps, g);

	// for (i = 0; i < ps.i; i ++){
	// 	p = &ps.ph[i];
	// 	printf("%ld %u %u %ld %d\n",i,  g->nd[g->eg[p->n[0]].in].name, g->nd[g->eg[p->n[p->i-1]].ou].name, p->sco, p->i);
	// }
	// exit(1);

	int each_num = 1000;
	int64_t prop_thres_lst_len = ps.i / each_num + 1;
	int64_t *prop_thres_lst = malloc(prop_thres_lst_len * sizeof(int64_t));
	for (i = 1, j = 0; i < ps.i / each_num; i ++){
		sco = ps.ph[i * each_num].sco;
		if (sco > perc && perc >= 0) break;
		prop_thres_lst[j++] = sco;
	}
	prop_thres_lst[j++] = perc;

	// for (i = 0; i < j; i ++){printf("%d %ld\n",i,prop_thres_lst[i]);}
	phs left_z_paths, ps_new;
	init_phs(&ps_new, 100, s + 1);
	init_phs(&left_z_paths, 100, s + 1);
	for (i = 0; i < j; i ++){
		// printf("1:%ld %d %ld\n",i, j, ps.ph[0].n[0]);
		z_clipping_by_updating(g, &ps, &ps_new, &left_z_paths, s, m, prop_thres_lst[i]);
	}
	free (prop_thres_lst);
	destroy_phs(&ps);
	destroy_phs(&ps_new);
	destroy_phs(&left_z_paths);
	// exit(1);
}

ph *get_dense_spot(graph *g, UINTL_T n, const int depth){

	int absent, dp;
	khash_t(set_) *r = kh_init(set_);
	khint_t k = kh_put(set_, r, n, &absent);
	kh_key(r, k) = n;
	
	UINTL_T j, v, w;

	ph *h, *h2, *h3;
	h = calloc(1, sizeof(ph));
	h2 = calloc(1, sizeof(ph));
	h3 = calloc(1, sizeof(ph));
	init_sn(h, 10);
	h->n[h->i++] = n;
	init_sn(h2, 100); 
	init_sn(h3, 100);

	while (h->i){
		h2->i = dp = 0;
		while (dp++ < depth && h->i){
			h3->i = 0;
			while (h->i){
				v = h->n[--h->i];
				for (j = 0; j < g->nd[v].odm; j++){
					if (g->eg[g->nd[v].oe[j]].l & MFLAG_FIT) continue;
					w = g->eg[g->nd[v].oe[j]].ou;
					k = kh_put(set_, r, w, &absent);
					if (!absent) continue;
					else kh_key(r, k) = w;
					if (g->nd[w].id > 1 || g->nd[w].od > 1 ){
						h2->n[h2->i++] = w;
						if (h2->i >= h2->im) reallocate_sn(h2, 1000);
						continue;
					}
					h3->n[h3->i++] = w;
					if (h3->i >= h3->im) reallocate_sn(h3, 1000);
				}
				for (j = 0; j < g->nd[v].idm; j++){
					if (g->eg[g->nd[v].ie[j]].l & MFLAG_FIT) continue;
					w = g->eg[g->nd[v].ie[j]].in;
					k = kh_put(set_, r, w, &absent);
					if (!absent) continue;
					else kh_key(r, k) = w;
					if (g->nd[w].id > 1 || g->nd[w].od > 1 ){
						h2->n[h2->i++] = w;
						if (h2->i >= h2->im) reallocate_sn(h2, 1000);
						continue;
					}
					h3->n[h3->i++] = w;
					if (h3->i >= h3->im) reallocate_sn(h3, 1000);
				}
			}
			SWAP(h, h3, ph*);
		}
		SWAP(h, h2, ph*);
	}
	destroy_sn(h2);
	free (h2);
	destroy_sn(h3);
	free (h3);

	h->i = 0;
	for (k = kh_begin(r) ; k < kh_end(r); k++){
		if (kh_exist(r, k)){
			h->n[h->i++] = v = kh_key(r, k);
			w = get_reversed_node(g, v);
			if (kh_get(set_, r, w) == kh_end(r)) h->n[h->i++] = w;
			if (h->i + 2 >= h->im) reallocate_sn(h, 1000);
		}
	}
	kh_destroy(set_, r);
	return h;
}

int get_all_spots(graph *g, ph **p, UINTL_T i, const int s){
	UINTL_T j;
	for (;i < g->ni; i++){
		if (g->nd[i].l || (g->nd[i].od <= 1 && g->nd[i].id <= 1)) continue;
		*p = get_dense_spot(g, i, s);
		for (j = 0; j < (*p)->i; j ++) g->nd[(*p)->n[j]].l = 1;
		return i;
	}
	return 0;
}

ph *get_start_end_nodes(graph *g, ph *n){
	UINTL_T i, j, v, w;
	ph *r = malloc(3 * sizeof(ph));
	for (i = 0; i < 3; i ++) init_sn(&r[i], 100);
	ph *s = &r[0], *e = &r[1], *x = &r[2];
	
	int l;
	
	// printf("\np:");
	// for (i = 0; i < n->i; i++) {printf("%u ", n->n[i]);}
	qsort(n->n, n->i, sizeof(UINTL_T), cmp_int);
	n->im = 0;

	// printf("\ns:");
	// for (i = 0; i < n->i; i++) {printf("%u ", n->n[i]);}
	
	for (i = 0; i < n->i; i++) {
		v = n->n[i];
		if (!g->nd[v].id) {
			s->n[s->i++] = v;
			if (s->i >= s->im) reallocate_sn(s, 100);
		}else{
			for (l = j = 0; j < g->nd[v].idm; j++){
				if (g->eg[g->nd[v].ie[j]].l & MFLAG_FIT) continue;
				w = g->eg[g->nd[v].ie[j]].in;
				// printf("(%u %u %d)",g->nd[v].name, g->nd[w].name, check_in_exclude(n, w));
				if (check_in_exclude(n, w)) continue;
				else{
					l = 1;
					x->n[x->i++] = w;
					if (x->i >= x->im) reallocate_sn(x, 100);
				}
			}
			if (l){
				s->n[s->i++] = v;
				if (s->i >= s->im) reallocate_sn(s, 100);
			}
		}

		if (!g->nd[v].od) {
			e->n[e->i++] = v;
			if (e->i >= e->im) reallocate_sn(e, 100);
		}else{
			for (l = j = 0; j < g->nd[v].odm; j++){
				if (g->eg[g->nd[v].oe[j]].l & MFLAG_FIT) continue;
				w = g->eg[g->nd[v].oe[j]].ou;
				if (check_in_exclude(n, w)) continue;
				else{
					l = 1;
					x->n[x->i++] = w;
					if (x->i >= x->im) reallocate_sn(x, 100);
				}
			}
			if (l){
				e->n[e->i++] = v;
				if (e->i >= e->im) reallocate_sn(e, 100);
			}
		}
	}
	return r;
}

void solve_spot_by_finding_paths_hete(graph *g, ph *n, const int m, const float perc){
	ph *r = get_start_end_nodes(g, n);
	ph *s = &r[0], *e = &r[1];

	UINTL_T i, j, v;
	for (i = j = 0; i < s->i && !j; i ++){
		v = s->n[i];
		if (check_in_exclude(e, v)) j = 1;
	}


	// printf("\nstart:");
	// for (i =0;i<s->i;i++) {printf("%u(%u) ", g->nd[s->n[i]].name,s->n[i]);}
	// printf("\nend:");
	// for (i =0;i<e->i;i++) {printf("%u(%u) ", g->nd[e->n[i]].name,e->n[i]);}
	// printf("\nex:");
	// for (i =0;i<r[2].i;i++) {printf("%u(%u) ", g->nd[r[2].n[i]].name,r[2].n[i] );}
	// printf("si:%d ei:%d %d\n",s->i,e->i,j );
	// exit(1);

	if (s->i && e->i && !j){
		bfs_r *bfs_ret = bfs_nodes_compound(g, s, &r[2], -1, -1, 0, perc);
		phs *ps = (phs *) bfs_ret->data2;
		ph *visited_edges = (ph *) bfs_ret->data1;

		// printf("\n#%d\n", ps->i);
		
		if (ps->i) rm_visited_edges(g, ps, visited_edges, m, 0);
		destroy_bfs_r(bfs_ret, 1);
	}
	for (i = 0; i < 3; i ++) destroy_sn(&r[i]);
	free (r);
}

void clean_complex_graph(graph *g, const int s, const int m, const float perc){
	clean_node_lable(g);
	ph *p = NULL;
	UINTL_T i = 0, j, k;
	while ((i = get_all_spots(g, &p, ++i, s)) > 0){
		// printf("%u\n", i);
		for (j = k = 0; j < p->i; j++){
			if (g->nd[p->n[j]].od > 1 || g->nd[p->n[j]].id > 1) k ++;
		}

		// printf("\n%u",k );
		// for (j = k = 0; j < p->i; j++){
		// 	if (g->nd[p->n[j]].od > 1 || g->nd[p->n[j]].id > 1) {printf(" %u",g->nd[p->n[j]].name);};
		// }
		
		if (k/2 < s * 200 && p->i/2 < s * 500) solve_spot_by_finding_paths_hete(g, p, m, perc);
		destroy_sn(p);
		free (p);
		p = NULL;
	}

	if (p) {
		destroy_sn(p);
		free (p);
	}
}

// phs *bfs_nodes_compound_get_path(khash_t(bfs) *h, const int s){
// 	ph_ *p;
// 	phs *ps = malloc(sizeof(phs));
// 	init_phs(ps, 20, s + 1);

// 	uint16_t l = 0;
// 	l |= (BFLAG_P1 | BFLAG_P2 | BFLAG_U);
	
// 	bfs *h_v;
// 	UINTL_T w, v;
// 	khiter_t k;
// 	for (k = kh_begin(h); k != kh_end(h); k++){
// 		if (kh_exist(h, k)){
// 			w = kh_key(h, k);
// 			h_v = &kh_val(h, k);
// 			if (h_v->l & l) {
// 				if (ps->i >= ps->im) reallocate_phs(ps, 100, s + 1);
// 				p = &ps->ph[ps->i++];
// 				p->i = 0;
// 				p->n[p->i++] = w;
// 				p->n[p->i++] = v = h_v->pnode;
// 				h_v->pnode = 0;
// 				while (!h_v->l & BFLAG_N){
// 					h_v = &kh_val(h, kh_get(bfs, h, v));
// 					p->n[p->i++] = v = h_v->pnode;
// 					assert (h_v->pnode);//TODO Check and del
// 					h_v->pnode = 0;
// 				}
// 			}
// 		}
// 	}
// 	return ps;
// }

void *bfs_nodes_compound_callback(UINTL_T v, UINTL_T w, khash_t(bfs) *h, void *opt){
	int *lastest = (int *)opt;
	bfs *bfs_w = &kh_val(h, kh_get(bfs, h, w));
	if (*lastest || !bfs_w->pnode) bfs_w->pnode = v;
	return NULL;
}


uint32_t find_v_index(bfs_node_info *n, uint32_t n_l, UINTL_T v){//n must be sort by start
	uint32_t low = 0, high = n_l, middle = 0;
	while (low < high) {
		middle = (low + high)/2;
		if (v == n[middle].start){
			return middle;
		}else if (v < n[middle].start) { 
			high = middle;
		}else if (v > n[middle].start) {
			low = middle + 1;
		}
	}
	assert (0);
	return -1;
}

// void get_z_paths_in_spot(graph *g, ph *n, const int max_len, const int min_len,  phs *ps){
// 	if (!ps->im) init_phs(ps, 100, max_len + 1);
// 	ph_ *p;
// 	ph **s_e_x = get_start_end_nodes(g, n);
// 	ph *x = s_e_x[2];

// 	UINTL_T i, j, v, t, k, w, m;
// 	for (i = 0; i < n->i; i++) {
// 		v = n->n[i];
// 		if (g->nd[v].od <= 1) continue;
// 		for (j = 0; j < g->nd[v].odm; j++){
// 			if (g->eg[g->nd[v].oe[j]].l & MFLAG_FIT) continue;
// 			w = g->eg[g->nd[v].oe[j]].ou;
// 			if (check_in_exclude(x, w)) continue;
// 			p = &ps->ph[ps->i];
// 			p->n[p->i++] = g->nd[v].oe[j];
// 			while (g->nd[w].id == 1 && p->i < max_len){
// 				for (m = k = 0; k < g->nd[w].odm; k++){
// 	        		if (g->eg[g->nd[w].oe[k]].l & MFLAG_FIT || check_in_exclude(x, g->eg[g->nd[w].oe[k]].ou)) continue;
// 	        		m ++;
// 	        		t = k;
// 	       		}
// 	       		if (m != 1) break;
// 	       		p->n[p->i++] = g->nd[w].oe[t];
// 	       		w = g->eg[g->nd[w].oe[t]].ou;
// 	       	}
// 			if (g->nd[w].id > 1 && p->i > min_thres) {
// 				ps->i++;
// 				if (ps.i >= ps.im) reallocate_phs(ps, 100, max_len + 1);
// 			}
// 		}
// 	}
// 	if (ps.i > 1){
// 		for (j = 0; j < ps.i; j ++){
// 			p = &ps.ph[j];
// 			n = g->eg[p->n[0]].in;
// 			p->sco = (uint64_t) g->eg[p->n[0]].sco * g->eg[p->n[0]].ide / 
// 						g->nd[n].oe[g->nd[n].odm] / 100;// /10000 * 100
// 			n = g->eg[p->n[p->i-1]].ou;
// 			p->sco += (uint64_t) g->eg[p->n[p->i-1]].sco * g->eg[p->n[p->i-1]].ide / 
// 			g->nd[n].ie[g->nd[n].idm] / 100;
// 			// p->sco -= g->nd[n].id * 5;
// 		}
// 		sort_phs(&ps);
// 	}

// 	for (i = 0; i < 3; i ++) destroy_sn(r[i]);
// 	free (s_e_x);
// }

static int sort_phs_by_len(const void *p1, const void *p2){
	return ((ph_ *) p2)->i - ((ph_ *) p1)->i;
}

void rm_long_bubble(graph *g, const int s, const int m, const float perc){
	phs ps;
	ph_ *p, *q;
	init_phs(&ps, 100, s + 1);

	UINTL_T i, j, k, t, n;
	int64_t max_sco;
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].od < 2) continue;
		ps.i = 0;
		max_sco = get_max_score(g, i, 1);
		for (j = 0; j < g->nd[i].odm; j++){
			if (g->eg[g->nd[i].oe[j]].l & MFLAG_FIT) continue;
			p = &ps.ph[ps.i];
			p->i = p->sco = 0;
			p->n[p->i++] = g->nd[i].oe[j];
			n = g->eg[g->nd[i].oe[j]].ou;
			
			while (g->nd[n].id == 1 && g->nd[n].od == 1 && p->i < s){
				t = get_validly_oe(g, n, 0);
				p->n[p->i++] = g->nd[n].oe[t];
				n = g->eg[g->nd[n].oe[t]].ou;
			}
			if (g->nd[n].id > 1){
				p->sco = (int64_t) g->eg[p->n[0]].sco * g->eg[p->n[0]].ide - max_sco * perc;
				for (k = 1; k < p->i; k++) p->sco += (int64_t) g->eg[p->n[k]].sco * g->eg[p->n[k]].ide * (1 - perc);
				ps.i++;
				if (ps.i >= ps.im) reallocate_phs(&ps, 100, s + 1);
			}
		}
		if (ps.i < 2) continue;

		//select the path with the largest ide of start & end edes as the final path.
		for (j = 0; j < ps.i; j ++){
			p = &ps.ph[j];
			if (p->sco == 0) continue;
			n = g->eg[p->n[p->i - 1]].ou;
			UINTL_T s = j, e = j;
			uint16_t s_ide = g->eg[p->n[0]].ide, e_ide = g->eg[p->n[p->i-1]].ide;

			for (k = 0; k < j; k ++){
				q = &ps.ph[k];
				if (g->eg[q->n[q->i - 1]].ou == n) break;
			}
			if (k < j) continue; //already checked

			for (k = j + 1; k < ps.i; k ++){
				q = &ps.ph[k];
				if (g->eg[q->n[q->i - 1]].ou == n){
					if (g->eg[q->n[0]].ide > s_ide) {
						s = k;
						s_ide = g->eg[q->n[0]].ide;
					}

					if (g->eg[q->n[q->i-1]].ide > e_ide){
						e = k;
						e_ide = g->eg[q->n[q->i-1]].ide;
					}	
				}
			}

			if (s == e) {
				for (k = j; k < ps.i; k ++){
					q = &ps.ph[k];
					if (k != s && g->eg[q->n[q->i - 1]].ou == n){
						q->sco = 0;
					}
				}
			}
		}

		msort_phs(&ps);
		for (j = 0; j < ps.i; j ++){
			p = &ps.ph[j];
			n = g->eg[p->n[p->i - 1]].ou;
			for (k = j + 1; k < ps.i; k ++){
				if (g->eg[ps.ph[k].n[ps.ph[k].i - 1]].ou == n){
					//don't rm bio loop
					// if (n == get_reversed_node(g, i) || (m && p->i >= m)){
						rm_edge(g, p->n[p->i-1]);
						rm_edge(g, get_reversed_edge(g, p->n[p->i-1]));
						rm_edge(g, p->n[0]);
						rm_edge(g, get_reversed_edge(g, p->n[0]));
					// }else{
					// 	while (--p->i >= 0){
					// 		rm_edge(g, p->n[p->i]);
					// 		rm_edge(g, get_reversed_edge(g, p->n[p->i]));
					// 	}
					// }
					break;
				}
			}
		}
	}
	destroy_phs(&ps);
}

// #if GENOME_SIZE == 2
//     KHASH_MAP_INIT_INT64(map_, UINTL_T)
// #else
//     KHASH_MAP_INIT_INT(map_, UINTL_T)
// #endif

// ph *bfs_nodes(graph *g, ph *n, ph *exclude, int depth, int max_child_num, int dire, 
// 		void *(*callback)(UINTL_T, khash_t(map_) *, void *), void *callback_opt){
// 	int absent, dp = 0;
// 	khint_t k;
// 	khash_t(map_) *r = kh_init(map_);

// 	UINTL_T i, j, v, w;
// 	for (i = 0; i < n->i; i++) kh_put(map_, r, n->n[i], &absent);//TODO CHECK kh_key
// 	// 	printf("%u %u\n", kh_key(r, kh_get(map_, r, n->n[i])), n->n[i]);

// 	ph *h, *h1;
// 	h = n;
// 	h1 = calloc(1, sizeof(ph));
// 	init_sn(h1, 100);

// 	while ((depth < 0 || dp++ < depth) && h->i){
// 		h1->i = 0;
// 		while (h->i){
// 			v = h->n[--h->i];
// 			if (dire & 0x1){
// 				for (j = 0; j < g->nd[v].odm; j++){
// 					if (g->eg[g->nd[v].oe[j]].l & MFLAG_FIT) continue;
// 					w = g->eg[g->nd[v].oe[j]].ou;
// 					if (check_in_exclude(exclude, w)) continue;
// 					if (max_child_num > 0 && (dire & 0x2 ? g->nd[w].id + g->nd[w].od : g->nd[w].od) >= max_child_num) continue;
// 					k = kh_put(map_, r, w, &absent);
// 					if (!absent) continue;
// 					else kh_value(r, k) = v;
// 					if (callback) {
// 						void *ret = callback(w, r, callback_opt);
// 						if (ret) goto END;
// 					}
// 					h1->n[h1->i++] = w;
// 					if (h1->i >= h1->im) reallocate_sn(h1, 1000);
// 				}
// 			}
// 			if (dire & 0x2){
// 				for (j = 0; j < g->nd[v].idm; j++){
// 					if (g->eg[g->nd[v].ie[j]].l & MFLAG_FIT) continue;
// 					w = g->eg[g->nd[v].ie[j]].in;
// 					if (check_in_exclude(exclude, w)) continue;
// 					if (max_child_num > 0 && (dire & 0x1 ? g->nd[w].id + g->nd[w].od : g->nd[w].id) >= max_child_num) continue;
// 					k = kh_put(map_, r, w, &absent);
// 					if (!absent) continue;
// 					else kh_value(r, k) = v;
// 					if (callback) {
// 						void *ret = callback(w, r, callback_opt);
// 						if (ret) goto END;
// 					}
// 					h1->n[h1->i++] = w;
// 					if (h1->i >= h1->im) reallocate_sn(h1, 1000);
// 				}
// 			}
// 		}
// 		SWAP(h, h1, ph*);
// 	}
// 	if (h != n) SWAP(h, h1, ph*);
// 	h1->i = 0;
// 	for (k = kh_begin(r) ; k < kh_end(r); k++){
// 		if (kh_exist(r, k)){
// 			h1->n[h1->i++] = kh_key(r, k);//TODO CHECK key
// 			if (h1->i >= h1->im) reallocate_sn(h1, 1000);
// 		}
// 	}
// 	kh_destroy(map_, r);
// 	return h1;
// END:
// 	kh_destroy(map_, r);
// 	if (h != n) SWAP(h, h1, ph*);
// 	destroy_sn(h1);
// 	free (h1);
// 	return NULL;
// }
