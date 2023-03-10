#include <assert.h>
#include "asg.h"
#include "opt.h"

void reallocate_edge(graph *s, const uint32_t len){
	s->em += len;
	s->eg = realloc(s->eg, sizeof(edge) * s->em);
	memset(s->eg + s->em - len, 0, len * sizeof(edge));
}

static int check_valid_edge1(edge *e, opt *p, uint32_t mide, uint32_t msco){
	if (e->sco == msco){
		return 2;
	}else if (mide >= p->min_ide){
		if (e->ide > mide * p->min_ide_ratio) return 1;
	}else{
		if (e->sco >= msco * p->min_sco_ratio) return 1;
	}
	return 0;
}

int check_valid_edge(edge *e, opt *p, ovlinfo *loli, ovlinfo *roli){
	uint32_t mide, msco;

	if(e->l & MFLAG_IL){
		mide = loli->lim; msco = loli->llm;
	}else{
		mide = loli->rim; msco = loli->rlm;
	}
	int v = check_valid_edge1(e, p, mide, msco);
	// printf(" %u %u %f %f %d",e->ide,e->sco, mide* p->min_ide_ratio,msco* p->min_sco_ratio, v);	
	if (v >= p->min_node_count) return v;

	if(e->l & MFLAG_OL){// should check reversed end
		mide = roli->rim; msco = roli->rlm;
	}else{
		mide = roli->lim; msco = roli->llm;
	}
	v += check_valid_edge1(e, p, mide, msco);
	// printf(" %f %f %d",mide* p->min_ide_ratio,msco* p->min_sco_ratio,v);	
	return v;
}

uint32_t check_exited_edge(graph *g, const UINTL_T in, const UINTL_T ou_){// only used in update_graph
	int j;
	UINT_T i;
	UINTL_T ou = ou_;
	for (j = 0; j < 2; j++){
		for (i = 0; i < g->nd[in].od; i ++){
			if (g->eg[g->nd[in].oe[i]].ou == ou) return g->eg[g->nd[in].oe[i]].sco;
		}
		for (i = 0; i < g->nd[in].id; i ++){
			if (g->eg[g->nd[in].ie[i]].in == ou) return g->eg[g->nd[in].ie[i]].sco;
		}
		ou = get_reversed_node(g, ou_);
	}
	return 0;
}

UINTL_T get_edge(graph *g, const UINTL_T in, const UINTL_T ou){
	UINT_T i;
	for (i = 0; i < g->nd[in].od; i ++){
		if (g->eg[g->nd[in].oe[i]].ou == ou) return g->nd[in].oe[i];
	}
	return -1;
}

UINTL_T rp_exited_edge(graph *g, const UINTL_T in_, const UINTL_T ou_){// replace existed edge
	int j;
	UINT_T i;
	UINTL_T e, in = in_, ou = ou_;
	// printf("rp_exited_edge: %u %u\n", in, ou);
	for (j = 0; j < 3; j++){
		// printf("%d %u %u\n",j, in, ou);
		for (i = 0; i < g->nd[in].od; i ++){
			if (g->eg[g->nd[in].oe[i]].ou == ou){
				e = g->nd[in].oe[i];
				rm_oe(g, in, e);
				rm_ie(g, ou, e);
				return e;
			}
		}
		for (i = 0; i < g->nd[in].id; i ++){
			if (g->eg[g->nd[in].ie[i]].in == ou) {
				e = g->nd[in].ie[i];
				rm_oe(g, ou, e);
				rm_ie(g, in, e);
				return e;
			}
		}
		if (j == 0) {
			in = in_;
			ou = get_reversed_node(g, ou_);
		}else if (j == 1){
			in = get_reversed_node(g, in_);
			ou = ou_;
		}
	}
	assert(0);
	return 0;
}

UINTL_T get_reversed_edge(graph *g, const UINTL_T n){
	UINT_T i;
	UINTL_T n1 = get_reversed_node(g, g->eg[n].in);
	UINTL_T n2 = get_reversed_node(g, g->eg[n].ou);
	for (i = 0; i < g->nd[n1].idm; i ++){
		if (g->eg[g->nd[n1].ie[i]].in == n2) return g->nd[n1].ie[i];
	}
	return n;//a node has been deleted and reversed edge did not existed
}

UINTL_T add_edge(graph *g, const UINTL_T in, const UINTL_T ou, const uint32_t ie, const uint32_t oe, \
		const uint32_t sco, const uint16_t ide, const uint32_t len, const uint32_t l) {
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

void rm_edge(graph *g, const UINTL_T n){
	if (!(g->eg[n].l & MFLAG_FIT)){
		g->eg[n].l |= MFLAG_FIT;
		g->nd[g->eg[n].in].od --;
		g->nd[g->eg[n].ou].id --;
	}
}


void rc_rm_edge(graph *g, const UINTL_T n){//recover rm_edge
	if (g->eg[n].l & MFLAG_FIT){
		g->eg[n].l ^= MFLAG_FIT;
		g->nd[g->eg[n].in].od ++;
		g->nd[g->eg[n].ou].id ++;
	}
}

uint32_t get_oreadlen(graph *g, UINTL_T n){
	return g->eg[n].len + g->eg[n].sco;
}

uint32_t get_ireadlen(graph *g, UINTL_T n){
	n = get_reversed_edge(g, n);
	return g->eg[n].len + g->eg[n].sco;
}

UINTL_T stat_valid_edge(graph *g){
    UINTL_T i, k;
    i = k = 0;
    for (; i < g->ei; i++){
        if (!(g->eg[i].l & MFLAG_FIT)) k++;
    }   
    return k;
}