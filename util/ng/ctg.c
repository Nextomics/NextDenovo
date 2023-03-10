#include <assert.h>
#include "ctg.h"
#include "asg.h"
#include "kit.h"
#include "node.h"
#include "edge.h"

static void reallocate_ctg(ctg_ *s, const uint32_t len){
	s->im += len;
	s->nd = realloc(s->nd, sizeof(ctg__) * s->im);
}

void init_ctgs(ctgs *s, const uint32_t len){
	uint32_t i;
	ctg_ *ctg;
	s->i = 0;
	s->im = len;
	s->ctg = malloc(s->im * sizeof(ctg_));
	for (i = 0; i < s->im; i++){
		ctg = &s->ctg[i];
		ctg->im = MAX_INIT_N;
		ctg->len = ctg->type = ctg->i = 0;
		ctg->nd = malloc(ctg->im * sizeof(ctg__));
	}
}

static void reallocate_ctgs(ctgs *s, const uint32_t len){
	uint32_t i;
	ctg_ *ctg;
	s->im += len;
	s->ctg = realloc(s->ctg, sizeof(ctg_) * s->im);
	for (i = s->im - len; i < s->im; i++){
		ctg = &s->ctg[i];
		ctg->im = MAX_INIT_N;
		ctg->len = ctg->type = ctg->i = 0;
		ctg->nd = malloc(ctg->im * sizeof(ctg__));
	}
}

void destroy_ctgs(ctgs *s){
	uint32_t i;
	for (i = 0; i < s->im; i++){
		free (s->ctg[i].nd);
	}
	free (s->ctg);
}

static uint64_t get_max_score(graph *g, UINTL_T n, int ou){ 

    UINTL_T i;
    uint64_t score, temp;
    if (ou){
        for (score = i = 0; i < g->nd[n].odm; i++){
            temp = (uint64_t) g->eg[g->nd[n].oe[i]].sco * g->eg[g->nd[n].oe[i]].ide;
            if (temp > score) score = temp;
        }    
    }else{
        for (score = i = 0; i < g->nd[n].idm; i++){
            temp = (uint64_t) g->eg[g->nd[n].ie[i]].sco * g->eg[g->nd[n].ie[i]].ide;
            if (temp > score) score = temp;
        }    
    }    
    return score;
}

void generate_ctg(graph *g, ctgs *s, const uint16_t l){
	UINTL_T i, j;
	ph p;
	init_sn(&p, 1000);
	for (i = 1; i < g->ni; i++){
		if (g->nd[i].id + g->nd[i].od == 0 || (g->nd[i].id == 1 && g->nd[i].od <= 1)) continue;
		p.n[p.i++] = i;
		if (p.i >= p.im) reallocate_sn(&p, 1000);
	}
	node *n = NULL;
	edge *e = NULL;
	ctg__ *nd = NULL;
	ctg_ *m = NULL;
	// init_ctg(&m, MAX_INIT_N);
	while (--p.i >= 0){
		j = p.n[p.i];
		n = &g->nd[j];
		if (n->l) continue;
		// printf("%d %d\n",s->i, s->im);
		if (s->i >= s->im) reallocate_ctgs(s, 100);
		m = &s->ctg[s->i++];
		// m.i = m.type = m.len = 0;
		if (!n->od){// single node
			nd = &m->nd[m->i++];
			// e = &g->eg[n->ie[0]];
			nd->name = n->name;
			nd->nidex = j;
			nd->ide = nd->ort = nd->irt = 1000;

			nd->s = nd->l = nd->lq = 0;
			nd->e = UINT32_MAX;
			m->type = 3;

			n->l = 1;
			n = &g->nd[get_reversed_node(g, j)];
			n->l = 1;
			for (i = 0; i < n->odm; i++){
				if (g->eg[n->oe[i]].l & MFLAG_FIT) continue;
				p.n[p.i++] = g->eg[n->oe[i]].ou;
				if (p.i >= p.im) reallocate_sn(&p, 1000);
			}
		}else if (n->id > 1 && n->od > 1){// single node

			nd = &m->nd[m->i++];
			// e = &g->eg[n->ie[0]];
			nd->name = n->name;
			nd->nidex = j;
			nd->ide = nd->ort = nd->irt = 1000;

			nd->s = nd->l = 0;
			nd->e = UINT32_MAX;
			nd->lq = 0;
			m->type = 4;

			n->l = 1;
			for (i = 0; i < n->odm; i++){
				if (g->eg[n->oe[i]].l & MFLAG_FIT) continue;
				p.n[p.i++] = g->eg[n->oe[i]].ou;
				if (p.i >= p.im) reallocate_sn(&p, 1000);
			}
			n = &g->nd[get_reversed_node(g, j)];
			n->l = 1;
			for (i = 0; i < n->odm; i++){
				if (g->eg[n->oe[i]].l & MFLAG_FIT) continue;
				p.n[p.i++] = g->eg[n->oe[i]].ou;
				if (p.i >= p.im) reallocate_sn(&p, 1000);
			}
		}else if (n->od > 1){//reversed node will be output
			for (i = 0; i < n->odm; i++){
				if (g->eg[n->oe[i]].l & MFLAG_FIT) continue;
				p.n[p.i++] = g->eg[n->oe[i]].ou;
				if (p.i >= p.im) reallocate_sn(&p, 1000);
			}
			s->i --;
			continue;
		}else{//od==1 && id>1 || id == 0
			nd = &m->nd[m->i++];
			// e = &g->eg[n->ie[0]];
			nd->name = n->name;
			nd->nidex = j;
			nd->ide = nd->ort = nd->irt = 1000;
			nd->s = 0;
			// nd->l = e->l & MFLAG_OL ? 1 : 0;
			m->type = 1;
			
			n->l = 1;
			n = &g->nd[get_reversed_node(g, j)];
			n->l = 1;			
			for (i = 0; i < n->odm; i++){
				if (g->eg[n->oe[i]].l & MFLAG_FIT) continue;
				p.n[p.i++] = g->eg[n->oe[i]].ou;
				if (p.i >= p.im) reallocate_sn(&p, 1000);
			}

			n = &g->nd[j];
			// printf("%d\n", p.n[p.i]);
			i = get_validly_oe(g, j, 0);
			e = &g->eg[n->oe[i]];
			i = e->ou;
			n = &g->nd[i];
			if (e->l & MFLAG_IL){
				nd->l = 1;
				nd->s = UINT32_MAX;
			}else{
				nd->l = 0;
			}
			// if (e->l & l) m->i--;
			// nd->l = e->l & MFLAG_IL ? 1 : 0;
			while (n->id == 1 && n->od == 1 && !n->l){		
				n->l = 1;
				g->nd[get_reversed_node(g, i)].l = 1;
				nd->e = e->ie;
				nd->lq = e->l & l ? 1 : 0;

				nd = &m->nd[m->i++];
				if (m->i >= m->im) reallocate_ctg(m, MAX_INIT_N);
				nd->name = n->name;
				nd->nidex = i;
				nd->ide = e->ide/10;
				nd->ort = (uint64_t) e->sco * e->ide * 1000 / get_max_score(g, e->in, 1);
				nd->irt = (uint64_t) e->sco * e->ide * 1000 / get_max_score(g, e->ou, 0);

				nd->s = e->oe;
				nd->l = e->l & MFLAG_OL ? 1 : 0;

				i = get_validly_oe(g, i, 0);
				e = &g->eg[n->oe[i]];
				i = e->ou;
				n = &g->nd[i];
			}

			nd->e = e->ie;
			nd->lq = e->l & l ? 1 : 0;
			if (n->id > 1 || n->l) {
				// printf("%u %u\n",i, m->i);
				// m->i --;
				continue;
			}
			// printf("%d\n", m.i);
			n->l = 1;
			g->nd[get_reversed_node(g, i)].l = 1;
			
			//for end node
			nd = &m->nd[m->i++];
			nd->name = n->name;
			nd->nidex = i;
			nd->ide = e->ide/10;
			nd->ort = (uint64_t) e->sco * e->ide * 1000 / get_max_score(g, e->in, 1);
			nd->irt = (uint64_t) e->sco * e->ide * 1000 / get_max_score(g, e->ou, 0);

			nd->s = e->oe;
			if (e->l & MFLAG_OL){
				nd->l = 1;
				nd->e = 0;
			}else{
				nd->l = 0;
				nd->e = UINT32_MAX;
			}
			nd->lq = e->l & l ? 1 : 0;
			// if (e->l & l) m->i--;
			// nd->l = e->l & MFLAG_OL ? 1 : 0;
			// nd->e = 0;

			for (i = 0; i < n->odm; i++){
				if (g->eg[n->oe[i]].l & MFLAG_FIT) continue;
				p.n[p.i++] = g->eg[n->oe[i]].ou;
				if (p.i >= p.im) reallocate_sn(&p, 1000);
			}
		}
		// out_ctg_path(&m);
	}
	destroy_sn(&p);
	// out_node(g, 195015);
	for (j = 1; j < g->ni; j++){

		if (g->nd[j].l || g->nd[j].id + g->nd[j].od == 0) continue;
		// printf("%d %d %d %d\n", j, g->nd[j].id, g->nd[j].od, p.i);
		assert(g->nd[j].id == 1 && g->nd[j].od == 1);
		
		if (s->i >= s->im) reallocate_ctgs(s, 100);
		m = &s->ctg[s->i++];
		m->type = 2;

		n = &g->nd[j];
		nd = &m->nd[m->i++];
		// e = &g->eg[n->ie[0]];
		nd->name = n->name;
		nd->nidex = j;
		nd->ide = nd->ort = nd->irt = 1000;
		nd->s = 0;
		// nd->l = e->l & MFLAG_OL ? 1 : 0;

		i = get_validly_oe(g, j, 0);
		e = &g->eg[n->oe[i]];
		i = e->ou;
		n = &g->nd[i];
		// nd->l = e->l & MFLAG_IL ? 1 : 0;
		if (e->l & MFLAG_IL){
			nd->l = 1;
			nd->s = UINT32_MAX;
		}else{
			nd->l = 0;
		}
		// if (e->l & l) m->i--;

		while (n->id == 1 && n->od == 1 && !n->l){
			n->l = 1;
			g->nd[get_reversed_node(g, i)].l = 1;
			nd->e = e->ie;
			nd->lq = e->l & l ? 1 : 0;

			nd = &m->nd[m->i++];
			if (m->i >= m->im) reallocate_ctg(m, MAX_INIT_N);
			nd->name = n->name;
			nd->nidex = i;

			nd->ide = e->ide/10;
			nd->ort = (uint64_t) e->sco * e->ide * 1000 / get_max_score(g, e->in, 1);
			nd->irt = (uint64_t) e->sco * e->ide * 1000 / get_max_score(g, e->ou, 0);
			nd->s = e->oe;
			nd->l = e->l & MFLAG_OL ? 1 : 0;			
			i = get_validly_oe(g, i, 0);
			e = &g->eg[n->oe[i]];
			i = e->ou;
			n = &g->nd[i];
		}
		m->i --;
		// if (e->l & l) m->i--;
	}
}
