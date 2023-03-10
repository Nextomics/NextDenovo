#include <assert.h>
#include "asg.h"

UINTL_T get_reversed_node(graph *g, const UINTL_T n) {
	node *nd = g->nd;
	if (n > 0 ) return nd[n-1].name == nd[n].name ? n - 1 : n + 1;
	else return n + 1;
}

static void reallocate_node(graph *s, const uint32_t len){
	s->nm += len;
	s->nd = realloc(s->nd, sizeof(node) * s->nm);
	memset(s->nd + s->nm - len, 0, len * sizeof(node));
}

void rm_oe(graph *g, const UINTL_T n, const UINTL_T e){
	UINT_T i;
	for (i = 0; i < g->nd[n].od; i ++){
		if (g->nd[n].oe[i] == e){
			while (i < g->nd[n].od - 1) {
				g->nd[n].oe[i] = g->nd[n].oe[i + 1];
				i ++;
			}
			break;
		}
	}
	g->nd[n].od --;
}
void rm_ie(graph *g, const UINTL_T n, const UINTL_T e){
	UINT_T i;
	for (i = 0; i < g->nd[n].id; i ++){
		if (g->nd[n].ie[i] == e){
			while (i < g->nd[n].id - 1) {
				g->nd[n].ie[i] = g->nd[n].ie[i + 1];
				i ++;
			}
			break;
		}
	}
	g->nd[n].id --;
}

UINTL_T add_node(graph *g, const UINTL_T name, const int len) {
	if (g->ni >= g->nm) reallocate_node(g, MAX_INIT_N);
	UINTL_T n = g->ni++;
	g->nd[n].name = name;
	g->nd[n].idm = len;
	g->nd[n].ie = malloc (g->nd[n].idm * sizeof(UINTL_T));
	g->nd[n].odm = len;
	g->nd[n].oe = malloc (g->nd[n].odm * sizeof(UINTL_T));
	return n;
}

void rm_node(graph *g, const UINTL_T n){
	node *no = &g->nd[n];
	if (no->odm + no->idm != 0){
		UINT_T i;
		for (i = 0; no->id; i++){
			if (!(g->eg[no->ie[i]].l & MFLAG_FIT)){
				g->eg[no->ie[i]].l |= MFLAG_FIT;
				g->nd[g->eg[no->ie[i]].in].od --;
				// printf("n %d in %d ou %d\n",n, g->eg[no->ie[i]].in, g->eg[no->ie[i]].ou);
				no->id --;
			}
			// printf("$%d %d\n", i, no->id);
		}
		// printf("%d %u %u\n", n,no->od, no->id);
		for (i = 0; no->od; i++){
			if (!(g->eg[no->oe[i]].l & MFLAG_FIT)){
				g->eg[no->oe[i]].l |= MFLAG_FIT;
				g->nd[g->eg[no->oe[i]].ou].id --;
				no->od --;
			}
		}
		assert(no->od == 0 && no->id == 0);
		free (no->ie);
		free (no->oe);
		no->ie = no->oe = NULL;
		no->odm = no->idm = 0;
	}
}

UINT_T get_validly_oe(graph *g, const UINTL_T n, const UINT_T i) {
	UINT_T j, k;
	for (j = k = 0; j < g->nd[n].odm; j++){
		if (g->eg[g->nd[n].oe[j]].l & MFLAG_FIT) continue;
		else if (k++ == i) return j;
	}
	assert (0);
}

int add_oe(graph *g, const UINTL_T n, const UINTL_T e){
	if (g->nd[n].od >= g->nd[n].odm){
		g->nd[n].odm += 8;
		g->nd[n].oe = realloc(g->nd[n].oe, g->nd[n].odm * sizeof(UINTL_T));
	}
	g->nd[n].oe[g->nd[n].od++] = e;
	return 0;
}

UINT_T get_validly_ie(graph *g, const UINTL_T n, const UINT_T i) {
	UINT_T j, k;
	for (j = k = 0; j < g->nd[n].idm; j++){
		if (g->eg[g->nd[n].ie[j]].l & MFLAG_FIT) continue;
		else if (k++ == i) return j;
	}
	assert (0);
}

int add_ie(graph *g, const UINTL_T n, const UINTL_T e){
	if (g->nd[n].id >= g->nd[n].idm){
		g->nd[n].idm += 8;
		g->nd[n].ie = realloc(g->nd[n].ie, g->nd[n].idm * sizeof(UINTL_T));
	}
	g->nd[n].ie[g->nd[n].id++] = e;
	return 0;
}

void clean_node_lable(graph *g){
	UINTL_T i;
	for (i = 0; i < g->ni; i++) g->nd[i].l = 0;
}

UINTL_T stat_valid_node(graph *g){
	UINTL_T i, k;
	i = k = 0;
	for (; i < g->ni; i++){
		if (g->nd[i].id == 0 && g->nd[i].od == 0) k++;
	}
	return i-k;
}

int check_node_lable(graph *g, const UINTL_T n, const uint16_t l){
	UINT_T i;
	for (i = 0; i < g->nd[n].idm; i++){
		if (g->eg[g->nd[n].ie[i]].l & MFLAG_FIT) continue;
		if (!(g->eg[g->nd[n].ie[i]].l & l)) return 0;
	}

	for (i = 0; i < g->nd[n].odm; i++){
		if (g->eg[g->nd[n].oe[i]].l & MFLAG_FIT) continue;
		if (!(g->eg[g->nd[n].oe[i]].l & l)) return 0;
	}

	return 1;
}