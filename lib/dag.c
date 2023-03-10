#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h> 
#include <assert.h>

#ifdef LGS_CORRECT
        #include "ctg_cns.h"
#else
        #include "nextcorrect.h"
#endif

// #include "kseq.h" 
// KSEQ_INIT(gzFile, gzread)

#define SEQ_MAX_COUNT 50
#define SCORE_MATCH 1
#define SCORE_MISMATCH -2
#define SCORE_GAP -2
#define GETMATCHSCORE(x, y) ((x) == (y) ? SCORE_MATCH : SCORE_MISMATCH)

typedef struct {
	uint8_t base;
	uint8_t indegree; //record inedge_index + 1
	uint8_t outdegree; //record outedge_index + 1
	uint32_t inedge[SEQ_MAX_COUNT];
	uint32_t outedge[SEQ_MAX_COUNT];
	uint16_t *alignedto;
	uint16_t alignedto_count;
	uint16_t alignedto_max;
	int32_t best_pnode; //used to consesus
	double best_score; //used to consesus
} node;

typedef struct {
	uint16_t innode_index;
	uint16_t outnode_index;
	uint8_t lable[SEQ_MAX_COUNT];
} edge;

typedef struct {
	int32_t x;
	int32_t y;
} matchroute;

typedef struct {
	matchroute *matchroute_;
	int64_t starty;
	int64_t endy; 
	uint32_t mroute_count;
} matchroutes;

typedef struct {
	uint16_t x;
	uint16_t y;
	long int s;
} score;

typedef struct {
	node *nodes;
	edge *edges;
	uint16_t *sorted_nodes;
	int32_t sorted_nodes_index;//need record -1
	uint16_t start_node[SEQ_MAX_COUNT]; // record the start node for each seq 
	uint16_t node_count;
	uint16_t node_max;
	uint32_t edge_count;
	uint32_t edge_max;
} graph;


// static uint8_t int_to_base[] = { 
// 	65, 84, 71, 67, 45 //A, T, G, C, -
// };

// static uint8_t base_to_int[] = { 
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 0, 4, 3,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  1, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 0, 4, 3,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
// 	4, 4, 4, 4,  1, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
// };

static score **score_init(const uint32_t x, const uint32_t y, graph *g, uint16_t *sorted_nodes_index){
	uint16_t i, k;
	// printf("%d %d\n",x, y );
	// for (i = 0; i < g->node_count; i++){
	// 	printf("%d ",  g->sorted_nodes[i]);
	// }
	// printf("\n");
	// score (*s)[y] = calloc(y, sizeof(score) * x);
	score **s = calloc(1, x * sizeof(score *) + x * y * sizeof(score));
	assert (s);
	score * const _ = (score *) (s + x);
	for (i = 0; i < x; i++){
		s[i] = _ + y * i;
	}

	for (i = 0; i < y; i++){
		s[0][i].s = i * SCORE_GAP;
	}

	long int bs_, s_;
	uint16_t node_index;
	for (i = 0; i < g->node_count; i++){
		node_index = g->sorted_nodes[i];
		sorted_nodes_index[node_index] = i;
		if (g->nodes[node_index].indegree == 0){
			bs_ = 0;
		}else{
			bs_ = s[sorted_nodes_index[g->edges[g->nodes[node_index].inedge[0]].innode_index] + 1][0].s;
			for (k = 1; k < g->nodes[node_index].indegree; k++){
				s_ = s[sorted_nodes_index[g->edges[g->nodes[node_index].inedge[k]].innode_index] + 1][0].s;
				if (s_ > bs_){
					bs_ = s_;
				}
			}
			// printf("%d %d\n", i+1, bs_);
		}
		s[i + 1][0].s = bs_ + SCORE_GAP;
	}
	// for (i=0; i< x; i++){
	// 	printf("\n%d ", i);
	// 	for (k=0; k<y;k++){
	// 		printf("%ld ",s[i][k].s );
	// 	}
	// }
	// printf("\n\n");
	return s;
}

static void graph_init(graph *g){
	g->node_count = g->edge_count = 0;
	g->node_max = 10000;
	g->edge_max = 2 * g->node_max;
	g->nodes = malloc(g->node_max * sizeof (node));
	assert (g->nodes);
	g->edges = calloc(g->edge_max, sizeof (edge));
	assert (g->edges);
	g->sorted_nodes = malloc(g->node_max * sizeof (uint16_t));
	for (uint32_t i = 0; i < g->node_max; i++){
		g->nodes[i].indegree = 0;
		g->nodes[i].outdegree = 0;
		g->nodes[i].alignedto_count = 0;
		g->nodes[i].alignedto_max = 100;
		g->nodes[i].alignedto = malloc(g->nodes[i].alignedto_max * sizeof (uint16_t));
	}
}

static void graph_destroy(graph *g){
	for (int i = 0; i <  g->node_max; i++){
		free (g->nodes[i].alignedto);
	}
	free (g->nodes);
	free (g->sorted_nodes);
	free (g->edges);
}

static void reallocate_nodes_mem(graph *g, const int number){
	// printf("1 %d %d\n", g->node_max, number);
	g->node_max += number;
	// printf("2 %d %d\n", g->node_max, number);

	g->nodes = realloc(g->nodes, g->node_max * sizeof (node));
	g->sorted_nodes = realloc(g->sorted_nodes, g->node_max * sizeof (uint16_t));
	for (uint32_t i = g->node_max - number; i < g->node_max; i++){
		g->nodes[i].indegree = 0;
		g->nodes[i].outdegree = 0;
		g->nodes[i].alignedto_count = 0;
		g->nodes[i].alignedto_max = 100;
		g->nodes[i].alignedto = malloc(g->nodes[i].alignedto_max * sizeof (uint16_t));
	}
}

static void reallocate_node_alignedto_mem(node *node_, const int number){
	node_->alignedto_max += number;
	node_->alignedto = realloc(node_->alignedto, node_->alignedto_max * sizeof (uint16_t));
}

static void reallocate_edges_mem(graph *g, const int number){
	g->edge_max += number;
	g->edges = realloc(g->edges, g->edge_max * sizeof (edge));
	memset(g->edges + g->edge_max - number, 0, number * sizeof (edge));
}

static inline void insert_node_alignedto(graph *g, uint16_t insert_node, uint16_t match_node){

	if (g->nodes[insert_node].alignedto_count + g->nodes[match_node].alignedto_count >= g->nodes[insert_node].alignedto_max - 10){//-10 was desigend
		reallocate_node_alignedto_mem(&g->nodes[insert_node], 100);
	}

	uint16_t i = 0;
	g->nodes[insert_node].alignedto[g->nodes[insert_node].alignedto_count++] = match_node;
	for (;i < g->nodes[match_node].alignedto_count; i++){
		g->nodes[insert_node].alignedto[g->nodes[insert_node].alignedto_count++] = g->nodes[match_node].alignedto[i];
	}

}

static inline int32_t insert_node(graph *g, const char base){
	g->nodes[g->node_count].base = base;
	if (++g->node_count >= g->node_max){
		reallocate_nodes_mem(g, 1000);
	}
	return g->node_count - 1;
}

static inline uint32_t insert_edge(graph *g, const uint16_t innode_index, const uint16_t outnode_index, const uint8_t lable){//need check edge is existed.
	g->edges[g->edge_count].innode_index = innode_index;
	g->edges[g->edge_count].outnode_index = outnode_index;
	g->edges[g->edge_count].lable[lable] = 1;
	// printf("edges: %d %d\n",lable,g->edge_count,g->edge_max );
	if (++g->edge_count >= g->edge_max){
		reallocate_edges_mem(g, 2000);
	}
	return g->edge_count - 1;
}

static inline int insert_lable_to_edge(graph *g, const uint16_t innode_index, const uint16_t outnode_index, const uint8_t lable){
	int not_existed = 1;//record edges existed
	for (uint16_t i = 0; i < g->nodes[innode_index].outdegree; i++){
		if (g->edges[g->nodes[innode_index].outedge[i]].outnode_index == outnode_index){
			g->edges[g->nodes[innode_index].outedge[i]].lable[lable] = 1;
			not_existed = 0;
			// printf("innode_index: %d outnode_index:%d seq:%d\n",innode_index, outnode_index, lable );
		}
	}
	return not_existed;
	// printf("\n");
}

static void insert_unmatched_nodes(const size_t seq_index, const char *seq, const size_t seq_len, graph *g, int32_t *firstnode, int32_t *headnode){
	uint16_t node_index;
	uint32_t edge_index;
	for (size_t i = 0; i < seq_len; i++){
		node_index = insert_node(g, seq[i]);
		if (*firstnode == -1) {
			*firstnode = node_index;
		}else{
			edge_index = insert_edge(g, *headnode, node_index, seq_index);
			g->nodes[*headnode].outedge[g->nodes[*headnode].outdegree ++ ] = edge_index;
			g->nodes[node_index].inedge[g->nodes[node_index].indegree ++ ] = edge_index;
		}
		*headnode = node_index;
	}
}

static inline uint16_t check_nodes_predecessors(graph *g, const uint16_t i){
	uint16_t predecessor_count = g->nodes[i].indegree;

	for (uint16_t j = 0; j < g->nodes[i].alignedto_count && predecessor_count == 0; j++){
		predecessor_count += g->nodes[g->nodes[i].alignedto[j]].indegree;
	}
	return predecessor_count;// return if predecessor_count != 0
}

static void align_seq_to_graph_updatescore(score **s, const uint16_t x, const uint16_t y, const char *seq, graph *g, uint16_t *sorted_nodes_index){
	int32_t pi;
	uint16_t i, j, k, node_index, bestx, besty;
	long int bests, bests1_, bests2_;
	for (g->sorted_nodes_index = 0; g->sorted_nodes_index < g->node_count; g->sorted_nodes_index++){
		node_index = g->sorted_nodes[g->sorted_nodes_index];

		// printf("#%d ",  node_index);

		i = sorted_nodes_index[node_index];
		for (j = 0; j < y; j++){
			bests = s[i+1][j].s + SCORE_GAP;bestx = i + 1;besty = j;
			// printf("\n%ld %d %d ins ", bests, i + 1, j);
			for (k = 0; k < g->nodes[node_index].indegree; k ++){
				pi = sorted_nodes_index[g->edges[g->nodes[node_index].inedge[k]].innode_index];
				bests1_ = s[pi+1][j+1].s +  SCORE_GAP;
				bests2_ = s[pi+1][j].s + GETMATCHSCORE(seq[j], g->nodes[node_index].base);

				
				// printf("%ld %d %d del ", bests1_, pi + 1, j+1);
				// printf("%ld %d %d match %ld %d %c %c", bests2_, pi + 1, j,  s[pi+1][j].s, GETMATCHSCORE(seq[j], g->nodes[node_index].base), seq[j], g->nodes[node_index].base);


				if (bests1_ > bests && bests1_ >= bests2_) {bests = bests1_;bestx = pi + 1;besty = j + 1;}
				else if (bests2_ > bests && bests2_ >= bests1_) {bests = bests2_;bestx = pi + 1;besty = j;}
			}
			if (g->nodes[node_index].indegree == 0){
				pi = -1;
				bests1_ = s[pi+1][j+1].s + SCORE_GAP;
				bests2_ = s[pi+1][j].s + GETMATCHSCORE(seq[j], g->nodes[node_index].base);
				if (bests1_ > bests && bests1_ >= bests2_) {bests = bests1_;bestx = pi + 1;besty = j + 1;}
				else if (bests2_ > bests && bests2_ >= bests1_) {bests = bests2_;bestx = pi + 1;besty = j;}
				// printf("%ld %d %d del1 ", bests1_, pi + 1, j+1);
				// printf("%ld %d %d match1 ", bests2_, pi + 1, j);
			}
			// printf("\n%d %d %ld %d %d\n",i+1, j+1, bests,bestx, besty );
			s[i+1][j+1].s = bests;s[i+1][j+1].x = bestx;s[i+1][j+1].y = besty;
		}
	}
}

static uint16_t align_seq_to_graph_getbestx(const uint16_t y, score **s, graph *g){
	uint16_t i, j, bestx;
	long int bests, bests_;
	i = j = bests = bestx = 0;
	for (; i < g->node_count; i ++){
		if (g->nodes[g->sorted_nodes[i]].outdegree == 0){
			bests_ = s[i + 1][y].s;
			if (j == 0 || bests_ > bests){ bestx = i + 1; bests = bests_; j = 1;}
		}
	}
	return bestx;
}

static void reverse_matchroute(matchroute *m, const uint32_t mroute_count){

	matchroute tmp;
	uint32_t left = 0;
	uint32_t right = mroute_count - 1;
	while (left < right){
		tmp = m[left];
		m[left++] = m[right];
		m[right--] = tmp;
	}
}

static void align_seq_to_graph_updatemroute(uint16_t bestx, uint16_t besty, score **s, uint16_t *sorted_nodes, matchroutes *m){
	uint16_t nextx, nexty;
	while (bestx != 0 || besty != 0){
		nextx = s[bestx][besty].x;
		nexty = s[bestx][besty].y;
		if (nextx != bestx) m->matchroute_[m->mroute_count].x = sorted_nodes[bestx - 1];
		if (nexty != besty) {
			m->matchroute_[m->mroute_count].y = m->starty = besty - 1;
			if (m->endy == -1) m->endy = m->matchroute_[m->mroute_count].y;
		}
		bestx = nextx;
		besty = nexty;
		m->mroute_count ++;
		// printf("bestx: %d besty:%d %d %d\n", bestx, besty, s[bestx][besty].x, s[bestx][besty].y);
	}
	reverse_matchroute(m->matchroute_, m->mroute_count);
}

static void align_seq_to_graph_updategraphy(const uint16_t y, const size_t seq_index, const char *seq, graph *g, matchroutes *m){
	int32_t i, j, firstnode, headnode, tailnode, node_index, edge_index, foundnode;
	firstnode = headnode = tailnode = node_index = foundnode = -1;
	char base;
	int updated_node, updated_headnode;// record whether the edge exixted, 0 = existed
	updated_node = updated_headnode = 1;
	if (m->starty > 0) {
		insert_unmatched_nodes(seq_index, seq, m->starty, g, &firstnode, &headnode);
	}
	if (m->endy < y - 1) insert_unmatched_nodes(seq_index, seq + m->endy + 1, y - m->endy, g, &tailnode, &node_index);
	for (i = 0; i < m->mroute_count; i++){
		// printf("#####%d %d %c\n",seq_index, i, seq[m->matchroute_[i].y]);
		if (m->matchroute_[i].y == -1) continue;
		updated_node = 0;
		base = seq[m->matchroute_[i].y];
		if (m->matchroute_[i].x == -1) updated_node = node_index = insert_node(g, base);
		else if (g->nodes[m->matchroute_[i].x].base == base) node_index = m->matchroute_[i].x;
		else{
			foundnode = -1;
			for (j = 0; j < g->nodes[m->matchroute_[i].x].alignedto_count; j++){
				if (g->nodes[g->nodes[m->matchroute_[i].x].alignedto[j]].base == base) node_index = foundnode = g->nodes[m->matchroute_[i].x].alignedto[j];
			}
			if (foundnode == -1){
				updated_node = node_index = insert_node(g, base);
				insert_node_alignedto(g, node_index, m->matchroute_[i].x);
				for (j = 0; j < g->nodes[node_index].alignedto_count; j++){
					g->nodes[g->nodes[node_index].alignedto[j]].alignedto[g->nodes[g->nodes[node_index].alignedto[j]].alignedto_count++] = node_index;
				}
			}
		}
		if (headnode != -1){
			// printf("%d ", headnode);
			if (updated_node || updated_headnode){
				// printf("needupdate\n");
				edge_index = insert_edge(g, headnode, node_index, seq_index);
				g->nodes[headnode].outedge[g->nodes[headnode].outdegree ++] = edge_index;
				g->nodes[node_index].inedge[g->nodes[node_index].indegree ++] = edge_index;
			}else{
				// printf("noupdate %d %d %d\n", headnode, node_index, seq_index);
				if (insert_lable_to_edge(g, headnode, node_index, seq_index)){
					edge_index = insert_edge(g, headnode, node_index, seq_index);
					g->nodes[headnode].outedge[g->nodes[headnode].outdegree ++] = edge_index;
					g->nodes[node_index].inedge[g->nodes[node_index].indegree ++] = edge_index;
				}
			}
		}
		headnode = node_index;
		updated_headnode = updated_node;
		if (firstnode == -1) firstnode = headnode;
	}
	if (tailnode != -1){
		edge_index = insert_edge(g, headnode, tailnode, seq_index);
		g->nodes[headnode].outedge[g->nodes[headnode].outdegree ++] = edge_index;
		g->nodes[tailnode].inedge[g->nodes[tailnode].indegree ++] = edge_index;
	}
	g->start_node[seq_index] = firstnode;
}

static void sort(const int32_t found, const uint16_t pnid_count, const uint16_t *pn_to_nodes, \
	const int32_t *node_to_pn, int8_t *completed, graph *g, uint32_t *stack_max, uint16_t *stack){
	uint16_t j, k, stack_count, started_count, pnid;
	int8_t started[pnid_count];
	memset(started, -1, sizeof(int8_t) * pnid_count);
	stack_count = started_count = 0;
	stack[stack_count ++] = found;
	// printf("#1 %d %d %d %d\n", stack_count, found, g->sorted_nodes_index, g->node_count);//1 20 25
	while (stack_count){
		pnid = stack[ -- stack_count];
		// printf("#2 %d %d %d %d %d %d\n", pnid, started[pnid], completed[pnid], pnid_count, stack_count, g->sorted_nodes_index);//20 1 40 46 0 3
		if (completed[pnid] == 1) continue;
		// printf("3 %d %d\n", g->nodes[pn_to_nodes[pnid]].outdegree, completed[pnid]);

		if (started[pnid] != -1){
			// printf("inloop\n");
			completed[pnid] = 1;
			// printf("sorted_nodes_index: %d \n", g->sorted_nodes_index);
			g->sorted_nodes[g->sorted_nodes_index--] = pn_to_nodes[pnid];
			for (j = 0; j < g->nodes[pn_to_nodes[pnid]].alignedto_count; j++){
				g->sorted_nodes[g->sorted_nodes_index--] = g->nodes[pn_to_nodes[pnid]].alignedto[j];
			}
			started[pnid] = -1;
			// printf("#3 %d %d %d\n",pnid, completed[0], completed[1]);
			continue;
		}

		started[pnid] = 1;


		stack[stack_count ++] = pnid;

		// printf("3 %d\n", g->nodes[pn_to_nodes[pnid]].outdegree);
		// printf("\n# %d %d ", stack_count, g->node_count);
		// for (j=0; j<stack_count; j++){printf(" %d", stack[j]);}

		if (*stack_max < stack_count + g->nodes[pn_to_nodes[pnid]].outdegree) { 
			*stack_max +=  g->node_count; 
			stack = realloc(stack, *stack_max * sizeof(uint16_t));
		}

		for (k = 0; k < g->nodes[pn_to_nodes[pnid]].outdegree; k++){
			stack[stack_count ++] = node_to_pn[g->edges[g->nodes[pn_to_nodes[pnid]].outedge[k]].outnode_index];
			// printf("\n%d %d\n", k, node_to_pn[g->edges[g->nodes[pn_to_nodes[pnid]].outedge[k]].outnode_index]);
		}
		
		// printf("\n%d %d\n", k, g->nodes[pn_to_nodes[pnid]].alignedto_count);
		// printf("\n$ %d %d %d", stack_count, g->node_count);
		
		// for (j=0; j<stack_count; j++){printf(" %d", stack[j]);}

		for (j = 0; j < g->nodes[pn_to_nodes[pnid]].alignedto_count; j++){
			node node_ = g->nodes[g->nodes[pn_to_nodes[pnid]].alignedto[j]];
			if (*stack_max < stack_count + node_.outdegree) { 
				*stack_max +=  g->node_count; 
				stack = realloc(stack, *stack_max * sizeof(uint16_t));
			}
			for (k = 0; k < node_.outdegree; k++){
				stack[stack_count ++] = node_to_pn[g->edges[node_.outedge[k]].outnode_index];
			}
		}
		// printf("\n## %d %d", stack_count, g->node_count);
		// for (j=0; j<stack_count; j++){printf(" %d", stack[j]);}
	}
}

static void toposort(graph *g){
	int32_t node_to_pn[g->node_count];
	memset(node_to_pn, -1, sizeof(int32_t) * g->node_count);
	uint16_t pn_to_nodes[g->node_count];

	uint16_t i = 0;
	uint16_t cur_pnid = 0;
	for (; i < g->node_count; i++){
		if (node_to_pn[i] == -1) {
			pn_to_nodes[cur_pnid] = i;
			node_to_pn[i] = cur_pnid;
			for (uint16_t j = 0; j < g->nodes[i].alignedto_count; j++){
				node_to_pn[g->nodes[i].alignedto[j]] = cur_pnid;
			}
			cur_pnid ++;
		}
	}

	int8_t completed[cur_pnid];
	memset(completed, -1, sizeof(int8_t) * cur_pnid);

	uint32_t stack_max = g->node_count * 2;
	uint16_t *stack = malloc(stack_max * sizeof(uint16_t));

	int32_t found;
	g->sorted_nodes_index = g->node_count - 1;
	while (g->sorted_nodes_index >= 0){
		found = -1;
		for (i = 0; i < cur_pnid; i++){
			if (completed[i] == -1 && check_nodes_predecessors(g, pn_to_nodes[i]) == 0){
				found = i;
				break;
			}
		}
		// printf("%d %d\n", found, g->sorted_nodes_index);
		assert (found != -1);
		sort(found, cur_pnid, pn_to_nodes, node_to_pn, completed, g, &stack_max, stack);
	}
	free (stack);
}

static void align_seq_to_graph_nw(score **s, const uint16_t x, const uint16_t y, const size_t seq_index, const char *seq, graph *g, uint16_t *sorted_nodes_index){

	align_seq_to_graph_updatescore(s, x, y, seq, g, sorted_nodes_index);
	// for (int i=0; i< x + 1; i++){
	// 	printf("\n%d ", i);
	// 	for (int k=0; k<y + 1;k++){
	// 		printf("%ld ",s[i][k].s );
	// 	}
	// }

	uint16_t bestx = align_seq_to_graph_getbestx(y, s, g);
	uint16_t besty = y;
	
	matchroutes mroute;
	mroute.starty = mroute.endy = -1;// here error
	mroute.mroute_count = 0;
	mroute.matchroute_ = malloc((x + y + 1) * sizeof (matchroute));
	memset(mroute.matchroute_, -1, (x + y + 1) * sizeof (matchroute));

	align_seq_to_graph_updatemroute(bestx, besty, s, g->sorted_nodes, &mroute);
	align_seq_to_graph_updategraphy(y, seq_index, seq, g, &mroute);
	toposort(g);
	free (mroute.matchroute_);
}

static inline int32_t get_nextnode(const uint16_t seq_index, const int32_t nodeid, const graph *g){
	int32_t next_nodeid = -1;
	for (uint16_t i = 0; i < g->nodes[nodeid].outdegree; i++){
		if (g->edges[g->nodes[nodeid].outedge[i]].lable[seq_index] == 1) {
			next_nodeid = g->edges[g->nodes[nodeid].outedge[i]].outnode_index;
			// printf("next_nodeid:%d seq_index:%d\n", next_nodeid, seq_index);
			break;
		}
	}
	return next_nodeid;
}

static inline int get_edge_lables(const uint8_t * lables, const int seq_count){
	int count = 0;
	for (int i = 0; i < seq_count; i++){
		count += lables[i];
	}
	return count;
}

static char *get_consensus_from_graph(const graph *g, const int seq_count){
	edge edge_;
	uint16_t i, nodeid, node_index;
	int32_t best_pnode, global_best_node;
	double score, best_score, global_best_score;
	global_best_node = best_score = global_best_score = -1;
	for (node_index = 0; node_index < g->node_count; node_index++){
		nodeid = g->sorted_nodes[node_index];
		best_pnode = -1;
		if (g->nodes[nodeid].indegree){
			for (i = 0; i < g->nodes[nodeid].indegree; i ++){
				edge_ = g->edges[g->nodes[nodeid].inedge[i]];
				score = g->nodes[edge_.innode_index].best_score + get_edge_lables(edge_.lable, seq_count) - 0.5 * g->nodes[nodeid].indegree;
				if (score > best_score || best_pnode == -1) {
					best_score = score;
					best_pnode = edge_.innode_index;
				}
			}
		}else{
			best_score = 0;
			best_pnode = -1;
		}
		g->nodes[nodeid].best_score = best_score;
		g->nodes[nodeid].best_pnode = best_pnode;
		if (best_score > global_best_score){
			global_best_score = best_score;
			global_best_node = nodeid;
		}

	}
	
	uint16_t corrected_seq_index = 0;
	char *corrected_seq = malloc((g->node_count + 1) * sizeof(char));
	while (global_best_node != -1){
		corrected_seq[corrected_seq_index++] = g->nodes[global_best_node].base;
		global_best_node = g->nodes[global_best_node].best_pnode;
	}
	corrected_seq[corrected_seq_index] = '\0';
	reverse_str(corrected_seq, corrected_seq_index);
	return corrected_seq;
}


// static char *get_alignstrings_from_graph(const graph *g, const int seq_count){
// 	int32_t found_idx_, nodeid, column_index[g->node_count];
// 	memset(column_index, -1, sizeof(int32_t) * g->node_count);
// 	uint16_t i, found_idx, node_index, current_column;
// 	current_column = node_index = 0;
// 	for (; node_index < g->node_count; node_index++){
// 		found_idx = g->node_count;
// 		nodeid = g->sorted_nodes[node_index];
// 		// printf("1## %d %d %d %d %d %d\n",nodeid, column_index[46], node_index, current_column, found_idx, g->node_count);
// 		for (i = 0; i < g->nodes[nodeid].alignedto_count; i++){
// 			found_idx_ = column_index[g->nodes[nodeid].alignedto[i]];
// 			if (found_idx_ != -1 && found_idx > found_idx_) found_idx = found_idx_;
// 			// printf("%d %d %d\n", i, found_idx, column_index[g->nodes[nodeid].alignedto[i]]);
// 		}
// 		// printf("2## %d %d %d %d %d %d\n",nodeid, column_index[46], node_index, current_column, found_idx, g->node_count);
// 		if (found_idx == g->node_count) found_idx = current_column ++;
// 		column_index[nodeid] = found_idx;
// 	}


// 	uint8_t (*consensus_matrix) [4] = calloc(current_column, sizeof (uint8_t[4]));
// 	char charlist[current_column + 1];
// 	charlist[current_column] = '\0';
// 	for (i = 0; i < seq_count; i++){
// 		memset(charlist, '-', sizeof(char) * current_column);
// 		nodeid = g->start_node[i];
// 		// printf("%d %d\n", nodeid, i);
// 		while (nodeid != -1){
// 			// printf("%d %d %d %d\n", current_column, column_index[nodeid] ,nodeid, g->node_count);
// 			charlist[column_index[nodeid]] = g->nodes[nodeid].base;
// 			consensus_matrix[column_index[nodeid]][base_to_int[g->nodes[nodeid].base]] ++;
// 			nodeid = get_nextnode(i, nodeid, g);
// 		}
// 		printf("%09d\t%s\n", i, charlist);
// 	}

// 	int j, totalbase, bestcount;
// 	char bestbase = 0;
// 	uint16_t corrected_seqs_index = 0;
// 	char *corrected_seq = malloc((current_column + 1) * sizeof(char));
// 	for (i = 0; i < current_column; i++){
// 		totalbase = seq_count;
// 		bestcount = 0;
// 		for (j = 0; j < 4; j++){
// 			totalbase -= consensus_matrix[i][j];
// 			if (consensus_matrix[i][j] > bestcount){
// 				bestcount = consensus_matrix[i][j];
// 				bestbase = int_to_base[j];
// 			}
// 		}
// 		// if (totalbase > bestcount) continue;
// 		if (totalbase > bestcount) bestbase = '-';
// 		corrected_seq[corrected_seqs_index++] = bestbase;
// 	}

// 	free (consensus_matrix);
// 	corrected_seq[corrected_seqs_index] = '\0';
// 	return corrected_seq;
// }

char *poa_to_consensus(const struct seq_ *seqs, const int seq_count){

	graph poa;
	graph_init(&poa);

	uint16_t x, y;
	size_t seq_index;
	size_t seq_len;
	assert (seq_count <= SEQ_MAX_COUNT);
	for (seq_index = 0; seq_index < seq_count; seq_index++)	{
		seq_len = seqs[seq_index].len;
		// assert (seq_len < 65536);
		if (seq_index == 0){
			int32_t firstnode, headnode;
			firstnode = headnode = -1;
			insert_unmatched_nodes(seq_index, seqs[seq_index].seq, seq_len, &poa, &firstnode, &headnode);
			for (x = 0; x < poa.node_count; x ++) poa.sorted_nodes[x] = x;
			poa.start_node[seq_index] = firstnode;
		}else{
			x = poa.node_count;
			y = seq_len;
			uint16_t sorted_nodes_index[poa.node_count];
			score **s = score_init(x + 1, y + 1, &poa, sorted_nodes_index);
			align_seq_to_graph_nw(s, x, y, seq_index, seqs[seq_index].seq, &poa, sorted_nodes_index);
			free(s);
		}
	}

	// char *corrected_seq1 = get_alignstrings_from_graph(&poa, seq_count);
	// printf("consensus1\t%s\n", corrected_seq1);
	// free (corrected_seq1);

	char *corrected_seq = get_consensus_from_graph(&poa, seq_count);
	// printf("consensus2\t%s\n",corrected_seq);
	graph_destroy(&poa);
	return corrected_seq;
}

// int main(int argc, char** argv) {
	
// 	int seq_count = 0;
// 	struct seq_ * *seqs = malloc(SEQ_MAX_COUNT * sizeof(struct seq_));

// 	gzFile fp = gzopen(argv[1], "r");
// 	kseq_t *seq = kseq_init(fp);
// 	while (kseq_read(seq) >= 0) {
// 		assert (seq->seq.l < lqseq_max_length);
// 		strcpy(seqs[seq_count].seqs, seq->seq.s);
// 		seqs[seq_count].len = seq->seq.l;
// 		seq_count += 1;
// 		assert (seq_count <= SEQ_MAX_COUNT);
// 	}
// 	if (seq_count > 2) {
// 		char *corrected_seq = poa_to_consensus((const lqseq*) seqs, seq_count);
// 		printf("consensus\t%s\n", corrected_seq);
// 		free (corrected_seq);
// 	}
// 	kseq_destroy(seq);
// 	gzclose(fp);
// 	free (seqs);	
// 	return (0);
// }
