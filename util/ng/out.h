#ifndef _OUT_H_
#define _OUT_H_
#include "common.h"
#include "ctg.h"
#include "../../lib/index.h"

void stat_ctg(ctgs *m, const uint64_t gs);
void out_node(graph *s, const UINTL_T j);
uint64_t out_ctg_path(ctgs *m, const idxs *idx, const int min_ctg_len);
uint64_t out_ctg_gfa(ctgs *m, graph *g, const idxs *idx, const int min_ctg_len);//TODO need more test
uint64_t out_ctg_fasta(ctgs *m, const idxs *idx, const int min_ctg_len, const int out_seq);
void out_graph_raw(graph *s);
void out_graph_graphml(graph *s);
void out_edges(graph *g);

#endif
