#ifndef _ASG_H_
#define _ASG_H_
#include "common.h"
#include "node.h"
#include "edge.h"
#include "opt.h"

void init_graph(graph *s, const uint32_t len);
void destroy_graph(graph *s);
void rm_node_con(graph *g, khash_t(ovlinfo_) *os, opt *p, list *l);
void sort_stat_oe(graph *g, det *d, opt *p);
void rm_edge_lq(graph *g, khash_t(ovlinfo_) *os, opt *p, list *l);
void rm_edge_chim(graph * g, const int s, const float sco, const uint32_t ide, const uint16_t l);
void rm_edge_ltc(graph *g, const int s, const float sco, const float idt, const int l);
void calc_edge_tc(graph *g, const int s);
void mark_edge_rep(graph *g, khash_t(ovlinfo_) *os, opt *p, const det *d, list *l);
void mark_edge_tr(graph *g, const int fuzz);
void rm_edge_spur(graph *g);
void mark_edge_hli(graph *g, const float max_ide_ratio);
void rm_edge_li(graph *g);
void mark_edge_hls(graph *g, const float max_sco_ratio);
void mark_node_chim(graph *g, const int s, const int m, const int md, int *kc);
//void mark_node_chim(graph *g);
void rm_edge_ls(graph *g);
void mark_edge_bs(graph *g);
void rm_sht_brh(graph *g, const int s, const int l);
void rm_z_clip_lable(graph *g, const int s, const int m, const uint16_t l);
void cal_node_io_bstsc(graph *g, const int s);
void rm_z_clip_score(graph *g, const int s, const int m);
void rm_sht_loop(graph *g, const int s);
void rm_sht_bubble(graph *g, const int s);
void rm_end_loop(graph *g, const int s);
void clean_complex_path(graph *g, const int s, const int r, const uint16_t l);
void update_graph(char *file, graph *g, khash_t(ovlinfo_) *os, opt *p, uint32_t *arr, list *l);


void clean_complex_single_path(graph *g, const int s, const int m, const float perc);
void clean_complex_multi_path(graph *g, const int s, const int m, const float perc);
void rm_z_clip_score3(graph *g, const int s, const int m, const int perc);
void clean_complex_graph(graph *g, const int s, const int m, const float perc);
void rm_z_clip_score2(graph *g, const int s, const int m, const int perc);
void rm_long_bubble(graph *g, const int s, const int m, const float perc);
#endif
