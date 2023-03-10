#ifndef _OPT_H_
#define _OPT_H_
#include <stdint.h>
#include <stdio.h>
#include <getopt.h>

typedef struct opt_ {
	FILE *fa; // input seq list, f
	// FILE *bl; // input black list, b
	int out_format;//out_seq; // output assembly sequences, a
	int debug; //debug mode, d
	int sort; // sort out edges by len
	int fuzz_len;//fuzz; // fuzz len for trans-reduction, F
	int chimera;//pre filter chimeric reads, c
	int keep_chimera_edge;//retain potential chimeric edges, G
	int bfs_depth;//ext_node_count;// depth of BFS to identify chimeric nodes, D
	int bfs_depth_multi;//ext_depth_multi;// max depth multiple of a node for BFS, P
	int keep_comp_path;// romove all compound path, k
	int out_spath_len;// min short branch len for output, q
	// int out_zpath_len;// min z path len for output, Q
	int refilt_con_read;//disable refilter contained reads, R
	int min_con_count; // min contained count to filter reads, u
	int out_ctg_len;// min contig len, E
	int max_aln_depth; // max aligned depth
	int min_aln_depth; // min aligned depth
	int sbranch_len;//short_branch_len; // max len of a short branch, l
	int zbranch_len; // max len of a z branch, z
	int sloop_len;//short_loop_len; // max len of a loop, L
	int bubble_len; // max len of a bubble, B
	int end_loop_len; // max len of a terminal loop
	int cpath_len;//comp_path_len; // max len of a compound path, C
	// int max_gap_len; // max gap len of two path
	// int max_hang_plen; // max over hang length at two end overlap to identity dovetails
	int max_hang_len;//max_hang_tlen; // max over hang length at two end overlap of dovetails, t
	// int thread;// T
	int min_node_count;//min valid node count of a read, N
	uint32_t min_ide; // min identity of alignments. i
	float min_ide_ratio; // min test-to-best identity ratio, I
	float max_ide_ratio; // max test-to-best identity ratio of a low quality edge, R
	float min_sco_ratio; // min test-to-best aligned length ratio, S
	float min_depth_multi; // min depth multiple of a repeat node, m
	float max_depth_multi; // max depth multiple of a node, n
	float max_sco_ratio; //min test-to-best score ratio of a node, r
	int min_edge_cov; //min depth of an edge, w
	int out_alt_ctg; // output alternative contigs, A
	float min_mat_ratio; //min test-to-best matches ratio, M
	float min_depth_ratio;//min_edepth_ratio; //min test-to-best depth ratio of an edge, T
	int median_aln_depth;
	int median_outdegree;
} opt;

void init_opt(opt *p);
int usage(opt *opts);

#endif
