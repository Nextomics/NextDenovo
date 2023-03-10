#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include "opt.h"
#include "kit.h"
#include "asg.h"
#include "ctg.h"
#include "out.h"

extern struct option long_options[];

static void init_os(khash_t(ovlinfo_) *os, khash_t(ovlh_) *os_aln, det *d, opt *p){
	UINTL_T n;
	khiter_t i, k;
	int absent;
	aln aln;
	ovlinfo *oli;
	size_t ols = sizeof(ovlinfo);
	ovlinfo_aln *oli_aln;

	int depth_i = 0, *depth = malloc(RANDOM_COUNT * sizeof(int));

	for (k = kh_begin(os_aln); k != kh_end(os_aln); ++k){
		if (!kh_exist(os_aln, k)) continue;
		oli_aln = &kh_val(os_aln, k);
		if (oli_aln->con >= p->min_con_count) {
			if (oli_aln->con < MAX_CON) free(oli_aln->alns);
			continue;
		}
		n = kh_key(os_aln, k);
		absent = find_alnse(oli_aln, &aln);
		if (p->chimera && absent > 1){
			free (oli_aln->alns);
			continue;
		}
		if (p->refilt_con_read && oli_aln->alnl.s < aln.s + p->max_hang_len \
				&& oli_aln->alnl.e > aln.e - p->max_hang_len) {
			if (++oli_aln->con >= p->min_con_count) {
				free (oli_aln->alns);
				continue;
			}
		}

		i = kh_put(ovlinfo_, os, n, &absent);
		oli = &kh_val(os, i);
		memset(oli, 0, ols);
		kh_key(os, i) = n;

		// if (aln.s < UINT16_MAX) oli->le = aln.s;
		// else oli->le = UINT16_MAX;
		// if (oli_aln->len - aln.e < UINT16_MAX) oli->re = oli_aln->len - aln.e;
		// else oli->re = UINT16_MAX;
		oli->le = aln.s;
		oli->re = oli_aln->len - aln.e;
		oli->lc = oli_aln->lc;
		oli->rc = oli_aln->rc;
		oli->con = oli_aln->con;
		oli->lim = oli_aln->lim;
		oli->rim = oli_aln->rim;
		oli->llm = oli_aln->llm;
		oli->rlm = oli_aln->rlm;

		if (depth_i < RANDOM_COUNT - 2){
			depth[depth_i++] = oli->lc;
			depth[depth_i++] = oli->rc;
		}
		if (oli->lc < p->min_aln_depth) d->l[0] ++;
		else if (oli->lc < p->max_aln_depth) d->n[0] ++;
		else d->h[0] ++;

		if (oli->rc < p->min_aln_depth) d->l[0] ++;
		else if (oli->rc < p->max_aln_depth) d->n[0] ++;
		else d->h[0] ++;
		free (oli_aln->alns);
	}
	p->median_aln_depth = quick_select(depth, 0, depth_i - 1, depth_i/2);
	free (depth);
}

void stat_graph(graph *g, char *fun, char *des){
	plog(2, "[DEBUG:FUN:%-26s:DES:%-33s] nodes: %"UINTL_FORMAT" edges: %"UINTL_FORMAT, 
		fun, des, stat_valid_node(g), stat_valid_edge(g));
}

int main(int argc, char *argv[])
{	
	double inittime = realtime();
	int option_index, c;
	opt opts;
	init_opt(&opts);
	while((c = getopt_long_only(argc, argv, "f:F:T:i:I:R:S:r:m:n:B:e:C:z:l:L:@cAdkGt:o:N:E:q:a:D:P:u:w:M:", \
		long_options, &option_index)) != -1) {
		switch(c) {
			case 'f': opts.fa = fopen(optarg, "r"); break;
			// case 'b': opts.bl = fopen(optarg, "r"); break;
			case 'F': opts.fuzz_len = atoi(optarg); break;
			case 'i': opts.min_ide = atof(optarg) * 10000; break;
			case 'I': opts.min_ide_ratio = atof(optarg); break;
			case 'R': opts.max_ide_ratio = atof(optarg); break;
			case 'S': opts.min_sco_ratio = atof(optarg); break;
			case 'r': opts.max_sco_ratio = atof(optarg); break;
			case 'm': opts.min_depth_multi = atof(optarg); break;
			case 'n': opts.max_depth_multi = atof(optarg); break;
			case 'B': opts.bubble_len = atoi(optarg); break;
			case 'e': opts.end_loop_len = atoi(optarg); break;
			case 'C': opts.cpath_len = atoi(optarg); break;
			// case 'g': opts.max_gap_len = atoi(optarg); break;
			case 'z': opts.zbranch_len = atoi(optarg); break;
			case 'l': opts.sbranch_len = atoi(optarg); break;
			case 'L': opts.sloop_len = atoi(optarg); break;
			// case 'p': opts.max_hang_plen = atoi(optarg); break;
			case 'N': opts.min_node_count = atoi(optarg); break;
			case 'u': opts.min_con_count = atoi(optarg); break;
			case 'c': opts.chimera = 0; break;
			// case 's': opts.sort = 0; break;
			case 'a': opts.out_format = atoi(optarg); break;
			// case 'R': opts.refilt_con_read = 0; break;
			case 'k': opts.keep_comp_path = 0; break;
			case 'G': opts.keep_chimera_edge = 1; break;
			case 'd': opts.debug = 1; break;
			case 't': opts.max_hang_len = atoi(optarg); break;
			case 'E': opts.out_ctg_len = atoi(optarg); break;
			case 'q': opts.out_spath_len = atoi(optarg); break;
			// case 'Q': opts.out_spath_len = opts.out_zpath_len = atoi(optarg); break;
			case 'D': opts.bfs_depth = atoi(optarg); break;
			case 'P': opts.bfs_depth_multi = atoi(optarg); break;
			case 'w': opts.min_edge_cov = atoi(optarg); break;
			case 'A': opts.out_alt_ctg = 1; break;
			case 'M': opts.min_mat_ratio = atof(optarg); break;
			case 'T': opts.min_depth_ratio = atof(optarg); break;
			case '@': return version();break;
			case 'o':
				if (freopen(optarg, "wb", stdout) == NULL){
					fprintf(stderr, "Failed to write output to file %s\n", optarg);
					exit(1);
				};
				break;
			default: return usage(&opts);
		}
	}

	if (optind + 1 > argc) return usage(&opts);
	if (opts.fa == NULL) {
		fprintf(stderr, "Failed open input seq list!\n");
		return usage(&opts);
	}
	FILE *fp = fopen(argv[optind], "r");
	if (fp == NULL) {
		fprintf(stderr, "Failed open input file %s!\n", argv[optind]);
		return usage(&opts);
	}
	if (opts.bubble_len == 500 && opts.out_alt_ctg) opts.bubble_len = 40;
	
	plog(2,"Initialize graph and reading...");
	char *file = NULL;
	size_t len = 0;
	ssize_t nread;
	list l;
	init_list(&l);
	khash_t(ovlh_) *os_aln = kh_init(ovlh_);
	while ((nread = getline(&file, &len, fp)) != -1) {
		if (nread != 1 && file[0] != '#'){
			file[nread - 1 ] = '\0';
			char file_[1024];
			sprintf(file_, "%s.bl", file);
			read_bl(os_aln, file_, &l);
		}
	}
	det dt = {{0}, {0}, {0}};
	khash_t(ovlinfo_) *os = kh_init(ovlinfo_);
	kh_resize(ovlinfo_, os, os_aln->n_buckets);
	init_os(os, os_aln, &dt, &opts);
	kh_destroy(ovlh_, os_aln);

	graph g;
	init_graph(&g, MAX_INIT_N);
	rewind(fp);
	uint32_t *arr = init_decode_table(NUM_LIMIT);
	while ((nread = getline(&file, &len, fp)) != -1) {
		if (nread != 1 && file[0] != '#'){
			file[nread - 1 ] = '\0';
			update_graph(file, &g, os, &opts, arr, &l);
		}
	}
	fclose(fp);
	destroy_table(arr);
	if (file) free(file);
	plog(2, "Initial Node(s): %"UINTL_FORMAT", Edge(s): %"UINTL_FORMAT, g.ni/2, g.ei/2);

	rm_node_con(&g, os, &opts, &l);//1
	rm_edge_lq(&g, os, &opts, &l);//2
	if (opts.debug) stat_graph(&g, "rm_edge_lq", "remove low quality edges");
	// out_graph_graphml(&g);exit(1);
	sort_stat_oe(&g, &dt, &opts);//3
	if (opts.debug) stat_graph(&g, "sort_stat_oe", "sort_stat_oe");
	mark_edge_rep(&g, os, &opts, &dt, &l);//4
	if (opts.debug) stat_graph(&g,"mark_edge_rep", "remove repeat edges");
	kh_destroy(ovlinfo_, os);//5
	destroy_list(&l);
	mark_node_chim(&g, opts.bfs_depth, opts.bfs_depth_multi, opts.median_outdegree, &opts.keep_chimera_edge);
	mark_edge_tr(&g, opts.fuzz_len);//6
	if (opts.debug) stat_graph(&g,"mark_edge_tr", "remove transitive edges");
	rm_edge_spur(&g);//7
	
	if (opts.max_ide_ratio != .0){
		mark_edge_hli(&g, opts.max_ide_ratio);
		rm_edge_li(&g);
	}

	mark_edge_hls(&g, opts.max_sco_ratio);//8
	rm_edge_ls(&g);//9
	mark_edge_bs(&g);//10
	if (opts.debug) stat_graph(&g,"mark_edge_bs", "remove low socre edges");
	rm_sht_brh(&g, opts.sbranch_len, opts.out_spath_len);//11
	if (opts.debug) stat_graph(&g, "rm_sht_brh1", "remove short branch");
	rm_z_clip_lable(&g, opts.zbranch_len, opts.out_spath_len, MFLAG_CN);//12
	if (opts.debug) stat_graph(&g, "rm_z_clip_lable1", "remove z branch by MFLAG_CN");
	rm_z_clip_lable(&g, opts.zbranch_len, opts.out_spath_len, MFLAG_CN);//12
	if (opts.debug) stat_graph(&g, "rm_z_clip_lable2", "remove z branch by MFLAG_CN");

	#ifdef P_VERSION
		rm_z_clip_lable(&g, opts.zbranch_len, opts.out_spath_len, MFLAG_LQ);//12
		if (opts.debug) stat_graph(&g, "rm_z_clip_lable3", "remove z branch by MFLAG_LQ");
	#else
		rm_z_clip_lable(&g, opts.zbranch_len, opts.out_spath_len, MFLAG_CC);//12
		if (opts.debug) stat_graph(&g, "rm_z_clip_lable3", "remove z branch by MFLAG_CC");
	#endif
	// out_edges(&g); exit(1);
	// out_graph_graphml(&g);exit(1);
	if (opts.out_alt_ctg){
		cal_node_io_bstsc(&g, 2);
		rm_z_clip_score(&g, opts.zbranch_len, opts.out_spath_len);//12
		rm_z_clip_score(&g, opts.zbranch_len, opts.out_spath_len);//12
	}else{
		clean_complex_single_path(&g, opts.bubble_len, 0, opts.min_mat_ratio);//solve_compound_paths_simple
		if (opts.debug) stat_graph(&g, "clean_complex_single_path", "remove single complex path");
		clean_complex_multi_path(&g, opts.bubble_len, opts.out_spath_len, opts.min_mat_ratio);//solve_compound_paths_common
		if (opts.debug) stat_graph(&g, "clean_complex_multi_path", "remove multi complex paths");
		calc_edge_tc(&g, 255);
		rm_z_clip_score3(&g, opts.zbranch_len, opts.out_spath_len, opts.min_depth_ratio * 100);
		if (opts.debug) stat_graph(&g, "rm_z_clip_score3", "remove z branch by score");
		clean_complex_graph(&g, 8, opts.out_spath_len, opts.min_mat_ratio);//solve_dense_spots_hete
		if (opts.debug) stat_graph(&g, "clean_complex_graph", "remove complex paths");
		rm_z_clip_score2(&g, 2 * opts.zbranch_len, opts.out_spath_len, -1);
		if (opts.debug) stat_graph(&g, "rm_z_clip_score2", "remove z branch by score");
	}

	rm_sht_brh(&g, opts.sbranch_len, opts.out_spath_len);//13
	if (opts.debug) stat_graph(&g, "rm_sht_brh2", "remove short branch");
	rm_sht_loop(&g, opts.sloop_len);//14
	if (opts.debug) stat_graph(&g, "rm_sht_loop", "remove short loop");
	if (opts.out_alt_ctg) rm_sht_bubble(&g, opts.bubble_len);//15
	else rm_long_bubble(&g, opts.bubble_len, opts.out_spath_len, opts.min_mat_ratio);//bubble_bursting2
	if (opts.debug) stat_graph(&g, "rm_bubble", "remove loop");
	rm_end_loop(&g, opts.end_loop_len);//16
	if (opts.debug) stat_graph(&g, "rm_end_loop", "remove end loop");
	clean_complex_path(&g, opts.cpath_len * 1.5, opts.keep_comp_path, 0);//17
	if (opts.debug) stat_graph(&g, "clean_complex_path", "remove all complex paths");
	clean_complex_path(&g, opts.cpath_len, opts.keep_comp_path, MFLAG_CN);//18
	if (opts.debug) stat_graph(&g, "clean_complex_path", "remove complex paths by MFLAG_CN");
	clean_complex_path(&g, opts.cpath_len, opts.keep_comp_path, MFLAG_LS);//19
	if (opts.debug) stat_graph(&g, "clean_complex_path", "remove complex paths by MFLAG_LS");
	rm_sht_brh(&g, opts.sbranch_len, opts.out_spath_len);
	if (opts.debug) stat_graph(&g, "rm_sht_brh3", "remove short branch");
	int bl = opts.keep_chimera_edge ? 10 : 50;
	if (opts.out_alt_ctg) calc_edge_tc(&g, opts.min_edge_cov);

	#ifdef P_VERSION
		rm_edge_chim(&g, bl, 0.5, 0.85 * 10000, MFLAG_LQ);
	#else
		rm_edge_chim(&g, bl, 0.5, 0.85 * 10000, MFLAG_CC);
	#endif
	if (opts.debug) stat_graph(&g, "rm_edge_chim", "remove potential chimeric edges");
	rm_edge_ltc(&g, 20, 0.33, 0.85, opts.min_edge_cov);
	if (opts.debug) stat_graph(&g, "rm_edge_ltc", "remove low sco & ide edges");
	rm_sht_brh(&g, 5, opts.out_spath_len);//TODO check
	if (opts.debug) stat_graph(&g, "rm_sht_brh4", "remove short branch");
	clean_node_lable(&g);
	plog(2,"Assembly done and outputting...");
	ctgs asmb;
	init_ctgs(&asmb, 1000);
	generate_ctg(&g, &asmb, MFLAG_LQ);
	idxs *idx  = NULL;
	uint64_t gs = 0;
	switch (opts.out_format){
		case 1:
			destroy_graph(&g);
			idx = init_index(opts.fa);
			gs = out_ctg_fasta(&asmb, idx, opts.out_ctg_len, opts.out_format);
			break;
		case 2:
			out_graph_graphml(&g);
			destroy_graph(&g);
			break;
		case 3:
			plog(3, "Output with GFA2 format has not been fully tested.");
			idx = init_index(opts.fa);
			gs = out_ctg_gfa(&asmb, &g, idx, opts.out_ctg_len);
			destroy_graph(&g);
			break;
		case 4:
			destroy_graph(&g);
			idx = init_index(opts.fa);
			gs = out_ctg_path(&asmb, idx, opts.out_ctg_len);
			break;
		case 5:
			out_graph_raw(&g);
			destroy_graph(&g);
			break;		
		default:
			destroy_graph(&g);
			idx = init_index(opts.fa);
			gs = out_ctg_fasta(&asmb, idx, opts.out_ctg_len, 0);
	}

	if (idx) destroy_index(idx);
	if (gs) stat_ctg(&asmb, gs);
	fclose (opts.fa);
	destroy_ctgs(&asmb);

	#ifdef RESTRICT
	if (gs > RESTRICT_GS && opts.out_format) plog(3, "Unfinished assembly, this is a limited version, currently" 
		" only supports assembly for genome size < %ld bp, please ask for help.", RESTRICT_GS);
	#endif
	
	plog(2, "CMD:");
	for (c = 0; c < argc; ++c) fprintf(stderr, " %s", argv[c]);
	fprintf(stderr, "\n");
	if (opts.out_format != 1) plog(3, "Please use '-a 1' to output assembled sequences");
	plog(2, "Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", realtime() - inittime, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
	
	return 0;
}
