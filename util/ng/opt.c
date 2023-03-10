#include "opt.h"

struct option long_options[] = {
	{"fa",               required_argument,  NULL,	'f'},
	// {"bl",	             required_argument,  NULL,	'b'},
	{"output",	         required_argument,  NULL,	'o'},
	{"out_format",	     required_argument,  NULL,	'a'},
	// {"sort",	         no_argument,        NULL,	's'},
	{"fuzz_len",	     required_argument,	 NULL,	'F'},
	// {"thread",	         required_argument,	 NULL,	'T'},
	{"chimera",	         no_argument,        NULL,	'c'},
	{"keep_comp_path",	 no_argument,        NULL,	'k'},
	// {"refilt_con_read",	 no_argument,        NULL,	'R'},
	{"keep_chimera_edge",no_argument,        NULL,	'G'},
	{"out_alt_ctg",      no_argument,        NULL,	'A'},

	{"min_ide",	         required_argument,	 NULL,	'i'},
	{"min_ide_ratio",	 required_argument,	 NULL,	'I'},
	{"max_ide_ratio",	 required_argument,	 NULL,	'R'},
	{"min_sco_ratio",	 required_argument,	 NULL,	'S'},
	{"max_sco_ratio",	 required_argument,	 NULL,	'r'},
	// {"max_aln_depth",	 required_argument,	 NULL,	'd'},
	{"min_depth_multi",	 required_argument,	 NULL,	'm'},
	{"max_depth_multi",	 required_argument,	 NULL,	'n'},
	{"bubble_len",       required_argument,	 NULL,	'B'},
	// {"end_loop_len",     required_argument,	 NULL,	'e'},
	{"cpath_len",        required_argument,	 NULL,	'C'},
	// {"max_gap_len",      required_argument,	 NULL,	'g'},
	{"zbranch_len",      required_argument,	 NULL,	'z'},
	{"sbranch_len",      required_argument,	 NULL,	'l'},
	{"sloop_len",        required_argument,	 NULL,	'L'},
	// {"max_hang_plen",    required_argument,	 NULL,	'p'},
	{"max_hang_len",     required_argument,	 NULL,	't'},
	{"min_node_count",   required_argument,	 NULL,	'N'},
	{"out_ctg_len",      required_argument,	 NULL,	'E'},
	{"out_spath_len",    required_argument,	 NULL,	'q'},
	// {"out_zpath_len",    required_argument,	 NULL,	'Q'},
	{"bfs_depth",        required_argument,	 NULL,	'D'},
	{"bfs_depth_multi",  required_argument,	 NULL,	'P'},
	{"min_con_count",    required_argument,  NULL,	'u'},
	{"min_edge_cov",     required_argument,  NULL,  'w'},
	{"min_mat_ratio",    required_argument,  NULL,  'M'},
	{"min_depth_ratio",  required_argument,  NULL,  'T'},

	{0,                  0,                  NULL,  0}
};

void init_opt(opt *p){
	p->fa = NULL;
	// p->bl = NULL;
	p->sort = 1;
	p->debug = 0;
	p->out_format = 1;
	p->keep_comp_path = 1;
	p->out_spath_len = 0;
	// p->out_zpath_len = 0;
	p->fuzz_len = 1000;
	// p->thread = 8;
	p->chimera = 1;
	p->refilt_con_read = 1;
	p->min_ide = 10;
	p->min_ide_ratio = 0.7;
	p->max_ide_ratio = 0;
	p->min_sco_ratio = 0.4;
	p->max_sco_ratio = 0.5;
	p->max_aln_depth = 500;
	p->min_aln_depth = 3;
	p->min_depth_multi = 1.5;
	p->max_depth_multi = 2000;
	p->bubble_len = 500;
	p->end_loop_len = 50;
	p->cpath_len = 20;
	// p->max_gap_len = 20;
	p->zbranch_len = 8;
	p->sbranch_len = 15;
	p->sloop_len = 5;
	// p->max_hang_plen = 3000;
	p->max_hang_len = 500;
	p->min_node_count = 2;
	p->out_ctg_len = 1000;
	p->bfs_depth = 2;
	p->bfs_depth_multi = 2;
	p->min_con_count = 2;
	p->keep_chimera_edge = 0;
	p->min_edge_cov = 3;
	p->out_alt_ctg = 0;
	p->min_mat_ratio = 0.9;
	p->min_depth_ratio = 0.6;
}

int usage(opt *opts)
{
	fprintf(stderr, "Usage: nextGraph [options] -f seq.fofn ovl.fofn\n\n");
	fprintf(stderr, "Options:\033[35m<TODO:UPDATE MORE DETAILS>\033[0m\n");
	fprintf(stderr, "  -f [FILE]                       input seq list [required]\n");
	// fprintf(stderr, "  -b [FILE]                     input black list [required]\n");
	fprintf(stderr, "  -o [FILE]                       output file [stdout]\n");
	// fprintf(stderr, "  -s                              disable sort out-edges by length \n");
	fprintf(stderr, "  -c                              disable pre-filter chimeric reads \n");
	// fprintf(stderr, "  -R                              disable re-filter contained reads \n");
	fprintf(stderr, "  -G                              retain potential chimeric edges \n");
	fprintf(stderr, "  -k                              delete complex bubble paths \n");
	fprintf(stderr, "  -A                              output alternative contigs, will increase assembly size \n");
	fprintf(stderr, "  -a --out_format [INT]           output format, 0=None, 1=fasta, 2=graphml, 3=gfa2, 4=path [%d]\n", opts->out_format);
	fprintf(stderr, "  -E --out_ctg_len [INT]          min contig length for output [%d]\n", opts->out_ctg_len);
	fprintf(stderr, "  -q --out_spath_len [INT]        min short branch len (0~16) for output, non-zero value will increase assembly size, 0=disable [%d]\n", opts->out_spath_len);
	// fprintf(stderr, "  -Q --out_zpath_len [INT]        min z branch length for output, 0=disable [%d]\n", opts->out_zpath_len);
	// fprintf(stderr, "  -T --thread [INT]             number of threads [%d]\n", opts->thread);
	fprintf(stderr, "  -i --min_ide [FLOAT]            min identity of alignments [%.2f]\n", (float) opts->min_ide/100);
	fprintf(stderr, "  -I --min_ide_ratio [FLOAT]      min test-to-best identity ratio [%.2f]\n", opts->min_ide_ratio);
	fprintf(stderr, "  -R --max_ide_ratio [FLOAT]      min test-to-best identity ratio of a low quality edge [%.2f]\n", opts->max_ide_ratio);
	fprintf(stderr, "  -S --min_sco_ratio [FLOAT]      min test-to-best aligned length ratio [%.2f]\n", opts->min_sco_ratio);
	fprintf(stderr, "  -r --max_sco_ratio [FLOAT]      max test-to-best score ratio of a low quality edge [%.2f]\n", opts->max_sco_ratio);
	fprintf(stderr, "  -M --min_mat_ratio [FLOAT]      min test-to-best aligned matches ratio [%.2f]\n", opts->min_mat_ratio);
	fprintf(stderr, "  -T --min_depth_ratio [FLOAT]    min test-to-best depth ratio of an edge [%.2f]\n", opts->min_depth_ratio);
	fprintf(stderr, "  -N --min_node_count [1,2]       min valid nodes of a read [%d]\n", opts->min_node_count);
	fprintf(stderr, "  -u --min_con_count [1,2]        min contained number to filter contained reads [%d]\n", opts->min_con_count);
	fprintf(stderr, "  -w --min_edge_cov [INT]         min depth of an edge [%d]\n", opts->min_edge_cov);
	//fprintf(stderr, "  -d --max_aln_depth [INT]        max aligned depth [%d]\n", opts->max_aln_depth);
	fprintf(stderr, "  -D --bfs_depth [INT]            depth of BFS to identify chimeric nodes [%d]\n", opts->bfs_depth);
	fprintf(stderr, "  -P --bfs_depth_multi [INT]      max depth multiple of a node for BFS [%d]\n", opts->bfs_depth_multi);
	fprintf(stderr, "  -m --min_depth_multi [FLOAT]    min depth multiple of a repeat node [%.2f]\n", opts->min_depth_multi);
	fprintf(stderr, "  -n --max_depth_multi [FLOAT]    max depth multiple of a node [%.2f]\n", opts->max_depth_multi);
	fprintf(stderr, "  -B --bubble_len [INT]           max len of a bubble [%d]\n", opts->bubble_len);
	// fprintf(stderr, "  -e --end_loop_len [INT]         max len of a terminal loop [%d]\n", opts->end_loop_len);
	fprintf(stderr, "  -C --cpath_len [INT]            max len of a compound path [%d]\n", opts->cpath_len);
	// fprintf(stderr, "  -g --max_gap_len [INT]          max gap len of two pathes [%d]\n", opts->max_gap_len);
	fprintf(stderr, "  -z --zbranch_len [INT]          max len of a z branch [%d]\n", opts->zbranch_len);
	fprintf(stderr, "  -l --sbranch_len [INT]          max len of a short branch [%d]\n", opts->sbranch_len);
	fprintf(stderr, "  -L --sloop_len [INT]            max len of a short loop [%d]\n", opts->sloop_len);
	// fprintf(stderr, "  -p --max_hang_plen [INT]        max over hang length of potential dovetails [%d]\n",opts->max_hang_plen);
	fprintf(stderr, "  -t --max_hang_len [INT]         max over hang length of dovetails [%d]\n", \
		opts->max_hang_len);
	fprintf(stderr, "  -F --fuzz_len [INT]             fuzz len for trans-reduction [%d]\n", opts->fuzz_len);
	return 1;
}
