#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "ketopt.h"
#include <libgen.h>
#include "../lib/ovl.h"

#define MM_VERSION "2.17-r941"

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

// char **bases_arr = NULL;
bytes_t *encode_tbl = NULL;
prev_t pid = {0, 0};

static ko_longopt_t long_options[] = {
	{ "bucket-bits",    ko_required_argument, 300 },
	{ "mb-size",        ko_required_argument, 'K' },
	{ "seed",           ko_required_argument, 302 },
	{ "no-kalloc",      ko_no_argument,       303 },
	{ "print-qname",    ko_no_argument,       304 },
	{ "no-self",        ko_no_argument,       'D' },
	{ "print-seeds",    ko_no_argument,       306 },
	{ "max-chain-skip", ko_required_argument, 307 },
	{ "min-dp-len",     ko_required_argument, 308 },
	{ "print-aln-seq",  ko_no_argument,       309 },
	{ "splice",         ko_no_argument,       310 },
	{ "cost-non-gt-ag", ko_required_argument, 'C' },
	{ "no-long-join",   ko_no_argument,       312 },
	{ "sr",             ko_no_argument,       313 },
	{ "frag",           ko_required_argument, 314 },
	{ "secondary",      ko_required_argument, 315 },
	{ "cs",             ko_optional_argument, 316 },
	{ "end-bonus",      ko_required_argument, 317 },
	{ "no-pairing",     ko_no_argument,       318 },
	{ "splice-flank",   ko_required_argument, 319 },
	{ "idx-no-seq",     ko_no_argument,       320 },
	{ "end-seed-pen",   ko_required_argument, 321 },
	{ "for-only",       ko_no_argument,       322 },
	{ "rev-only",       ko_no_argument,       323 },
	{ "heap-sort",      ko_required_argument, 324 },
	{ "all-chain",      ko_no_argument,       'P' },
	{ "dual",           ko_required_argument, 326 },
	{ "max-clip-ratio", ko_required_argument, 327 },
	{ "min-occ-floor",  ko_required_argument, 328 },
	{ "MD",             ko_no_argument,       329 },
	{ "lj-min-ratio",   ko_required_argument, 330 },
	{ "score-N",        ko_required_argument, 331 },
	{ "eqx",            ko_no_argument,       332 },
	{ "paf-no-hit",     ko_no_argument,       333 },
	{ "split-prefix",   ko_required_argument, 334 },
	{ "no-end-flt",     ko_no_argument,       335 },
	{ "hard-mask-level",ko_no_argument,       336 },
	{ "cap-sw-mem",     ko_required_argument, 337 },
	{ "max-qlen",       ko_required_argument, 338 },
	{ "max-chain-iter", ko_required_argument, 339 },
	{ "junc-bed",       ko_required_argument, 340 },
	{ "junc-bonus",     ko_required_argument, 341 },
	{ "sam-hit-only",   ko_no_argument,       342 },
	{ "help",           ko_no_argument,       'h' },
	{ "max-intron-len", ko_required_argument, 'G' },
	{ "version",        ko_no_argument,       'V' },
	{ "min-count",      ko_required_argument, 'n' },
	{ "min-chain-score",ko_required_argument, 'm' },
	{ "mask-level",     ko_required_argument, 'M' },
	{ "min-dp-score",   ko_required_argument, 's' },
	{ "sam",            ko_no_argument,       'a' },
//  NextDenovo options
	{"step",           ko_required_argument,  400 },
	{"minlen",         ko_required_argument,  401 },
	{"minide",         ko_required_argument,  402 },
	{"maxhan1",        ko_required_argument,  403 },
	{"maxhan2",        ko_required_argument,  404 },
	{"outraw",         ko_no_argument,        405 },
	{"minmatch",       ko_required_argument,  406 },
	{"mode",           ko_required_argument,  407 },
	{"kn",             ko_required_argument,  408 },
	{"wn",             ko_required_argument,  409 },
	{"cn",             ko_required_argument,  410 },
	{"outctn",         ko_no_argument,        411 },
	{"df",             ko_required_argument,  412 },
	{"dvt",            ko_no_argument,        413 },

	{ 0, 0, 0 }
};

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

static inline void yes_or_no(mm_mapopt_t *opt, int flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

void get_idx_offset(const char *fn, idx_t *ptr)
{
	ptr->is_idx = 0;
	ptr->name = strdup(fn);
	ptr->h = kh_init(idx_offset);

	int fn_len = strlen(fn);
	char *idx_name = (char *)calloc(fn_len + 1, sizeof(char));;
	char *file_name = basename((char *)fn);
	int fname_len = strlen(file_name);
	memcpy(idx_name, fn, fn_len - fname_len);
	memcpy(idx_name + fn_len - fname_len, ".", 1); 
	memcpy(idx_name + fn_len - fname_len + 1, file_name, fname_len - 5); 
	memcpy(idx_name + fn_len - fname_len + 1 + fname_len - 5, ".idx", 4); 

	FILE *stream;
	char *line = NULL;
	size_t len = 0;
	stream = fopen(idx_name, "r");
	if (stream == NULL) {
		fprintf(stderr, "Failed open idx file\n");
		exit(1);
	}
	uint32_t i, offset, length;
	int absent; khint_t k;
	while(getline(&line, &len, stream) != -1) {
		sscanf(line, "%u\t%u\t%u", &i, &offset, &length);
		k = kh_put(idx_offset, ptr->h, i, &absent);
		kh_val(ptr->h, k) = offset;
	}
	free(line);
	free(idx_name); fclose(stream);
}

int main(int argc, char *argv[])
{	
	const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:hF:LC:yYPo:";
	ketopt_t o = KETOPT_INIT;
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	int i, c, n_threads = 3, n_parts, old_best_n = -1;
	char *fnw = 0, *rg = 0, *junc_bed = 0, *s;
	FILE *fp_help = stderr;
	mm_idx_reader_t *idx_rdr;
	mm_idx_t *mi;

	mm_verbose = 3;
	liftrlimit();
	mm_realtime0 = realtime();
	mm_set_opt(0, &ipt, &opt);

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) { // test command line options and apply option -x/preset first
		if (c == 'x') {
			if (mm_set_opt(o.arg, &ipt, &opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", o.arg);
				return 1;
			}
		} else if (c == ':') {
			fprintf(stderr, "[ERROR] missing option argument\n");
			return 1;
		} else if (c == '?') {
			fprintf(stderr, "[ERROR] unknown option in \"%s\"\n", argv[o.i - 1]);
			return 1;
		}else if (c == 400) {
			opt.step = atoi(o.arg);
			opt.outraw = 1;
			if (opt.step == 1){ 
				opt.minlen = 500; opt.outide = 0;
			}else if (opt.step == 2) {
				opt.minide = 0.05; opt.outide = 1; opt.minlen = 2000;
				opt.maxhan1 = 5000; opt.maxhan2 = 500;opt.minmatch = 100;
				opt.kn = 17; opt.wn = 10; opt.cn = 20;
			}else if (opt.step == 3) {
				opt.minlen = 2000;
			}else{
				fprintf(stderr, "[ERROR] unknown preset step '%s'\n", o.arg);
				return 1;
			}
		}
	}
	// assert (opt.step);
	o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, opt_str, long_options)) >= 0) {
		if (c == 'w') ipt.w = atoi(o.arg);

		else if (c == 401) opt.minlen = (int32_t)mm_parse_num(o.arg);
		else if (c == 402) opt.minide = atof(o.arg);
		else if (c == 403) opt.maxhan1 = (int32_t)mm_parse_num(o.arg);
		else if (c == 404) opt.maxhan2 = (int32_t)mm_parse_num(o.arg);
		else if (c == 405) opt.outraw = 1;
		else if (c == 406) opt.minmatch = (int32_t)mm_parse_num(o.arg);
		else if (c == 407) opt.mode = atoi(o.arg);
		else if (c == 408) opt.kn = atoi(o.arg);
		else if (c == 409) opt.wn = atoi(o.arg);
		else if (c == 410) opt.cn = atoi(o.arg);
		else if (c == 411) opt.outctn = 1;
		else if (c == 412) opt.d_factor = atof(o.arg);
		else if (c == 413) opt.dvt = 1;

		else if (c == 'k') ipt.k = atoi(o.arg);
		else if (c == 'H') ipt.flag |= MM_I_HPC;
		else if (c == 'd') fnw = o.arg; // the above are indexing related options, except -I
		else if (c == 'r') opt.bw = (int)mm_parse_num(o.arg);
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'v') mm_verbose = atoi(o.arg);
		else if (c == 'g') opt.max_gap = (int)mm_parse_num(o.arg);
		else if (c == 'G') mm_mapopt_max_intron_len(&opt, (int)mm_parse_num(o.arg));
		else if (c == 'F') opt.max_frag_len = (int)mm_parse_num(o.arg);
		else if (c == 'N') old_best_n = opt.best_n, opt.best_n = atoi(o.arg);
		else if (c == 'p') opt.pri_ratio = atof(o.arg);
		else if (c == 'M') opt.mask_level = atof(o.arg);
		else if (c == 'c') opt.flag |= MM_F_OUT_CG | MM_F_CIGAR;
		else if (c == 'D') opt.flag |= MM_F_NO_DIAG;
		else if (c == 'P') opt.flag |= MM_F_ALL_CHAINS;
		else if (c == 'X') opt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no
		else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
		else if (c == 'Y') opt.flag |= MM_F_SOFTCLIP;
		else if (c == 'L') opt.flag |= MM_F_LONG_CIGAR;
		else if (c == 'y') opt.flag |= MM_F_COPY_COMMENT;
		else if (c == 'T') opt.sdust_thres = atoi(o.arg);
		else if (c == 'n') opt.min_cnt = atoi(o.arg);
		else if (c == 'm') opt.min_chain_score = atoi(o.arg);
		else if (c == 'A') opt.a = atoi(o.arg);
		else if (c == 'B') opt.b = atoi(o.arg);
		else if (c == 's') opt.min_dp_max = atoi(o.arg);
		else if (c == 'C') opt.noncan = atoi(o.arg);
		else if (c == 'I') ipt.batch_size = mm_parse_num(o.arg);
		else if (c == 'K') opt.mini_batch_size = (int)mm_parse_num(o.arg);
		else if (c == 'R') rg = o.arg;
		else if (c == 'h') fp_help = stdout;
		else if (c == '2') opt.flag |= MM_F_2_IO_THREADS;
		else if (c == 'o') {
			if (strcmp(o.arg, "-") != 0) {
				if (freopen(o.arg, "wb", stdout) == NULL) {
					fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m\n", o.arg);
					exit(1);
				}
				if (opt.step == 2){
					char outfile[1024] = {0};
					assert (sprintf(outfile, "%s.bl", o.arg) != -1);
					if ((opt.bl = fopen(outfile, "w")) == NULL) {
						fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m\n", outfile);
						exit(1);
					}
					opt.os = kh_init(ovlh_);
				}
				if (opt.step) opt.outraw = 0;
			}
		}
		else if (c == 300) ipt.bucket_bits = atoi(o.arg); // --bucket-bits
		else if (c == 302) opt.seed = atoi(o.arg); // --seed
		else if (c == 303) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 304) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 306) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED, n_threads = 1; // --print-seed
		else if (c == 307) opt.max_chain_skip = atoi(o.arg); // --max-chain-skip
		else if (c == 339) opt.max_chain_iter = atoi(o.arg); // --max-chain-iter
		else if (c == 308) opt.min_ksw_len = atoi(o.arg); // --min-dp-len
		else if (c == 309) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ, n_threads = 1; // --print-aln-seq
		else if (c == 310) opt.flag |= MM_F_SPLICE; // --splice
		else if (c == 312) opt.flag |= MM_F_NO_LJOIN; // --no-long-join
		else if (c == 313) opt.flag |= MM_F_SR; // --sr
		else if (c == 317) opt.end_bonus = atoi(o.arg); // --end-bonus
		else if (c == 318) opt.flag |= MM_F_INDEPEND_SEG; // --no-pairing
		else if (c == 320) ipt.flag |= MM_I_NO_SEQ; // --idx-no-seq
		else if (c == 321) opt.anchor_ext_shift = atoi(o.arg); // --end-seed-pen
		else if (c == 322) opt.flag |= MM_F_FOR_ONLY; // --for-only
		else if (c == 323) opt.flag |= MM_F_REV_ONLY; // --rev-only
		else if (c == 327) opt.max_clip_ratio = atof(o.arg); // --max-clip-ratio
		else if (c == 328) opt.min_mid_occ = atoi(o.arg); // --min-occ-floor
		else if (c == 329) opt.flag |= MM_F_OUT_MD; // --MD
		else if (c == 330) opt.min_join_flank_ratio = atof(o.arg); // --lj-min-ratio
		else if (c == 331) opt.sc_ambi = atoi(o.arg); // --score-N
		else if (c == 332) opt.flag |= MM_F_EQX; // --eqx
		else if (c == 333) opt.flag |= MM_F_PAF_NO_HIT; // --paf-no-hit
		else if (c == 334) opt.split_prefix = o.arg; // --split-prefix
		else if (c == 335) opt.flag |= MM_F_NO_END_FLT; // --no-end-flt
		else if (c == 336) opt.flag |= MM_F_HARD_MLEVEL; // --hard-mask-level
		else if (c == 337) opt.max_sw_mat = mm_parse_num(o.arg); // --cap-sw-mat
		else if (c == 338) opt.max_qlen = mm_parse_num(o.arg); // --max-qlen
		else if (c == 340) junc_bed = o.arg; // --junc-bed
		else if (c == 342) opt.flag |= MM_F_SAM_HIT_ONLY; // --sam-hit-only
		else if (c == 314) { // --frag
			yes_or_no(&opt, MM_F_FRAG_MODE, o.longidx, o.arg, 1);
		} else if (c == 315) { // --secondary
			yes_or_no(&opt, MM_F_NO_PRINT_2ND, o.longidx, o.arg, 0);
		} else if (c == 316) { // --cs
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR;
			if (o.arg == 0 || strcmp(o.arg, "short") == 0) {
				opt.flag &= ~MM_F_OUT_CS_LONG;
			} else if (strcmp(o.arg, "long") == 0) {
				opt.flag |= MM_F_OUT_CS_LONG;
			} else if (strcmp(o.arg, "none") == 0) {
				opt.flag &= ~MM_F_OUT_CS;
			} else if (mm_verbose >= 2) {
				fprintf(stderr, "[WARNING]\033[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\033[0m\n");
			}
		} else if (c == 319) { // --splice-flank
			yes_or_no(&opt, MM_F_SPLICE_FLANK, o.longidx, o.arg, 1);
		} else if (c == 324) { // --heap-sort
			yes_or_no(&opt, MM_F_HEAP_SORT, o.longidx, o.arg, 1);
		} else if (c == 326) { // --dual
			yes_or_no(&opt, MM_F_NO_DUAL, o.longidx, o.arg, 0);
		} else if (c == 'S') {
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR | MM_F_OUT_CS_LONG;
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING]\033[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\033[0m\n");
		} else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'f') {
			double x;
			char *p;
			x = strtod(o.arg, &p);
			if (x < 1.0) opt.mid_occ_frac = x, opt.mid_occ = 0;
			else opt.mid_occ = (int)(x + .499);
			if (*p == ',') opt.max_occ = (int)(strtod(p+1, &p) + .499);
		} else if (c == 'u') {
			if (*o.arg == 'b') opt.flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV; // both strands
			else if (*o.arg == 'f') opt.flag |= MM_F_SPLICE_FOR, opt.flag &= ~MM_F_SPLICE_REV; // match GT-AG
			else if (*o.arg == 'r') opt.flag |= MM_F_SPLICE_REV, opt.flag &= ~MM_F_SPLICE_FOR; // match CT-AC (reverse complement of GT-AG)
			else if (*o.arg == 'n') opt.flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV); // don't try to match the GT-AG signal
			else {
				fprintf(stderr, "[ERROR]\033[1;31m unrecognized cDNA direction\033[0m\n");
				return 1;
			}
		} else if (c == 'z') {
			opt.zdrop = opt.zdrop_inv = strtol(o.arg, &s, 10);
			if (*s == ',') opt.zdrop_inv = strtol(s + 1, &s, 10);
		} else if (c == 'O') {
			opt.q = opt.q2 = strtol(o.arg, &s, 10);
			if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt.e = opt.e2 = strtol(o.arg, &s, 10);
			if (*s == ',') opt.e2 = strtol(s + 1, &s, 10);
		}
	}
	if ((opt.flag & MM_F_SPLICE) && (opt.flag & MM_F_FRAG_MODE)) {
		fprintf(stderr, "[ERROR]\033[1;31m --splice and --frag should not be specified at the same time.\033[0m\n");
		return 1;
	}
	if (!fnw && !(opt.flag&MM_F_CIGAR) && opt.step != 2 && opt.mode != 3)
		ipt.flag |= MM_I_NO_SEQ;
	if (mm_check_opt(&ipt, &opt) < 0)
		return 1;
	if (opt.best_n == 0) {
		fprintf(stderr, "[WARNING]\033[1;31m changed '-N 0' to '-N %d --secondary=no'.\033[0m\n", old_best_n);
		opt.best_n = old_best_n, opt.flag |= MM_F_NO_PRINT_2ND;
	}

	if (argc == o.ind || fp_help == stdout) {
		fprintf(fp_help, "Usage: minimap2-nd [options] <target.fa>|<target.idx> [query.fa] [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  NextDenovo:\n");
		fprintf(fp_help, "    --step [1|2|3] preset options for NextDenovo [required]\n");
		// fprintf(fp_help, "    --outraw       do not compresse alignments\n");
		fprintf(fp_help, "    --minlen INT   min overlap length [%d]\n", opt.minlen);
		fprintf(fp_help, "    --minmatch INT min match length [%d]\n", opt.minmatch);
		fprintf(fp_help, "    --minide FLOAT min identity [%.2f]\n", opt.minide);
		fprintf(fp_help, "    --kn INT       k-mer size (no larger than 28), used to re-align [%d]\n", opt.kn);
		fprintf(fp_help, "    --wn INT       minizer window size, used to re-align [%d]\n",  opt.wn);
		fprintf(fp_help, "    --df FLOAT     f factor [%f]\n",  opt.d_factor);
		fprintf(fp_help, "    --mode [0,1,2] re-align mode, 0:disable 1:fast mode, low accuracy 2:slow mode, high accuracy 3:hifi mode [%d]\n", opt.mode);
		fprintf(fp_help, "    --cn INT       do re-align for every INT reads, larger is faster [%d]\n", opt.cn);
		fprintf(fp_help, "    --maxhan1 INT  max over hang length, used to re-align [%d]\n", opt.maxhan1);
		fprintf(fp_help, "    --maxhan2 INT  max over hang length, used to filter contained reads [%d]\n",opt.maxhan2);
		fprintf(fp_help, "    --outctn       output contained alignments [%d]\n", opt.outctn);
		fprintf(fp_help, "    --dvt          only output dovetail or contained alignments [%d]\n", opt.dvt);
		fprintf(fp_help, "    -x ava-hifi    Hifi read overlap\n");

		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -H           use homopolymer-compressed k-mer (preferrable for PacBio)\n");
		fprintf(fp_help, "    -k INT       k-mer size (no larger than 128) [%d]\n", ipt.k);
		fprintf(fp_help, "    -w INT       minimizer window size [%d]\n", ipt.w);
		fprintf(fp_help, "    -I NUM       split index for every ~NUM input bases [4G]\n");
		fprintf(fp_help, "    -d FILE      dump index to FILE []\n");
		fprintf(fp_help, "  Mapping:\n");
		fprintf(fp_help, "    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]\n");
		fprintf(fp_help, "    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]\n");
		fprintf(fp_help, "    -r NUM       bandwidth used in chaining and DP-based alignment [%d]\n", opt.bw);
		fprintf(fp_help, "    -n INT       minimal number of minimizers on a chain [%d]\n", opt.min_cnt);
		fprintf(fp_help, "    -m INT       minimal chaining score (matching bases minus log gap penalty) [%d]\n", opt.min_chain_score);
//		fprintf(fp_help, "    -T INT       SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres); // TODO: this option is never used; might be buggy
		fprintf(fp_help, "    -X           skip self and dual mappings (for the all-vs-all mode)\n");
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
		fprintf(fp_help, "  Alignment:\n");
		fprintf(fp_help, "    -A INT       matching score [%d]\n", opt.a);
		fprintf(fp_help, "    -B INT       mismatch penalty [%d]\n", opt.b);
		fprintf(fp_help, "    -O INT[,INT] gap open penalty [%d,%d]\n", opt.q, opt.q2);
		fprintf(fp_help, "    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [%d,%d]\n", opt.e, opt.e2);
		fprintf(fp_help, "    -z INT[,INT] Z-drop score and inversion Z-drop score [%d,%d]\n", opt.zdrop, opt.zdrop_inv);
		fprintf(fp_help, "    -s INT       minimal peak DP alignment score [%d]\n", opt.min_dp_max);
		fprintf(fp_help, "    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
		fprintf(fp_help, "  Input/Output:\n");
		fprintf(fp_help, "    -a           output in the SAM format (PAF by default)\n");
		fprintf(fp_help, "    -o FILE      output alignments to FILE [stdout]\n");
		fprintf(fp_help, "    -L           write CIGAR with >65535 ops at the CG tag\n");
		fprintf(fp_help, "    -R STR       SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
		fprintf(fp_help, "    -c           output CIGAR in PAF\n");
		fprintf(fp_help, "    --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]\n");
		fprintf(fp_help, "    --MD         output the MD tag\n");
		fprintf(fp_help, "    --eqx        write =/X CIGAR operators\n");
		fprintf(fp_help, "    -Y           use soft clipping for supplementary alignments\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]\n");
//		fprintf(fp_help, "    -v INT       verbose level [%d]\n", mm_verbose);
		fprintf(fp_help, "    --version    show version number\n");
		fprintf(fp_help, "  Preset:\n");
		fprintf(fp_help, "    -x STR       preset (always applied before other options; see minimap2.1 for details) []\n");
		fprintf(fp_help, "                 - map-pb/map-ont: PacBio/Nanopore vs reference mapping\n");
		fprintf(fp_help, "                 - ava-pb/ava-ont: PacBio/Nanopore read overlap\n");
		fprintf(fp_help, "                 - asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5%% sequence divergence\n");
		fprintf(fp_help, "                 - splice: long-read spliced alignment\n");
		fprintf(fp_help, "                 - sr: genomic short-read mapping\n");
		fprintf(fp_help, "\nSee `man ./minimap2.1' for detailed description of these and other advanced command-line options.\n");
		return fp_help == stdout? 0 : 1;
	}

	// if (!opt.step){
	// 	fprintf(stderr, "[ERROR] failed to find option: --step.\n");
	// 	return 1;
	// }else{
	if (opt.step == 1 || opt.step == 2){
		encode_tbl = init_encode_table(NUM_LIMIT);
		if (opt.step == 2 && !opt.outraw) init_ovl_mode(stdout, 10);
		if (opt.step == 2 && opt.mode == 1 && opt.minide < 0.01) opt.minide = 0.01;
		// if (opt.cn < 20) opt.cn = 20;
		if (opt.mode == 1 && opt.cn == 20) opt.cn = 50;
	}

	if ((opt.flag & MM_F_SR) && argc - o.ind > 3) {
		fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
		return 1;
	}
	idx_rdr = mm_idx_reader_open(argv[o.ind], &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s'\n", argv[o.ind]);
		return 1;
	}
	if (!idx_rdr->is_idx && fnw == 0 && argc - o.ind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
		mm_idx_reader_close(idx_rdr);
		return 1;
	}

	if (opt.best_n == 0 && (opt.flag&MM_F_CIGAR) && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m `-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments.\033[0m\n");
	while ((mi = mm_idx_reader_read(idx_rdr, n_threads)) != 0) {
		if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
			fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
			mm_idx_destroy(mi);
			mm_idx_reader_close(idx_rdr);
			return 1;
		}
		if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
			if (mm_idx_reader_eof(idx_rdr)) {
				mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
			} else {
				mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
				if (opt.split_prefix == 0 && mm_verbose >= 2)
					opt.split_prefix = "nd";
					fprintf(stderr, "[M::%s::%.3f*%.2f] auto set --split-prefix to %s\n",
						__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), opt.split_prefix);
					// fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted. Please use --split-prefix.\033[0m\n");
			}
		}
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
		if (argc != o.ind + 1) mm_mapopt_update(&opt, mi);
		if (mm_verbose >= 3) mm_idx_stat(mi);
		if (junc_bed) mm_idx_bed_read(mi, junc_bed, 1);
		if (!(opt.flag & MM_F_FRAG_MODE)) {
			for (i = o.ind + 1; i < argc; ++i)
				mm_map_file(mi, argv[i], &opt, n_threads);
		} else {
			mm_map_file_frag(mi, argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_threads);
		}
		mm_idx_destroy(mi);
	}
	n_parts = idx_rdr->n_parts;
	mm_idx_reader_close(idx_rdr);

	if (opt.split_prefix)
		mm_split_merge(argc - (o.ind + 1), (const char**)&argv[o.ind + 1], &opt, n_parts);

	if (fflush(stdout) == EOF) {
		fprintf(stderr, "[ERROR] failed to write the results\n");
		exit(EXIT_FAILURE);
	}

	if (encode_tbl) destroy_table(encode_tbl);
	if (opt.bl) {
		out_bl((khash_t(ovlh_) *)opt.os, opt.bl);
		fclose (opt.bl);
		kh_destroy(ovlh_, (khash_t(ovlh_) *)opt.os);
	}

	if (mm_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - mm_realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
