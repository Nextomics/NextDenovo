#include <stdio.h>
#include <stdint.h>
#include "../lib/bseq.h"

static int usage()
{
	fprintf(stderr, "Usage: seq_bit [options] input.[fa|bit]\n");
	fprintf(stderr, "Compress or uncompress fa/bit files, seq name must be number\n\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (optind + 1 > argc) return usage();
	gzFile fp = gzopen(argv[optind], "r");
	if (fp == NULL) {
		fprintf(stderr, "Failed open input file!\n");
		return usage();
	}
	
	int lable = 1;
	kseq_t *seq = kseq_init(fp);
	while ((kseq_r(seq)) >= 0) {
		if (seq->fm == 1){
			if (lable) {init_seq_mode(stdout); lable = 0;}
			seq2bit(stdout, strtoul(seq->name.s, NULL, 10), seq->seq.l, seq->seq.s);
		}else{
			printf(">%s\n%s\n", seq->name.s, seq->seq.s);
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
