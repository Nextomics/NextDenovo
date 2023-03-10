#include <zlib.h>
#include "mseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	char *seq;
	char *buffer;
	char **bases_arr;
	uint32_t *tmp;
	uint32_t m;
} sbbuf;

int init_seq_mode(FILE *fp);
uint32_t seq2bit(FILE *fp, uint32_t name, uint32_t len, char *seq);
void seq2bit1(uint32_t *s, uint32_t len, char *seq);
void bit2seq1(uint32_t *s, uint32_t len, char *seq);
sbbuf *init_sbbuf(uint32_t len);
void destroy_sbbbuf(sbbuf *buf);
void subbit(sbbuf *buf, FILE *fp, int rev, uint64_t offset, uint32_t start, uint32_t end);
void subbit_(sbbuf *buf, uint32_t *fp, int rev, uint64_t offset, uint32_t start, uint32_t end);
void subfa(sbbuf *buf, FILE *fp, int rev, uint64_t offset, uint32_t start, uint32_t end);
int kseq_r(kseq_t *seq);
