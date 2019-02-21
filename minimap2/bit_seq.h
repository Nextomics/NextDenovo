#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)

typedef struct __kstring_t {
  size_t l, m;
  char *s;
} kstring_t;

typedef struct __kstream_t {
  unsigned char *buf;
  int begin, end, is_eof;
  FILE *f;
} kstream_t;

typedef struct {
  kstring_t name, comment, seq, qual;
  int last_char;
  kstream_t *f;
  uint64_t offset;
} kseq_t;

extern char **bases_arr;

char **init_bases(int n)
{
  int i, j;
  char **arr = (char **)malloc(sizeof(char *) * n);
  for(i = 0; i < n; i++) {
    arr[i] = (char *)malloc(sizeof(char) * 9);
    for(j = 0; j < 8; j++)
      arr[i][j] = "ACGT"[i >> (14 - (j << 1)) & 3];
    arr[i][8] = '\0';
  }
  return arr;
}

kstream_t *ks_init(FILE *fp)
{
  kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));
  ks->f = fp;
  ks->buf = NULL;
  ks->begin = ks->end = ks->is_eof = 0;
  return ks;
}

kseq_t *kseq_init(FILE *fp)
{
  kseq_t *s = (kseq_t *)calloc(1, sizeof(kseq_t));
  s->f = ks_init(fp);
  s->name.s = s->comment.s = s->seq.s = s->qual.s = NULL;
  s->name.l = s->name.m = s->comment.l = s->comment.m = 0;
  s->seq.l = s->seq.m = s->qual.l = s->qual.m = 0;
  if (!bases_arr) bases_arr = init_bases(65536);
  return s;
}

void ks_destroy(kstream_t *ks)
{
  if (ks) {
    free(ks->buf);
    free(ks);
  }
}

void kseq_destroy(kseq_t *ks)
{
  if (!ks) return;
  free(ks->comment.s); free(ks->qual.s);
  ks_destroy(ks->f);
  free(ks);
}

int kseq_read(kseq_t *seq)
{
  if (seq->seq.l) {
    free(seq->name.s);
    free(seq->seq.s);
    seq->name.s = seq->seq.s = NULL;
  }

  uint32_t name, length;
  //Parse sequence name
  fseeko(seq->f->f, seq->offset, SEEK_SET);
  fread(&name, sizeof(uint32_t), 1, seq->f->f);
  //Check end of file!
  if(feof(seq->f->f)) {
    seq->f->is_eof = 1;
    return EOF;
  }
  seq->offset += 4;
  //Parse sequence length
  fseeko(seq->f->f, seq->offset, SEEK_SET);
  fread(&length, sizeof(uint32_t), 1, seq->f->f);
  seq->offset += 4;
  //Parse sequence, convert 2bit to bases
  fseeko(seq->f->f, seq->offset, SEEK_SET);
  uint32_t end = (((length - 1) >> 4) + 1) << 4;
  uint32_t i, cnt;
  cnt = end >> 4;
  uint32_t *tmp = (uint32_t *)malloc(sizeof(uint32_t) * cnt);
  seq->seq.s = (char *)malloc(sizeof(char) * (end + 1));
  fread(tmp, sizeof(uint32_t), cnt, seq->f->f);
  for(i = 0; i < cnt; i++) {
    memcpy(seq->seq.s + (i << 4), bases_arr[tmp[i] >> 16 & 0xFFFF], 8);
    memcpy(seq->seq.s + (i << 4) + 8, bases_arr[tmp[i] & 0xFFFF], 8);
  }
  seq->seq.s[length] = '\0';
  seq->seq.l = length;
  seq->name.s = (char *)malloc(sizeof(char) * 128);
  sprintf(seq->name.s, "%d", name);
  seq->name.l = strlen(seq->name.s);
  seq->offset += (end >> 2);

  free(tmp);
  return seq->seq.l;
}
