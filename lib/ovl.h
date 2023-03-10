#ifndef OVL_H
#define OVL_H
#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>

#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x ,y) ((x) > (y) ? (y) : (x))
#define NUM_LIMIT 100000000UL
#define ECBUFSIZE 104857600
#define DCBUFSIZE 1048576
#define EDGEBACKLEN 10

typedef struct {
  uint8_t *buffer;
  uint32_t buffer_i;
  size_t read_n;
} buffer_t;

typedef struct {
  uint8_t rev;
  uint32_t qname, qs, qe;
  uint32_t tname, ts, te;
  uint32_t match;
} overlap;

typedef struct {
    uint8_t rev;
    uint32_t qname, qs, qe, qlen; 
    uint32_t tname, ts, te, tlen; 
    uint32_t identity;
} overlap_i;

typedef struct {
  uint8_t bytes[5];
  int n;
} bytes_t;

typedef struct {
  uint32_t prev_qname;
  uint32_t prev_tname;
} prev_t;

bytes_t *init_encode_table(uint32_t n);
uint32_t *init_decode_table(uint64_t n);
void destroy_table(void *arr);
buffer_t *init_buffer(uint32_t buffersize);
void flush_buffer(FILE *fp, buffer_t *buf);
void init_ovl_mode(FILE *fp, int mode);
int find_ovlb_mode(FILE *fp);
int find_ovlt_mode(FILE *fp);
void encode_ovl(FILE *fp, bytes_t *arr, prev_t *id, overlap *ovl, buffer_t *buf);
int decode_ovl(FILE *fp, uint32_t *arr, prev_t *id, uint32_t *ovl, buffer_t *buf, int n);
void encode_ovl_i(FILE *fp, bytes_t *arr, prev_t *id, overlap_i *ovl, buffer_t *buf);

#ifndef LIST32
#define LIST32
typedef struct {
  uint32_t *data;
  uint64_t i;
  uint64_t im;
} list;
int init_list(list *l);
int reallocate_list(list *l, uint64_t s);
void destroy_list(list *l);
#endif

#include "../util/khash.h"
#define INIT_ALNM 5
#define MAX_CON 2 //[1,2]

typedef struct {
  uint32_t s, e;
} aln;

typedef struct {//size 40
  uint16_t lc; // 5' depth
  uint16_t rc; // 3' depth
  uint16_t alni; // index of aln
  uint16_t alnm; // max size of aln
  uint32_t con:2; //contained
  uint32_t lim:15; // max 5' identity
  uint32_t rim:15; // max 3' identity
  uint32_t llm; // max 5' aln len
  uint32_t rlm; // max 3' aln len
  uint32_t len;
  aln alnl;
  aln *alns; // [s, e]
} ovlinfo_aln;

// typedef struct {//size 28
//   uint16_t lc; // 5' depth
//   uint16_t rc; // 3' depth
//   uint16_t le; // min align start
//   uint16_t re; // max align end
//   uint32_t con:2; //contained
//   uint32_t lim:15; // max 5' identity
//   uint32_t rim:15; // max 3' identity
//   uint32_t llm; // max 5' aln len
//   uint32_t rlm; // max 3' aln len
//   uint32_t lid; // 5' node index
//   uint32_t rid; // 3' node index
// } ovlinfo;

typedef struct {//size 28
  uint32_t lc; // 5' depth
  uint32_t rc; // 3' depth
  uint32_t le; // min align start
  uint32_t re; // max align end
  uint32_t con:2; //contained
  uint32_t lim:15; // max 5' identity
  uint32_t rim:15; // max 3' identity
  uint32_t llm; // max 5' aln len
  uint32_t rlm; // max 3' aln len
#if GENOME_SIZE == 2
  uint64_t lid; // 5' node index
  uint64_t rid; // 3' node index
#else
  uint32_t lid; // 5' node index
  uint32_t rid; // 3' node index
#endif
} ovlinfo;

KHASH_MAP_INIT_INT(ovlh_, ovlinfo_aln);
KHASH_MAP_INIT_INT(ovlinfo_, ovlinfo);

uint32_t find_alnse(ovlinfo_aln *s, aln *aln);
void out_bl(khash_t(ovlh_) *os, FILE *f);
void read_bl(khash_t(ovlh_) *os, char *file, list *l);
int filter_ovl(overlap_i *ovl, khash_t(ovlh_) *os, int32_t maxhan1, int32_t maxhan2);

#endif
