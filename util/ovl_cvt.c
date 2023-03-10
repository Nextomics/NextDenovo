#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include "../lib/ovl.h"

void encode(FILE *fp, FILE *fq, bytes_t *arr, prev_t *id, int mode)
{
  buffer_t *buf = init_buffer(ECBUFSIZE);
  char *line = NULL;
  size_t len = 0;
  init_ovl_mode(fq, mode);
  if (mode == 8){
    overlap ovl;
    while(getline(&line, &len, fp) != -1) {
      sscanf(line, "%u\t%hhu\t%u\t%u\t%u\t%u\t%u\t%u\n", &ovl.qname, &ovl.rev,\
        &ovl.qs, &ovl.qe, &ovl.tname, &ovl.ts, &ovl.te, &ovl.match);
      encode_ovl(fq, arr, id, &ovl, buf);
    }
  }else{
    overlap_i ovl;
    while(getline(&line, &len, fp) != -1) {
      sscanf(line, "%u\t%hhu\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", &ovl.qname, &ovl.rev, &ovl.qs, \
        &ovl.qe, &ovl.tname, &ovl.ts, &ovl.te, &ovl.qlen, &ovl.tlen, &ovl.identity);
      encode_ovl_i(fq, arr, id, &ovl, buf);
    }
  }
  flush_buffer(fq, buf);
  if (line) free(line);
}

void decode(FILE *fp, FILE *fq, uint32_t *arr, prev_t *id, int mode)
{
  buffer_t *buf = init_buffer(DCBUFSIZE);
  uint32_t ovl[10];
  uint32_t qlen = 0, tlen = 0;
  while(decode_ovl(fp, arr, id, ovl, buf, mode) >= 0) {
    if (mode == 8){
      fprintf(fq, "%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", ovl[0], ovl[1], \
        ovl[2], ovl[3], ovl[4], ovl[5], ovl[6], ovl[7]);
    }else{
      if (ovl[7]) qlen = ovl[7];
      if (ovl[8]) tlen = ovl[8];
      fprintf(fq, "%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", ovl[0], ovl[1], \
        ovl[2], ovl[3], ovl[4], ovl[5], ovl[6], qlen, tlen, ovl[9]);
    }
  }
  flush_buffer(NULL, buf);
}

static int usage()
{
  fprintf(stderr, "Usage: ovl_cvt [options] input.ovl\n");
  fprintf(stderr, "Compress or uncompress ovl files\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, " -m INT    conversion mode (0 for compress, 1 for uncompress)\n");
  return 1;
}

int main(int argc, char *argv[])
{
  int c, mode = 0;
  while((c = getopt(argc, argv, "m:")) != -1) {
    switch (c) {
      case 'm':
        mode = atoi(optarg);
        break;
      default:
        return usage();
    }
  }
  if (optind + 1 > argc) return usage();
  FILE *fp;
  fp = mode ? fopen(argv[optind], "rb"): fopen(argv[optind], "r");
  if (fp == NULL) {
    fprintf(stderr, "Failed open input file!\n");
    return usage();
  }
  prev_t id;
  id.prev_qname = id.prev_tname = 0;
  if (mode) {
    uint32_t *arr = init_decode_table(NUM_LIMIT);
    mode = find_ovlb_mode(fp);
    decode(fp, stdout, arr, &id, mode);
    destroy_table(arr);
  } else {
    bytes_t *arr = init_encode_table(NUM_LIMIT);
    mode = find_ovlt_mode(fp);
    encode(fp, stdout, arr, &id, mode);
    destroy_table(arr);
  }
  fclose(fp);
  return 0;
}