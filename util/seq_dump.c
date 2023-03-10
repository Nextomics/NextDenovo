#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <zlib.h>
#include <libgen.h>
#include <sys/stat.h>
#include <sys/time.h>
#include "../lib/bseq.h"
// #include "kseq.h"
// KSEQ_INIT(gzFile, gzread)

#define LEN_LIMIT 1000000UL
#define FILE_MAX_SIZE (uint64_t)-1

static int part_cnt = 1;
static uint32_t part_n = 0;
static int seed_cnt = 1;

#define CHAR_NUM 0x400

typedef struct {
  char part_pre[CHAR_NUM], seed_pre[CHAR_NUM];
  char part_idx_pre[CHAR_NUM], seed_idx_pre[CHAR_NUM];
  uint32_t filter_length, seed_filter_length;
  uint64_t block_size;
  int seed_n;
} opt_t;

typedef struct {
  FILE *fn, *fidx;
  uint64_t offset;
  uint64_t length;
} fout_t;

void convert_2bit(kseq_t *seq, FILE *output, FILE *idx, uint64_t *start, uint32_t n_cnt)
{  
  uint32_t n = seq->seq.l < LEN_LIMIT ? seq->seq.l : LEN_LIMIT;
  fprintf(idx, "%u\t%lu\t%u\n", n_cnt, *start + 8, n);
  *start += seq2bit(output, n_cnt, n, seq->seq.s);
}

static void open_part(fout_t *ptr, char *part_pre, char *idx_pre, int cnt)
{
  char fn[CHAR_NUM], fidx[CHAR_NUM];
  sprintf(fn, "%s%03d%s", part_pre, cnt, ".2bit");
  sprintf(fidx, "%s%03d%s", idx_pre, cnt, ".idx");
  ptr->fn = fopen(fn, "w");
  ptr->fidx = fopen(fidx, "w");
  if (ptr->fn == NULL || ptr->fidx == NULL) {
    fprintf(stderr, "Failed create split file or index!\n");
    fprintf(stderr, "2bit file name: %s\n", fn);
    fprintf(stderr, "idx file name: %s\n", fidx);
    exit(1);
  }
  ptr->offset = init_seq_mode(ptr->fn);
}

static void split_data(char *line, opt_t opt, fout_t *ptr, fout_t *part_ptr, fout_t *bonus_ptr)
{
  static int flag = 0;
  static int seed_num;

  gzFile fp;
  fp = gzopen(line, "r");
  gzbuffer(fp, 131072);
  if (fp == NULL) {
    fprintf(stderr, "Error! %s does not exist!", line);
    exit(1);
  }
  kseq_t *seq;
  seq = kseq_init(fp);
  while(kseq_read(seq) >= 0) {
    if (seq->seq.l >= opt.filter_length && seq->seq.l < opt.seed_filter_length) {
      part_ptr->length += seq->seq.l;
      if (part_ptr->length > opt.block_size) {
        fclose(part_ptr->fn); fclose(part_ptr->fidx);
        // part_ptr->offset = 0;
        part_ptr->length = seq->seq.l;
        part_cnt++;
        open_part(part_ptr, opt.part_pre, opt.part_idx_pre, part_cnt);
      }
      convert_2bit(seq, part_ptr->fn, part_ptr->fidx, &part_ptr->offset, part_n);
      ++part_n;
    } else if (seq->seq.l >= opt.seed_filter_length && seq->seq.l < LEN_LIMIT) {
      if (!flag) {
        if (seed_cnt > opt.seed_n) seed_cnt = 1;
        if (ptr[seed_cnt - 1].length + seq->seq.l < FILE_MAX_SIZE) {
          convert_2bit(seq, ptr[seed_cnt - 1].fn, ptr[seed_cnt - 1].fidx, \
          &ptr[seed_cnt - 1].offset, part_n);
          ptr[seed_cnt - 1].length += seq->seq.l;
          ++seed_cnt;
        } else {
          flag = 1;
          seed_num = opt.seed_n + 1;
          open_part(bonus_ptr, opt.seed_pre, opt.seed_idx_pre, seed_num);
          convert_2bit(seq, bonus_ptr->fn, bonus_ptr->fidx, &bonus_ptr->offset, part_n);
          bonus_ptr->length += seq->seq.l;
          ++seed_num;
        }
      } else {
        if (bonus_ptr->length + seq->seq.l < FILE_MAX_SIZE) {
          bonus_ptr->length += seq->seq.l;
        } else {
          fclose(bonus_ptr->fn); fclose(bonus_ptr->fidx);
          // bonus_ptr->offset = 0;
          bonus_ptr->length = seq->seq.l;
          open_part(bonus_ptr, opt.seed_pre, opt.seed_idx_pre, seed_num);
          ++seed_num;
        }
        convert_2bit(seq, bonus_ptr->fn, bonus_ptr->fidx, &bonus_ptr->offset, part_n);
      }
      ++part_n;
    }
  }
  kseq_destroy(seq);
  gzclose(fp);
}

fout_t *init_seed_file(opt_t opt)
{
  int i;
  char fname[CHAR_NUM], fidx[CHAR_NUM];
  fout_t *ptr = (fout_t *)malloc(sizeof(fout_t) * opt.seed_n);
  for(i = 0; i < opt.seed_n; i++) {
    sprintf(fname, "%s%03d%s", opt.seed_pre, i + 1, ".2bit");
    sprintf(fidx, "%s%03d%s", opt.seed_idx_pre, i + 1, ".idx");
    ptr[i].fn = fopen(fname, "w");
    ptr[i].fidx = fopen(fidx, "w");
    ptr[i].offset = init_seq_mode(ptr[i].fn);
    ptr[i].length = 0UL;
  }
  return ptr;
}

void close_seed_file(opt_t opt, fout_t *ptr)
{
  int i;
  for(i = 0; i < opt.seed_n; i++) {
    fclose(ptr[i].fn);
    fclose(ptr[i].fidx);
  }
}

static inline uint64_t mm_parse_num(const char *str)
{
  double x;
  char *p;
  x = strtod(str, &p);
  if (*p == 'G' || *p == 'g') x *= 1e9;
  else if (*p == 'M' || *p == 'm')  x *= 1e6;
  else if (*p == 'K' || *p == 'k')  x *= 1e3;
  return (uint64_t)(x + .499);
}

static int usage()
{
  fprintf(stderr, "Usage: seq_dump [options] input.fofn\n");
  fprintf(stderr, "Dump fasta/fastq[.gz] files\n\n");

  fprintf(stderr, "Options:\n");
  fprintf(stderr, " -f minimum read length\n");
  fprintf(stderr, " -s minimum seed length\n");
  fprintf(stderr, " -b block size (Mb or Gb), 0=inf\n");
  fprintf(stderr, " -n number of seed subfiles in total\n");
  fprintf(stderr, " -d output directory\n");
  return 1;
}

int main(int argc, char *argv[])
{
  int c;
  opt_t opt;
  while((c = getopt(argc, argv, "f:s:b:n:d:")) != -1) {
    switch (c) {
      case 'f':
        opt.filter_length = mm_parse_num(optarg);
        break;
      case 's':
        if (mm_parse_num(optarg) <= opt.filter_length) {
          fprintf(stderr, "\033[1;31mError! Seed filter length should be larger than filter length!\033[0m\n");
          return usage();
        }
        opt.seed_filter_length = mm_parse_num(optarg);
        break;
      case 'b':
        // if (mm_parse_num(optarg) > 16e9) {
        //   fprintf(stderr, "\033[1;31mError! Block size should be less than 16Gb!\033[0m\n");
        //   return usage();
        // }
        opt.block_size = mm_parse_num(optarg);
        if (opt.block_size == 0) opt.block_size = FILE_MAX_SIZE;
        break;
      case 'n':
        opt.seed_n = atoi(optarg);
        break;
      case 'd':
        if (access(optarg, 0) != 0) {
          if (mkdir(optarg, 0755) == -1) {
            fprintf(stderr, "Failed create output directory!\n");
            exit(1);
          }
        }
        sprintf(opt.part_pre, "%s%s", optarg, "/input.part.");
        sprintf(opt.part_idx_pre, "%s%s", optarg, "/.input.part.");
        sprintf(opt.seed_pre, "%s%s", optarg, "/input.seed.");
        sprintf(opt.seed_idx_pre, "%s%s", optarg, "/.input.seed.");
        break;
      default:
        return usage();
    }
  }
  if (optind + 1 > argc) return usage();
  char *fofn = strdup(argv[optind]);
  char *path = dirname(fofn);
  FILE *stream = fopen(argv[optind], "r");
  if (stream == NULL) {
    fprintf(stderr, "Failed open input file list!\n");
    return usage();
  }

  fout_t *seed_ptr = init_seed_file(opt);
  fout_t *part_ptr = (fout_t *)calloc(1, sizeof(fout_t));
  fout_t *seed_bonus_ptr = (fout_t *)calloc(1, sizeof(fout_t));
  open_part(part_ptr, opt.part_pre, opt.part_idx_pre, part_cnt);

  char *line = NULL;
  size_t len = 0;
  ssize_t nread;
  while((nread = getline(&line, &len, stream)) != -1) {
    if (nread != 1 && line[0] != '#') {
      line[strlen(line) - 1] = '\0';
      if (line[0] == '/')
        split_data(line, opt, seed_ptr, part_ptr, seed_bonus_ptr);
      else {
        char file_path[1024];
        sprintf(file_path, "%s/%s", path, line);
        split_data(file_path, opt, seed_ptr, part_ptr, seed_bonus_ptr);
      }
    }
  }
  free(line);
  close_seed_file(opt, seed_ptr);

  fclose(part_ptr->fn), fclose(part_ptr->fidx);
  if (seed_bonus_ptr->fn) fclose(seed_bonus_ptr->fn), fclose(seed_bonus_ptr->fidx);
  free(part_ptr); free(seed_bonus_ptr);

  free(fofn); free(seed_ptr);
  fclose(stream);
  return 0;
}
