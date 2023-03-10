#include <inttypes.h>

#define DAG_MAX_LENGTH 10000
#define DAG_MAX_RATIO 0.8//0.5 TODO pls test
#define DAG_MIN_QV 40
#define LQBASE_MIN_QV 20 //40
#define LQREG_MAX_GAP 10 //50
#define LQREG_MAX_LEN 100 //500
#define LQREG_MAX_COUNT 10

#define LQSEQ_MAX_CAN_COUNT 40
#define LQSEQ_MAX_COUNT 30
#define LQSEQ_MAX_REV_LEN 2000

#define KMER_RANGE 40
#define KMER_MAX_SEQ 10
#define KMER_LEN 8 // do not change
#define KMER_LEN_COUNT 65536 // do not change

#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x ,y) ((x) > (y) ? (y) : (x))
#define mabs(x, y) ((x) > (y) ? (x) - (y) : (y) - (x))
#define SWAP(x, y, T) do { T SWAP = x; x = y; y = SWAP; } while (0)


extern unsigned int lqseq_max_length;

typedef struct {
    uint8_t q_base;
    uint16_t delta;
    int t_pos;
} align_tag;

typedef struct {
    align_tag p; // current position
    align_tag pp; // pre pf currect position
    align_tag ppp; // pre of pre of current position
} align_tag_t;

typedef struct {
    unsigned int len;
    unsigned int aln_t_s;
    align_tag *align_tags;
} align_tags_t;

typedef struct {
    unsigned int start;
    unsigned int end;
    unsigned int lqlen;
    unsigned int lq_total_len;
} lqreg;

typedef struct {
    uint8_t indexs;// selected index of start reads
    uint8_t indexe;//selected index of end reads
    int len;
    unsigned int lqcount;
    unsigned int start;
    unsigned int end;
    unsigned int sudoseed_len;
    char *sudoseed;
    struct seq_ {
        uint16_t order;
        uint16_t kscore;
        uint16_t len;
        char seq[DAG_MAX_LENGTH];
    } *seqs;
} lqseq;

typedef struct {
    unsigned int len;
    float identity;
    char *seq;
} consensus_trimed;

typedef struct {
    unsigned int pos;
    char base;
} consensus_base;

typedef struct {
    unsigned int max_size;
    unsigned int len;
    unsigned int uncorrected_len;
    unsigned int lstrip;
    unsigned int rstrip;
    consensus_base *cns_bases;
} consensus_data;

typedef struct {
    align_tag pp;
    align_tag ppp;
    int64_t score:48;
    // int64_t score;
    uint16_t link_count;
} msa_p_d_b_pp_ppp;//pp_ppp

typedef struct {
    uint64_t max_size;
    uint64_t len;
    msa_p_d_b_pp_ppp *pp_ppp;
} msa_p_d_b_pp_ppps;

typedef struct {
    unsigned int max_size;
    unsigned int len;
    align_tag best_pp;
    msa_p_d_b_pp_ppp *pp_ppp;
    int64_t best_score:48;
    // int64_t best_score;
    uint16_t best_link_count;
} msa_p_d_b;

typedef struct {
    uint16_t max_size;
    uint16_t coverage;
    msa_p_d_b **d_b;//p->d
} msa_p;

typedef struct {
    unsigned int shift;
    unsigned int aln_len;
    unsigned int aln_t_s;
    unsigned int aln_t_e;
    unsigned int aln_t_len;
    unsigned int aln_q_len;
    char * q_aln_str;
    char * t_aln_str;
} alignment;

// typedef struct {
//     unsigned int min_len_aln;
//     unsigned int max_cov_aln;
//     unsigned int min_len_seed;
//     unsigned int min_cov_aln;
//     unsigned int min_cov_raw;
//     // int split = 0;
//     unsigned int max_count;
//     unsigned int max_thread;
//     unsigned int max_task;
//     float min_error_corrected_ratio;
// } attr;

// typedef struct {
//     const attr *args;
//     unsigned int count;
//     uint32_t seed_name;
//     uint32_t *aln_start;
//     uint32_t *aln_end;
//     char **seqs;
//     consensus_trimed_data consensus_trimed;
// } attrs;

void malloc_vd(int **, uint8_t ***, uint64_t);
void align(char *, int, char *, int, alignment *, int *, uint8_t **);
void align_hq(char *, int, char *, int, alignment *, int *, uint8_t **);
void align_nd(const char *s1, const uint32_t s1_l, const char *s2, const uint32_t s2_l, alignment *aln);

void reverse_str(char *, int);

consensus_trimed *nextCorrect(char **, unsigned int *, unsigned int *, unsigned int, \
    unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float, unsigned int, unsigned int, int);

void free_consensus_trimed(consensus_trimed *);

char *poa_to_consensus(const struct seq_ *,  const int seq_count);
