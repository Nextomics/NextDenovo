#include <inttypes.h>

#define DAG_MAX_LENGTH 10000
#define LQBASE_MIN_QV 20
#define LQSEQ_MAX_COUNT 30
#define LQSEQ_MAX_REV_LEN 2000
#define LQSEQ_INIT_LEN 500

#define KMER_RANGE 40
#define KMER_MAX_SEQ 10
#define KMER_LEN 8 // do not change
#define KMER_LEN_COUNT 65536 // do not change

// #define GAP_MIN_LEN 3
// #define GAP_MIN_RATIO1 0.01
#define GAP_FLANK_LEN 10
#define GAP_BETWEEN_LEN 30
#define GAP_MIN_RATIO2 0.1
#define GAP_MIN_RATIO3 0.6
#define DEL_MIN_LEN 20 //TODO CHECK
#define DEL_MIN_DEPTH_RATIO 0.3

#define LQSEQ_MIN_LEN 0 //8
#define HQSEQ_MIN_LEN 4
#define HQ_MIN_QV 60
#define LQSEQ_MAX_CAN_COUNT 60

#define INS_MIN_CHECK_LEN 100000
#define INS_RADOM_COUNT 50000
#define INS_RADOM_LEN 15000000
#define INS_WIN_STEP 10
#define INS_WIN_DIV 20
#define INS_MIN_DEPTH_RATIO 0.1
#define INS_MIN_DEPTH_RATIO_REFQV 0.3
#define INS_WIN_MIN_SIZE 500
#define INS_CLUSTER_SIZE 1000

#define CLUSTER_MIN_DEPTH_RATIO 0.2

#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x ,y) ((x) > (y) ? (y) : (x))
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
    // unsigned int len;
    unsigned int aln_t_s;
    unsigned int aln_t_e;
    uint8_t *align_tags;
    // align_tag *align_tags;
} align_tags_t;

typedef struct {
    unsigned int start;
    unsigned int end;
    unsigned int lqlen;
    unsigned int lq_total_len;
} lqreg;

typedef struct {
    uint8_t l;
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
        uint32_t len;
        // char seq[DAG_MAX_LENGTH];
        char *seq;
    } *seqs;
} lqseq;

typedef struct {
    unsigned int len;
    float identity;
    char *seq;
} consensus_trimed;

typedef struct {
    consensus_trimed *data;
    int i_m;
} consensus_trimed_data;

typedef struct {
    unsigned int pos;
    char qv;
    char base;
} consensus_base;

typedef struct {
    unsigned int max_size;
    unsigned int len;
    unsigned int uncorrected_len;
    unsigned int lstrip;//TODO REMOVE
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
    // align_tag best_pp;
    msa_p_d_b_pp_ppp *pp_ppp;
    // int64_t best_score:48;
    // int64_t best_score;
    // uint16_t best_link_count;
} msa_p_d_b;

typedef struct {
    uint16_t max_size;
    uint16_t coverage;
    uint16_t l_del;
    uint16_t l_ins;
    msa_p_d_b **d_b;//p->d
} msa_p;

typedef struct {
    unsigned int shift;
    unsigned int aln_len;
	unsigned int max_aln_len;
    unsigned int aln_t_s;
    unsigned int aln_t_e;//not include
    unsigned int aln_t_len;
    unsigned int aln_q_s;
    unsigned int aln_q_e;//not include
    unsigned int aln_q_len;
    char * q_aln_str;
    char * t_aln_str;
} alignment;


typedef struct ref_qv{
    uint32_t ide:12;
    uint32_t ort:10;
    uint32_t irt:10;
    uint32_t p;
} ref_qv;

typedef struct {
    char *n;
    uint32_t *s;
    ref_qv *qv;
    uint32_t qv_l;
    uint32_t length;
} ref_;

typedef struct {
    ref_ *ref;
    uint32_t i;
    uint32_t i_m;
} refs_;

typedef struct {
    uint32_t s, e;//start, end, 0-based
} pos;

typedef struct {
    pos gap;
    uint32_t fs, ds;
    uint32_t score;
} gap;

typedef struct {
    pos gap;
    uint32_t p_id, s_id;
    uint32_t p_s, s_s;
    uint32_t dl_m, l;
    uint8_t *dseq; //TODO shorten
} gap_;

typedef struct {
    uint32_t i, i_m;
    gap_ *gap;
} gaps;

typedef struct {
    pos r;
	uint32_t median;
    uint32_t i_m;
    gap_ *gap[LQSEQ_MAX_CAN_COUNT<<1];
} gap_cluster;

typedef struct {
    uint32_t i, i_m;
    gap_cluster *clusters;
} gap_clusters;

typedef struct {
    uint32_t fs, ds;
    uint32_t i, i_m;
    uint32_t *cigar;
} sup_aln;

typedef struct {
    uint32_t i, i_m;
    sup_aln *aln;
} sup_alns;

typedef struct {
    uint32_t i, i_m;
    pos *reg;
} ld_regs;

typedef struct {
    int8_t strand;
    uint32_t pos;
    char *rname;
    char *cigar;
} satag;

typedef struct {
    satag *sa;
    int i;
    int i_m;
} satags;

typedef struct {
    consensus_data **consensus;
    int w, s;
    int i, i_m;
} consensuss_data;


void malloc_vd(int **, uint8_t ***, uint64_t);
void align(char *, int, char *, int, alignment *, int *, uint8_t **);

void reverse_str(char *, int);

consensus_trimed *nextCorrect(char **, unsigned int *, unsigned int *, unsigned int, \
    unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, float, unsigned int, unsigned int);

void free_consensus_trimed(consensus_trimed *);

char *poa_to_consensus(const struct seq_ *,  const int seq_count);
