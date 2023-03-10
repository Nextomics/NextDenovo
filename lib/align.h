#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>


#ifdef LGS_CORRECT
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
#else
    typedef struct {
        unsigned int shift;
        unsigned int aln_len;
        unsigned int aln_t_s;
        unsigned int aln_t_e;//not include
        unsigned int aln_t_len;
        unsigned int aln_q_len;//not include
        char * q_aln_str;
        char * t_aln_str;
    } alignment;
#endif


typedef struct {
    unsigned int aln_len; //align len
    unsigned int aln_mlen;//match len
    unsigned int aln_t_s;
    unsigned int aln_t_e;//not include
    unsigned int aln_q_s;
    unsigned int aln_q_e;//not include
} alignpos;

void malloc_vd(int **V, uint8_t ***D, uint64_t max_mem_d);
void clean_V(int *V, int max_mem_d);
void destory_vd(int *V, uint8_t **D);
void revcomp_bseq(char *str, int len);
void str_tolower(char *p);
void str_toupper(char *p);
void ide(const char *query_seq, int q_len, const char *target_seq, int t_len,
	int *V, uint8_t **D, int max_d, int band_size, int *mlen, int *blen);
void alnpos(const char *query_seq, int q_len, const char *target_seq, int t_len,
    int *V, uint8_t **D, int max_d, int band_size, alignpos *aln);
void extend_fwd(const char *query_seq, int q_len, const char *target_seq, int t_len,
	int *V, uint8_t **D, int max_d, int band_size, float d_factor, int *bstx, int *bsty);
void extend_rev(const char *query_seq, int q_len, const char *target_seq, int t_len,
	int *V, uint8_t **D, int max_d, int band_size, float d_factor, int *bstx, int *bsty);
void align_hq(char *query_seq, int q_len, char *target_seq, int t_len, 
	alignment *align_rtn, int *V, uint8_t **D);
void align(char *query_seq, int q_len, char *target_seq, int t_len, 
	alignment *align_rtn, int *V, uint8_t **D);
