#ifndef _CTG_H_
#define _CTG_H_
#include "common.h"

typedef struct {
	uint32_t lq:2;
	uint32_t ide:10; // identity
	uint32_t ort:10; // matched out_len/Max ratio
	uint32_t irt:10; // matched in_len/Max ratio
	uint32_t l; // strand, 0=+, 1=-
	uint32_t s; // start, UINT32_MAX == reads end
	uint32_t e; // end
	UINTL_T name; // reads name
	UINTL_T nidex; // node index
} ctg__;

typedef struct {
	ctg__ *nd;
	// 0:unknown; 1:linear; 2:loop; 3:od==0; 4:id>1 && od>1;
	// {"unknown", "linear", "loop", "breakpoint", "junction"}
	uint32_t type; 
	uint32_t len;
	uint32_t i;
	uint32_t im;
} ctg_;

typedef struct {
	ctg_ *ctg;
	uint32_t i;
	uint32_t im;
} ctgs;

void init_ctgs(ctgs *s, const uint32_t len);
void destroy_ctgs(ctgs *s);
void generate_ctg(graph *g, ctgs *s, const uint16_t l);

#endif
