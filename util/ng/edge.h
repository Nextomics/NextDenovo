#ifndef _EDGE_H_
#define _EDGE_H_
#include "common.h"
#include "opt.h"
#include "../../lib/ovl.h"

#define MFLAG_FIT	0x1
#define MFLAG_REP1	0x2
#define MFLAG_REP2	0x4
#define MFLAG_TR	0x8
#define MFLAG_LS	0x10
#define MFLAG_HS	0x20
#define MFLAG_BS	0x40
#define MFLAG_CN	0x80
#define MFLAG_CC	0x100
#define MFLAG_LQ	0x200
#define MFLAG_HI	0x400
#define MFLAG_LI	0x800

#define MFLAG_IL	0x2000
#define MFLAG_OL	0x4000
#define MFLAG_TT	0x8000 //temp flag

void reallocate_edge(graph *s, const uint32_t len);
int check_valid_edge(edge *e, opt *p, ovlinfo *loli, ovlinfo *roli); //TODO: del contained edge
uint32_t check_exited_edge(graph *g, const UINTL_T in, const UINTL_T ou);// only used in update_graph
UINTL_T rp_exited_edge(graph *g, const UINTL_T in, const UINTL_T ou);// replace existed edge
UINTL_T get_reversed_edge(graph *g, const UINTL_T n);
UINTL_T get_edge(graph *g, const UINTL_T in, const UINTL_T ou);
uint32_t get_oreadlen(graph *g, UINTL_T n);
uint32_t get_ireadlen(graph *g, UINTL_T n);
UINTL_T stat_valid_edge(graph *g);
/*uint32_t add_edge(graph *g, const uint32_t in, const uint32_t ou, const uint32_t ie,\
	const uint32_t oe, const uint32_t sco, const uint16_t ide, const uint32_t len, const uint32_t l);
void rm_edge(graph *g, const uint32_t n);*/

#endif
