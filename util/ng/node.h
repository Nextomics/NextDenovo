#ifndef _NODE_H_
#define _NODE_H_
#include "common.h"

UINTL_T get_reversed_node(graph *g, const UINTL_T n);
UINTL_T add_node(graph *g, const UINTL_T name, const int len);
UINT_T get_validly_oe(graph *g, const UINTL_T n, const UINT_T i);
UINT_T get_validly_ie(graph *g, const UINTL_T n, const UINT_T i);
void rm_oe(graph *g, const UINTL_T n, const UINTL_T e);
void rm_ie(graph *g, const UINTL_T n, const UINTL_T e);
void rm_node(graph *g, const UINTL_T n);
int add_oe(graph *g, const UINTL_T n, const UINTL_T e);
int add_ie(graph *g, const UINTL_T n, const UINTL_T e);
void clean_node_lable(graph *g);
UINTL_T stat_valid_node(graph *g);
int check_node_lable(graph *g, const UINTL_T n, const uint16_t l);

#endif
