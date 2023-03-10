#ifndef _COM_H_
#define _COM_H_
#include <stdint.h>
#include "kit.h"

typedef struct node_ {// size 32
	UINT_T l; // lable, >= 2 bit
	UINT_T id; // in degree
	UINT_T idm; // max in degree
	UINT_T od; // out degree
	UINT_T odm; // max out degree
	UINTL_T name; // node name
	UINTL_T *ie; // index of in edge
	UINTL_T *oe; // index of out edge
} node;

typedef struct edge_ {//size 24
	uint16_t l; // lable
	uint16_t ide; // aln identity
	UINTL_T in; // index of in node
	UINTL_T ou; // index of out node
	uint32_t tc:8; // triangle count
	uint32_t ie:24; // in end
	uint32_t oe; // out end
	uint32_t len; // edge length
	uint32_t sco; // aln length
} edge;

typedef struct graph_ {
	node *nd; // node
	edge *eg; // edge
	UINTL_T ni; // node index
	UINTL_T nm; // node max
	UINTL_T ei; // edge index
	UINTL_T em; // edge max
} graph;

#endif
