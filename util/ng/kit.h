#ifndef _KIT_H_
#define _KIT_H_
#include <inttypes.h>

// #define GENOME_SIZE 1

#if GENOME_SIZE == 1
	#define UINT_T uint32_t
	#define UINT_FORMAT PRIu32
	#define UINTL_T uint32_t
	#define UINTL_FORMAT PRIu32
#elif GENOME_SIZE == 2
	#define UINT_T uint32_t
	#define UINT_FORMAT PRIu32
	#define UINTL_T uint64_t
	#define UINTL_FORMAT PRIu64
#else
	#define UINT_T uint16_t
	#define UINT_FORMAT PRIu16
	#define UINTL_T uint32_t
	#define UINTL_FORMAT PRIu32
#endif

#define MAX_INIT_N 65536
#define RANDOM_COUNT 10000000
#define SWAP(x, y, T) do { T SWAP = x; x = y; y = SWAP; } while (0)

//#define RESTRICT
#define RESTRICT_GS 3500000000L

//#define P_VERSION

typedef struct {
	UINTL_T *n;
	int64_t i;
	int64_t im;
} ph;

typedef struct {
	//uint64_t d[2]; // raw depth, dovetail depth
	UINTL_T l[2]; // low count
	UINTL_T n[2]; // record count
	UINTL_T h[2]; // high record count
} det;

double cputime();
long peakrss();
double realtime();
void plog(const int level, char *format, ...);
void init_sn(ph *s, const int len);
void reallocate_sn(ph *s, const int len);
void destroy_sn(ph *s);
int version();
int quick_select(int *r, int s, int e, int k);

#endif
