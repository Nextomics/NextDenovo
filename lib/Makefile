CC=gcc
CCC=c++
CFLAGS = -std=c99 -O3 -Wall -s
FPIC = -fPIC -shared
HTSDIR = ./htslib/
EXE=nextcorrect
LIB_OBJECT = nextcorrect.so ovlseq.so ctg_cns.so ckit.so
CTG_CNS_CFLAGS = -D_SVID_SOURCE -DLGS_CORRECT -I $(HTSDIR) -lm -lbz2 -llzma -lpthread -lz

all: $(LIB_OBJECT)

$(HTSDIR)/libhts.a:
	@+make -C $(HTSDIR) CPPFLAGS+=-fPIC

nextcorrect.so: align.c dag.c nextcorrect.c
	$(CC) $(CFLAGS) $(FPIC) -o $@ $^

ovlseq.so: ovlseq.c bseq.c ovl.c index.c
	$(CC) $(CFLAGS) $(FPIC) -Wno-unused-function -o $@ $^ -lz

nextcorrect: align.c dag.c nextcorrect.c
	$(CC) $(CFLAGS) -o $@ $^

ctg_cns.so: bseq.c align.c bsort.c dag.c ctg_cns.c $(HTSDIR)/libhts.a
	$(CC) $(CFLAGS) $(FPIC) -Wno-unused-function -o $@ $^ $(CTG_CNS_CFLAGS)

ckit.so: kit.c bseq.c
	$(CC) $(CFLAGS) $(FPIC) -Wno-unused-function -o $@ $^ -lz

ckit: kit.c bseq.c
	$(CC) $(CFLAGS) -Wno-unused-function -o $@ $^ -lz

ctg_cns: bseq.c align.c dag.c bsort.c ctg_cns.c $(HTSDIR)/libhts.a
	$(CC) $(CFLAGS) -Wno-unused-function -o $@ $^ $(CTG_CNS_CFLAGS)

ovlseq: ovlseq.c bseq.c ovl.c index.c
	$(CC) $(CFLAGS) -lz -Wno-unused-function -o $@ $^

clean:
	rm -rf *.o *.pyc $(LIB_OBJECT) $(EXE) ovlseq nextcorrect ctg_cns kit __pycache__
	@+make clean -C $(HTSDIR)
