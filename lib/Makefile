CC=gcc
INCLUDES =
CFLAGS = -std=c99 -O3 -Wall -Wno-unused-function
FPIC = -fPIC -shared
HTSDIR = ./htslib/
LIB_OBJECT = nextcorrect.so ovlseq.so ctg_cns.so ckit.so
CTG_CNS_CFLAGS = -D_SVID_SOURCE -DLGS_CORRECT -I $(HTSDIR) -lm -lpthread -lz $(HTSLIB_CFLAGS)

all: $(LIB_OBJECT)

$(HTSDIR)/libhts.a:
	cd  $(HTSDIR) && ./configure  --disable-bz2 --disable-lzma
	@+make -C $(HTSDIR) CPPFLAGS+=-fPIC $(INCLUDES)

nextcorrect.so: align.c dag.c nextcorrect.c
	$(CC) $(CFLAGS) $(FPIC) $(INCLUDES) -o $@ $^ 

ovlseq.so: ovlseq.c bseq.c ovl.c index.c
	$(CC) $(CFLAGS) $(FPIC) $(INCLUDES) -o $@ $^ -lz

nextcorrect: align.c dag.c nextcorrect.c
	$(CC) $(CFLAGS) -o $@ $^

ctg_cns.so: bseq.c align.c bsort.c dag.c ctg_cns.c $(HTSDIR)/libhts.a
	$(CC) $(CFLAGS) $(FPIC) $(INCLUDES) -o $@ $^ $(CTG_CNS_CFLAGS) 

ckit.so: kit.c bseq.c
	$(CC) $(CFLAGS) $(FPIC) $(INCLUDES) -o $@ $^ -lz

ckit: kit.c bseq.c
	$(CC) $(CFLAGS) -o $@ $^ -lz

ctg_cns: bseq.c align.c dag.c bsort.c ctg_cns.c $(HTSDIR)/libhts.a
	$(CC) $(CFLAGS) -o $@ $^ $(CTG_CNS_CFLAGS)

ovlseq: ovlseq.c bseq.c ovl.c index.c
	$(CC) $(CFLAGS) -lz -Wno-unused-function -o $@ $^

clean:
	rm -rf *.o *.pyc $(LIB_OBJECT) nextcorrect ovlseq nextcorrect ctg_cns kit __pycache__
	@+make clean -C $(HTSDIR)
