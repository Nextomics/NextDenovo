CC = gcc
CFLAGS = -Wall -O3 -Wno-unused-function 
INCLUDES=
PROGS = seq_dump ovl_sort seq_stat ovl_cvt seq_bit bam_sort nextgraph
HTSDIR = ../lib/htslib/
NEXTGRPHDIR = ./ng

all:$(PROGS)

$(HTSDIR)/libhts.a:
	cd  $(HTSDIR) && ./configure --without-curses --disable-bz2 --disable-lzma
	@+make -C $(HTSDIR) CPPFLAGS+=-fPIC $(INCLUDES)

seq_dump: seq_dump.c ../lib/bseq.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ -lz

ovl_sort: ovl_sort.c thpool.c ../lib/ovl.c
	$(CC) $(CFLAGS) -std=c99 -pthread $(INCLUDES) -o $@ $^

ovl_cvt: ovl_cvt.c ../lib/ovl.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

seq_bit: seq_bit.c ../lib/bseq.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ -lz

seq_stat: seq_stat.c
	$(CC) $(CFLAGS) -std=c99 $(INCLUDES) -o $@ $^ -lz

bam_sort: bam_sort.c ../lib/bsort.c $(HTSDIR)/libhts.a
		$(CC) $(CFLAGS) -o $@ $^ -I $(HTSDIR) $(INCLUDES) -lpthread -lm -lz

# align: realign.c ../lib/ovl.c ../lib/bseq.c ../lib/align.c
# 	$(CC) $(CFLAGS) -Wno-unused-function -o $@ $^ -lz

nextgraph:
	make -C $(NEXTGRPHDIR)
	@+cp $(NEXTGRPHDIR)/$@ ./

clean:
	rm -f *.o $(PROGS)
	make clean -C $(NEXTGRPHDIR)
	@+make clean -C $(HTSDIR)
