CC = gcc
CFLAGS = -Wall -O3
PROGS = seq_dump ovl_sort seq_stat ovl_cvt seq_bit bam_sort nextgraph
HTSDIR = ../lib/htslib/
NEXTGRPHDIR = ./ng

all:$(PROGS)

$(HTSDIR)/libhts.a:
	@+make -C $(HTSDIR) CPPFLAGS+=-fPIC

seq_dump: seq_dump.c ../lib/bseq.c
	$(CC) $(CFLAGS) -Wno-unused-function -o $@ $^ -lz

ovl_sort: ovl_sort.c thpool.c ../lib/ovl.c
	$(CC) $(CFLAGS) -std=c99 -pthread -o $@ $^

ovl_cvt: ovl_cvt.c ../lib/ovl.c
	$(CC) $(CFLAGS) -o $@ $^

seq_bit: seq_bit.c ../lib/bseq.c
	$(CC) $(CFLAGS) -Wno-unused-function -o $@ $^ -lz

seq_stat: seq_stat.c
	$(CC) $(CFLAGS) -std=c99 -o $@ $^ -lz

bam_sort: bam_sort.c ../lib/bsort.c $(HTSDIR)/libhts.a
		$(CC) $(CFLAGS) -o $@ $^ -I $(HTSDIR) -lpthread -lm -llzma -lbz2 -lz

# align: realign.c ../lib/ovl.c ../lib/bseq.c ../lib/align.c
# 	$(CC) $(CFLAGS) -Wno-unused-function -o $@ $^ -lz

nextgraph:
	make -C $(NEXTGRPHDIR)
	@+cp $(NEXTGRPHDIR)/$@ ./

clean:
	rm -f *.o $(PROGS)
	make clean -C $(NEXTGRPHDIR)
	@+make clean -C $(HTSDIR)
