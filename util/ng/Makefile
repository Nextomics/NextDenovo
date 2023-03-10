CC=			gcc
#GENOME_SIZE=0, 1, 2 for small, middle, large genome size
CPPFLAGS=   -DGENOME_SIZE=0
CFLAGS=		-Wall -O3  -Wno-unused-function $(CPPFLAGS)
INCLUDES=	-I.
OBJS=		asg.o ctg.o edge.o kit.o main.o node.o opt.o out.o
PROG=		nextgraph
LIBS=		../../lib/bseq.c ../../lib/index.c ../../lib/ovl.c -lm -lz

.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

nextgraph:$(OBJS) main.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)
clean:
	rm -fr *.o a.out $(PROG) *~
