.. _utilities:

Utilities
=========

.. _seq_stat:

seq_stat
---------
  
  seq_stat can be used to perform some simple statistics (such as length distribution, total amount of data and sequencing depth) on the input data, and give the recommended minimum seed length.

  .. object:: INPUT

    .. :type: file
    .. :default: ``stdout``

    - read files list, one line one file
  .. object:: OUTPUT (stdout)

    - Read length histogram
    - Read length info.
    - Total Bases info.
    - Recommended minimum seed length
  .. object:: OPTIONS

    -f  skip reads with length shorter than this value [1kb].
    -g  estimated genome size [5Mb].
    -d  expected seed depth, used to be corrected [45].
    -a  disable automatic adjustment.

.. _seq_dump:

seq_dump
---------

  sql_dump is used to classify reads based on a given seed length threshold, and split and compress different categories to subfiles (bit format).
  
  .. object:: INPUT

    - read files list, one line one file
  .. object:: OUTPUT
   
    The output consists of four parts::

      - input.part*2bit (non-seed reads)
      - .input.part*idx (index of non-seed reads)
      - input.seed*2bit (seed reads)
      - .input.seed*idx (index of seed reads)
  .. object:: OPTIONS

    -f  minimum read length.
    -s  minimum seed length.
    -b  block size (Mb or Gb, < 16Gb).
    -n  number of seed subfiles in total.
    -d  output directory.

.. _seq_bit:

seq_bit
--------

  seq_bit can be used to compress fasta files to bit files or uncompress bit files to fasta files.
  
  .. object:: INPUT

    - one seq file.
  .. object:: OUTPUT (stdout)

    - sequences with fasta or bit format.

.. _minimap2-nd:

minimap2-nd
-----------

  minimap2-nd is a modified version of `minimap2 <https://github.com/lh3/minimap2>`__, which is used to find all overlaps between raw reads and dovetail overlaps between corrected seeds. Compared to minimap2, minimap-nd has four minor modifications::

    1. Add support for input files in bit format.
    2. Add a filter step for output.
    3. Compress output when output to a file.
    4. Add a re-align step for potential dovetail overlaps.

  .. object:: EXTRA OPTIONS

    --step <1,2,3>  preset options for NextDenovo, **[required]**.
    --minlen INT    min overlap length [500]
    --minmatch INT  min match length [100]
    --minide FLOAT  min identity [0.05]
    --mode <1,2>    re-align mode, 1:fast mode, low accuracy 2:slow mode, high accuracy [2]
    --kn INT        k-mer size (no larger than 28), used to re-align [17]
    --wn INT        minizer window size, used to re-align [10]
    --cn INT        do re-align for every INT reads [20]
    --maxhan1 INT   max over hang length, used to re-align [5000]
    --maxhan2 INT   max over hang length, used to filter contained reads [500]

.. _ovl_sort:

ovl_sort
---------

  ovl_sort is used to sort and remove redundancy overlaps by number of matches for a given seed.
  
  .. object:: INPUT

    - overlap files, one line one file. 
    - index file of seeds need to be sorted.
  .. object:: OUTPUT
  
    - sorted overlap file.
  .. object:: OPTIONS

    -i    index file of seeds need to be sorted **[required]**
    -m    set max total available buffer size, suffix K/M/G [40G]
    -t    number of threads to use [8]
    -k    max depth of each overlap, should <= average sequencing depth [40]
    -l    max over hang length to filter [300]
    -o    output file name **[required]**
    -d    temporary directory [$CWD]

.. _ovl_cvt:

ovl_cvt
--------

  ovl_cvt can be used to compress or uncompress overlap files.

  .. object:: INPUT

    - one overlap file
  .. object:: OUTPUT (stdout)

    - compressed or uncompressed overlaps
  .. object:: OPTIONS

    -m INT    conversion mode (0 for compress, 1 for uncompress)

.. _nextgraph:

nextgraph
---------

  NextGraph is used to construct a string graph with corrected reads. The main algorithms are similar to other mainstream assemblers except using a graph-based algorithm to identify chimeric nodes and a scoring-based strategy to identify incorrect edges. It can output an assembly in `Fasta <https://en.wikipedia.org/wiki/FASTA_format>`__, `GFA2 <https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md>`__, `GraphML <https://en.wikipedia.org/wiki/GraphML>`__, Path formats, or only statistical information (for quickly optimize parameters).

  .. object:: INPUT

    - read files list, one line one file.
    - overlap files list, one line one file.
  .. object:: OUTPUT

    - assembly statistical information.
    - assembly sequences.

  .. object:: OPTIONS

    -f FILE                        input seq list [required]
    -o FILE                        output file [stdout]
    -c                             disable pre-filter chimeric reads 
    -G                             retain potential chimeric edges 
    -k                             delete complex bubble paths 
    -A                             output alternative contigs 
    -a, --out_format INT           output format, 0=None, 1=fasta, 2=graphml, 3=gfa2, 4=path [1]
    -E, --out_ctg_len INT          min contig length for output [1000]
    -q, --out_spath_len INT        min short branch len for output, 0=disable [0]
    -i, --min_ide FLOAT            min identity of alignments [0.10]
    -I, --min_ide_ratio FLOAT      min test-to-best identity ratio [0.70]
    -S, --min_sco_ratio FLOAT      min test-to-best aligned length ratio [0.40]
    -r, --max_sco_ratio FLOAT      max test-to-best score ratio of a low quality edge [0.50]
    -M, --min_mat_ratio FLOAT      min test-to-best aligned matches ratio [0.90]
    -T, --min_depth_ratio FLOAT    min test-to-best depth ratio of an edge [0.60]
    -N, --min_node_count <1,2>     min valid nodes of a read [2]
    -u, --min_con_count <1,2>      min contained number to filter contained reads [2]
    -w, --min_edge_cov INT         min depth of an edge [3]
    -D, --bfs_depth INT            depth of BFS to identify chimeric nodes [2]
    -P, --bfs_depth_multi INT      max depth multiple of a node for BFS [2]
    -m, --min_depth_multi FLOAT    min depth multiple of a repeat node [1.50]
    -n, --max_depth_multi FLOAT    max depth multiple of a node [2000.00]
    -B, --bubble_len INT           max len of a bubble [500]
    -C, --cpath_len INT            max len of a compound path [20]
    -z, --zbranch_len INT          max len of a z branch [8]
    -l, --sbranch_len INT          max len of a short branch [15]
    -L, --sloop_len INT            max len of a short loop [10]
    -t, --max_hang_len INT         max over hang length of dovetails [500]
    -F, --fuzz_len INT             fuzz len for trans-reduction [1000]

.. _bam_sort:

bam_sort
---------
  
  bam_sort is used to sort bam files.
  
  .. object:: INPUT

     - bam file need to be sorted.
  .. object:: OUTPUT

    - sorted bam file.
    - index file.
  .. object:: OPTIONS

    -i         Write index file.
    -m INT     Set maximum memory per thread; suffix K/M/G recognized [1024M]
    -o FILE    Write final output to FILE rather than standard output.
    -T PREFIX  Write temporary files to PREFIX.nnnn.bam.
    
    -@ INT
      Number of additional threads to use [0]
