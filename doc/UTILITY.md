## seq_stat
* **DESCRIPTION**    
seq_stat can be used to perform some simple statistics (such as length distribution, total amount of data and sequencing depth) on the input data, and give the recommended minimum seed length.  
* **INPUT**  
    `- read files list, one line one file`  
* **OUTPUT** (stdout)
<pre>
    - Read length histogram
    - Read length info.
    - Total Bases info.
    - Recommended minimum seed length
</pre>
* **OPTIONS**
<pre>
    -f skip reads with length shorter than this value, default 1kb
    -g estimated genome size, default: 5Mb
    -d expected seed depth, used to be corrected, default: 45
</pre>

## seq_dump  
* **DESCRIPTION**   
sql_dump is used to classify reads based on a given seed length threshold, and split and compress different categories to subfiles (bit format).    
* **INPUT**  
    `- read files list, one line one file`  
* **OUTPUT**   
The output consists of four parts:
<pre>
    - input.part*2bit (non-seed reads)
    - .input.part*idx (index of non-seed reads)
    - input.seed*2bit (seed reads)
    - .input.seed*idx (index of seed reads)
</pre>
* **OPTIONS**
<pre>
    -f minimum read length
    -s minimum seed length
    -b block size (Mb or Gb, < 16Gb)
    -n number of seed subfiles in total
    -d output directory
</pre>

## seq_bit
* **DESCRIPTION**  
seq_stat can be used to compress fasta files to bit files or uncompress bit files to fasta files. 
* **INPUT**   
    `- one seq file`  
* **OUTPUT** (stdout)  
    `- sequences with fasta or bit format`  

## minimap2-nd
* **DESCRIPTION**  
minimap2-nd is a modified version of [minimap2](https://github.com/lh3/minimap2), which is used to find all overlaps between raw reads and dovetail overlaps between corrected seeds. Compared to minimap2, minimap-nd has four minor modifications: 
    1. Add support for input files in bit format.
    2. Add a filter step for output.
    3. Compress output when output to a file.
    4. Add a re-align step for potential dovetail overlaps.
* **EXTRA OPTIONS**  
<pre>
    --step [1|2]   preset options for NextDenovo [required]
    --minlen INT   min overlap length [500]
    --minmatch INT min match length [100]
    --minide FLOAT min identity [0.05]
    --norealign    output the initial alignments, much faster, less accurate
    --kn INT       k-mer size (no larger than 28), used to re-align [17]
    --wn INT       minizer window size, used to re-align [10]
    --maxhan1 INT  max over hang length, used to re-align [5000]
    --maxhan2 INT  max over hang length, used to filter contained reads [500]
</pre>

## ovl_sort
* **DESCRIPTION**    
ovl_sort is used to sort and remove redundancy overlaps by number of matches for a given seed.
* **INPUT** 
<pre>
    - overlap files, one line one file. 
    - index file of seeds need to be sorted.
</pre>
* **OUTPUT**   
    `- sorted overlap file`
* **OPTIONS**   
<pre>
    -i  index file of seeds need to be sorted [required]
    -m  set max available buffer size, suffix K/M/G [40G]
    -t  number of threads to use [8]
    -k  max depth of each overlap [40]
    -l  max over hang length to filter [300]
    -o  output file name [required]
    -d  temporary directory [$CWD]
</pre>

## ovl_cvt
* **DESCRIPTION**    
ovl_cvt can be used to compress or uncompress overlap files.
* **INPUT**    
    `- one overlap file`
* **OUTPUT** (stdout)  
    `- compressed or uncompressed overlaps`
* **OPTIONS**   
    `-m INT    conversion mode (0 for compress, 1 for uncompress)`
## nextgraph
* **DESCRIPTION**  
NextGraph is used to construct a string graph with corrected reads. The main algorithms are similar to other mainstream assemblers except using a graph-based algorithm to identify chimeric nodes and a scoring-based strategy to identify incorrect edges. It can output an assembly in [Fasta](https://en.wikipedia.org/wiki/FASTA_format), [GFA2](https://github.com/GFA-spec/GFA-spec/blob/master/GFA2.md), [GraphML](https://en.wikipedia.org/wiki/GraphML), Path formats, or only statistical information (for quickly optimize parameters).
* **INPUT** 
<pre>
    - read files list, one line one file
    - overlap files list, one line one file
</pre>
* **OUTPUT**
<pre>
    - assembly statistical information
    - assembly sequences
</pre>
* **OPTIONS**  
<pre>
    -f [FILE]                       input seq list [required]
    -o [FILE]                       output file [stdout]
    -s                              disable sort out-edges by length 
    -c                              disable pre-filter chimeric reads 
    -k                              delete compound branch pathes 
    -R                              disable re-filter contained reads 
    -a --out_seq [INT]              output format, 0=None, 1=fasta, 2=graphml, 3=gfa2, 4=path [1]
    -q --keep_short_path_len [INT]  min short branch length for output, 0=disable [0]
    -Q --keep_z_path_len [INT]      min z branch length for output, 0=disable [0]
    -F --fuzz [INT]                 fuzz len for trans-reduction [1000]
    -D --ext_node_count [INT]       depth of BFS to identify chimeric nodes [2]
    -P --ext_depth_multi [INT]      max depth multiple of a node for BFS [2]
    -E --min_ctg_len [INT]          min contig len for output [1000]
    -i --min_ide [FLOAT]            min identity of alignments [0.10]
    -I --min_ide_ratio [FLOAT]      min test-to-best identity ratio [0.70]
    -S --min_sco_ratio [FLOAT]      min test-to-best aligned len ratio [0.40]
    -N --min_node_count [1,2]       min valid ends of a read [2]
    -u --min_con_count [1,2]        min contained number to filter reads [2]
    -r --max_sco_ratio [FLOAT]      max high score ratio [0.50]
    -d --max_aln_depth [INT]        max aligned depth [500]
    -m --min_depth_multi [FLOAT]    min depth multiple of a repeat node [1.50]
    -n --max_depth_multi [FLOAT]    max depth multiple of a node [2000.00]
    -B --short_bubble_len [INT]     max len of a bubble [40]
    -e --end_loop_len [INT]         max len of a terminal loop [50]
    -C --comp_path_len [INT]        max len of a compound path [20]
    -z --zclip_len [INT]            max len of a z branch [8]
    -l --short_branch_len [INT]     max len of a short branch [15]
    -L --short_loop_len [INT]       max len of a short loop [10]
    -p --max_hang_plen [INT]        max over hang length of potential dovetails [3000]
    -t --max_hang_tlen [INT]        max over hang length of dovetails [500]
</pre>




