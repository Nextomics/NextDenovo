.. _chm13_120x_ont:

.. title:: CHM13hTERT human cell line with 120x Oxford Nanopore data

Assessment of the `CHM13 <https://github.com/nanopore-wgs-consortium/CHM13>`__ genome (120X NanoPore data) assemblies using NextDenovo, Canu, Flye, Shasta
----------------------------------------------------------------------------------------------------------------------------------------------------------

1. **Download reads**
  
  .. code:: console

    wget https://s3.amazonaws.com/nanopore-human-wgs/chm13/nanopore/rel3/rel3.fastq.gz

2. **Prepare input file (input.fofn)**

  .. code:: console

    ls rel3.fastq.gz > input.fofn

3. **Calculate the recommended minimum seed length**
   
  .. code:: console

    bin/seq_stat -f 20k -g 3.1g input.fofn > input.fofn.stat

  The following is the partial content of file input.fofn.stat, and it shows the recommended minimum seed length is ``71477`` bp at the last line.

  ::

    [Read length stat]
    Types            Count (#) Length (bp)
    N10                 110235  202350
    N20                 295841  140527
    N30                 550659  105149
    N40                 885683   81094
    N50                1313448   64445
    N60                1843885   52578
    N70                2489549   43370
    N80                3274007   35502
    N90                4246284   28093

    Types               Count (#)           Bases (bp)  Depth (X)
    Raw                  28449385         367231282800     118.46
    Filtered             22918160          59455023122      19.18
    Clean                 5531225         307776259678      99.28

    *Suggested seed_cutoff (genome size: 3100000000, expected seed depth: 45) : 71477 bp

4. **Prepare config file (run.cfg)**
   
   .. code:: console

    [General]
    job_type = sge
    job_prefix = nextDenovo
    task = all # 'all', 'correct', 'assemble'
    rewrite = yes # yes/no
    deltmp = yes
    rerun = 3
    parallel_jobs = 25
    input_type = raw
    input_fofn = ./input.fofn
    workdir = ./chm13_asm

    [correct_option]
    read_cutoff = 20k
    seed_cutoff = 71477
    blocksize = 5g
    pa_correction = 5
    seed_cutfiles = 5
    sort_options = -m 150g -t 30 -k 50
    minimap2_options_raw = -x ava-ont -t 8
    correction_options = -p 30

    [assemble_option]
    minimap2_options_cns = -x ava-ont -t 8 -k17 -w17
    nextgraph_options = -a 1

5. **Run**

  .. code:: console

    nohup nextDenovo run.cfg &

6. **Get result**
   
  - Final corrected reads file (use the ``-b`` parameter to get more corrected reads)::
  
      chm13_asm/02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta
  
  - Final assembly result::  
  
      chm13_asm/03.ctg_graph/nd.asm.fasta

  The folowing is the assembly statistics::

    Type           Length (bp)            Count (#)
    N10            179297054                   2
    N20            169128386                   3
    N30            131652719                   6
    N40            120761272                   8
    N50            106090521                  10
    N60             95206689                  13
    N70             80513393                  16
    N80             59725892                  21
    N90             39058727                  27

    Min.               84432                   -
    Max.           237405279                   -
    Ave.            35344197                   -
    Total         2898224197                  82

7. **Download reference**   
  
  .. code:: console
  
    wget https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v0.7.fasta.gz
    gzip -d chm13.draft_v0.7.fasta.gz

8. **Run Quast v5.0.2**
  
  .. code:: console

    quast.py --eukaryote --large --min-identity 80 --threads 30 -r ./chm13.draft_v0.7.fasta --fragmented nd.asm.fasta

  .. object:: Quast result

  +--------------------------------+---------------+------------------+----------------+---------------+
  |                                | NextDenovo    | Canu             | Flye           | Shasta        |
  +================================+===============+==================+================+===============+
  | # contigs                      | 82            | 1223             | 472            | 297           |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | Largest contig                 | 237405279     | 139909728        | 132009996      | 130803838     |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | Total length                   | 2898224197    | 2991947723       | 2920201070     | 2823384269    |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | # **misassemblies**            | 1227          | 6396             | 3230           | 187           |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | # misassembled contigs         | 61            | 875              | 193            | 78            |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | Misassembled contigs length    | 2740877545    | 2458710426       | 2440399207     | 1351075153    |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | # **local misassemblies**      | 433           | 1164             | 981            | 129           |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | # possible TEs                 | 42            | 160              | 96             | 14            |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | # unaligned mis. contigs       | 11            | 73               | 17             | 0             |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | # unaligned contigs            | 0 + 64 part   | 168 + 248 part   | 8 + 135 part   | 0 + 37 part   |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | Unaligned length               | 22021119      | 30076945         | 14583673       | 393547        |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | Genome fraction (%)            | 97.421        | 98.391           | 97.392         | 96.149        |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | Duplication ratio              | 1.007         | 1.027            | 1.018          | 1.002         |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | # **mismatches per 100 kbp**   | 29.43         | 77.26            | 74.04          | 15.56         |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | # **indels per 100 kbp**       | 170.98        | 327.08           | 447.97         | 141.25        |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | Largest alignment              | 111497488     | 104447985        | 111814657      | 111679369     |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | Total aligned length           | 2865321418    | 2943726417       | 2894073152     | 2821352191    |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | **N50**                        | 106090521     | 77964612         | 70319350       | 58111632      |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | NG50                           | 106090521     | 77964612         | 70319350       | 58088067      |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | L50                            | 10            | 15               | 16             | 17            |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | LG50                           | 10            | 15               | 16             | 18            |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | **NA50**                       | 57779597      | 47440498         | 46858921       | 47392260      |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | NGA50                          | 57779597      | 47440498         | 46546094       | 44539326      |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | LA50                           | 18            | 21               | 19             | 19            |
  +--------------------------------+---------------+------------------+----------------+---------------+
  | LGA50                          | 18            | 21               | 20             | 20            |
  +--------------------------------+---------------+------------------+----------------+---------------+

  .. note:: The results of Canu, Flye and Shasta are copied from `here <https://github.com/human-pangenomics/assembly-analysis>`__, the complete result of NextDenovo can be seen from :download:`here <./TEST2.pdf>`.