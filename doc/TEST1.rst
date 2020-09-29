.. title:: Tutorial

Assemble the genome of HG002_NA24385_son using NextDenovo
-----------------------------------------------------------

1. **Download reads**

  .. code:: console

    wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/final/ultra-long-ont.fastq.gz

2. **Prepare input file (input.fofn)**

   .. code:: console

    ls ultra-long-ont.fastq.gz > input.fofn

3. **Calculate the recommended minimum seed length**
  
  .. code:: console
  
    bin/seq_stat -f 1k -g 3g -d 45 input.fofn > input.fofn.stat

  The following is the partial content of file input.fofn.stat, and it shows the recommended minimum seed length is ``25799`` bp at the last line.

  ::

    [Read length stat]
    Types            Count (#) Length (bp)
    N10                  56450  231137
    N20                 161886  143469
    N30                 324002   95895
    N40                 558854   68596
    N50                 876274   51703
    N60                1297503   38714
    N70                1872389   27609
    N80                2723656   17503
    N90                4236082    8503

    Types               Count (#)           Bases (bp)  Depth (X)
    Raw                  14063218         190233179420      63.41
    Filtered              4221019           2025702592       0.68
    Clean                 9842199         188207476828      62.74

    *Suggested seed_cutoff (genome size: 3000000000, expected seed depth: 45) : 25799 bp

4. **Prepare config file (run.cfg)** 

  .. code:: console

    [General]
    job_type = sge # here we use SGE to manage jobs
    job_prefix = nextDenovo
    task = all # 'all', 'correct', 'assemble'
    rewrite = yes # yes/no
    deltmp = yes
    rerun = 3
    parallel_jobs = 22
    input_type = raw
    input_fofn = ./input.fofn # input file
    workdir = ./HG002_NA24385_son_assemble

    [correct_option]
    read_cuoff = 1k
    seed_cutoff = 25799 # the recommended minimum seed length
    blocksize = 5g
    pa_correction = 5
    seed_cutfiles = 5
    sort_options = -m 50g -t 30 -k 50
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
  
      HG002_NA24385_son_assemble/02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta

  - Final assembly result::
  
      HG002_NA24385_son_assemble/03.ctg_graph/nd.asm.fasta

  you can get some basic statistical information from file ``HG002_NA24385_son_assemble/03.ctg_graph/nd.asm.fasta.stat``, the folowing is the assembly statistics with default parameters::   

    Type           Length (bp)            Count (#)
    N10            168924870                   2
    N20            127260938                   4
    N30             94622851                   7
    N40             85456034                  10
    N50             79737202                  13
    N60             69943198                  17
    N70             58504138                  21
    N80             40548231                  27
    N90             19732879                  36

    Min.               82439                   -
    Max.           220056807                   -
    Ave.            24389616                   -
    Total         2877974703                 118

  .. note:: This result will have some minor changes with the version upgrade.