.. title:: Tutorial

Assemble the genome of HG002_NA24385_son using NextDenovo
-----------------------------------------------------------

1. **Download reads**

  .. code-block:: shell

    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V2.3.4_2019-06-26/ultra-long-ont.fastq.gz

2. **Prepare input file (input.fofn)**

   .. code-block:: shell

    ls ultra-long-ont.fastq.gz > input.fofn

3. **Prepare config file (run.cfg)** 

  .. code-block:: shell

    [General]
    job_type = sge # here we use SGE to manage jobs
    job_prefix = nextDenovo
    task = all
    rewrite = yes
    deltmp = yes 
    parallel_jobs = 22
    input_type = raw
    read_type = ont # clr, ont, hifi
    input_fofn = input.fofn
    workdir = HG002_NA24385_son_assemble

    [correct_option]
    read_cutoff = 1k
    genome_size = 3g # estimated genome size
    sort_options = -m 50g -t 30
    minimap2_options_raw = -t 8
    pa_correction = 5
    correction_options = -p 30

    [assemble_option]
    minimap2_options_cns = -t 8 
    nextgraph_options = -a 1 

4. **Run**

  .. code-block:: shell
    
    nohup nextDenovo run.cfg &

5. **Get result**

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