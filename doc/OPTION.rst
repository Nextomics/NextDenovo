.. _parameterreference:

NextDenovo Parameter Reference
==============================

NextDenovo requires at least one read file (option: ``input_fofn``) as input, it works with gzip'd FASTA and FASTQ formats and uses a ``config file`` to pass options.

Input
-----

- ``input_fofn`` (one file one line)

  .. code-block:: shell

    ls reads1.fasta reads2.fastq reads3.fasta.gz reads4.fastq.gz ... > input.fofn
- ``config file``

  A config file is a text file that contains a set of parameters (key=value pairs) to set runtime parameters for NextDenovo. The following is a typical config file, which is also located in ``doc/run.cfg``.
  
  .. code-block:: bash

    [General]
    job_type = local
    job_prefix = nextDenovo
    task = all
    rewrite = yes
    deltmp = yes 
    parallel_jobs = 20
    input_type = raw
    read_type = clr # clr, ont, hifi
    input_fofn = input.fofn
    workdir = 01_rundir

    [correct_option]
    read_cutoff = 1k
    genome_size = 1g # estimated genome size
    sort_options = -m 20g -t 15
    minimap2_options_raw = -t 8
    pa_correction = 3
    correction_options = -p 15

    [assemble_option]
    minimap2_options_cns = -t 8 
    nextgraph_options = -a 1  

Output
------

- ``workdir/03.ctg_graph/nd.asm.fasta``

  Contigs with fasta format, the fasta header includes ID, type, length, node count, a consecutive lowercase region in the sequence implies a weak connection, and a low quality base is marked with a single lowercase base.
- ``workdir/03.ctg_graph/nd.asm.fasta.stat``

  Some basic statistical information (N10-N90, Total size et al.).

.. _options:

Options
-------

Global options
##############

  .. option:: job_type = sge           
    
    local, sge, pbs, lsf, slurm... (default: sge)

  .. option:: job_prefix = nextDenovo  

    prefix tag for jobs. (default: nextDenovo)
  .. option:: task = <all, correct, assemble>     

    task need to run, correct = only do the correction step, assemble = only do the assembly step (only work if ``input_type`` = corrected or ``read_type`` = hifi), all = correct + assemble. (default: all)
  .. option::  rewrite = no  

    overwrite existed directory [yes, no]. (default: no)
  .. option::  deltmp = yes      

    delete intermediate results. (default: yes)
  .. option::  rerun = 3         

    re-run unfinished jobs untill finished or reached ``rerun`` loops, 0=no. (default: 3)
  .. option::  parallel_jobs = 10       

    number of tasks used to run in parallel. (default: 10)
  .. option::  input_type = raw         

    input reads type [raw, corrected]. (default: raw)
  .. option::  input_fofn = input.fofn  

    input file, one line one file. (**required**)

.. _read_type:

  .. option::  read_type = {clr, hifi, ont}  

    reads type, clr=PacBio continuous long read, hifi=PacBio highly accurate long reads, ont=NanoPore 1D reads. (**required**)
  .. option::  workdir = 01.workdir     

    work directory. (default: ./)
  .. option::  usetempdir = /tmp/test   

    temporary directory in compute nodes to avoid high IO wait. (default: None)
  .. option::  nodelist = avanode.list.fofn

    a list of hostnames of available nodes, one node one line, used with usetempdir for non-sge job_type.
  .. option:: submit = auto   

    command to submit a job, auto = automatically set by `Paralleltask <https://github.com/moold/ParallelTask>`__.
  .. option:: kill = auto   

    command to kill a job, auto = automatically set by `Paralleltask <https://github.com/moold/ParallelTask>`__.
  .. option:: check_alive = auto   

    command to check a job status, auto = automatically set by `Paralleltask <https://github.com/moold/ParallelTask>`__.
  .. option:: job_id_regex = auto   

    the job-id-regex to parse the job id from the out of ``submit``, auto = automatically set by `Paralleltask <https://github.com/moold/ParallelTask>`__.
  .. option:: use_drmaa = no   

    use drmaa to submit and control jobs.

Correction options
##################

  .. option::  read_cutoff = 1k   

    filter reads with length < ``read_cutoff``. (default: 1k)

.. _genome_size:

  .. option::  genome_size = 1g   

    estimated genome size, suffix K/M/G recognized, used to calculate ``seed_cutoff``/``seed_cutfiles``/``blocksize`` and average depth, it can be omitted when manually setting ``seed_cutoff``.
  .. option::  seed_depth = 45   

    expected seed depth, used to calculate ``seed_cutoff``, co-use with ``genome_size``, you can try to set it 30-45 to get a better assembly result. (default: 45)
  .. option::  seed_cutoff = 0   

    minimum seed length, <=0 means calculate it automatically using :ref:`bin/seq_stat <seq_stat>`.
  .. option::  seed_cutfiles = 5    

    split seed reads into ``seed_cutfiles`` subfiles. (default: ``pa_correction``)
  .. option::  blocksize = 10g      

    block size for parallel running, split non-seed reads into small files, the maximum size of each file is ``blocksize``. (default: 10g)
  .. option::  pa_correction = 3        

    number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage, overwrite ``parallel_jobs`` only for this step. (default: 3)
  .. option::  minimap2_options_raw = -t 10  

    minimap2 options, used to find overlaps between raw reads, see :ref:`minimap2-nd <minimap2-nd>` for details.
  .. option::  sort_options = -m 40g -t 10 

    sort options, see :ref:`ovl_sort <ovl_sort>` for details.  
  .. option::  correction_options = -p 10 

    correction options, see following::

      -p, --process, set the number of processes used for correcting. (default: 10)
      -b, --blacklist, disable the filter step and increase more corrected data.
      -s, --split, split the corrected seed with un-corrected regions. (default: False)
      -fast, 0.5-1 times faster mode with a little lower accuracy. (default: False)
      -dbuf, disable caching 2bit files and reduce ~TOTAL_INPUT_BASES/4 bytes of memory usage. (default:False)
      -max_lq_length, maximum length of a continuous low quality region in a corrected seed, larger max_lq_length will produce more corrected data with lower accuracy. (default: auto [pb/1k, ont/10k])

Assembly options
##################

  .. option::  minimap2_options_cns = -t 8 -k17 -w17 

    minimap2 options, used to find overlaps between corrected reads.
  .. option::  minimap2_options_map = -t 10

    minimap2 options, used to map reads back to the assembly.
  .. option::  nextgraph_options = -a 1

    nextgraph options, see :ref:`nextgraph <nextgraph>` for details.  
