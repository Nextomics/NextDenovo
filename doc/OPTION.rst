.. _parameterreference:

NextDenovo Parameter Reference
==============================

NextDenovo requires at least one read file (option: ``input_fofn``) as input, it works with gzip'd FASTA and FASTQ formats and uses a ``config file`` to pass options.

Input
-----

- ``input_fofn`` (one file one line)

  .. code:: console

    ls reads1.fasta reads2.fastq reads3.fasta.gz reads4.fastq.gz ... > input.fofn
- ``config file``

  A config file is a text file that contains a set of parameters (key=value pairs) to set runtime parameters for NextDenovo. The following is a typical config file, which is also located in ``doc/run.cfg``.
  
  .. code:: console

    [General]
    job_type = local
    job_prefix = nextDenovo
    task = all
    rewrite = yes # yes/no
    deltmp = yes 
    rerun = 3
    parallel_jobs = 20
    input_type = raw
    input_fofn = input.fofn
    workdir = 01_rundir

    [correct_option]
    read_cutoff = 1k
    seed_cutoff = 29999 
    blocksize = 2g
    pa_correction = 3
    seed_cutfiles = 10
    sort_options = -m 20g -t 10 -k 40 
    minimap2_options_raw = -x ava-ont -t 8 
    correction_options = -p 15

    [assemble_option]
    minimap2_options_cns = -x ava-ont -t 8 -k17 -w17 
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
    
    local, sge, pbs... (default: sge)

  .. option:: job_prefix = nextDenovo  

    prefix tag for jobs. (default: nextDenovo)
  .. option:: task = <all, correct, assemble>     

    task need to run. (default: all)
  .. option::  rewrite = no  

    overwrite existed directory [yes, no]. (default: no)
  .. option::  deltmp = yes      

    delete intermediate results. (default: yes)
  .. option::  rerun = 3         

    re-run unfinished jobs untill finished or reached ${rerun} loops, 0=no. (default: 3)
  .. option::  parallel_jobs = 10       

    number of tasks used to run in parallel. (default: 10)
  .. option::  input_type = raw         

    input reads type [raw, corrected]. (default: raw)
  .. option::  input_fofn = input.fofn  

    input file, one line one file. (**required**)
  .. option::  workdir = 01.workdir     

    work directory. (default: ./)
  .. option::  usetempdir = /tmp/test   

    temporary directory in compute nodes to avoid high IO wait. (default: None)
  .. option::  nodelist = avanode.list.fofn

    a list of hostnames of available nodes, one node one line, used with usetempdir for non-sge job_type.
  .. option:: cluster_options = auto

    a template to define the resource requirements for each job, which will pass to `DRMAA <https://github.com/pygridtools/drmaa-python/wiki/FAQ>`__ as the nativeSpecification field.

Correction options
##################

  .. option::  read_cutoff = 1k   

    filter reads with length < read_cutoff. (default: 1k)
  .. option::  seed_cutoff = 25k   

    minimum seed length. (**required**)
  .. option::  seed_cutfiles = 5    

    split seed reads into ${seed_cutfiles} subfiles. (default: ${pa_correction})
  .. option::  blocksize = 10g      

    block size for parallel running. (default: 10g)
  .. option::  pa_correction = 5        

    number of corrected tasks used to run in parallel, overwrite ${parallel_jobs} only for this step. (default: 15)
  .. option::  minimap2_options_raw = -x ava-ont -t 10  

    minimap2 options, used to find overlaps between raw reads and set PacBio/Nanopore read overlap, see :ref:`minimap2-nd <minimap2-nd>` for details. (**required**)
  .. option::  sort_options = -m 40g -t 10 -k 50 

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

  .. option::  minimap2_options_cns = -x ava-ont -t 8 -k17 -w17 

    minimap2 options, used to find overlaps between corrected reads. (default: -k17 -w17)
  .. option::  minimap2_options_map = -x map-ont

    minimap2 options, used to map reads back to the assembly.
  .. option::  nextgraph_options = -a 1

    nextgraph options, see :ref:`nextgraph <nextgraph>` for details.  