.. _faq:

Frequently Asked Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
  :local:

.. _how-to-optimize-parallel-computing-parameters:

How to optimize parallel computing parameters?
----------------------------------------------
The main parallel computing parameters include ``parallel_jobs``, ``pa_correction``, ``-t`` in `minimap2_options_raw`, `minimap2_options_cns`, `minimap2_options_map` and ``-p`` in `correction_options`. Each ``parallel_jobs`` subjob requires 32~64 gb RAM depending on the read length, each ``pa_correction`` subjob requires ~TOTAL_INPUT_BASES/4 bytes RAM.  

1. For an assembly on a local computer with ``P`` cores and ``M`` gb memory. A typical configuration file can be set like this (not the best, but better than the default):
	
  .. code-block:: shell

    [General]
    job_type = local 
    parallel_jobs = M/64 #here, 64 can optimize to 32~64
    ...

    [correct_option]
    pa_correction = M/(TOTAL_INPUT_BASES * 1.2/4)
    sort_options = -m TOTAL_INPUT_BASES * 1.2/4g -t P/pa_correction
    correction_options = -p P/pa_correction
    minimap2_options_raw = -t P/parallel_jobs
    ...

    [assemble_option]
    minimap2_options_cns = -t P/parallel_jobs 
    ...

2. For an assembly on a computer cluster with ``N`` computer nodes and each computer node has ``P`` cores and ``M`` gb memory. A typical configuration file can be set like this (not the best, but better than the default):

  .. code-block:: shell

  	let parallel_jobs_local = M/64 #here, 64 can optimize to 32~64
  	let pa_correction_local = M/(TOTAL_INPUT_BASES * 1.2/4)

  .. code-block:: shell

    [General]
    job_type = sge 
    parallel_jobs = parallel_jobs_local * N 
    ...

    [correct_option]
    pa_correction = pa_correction_local * N
    sort_options = -m TOTAL_INPUT_BASES * 1.2/4g -t P/pa_correction_local
    correction_options = -p P/pa_correction_local
    minimap2_options_raw = -t P/parallel_jobs_local
    ...

    [assemble_option]
    minimap2_options_cns = -t P/parallel_jobs_local 
    ...

What's the difference between ``nd.asm.p.fasta`` and the final assembly result ``nd.asm.fasta``?
------------------------------------------------------------------------------------------------
In theroy, ``nd.asm.p.fasta`` contains more structural & base errors than ``nd.asm.fasta``, you can chose ``nd.asm.p.fasta`` as the final assembly result, but validate the assembly quality first.

How to adjust parameters if the assembly size is smaller than the expected genome size?
-------------------------------------------------------------------------------------------
For highly heterozygous genomes, try to set ``nextgraph_options = -a 1 -A``, otherwise you can set ``-q`` from 5 to 16 in ``nextgraph_options``, our tests show that setting ``nextgraph_options = -a 1 -q 10`` can usually get the best result.

Which job scheduling systems are supported by NextDenovo?
----------------------------------------------------------

NextDenovo uses `DRMAA <https://en.wikipedia.org/wiki/DRMAA>`__ to submit, control, and monitor jobs, so theoretically it supports all DRMAA-compliant systems, such as LOCAL, SGE, PBS, SLURM. See `ParallelTask <https://github.com/moold/ParallelTask>`_ to configure drmaa.

How to continue running unfinished tasks?
----------------------------------------------------------

No need to make any changes, simply run the same command again.

How to reduce the total number of subtasks?
----------------------------------------------------------

Please increase blocksize and reduce seed\_cutfiles.

How to speed up NextDenovo?
----------------------------------------------------------

Currently, the bottlenecks of NextDenovo are minimap2 and IO. For minimap2, please see `here <https://github.com/lh3/minimap2/issues/322>`__ to accelerate minimap2, besides, you can increase ``-l`` to reduce result size and disk consumption. For IO, you can check how many activated subtasks using top/htop, in theory, it should be equal to the ``-p`` parameter defined in correction\_options. Use usetempdir will reduce IO wait, especially if usetempdir is on a SSD driver.

How to specify the queue/cpu/memory/bash to submit jobs?
----------------------------------------------------------

Please use cluster\_options, NextDenovo will replace ``{vf}``, ``{cpu}``, ``{bash}`` with specific values needed for each jobs.

RuntimeError: Could not find drmaa library. Please specify its full path using the environment variable DRMAA\_LIBRARY\_PATH.
-------------------------------------------------------------------------------------------------------------------------------------------------
   
Please setup the environment variable: DRMAA\_LIBRARY\_PATH, see `here <https://github.com/pygridtools/drmaa-python>`__ for more details.

ERROR: drmaa.errors.DeniedByDrmException: code 17: error: no suitable queues.
---------------------------------------------------------------------------------------

This is usually caused by a wrong setting of cluster\_options, please check cluster\_options first. If you use SGE, you also can add ``-w n`` to cluster\_options, it will switch off validation for invalid resource requests. Please add a similar option for other job scheduling systems.
