.. _faq:

Frequently Asked Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
  :local:

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
