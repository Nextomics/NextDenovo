.. _qstart:

.. image:: https://img.shields.io/github/downloads/Nextomics/NextDenovo/total?logo=github
   :target: https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz
   :alt: Download
.. image:: https://img.shields.io/github/release/Nextomics/NextDenovo.svg
   :target: https://github.com/Nextomics/NextDenovo/releases
   :alt: Version
.. .. image:: https://img.shields.io/github/issues/Nextomics/NextDenovo.svg
..    :target: https://github.com/Nextomics/NextDenovo/issues
..    :alt: Issues
.. .. image:: https://img.shields.io/badge/切换-中文版本-9cf
..    :target: https://github.com/Nextomics/NextDenovo/issues
..    :alt: 中文版本

==========
NextDenovo
==========

NextDenovo is a string graph-based *de novo* assembler for long reads. It uses a "correct-then-assemble" strategy similar to canu, but requires significantly less computing resources and storages. After assembly, the per-base accuracy is about 98-99.8%, to further improve single base accuracy, please use `NextPolish <https://github.com/Nextomics/NextPolish>`_.

NextDenovo contains two core modules: NextCorrect and NextGraph. NextCorrect can be used to correct long noisy reads with approximately 15% sequencing errors, and NextGraph can be used to construct a string graph with corrected reads. It also contains a modified version of `minimap2 <https://github.com/lh3/minimap2>`_ and some useful utilities (see :ref:`utilities <utilities>` for more details).

We benchmarked NextDenovo against other assemblers using Oxford Nanopore long reads from :ref:`human <chm13_120x_ont>` and :ref:`Drosophila melanogaster <dmel_69x_ont>`, and PacBio continuous long reads (CLR) from :ref:`Arabidopsis thaliana <alta_192x_clr>`. NextDenovo produces more contiguous assemblies with fewer contigs compared to the other tools. NextDenovo also shows a high assembly accurate level in terms of assembly consistency and single-base accuracy.

.. Table of Contents
.. -----------------

.. -  `Installation <#install>`_
.. -  `Quick start <#start>`_
.. -  `Tutorial <./doc/TEST1.md>`_
.. -  `Parameters <./doc/OPTION.md>`_
.. -  `Benchmark <#benchmark>`_
.. -  `Utilities <./doc/UTILITY.md>`_
.. -  `Getting help <#help>`_
.. -  `Copyright <#copyright>`_
.. -  `Cite <#cite>`_
.. -  `Limitations <#limit>`_
.. -  `FAQ <#faq>`_
.. -  `Star <#star>`_

Installation
~~~~~~~~~~~~

-  **DOWNLOAD**  

   click `here <https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz>`__ or use the following command:

   .. code:: console

      wget https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz

   .. note:: If you get an error like ``version 'GLIBC_2.14' not found`` or ``liblzma.so.0: cannot open shared object file``, Please download `this version <https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo-CentOS6.9.tgz>`_.

-  **REQUIREMENT**

   -  `Python <https://www.python.org/download/releases/>`_ (Support python 2 and 3):
   
      -  `Psutil <https://psutil.readthedocs.io/en/latest/>`_
      -  `Drmaa <https://github.com/pygridtools/drmaa-python>`_ (Only required by running under non-local system)

-  **INSTALL**

      |ss| ``tar -vxzf NextDenovo.tgz && cd NextDenovo && make`` |se|

-  **UNINSTALL**
   
      |ss| ``cd NextDenovo && make clean`` |se|

-  **TEST**
   
   .. code:: console

      nextDenovo test_data/run.cfg 


Quick Start
~~~~~~~~~~~

#. Prepare input.fofn

   .. code:: console

      ls reads1.fasta reads2.fastq reads3.fasta.gz reads4.fastq.gz ... > input.fofn
#. Create run.cfg

   .. code:: console

      cp doc/run.cfg ./
   
   .. note:: Please change the value of seed\_cutoff using :ref:`bin/seq\_stat <seq_stat>` and refer to :ref:`doc/OPTION <options>` to optimize parallel computing parameters.

#. Run

   .. code:: console

      nextDenovo run.cfg

#. Result

   -  Sequence: ``01_rundir/03.ctg_graph/nd.asm.fasta``
   -  Statistics: ``01_rundir/03.ctg_graph/nd.asm.fasta.stat``

Getting Help
~~~~~~~~~~~~

-  **HELP**

   Feel free to raise an issue at the `issue page <https://github.com/Nextomics/NextDenovo/issues/new/choose>`_. They would also be helpful to other users.

   .. important:: Please ask questions on the issue page first. They are also helpful to other users and avoid answering the same questions again and again.
-  **CONTACT**
   
   For additional help, please send an email to huj\_at\_grandomics\_dot\_com.

Copyright
~~~~~~~~~

NextDenovo is freely available for academic use and other non-commercial use. For commercial use, please contact `NextOmics <https://www.nextomics.cn/en/>`_.

Cite
~~~~

We are now preparing the manuscript of NextDenovo, so if you use NextDenovo now, please cite the official website (https://github.com/Nextomics/NextDenovo)

Limitations
~~~~~~~~~~~

#. The current version of NextDenovo is not suitable for assembly with PacBio HiFi reads, becasue Minimap2 does not optimize for HiFi reads overlapping.
#. NextDenovo is optimized for assembly with seed\_cutoff >= 10kb. This should not be a big problem because it only requires the longest 30x-45x seeds length >= 10kb. For shorter seeds, it may produce unexpected results for some complex genomes and need be careful to check the quality.

Star
~~~~

You can track updates by tab the "Star" button on the upper-right corner at the `github page <https://github.com/Nextomics/NextDenovo>`_.

.. |ss| raw:: html

   <strike>

.. |se| raw:: html

   </strike>
