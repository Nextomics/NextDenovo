.. _qstart:

.. image:: https://img.shields.io/github/downloads/Nextomics/NextDenovo/total?logo=github
   :target: https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz
   :alt: Download
.. image:: https://img.shields.io/github/release/Nextomics/NextDenovo.svg
   :target: https://github.com/Nextomics/NextDenovo/releases
   :alt: Version
.. image:: https://img.shields.io/github/issues/Nextomics/NextDenovo.svg
   :target: https://github.com/Nextomics/NextDenovo/issues
   :alt: Issues
.. image:: https://readthedocs.org/projects/nextdenovo/badge/?version=latest
   :target: https://nextdenovo.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. .. image:: https://img.shields.io/badge/切换-中文版本-9cf
..    :target: https://github.com/Nextomics/NextDenovo/issues
..    :alt: 中文版本

==========
NextDenovo
==========

NextDenovo is a string graph-based *de novo* assembler for long reads (CLR, HiFi and ONT). It uses a "correct-then-assemble" strategy similar to canu (no correction step for PacBio HiFi reads), but requires significantly less computing resources and storages. After assembly, the per-base accuracy is about 98-99.8%, to further improve single base accuracy, try `NextPolish <https://github.com/Nextomics/NextPolish>`_.

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

-  **REQUIREMENT**

   -  `Python <https://www.python.org/download/releases/>`__ (Support python 2 and 3):
   
      -  `Paralleltask <https://github.com/moold/ParallelTask>`__ 

      .. code-block:: shell
         
         pip install paralleltask

-  **INSTALL**  

   click `here <https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz>`__ or use the following command:

   .. code-block:: shell

      wget https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz
      tar -vxzf NextDenovo.tgz && cd NextDenovo

   If you want to compile from the source, run:

   .. code-block:: shell

      git clone git@github.com:Nextomics/NextDenovo.git
      cd NextDenovo && make

-  **TEST**
   
   .. code-block:: shell

      nextDenovo test_data/run.cfg 


Quick Start
~~~~~~~~~~~

#. Prepare input.fofn

   .. code-block:: shell

      ls reads1.fasta reads2.fastq reads3.fasta.gz reads4.fastq.gz ... > input.fofn
#. Create run.cfg

   .. code-block:: shell

      cp doc/run.cfg ./
   
   .. note:: Please set :ref:`read_type <read_type>` and :ref:`genome_size <genome_size>`, and refer to :ref:`doc/FAQ <how-to-optimize-parallel-computing-parameters>` and :ref:`doc/OPTION <options>` to optimize parallel computing parameters.

#. Run

   .. code-block:: shell

      nextDenovo run.cfg

#. Result

   -  Sequence: ``01_rundir/03.ctg_graph/nd.asm.fasta``
   -  Statistics: ``01_rundir/03.ctg_graph/nd.asm.fasta.stat``

Getting Help
~~~~~~~~~~~~

-  **HELP**

   Feel free to raise an issue at the `issue page <https://github.com/Nextomics/NextDenovo/issues/new/choose>`_.

   .. important:: Please ask questions on the issue page first. They are also helpful to other users.
-  **CONTACT**
   
   For additional help, please send an email to huj\_at\_grandomics\_dot\_com.

Copyright
~~~~~~~~~

NextDenovo is only freely available for academic use and other non-commercial use. For commercial use, please contact `GrandOmics <https://www.nextomics.cn/en/>`_.

Cite
~~~~

Hu, J. et al. An efficient error correction and accurate assembly tool for noisy long reads. bioRxiv 2023.03.09.531669 (2023) doi:10.1101/2023.03.09.531669.

Limitations
~~~~~~~~~~~

#. NextDenovo is optimized for assembly with seed\_cutoff >= 10kb. This should not be a big problem because it only requires the longest 30x-45x seeds length >= 10kb. For shorter seeds, it may produce unexpected results for some complex genomes and need be careful to check the quality.

Star
~~~~

You can track updates by tab the ``Star`` button on the upper-right corner at the `github page <https://github.com/Nextomics/NextDenovo>`_.

.. |ss| raw:: html

   <strike>

.. |se| raw:: html

   </strike>
