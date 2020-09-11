[![Downloads](https://img.shields.io/github/downloads/Nextomics/NextDenovo/total?logo=github)](https://github.com/Nextomics/NextDenovo/releases/download/v2.3.1/NextDenovo.tgz)
[![Release](https://img.shields.io/github/release/Nextomics/NextDenovo.svg)](https://github.com/Nextomics/NextDenovo/releases)
[![Issues](https://img.shields.io/github/issues/Nextomics/NextDenovo.svg)](https://github.com/Nextomics/NextDenovo/issues)

# NextDenovo
NextDenovo is a string graph-based *de novo* assembler for TGS long reads. It uses a "correct-then-assemble" strategy similar to canu, but requires significantly less computing resources and storages. After assembly, the per-base error rate is about 98-99.8%, to further improve single base accuracy, please use [NextPolish](https://github.com/Nextomics/NextPolish).

NextDenovo contains two core modules: NextCorrect and NextGraph. NextCorrect can be used to correct TGS long reads with approximately 15% sequencing errors, and NextGraph can be used to construct a string graph with corrected reads. It also contains a modified version of [minimap2](https://github.com/lh3/minimap2) and some useful utilities (see [here](./doc/UTILITY.md) for more details).

We benchmarked NextDenovo against other assemblers using Oxford Nanopore long reads from human and drosophila melanogaster, and PacBio continuous long reads (CLR) from arabidopsis thaliana. NextDenovo produces more contiguous assemblies with fewer contigs compared to the other tools. NextDenovo also shows a high assembly accurate level in terms of assembly consistency and single-base accuracy, see [here](#benchmark) for details.

## Table of Contents

- [Installation](#install)
- [Quick start](#start)
- [Tutorial](./doc/TEST1.md)
- [Parameters](./doc/OPTION.md)
- [Benchmark](#benchmark)
- [Utilities](./doc/UTILITY.md)
- [Getting help](#help)
- [Copyright](#copyright)
- [Cite](#cite)
- [Limitations](#limit)
- [FAQ](#faq)
- [Star](#star)

### <a name="install"></a>Installation
* **DOWNLOAD**   
click [here](https://github.com/Nextomics/NextDenovo/releases/download/v2.3.1/NextDenovo.tgz) or use the following command:   
`wget https://github.com/Nextomics/NextDenovo/releases/download/v2.3.1/NextDenovo.tgz`   
***Note:*** If you get an error like `version 'GLIBC_2.14' not found` or `liblzma.so.0: cannot open shared object file`, Please download [this version](https://github.com/Nextomics/NextDenovo/releases/download/v2.3.1/NextDenovo-CentOS6.9.tgz).

* **REQUIREMENT**
	* [Python](https://www.python.org/download/releases/) (Support python 2 and 3):
		* [Psutil](https://psutil.readthedocs.io/en/latest/)
		* [Drmaa](https://github.com/pygridtools/drmaa-python) (Only required by running under non-local system)

* **INSTALL**  
~~`tar -vxzf NextDenovo.tgz && cd NextDenovo && make`~~

* **UNINSTALL**  
~~`cd NextDenovo && make clean`~~

* **TEST**  
`nextDenovo test_data/run.cfg`

### <a name="start"></a>Quick Start
1. Prepare input.fofn  
	`ls reads1.fasta reads2.fastq reads3.fasta.gz reads4.fastq.gz ... > input.fofn`
2. Create run.cfg  
    `cp NextDenovo/doc/run.cfg ./`    
    ***Note:*** Please change the value of seed_cutoff using [bin/seq_stat](./doc/UTILITY.md#seq_stat) and refer to [doc/OPTION](doc/OPTION.md) to optimize parallel computing parameters. 
3. Run  
    `nextDenovo run.cfg`
4. Result
    - sequence: `01_rundir/03.ctg_graph/nd.asm.fasta`
    - statistics: `01_rundir/03.ctg_graph/nd.asm.fasta.stat`

### <a name="benchmark"></a>Benchmark
+ [CHM13hTERT human cell line with 120x Oxford Nanopore data](./doc/TEST2.md) 
+ [Arabidopsis thaliana with 192X PacBio CLR data (~1% heterozygosity)](./doc/TEST3.md) 
+ [Drosophila melanogaster with 69X Oxford Nanopore data](./doc/TEST4.md)

### <a name="help"></a>Getting Help
* **HELP**   
Please raise an issue at the [issue page](https://github.com/Nextomics/NextDenovo/issues/new/choose). They would also be helpful to other users.

* **CONTACT**    
For additional help, please send an email to huj_at_grandomics_dot_com.

### <a name="copyright"></a>Copyright
NextDenovo is freely available for academic use and other non-commercial use. For commercial use, please contact [NextOmics](https://www.nextomics.cn/en/).

### <a name="cite"></a>Cite
We are now preparing the manuscript of NextDenovo, so if you use NextDenovo now, please cite the official website (https://github.com/Nextomics/NextDenovo)

### <a name="limit"></a>Limitations
1. The current version of NextDenovo is not suitable for assembly with PacBio HiFi reads, becasue Minimap2 does not optimize for HiFi reads overlapping.
2. NextDenovo is optimized for assembly with seed_cutoff >= 10kb. This should not be a big problem because it only requires the longest 30x-45x seeds length >= 10kb. For shorter seeds, it may produce unexpected results for some complex genomes and need be careful to check the quality.

### <a name="faq"></a>Frequently Asked Questions
1. Which job scheduling systems are supported by NextDenovo?  
NextDenovo uses [DRMAA](https://en.wikipedia.org/wiki/DRMAA) to submit, control, and monitor jobs, so theoretically it supports all DRMAA-compliant systems, such as LOCAL, SGE, PBS, SLURM. See [ParallelTask](https://github.com/moold/ParallelTask) to configure drmaa.
2. **How to continue running unfinished tasks?**  
No need to make any changes, simply run the same command again.
3. How to reduce the total number of subtasks?  
Please increase blocksize and reduce seed_cutfiles.
4. How to speed up NextDenovo?  
Currently, the bottlenecks of NextDenovo are minimap2 and IO. For minimap2, please see [here](https://github.com/lh3/minimap2/issues/322) to accelerate minimap2, besides, you can increase -l to reduce result size and disk consumption. For IO, you can check how many activated subtasks using top/htop, in theory, it should be equal to the -p parameter defined in correction_options. Use usetempdir will reduce IO wait, especially if usetempdir is on a SSD driver.
5. How to specify the queue cpu/memory/bash to submit jobs?  
Please use cluster_options, NextDenovo will replace {vf}, {cpu}, {bash} with specific values needed for each jobs.
6. RuntimeError: Could not find drmaa library.  Please specify its full path using the environment variable DRMAA_LIBRARY_PATH.   
Please setup the environment variable: DRMAA_LIBRARY_PATH, see [here](https://github.com/pygridtools/drmaa-python) for more details.
7. ERROR: drmaa.errors.DeniedByDrmException: code 17: error: no suitable queues.  
This is usually caused by a wrong setting of cluster_options, please check cluster_options first. If you use SGE, you also can add '-w n' to cluster_options, it will switch off validation for invalid resource requests. Please add a similar option for other job scheduling systems. 

### <a name="star"></a>Star
You can track updates by tab the "Star" button on the upper-right corner at this page.
