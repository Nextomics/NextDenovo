[![Downloads](https://img.shields.io/github/downloads/Nextomics/NextDenovo/total?logo=github)](https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz)
[![Release](https://img.shields.io/github/release/Nextomics/NextDenovo.svg)](https://github.com/Nextomics/NextDenovo/releases)
[![Documentation Status](https://readthedocs.org/projects/nextdenovo/badge/?version=latest)](https://nextdenovo.readthedocs.io/en/latest/?badge=latest)

# NextDenovo
NextDenovo is a string graph-based *de novo* assembler for long reads (CLR, ~~HiFi~~ and ONT). It uses a "correct-then-assemble" strategy similar to canu (no correction step for PacBio HiFi reads), but requires significantly less computing resources and storages. After assembly, the per-base accuracy is about 98-99.8%, to further improve single base accuracy, try [NextPolish](https://github.com/Nextomics/NextPolish).

We benchmarked NextDenovo against other assemblers using Oxford Nanopore long reads from [human](https://nextdenovo.readthedocs.io/en/latest/TEST2.html) and [Drosophila melanogaster](https://nextdenovo.readthedocs.io/en/latest/TEST4.html), and PacBio continuous long reads (CLR) from [Arabidopsis thaliana](https://nextdenovo.readthedocs.io/en/latest/TEST3.html). NextDenovo produces more contiguous assemblies with fewer contigs compared to the other tools. NextDenovo also shows a high assembly accurate level in terms of assembly consistency and single-base accuracy.

## Installation

* **REQUIREMENT**
	* [Python](https://www.python.org/download/releases/) (Support python 2 and 3):
		* [Paralleltask](https://github.com/moold/ParallelTask) `pip install paralleltask`

* **INSTALL**  

	click [here](https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz) or use the following command:   
	```sh
	wget https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz
	tar -vxzf NextDenovo.tgz && cd NextDenovo
	```   

	If you want to compile from the source, run:

	```sh
	git clone git@github.com:Nextomics/NextDenovo.git
	cd NextDenovo && make
	```

* **TEST**  
`nextDenovo test_data/run.cfg`

## Learn

* [Quick Start](https://nextdenovo.readthedocs.io/en/latest/QSTART.html#quick-start) - no experience required, download and assemble now
* [Tutorial](https://nextdenovo.readthedocs.io/en/latest/TEST1.html) - step by step introduction to assemble the HG002 genome
* [FAQ](https://nextdenovo.readthedocs.io/en/latest/FAQ.html) - frequently asked questions
* [Parameter Reference](https://nextdenovo.readthedocs.io/en/latest/OPTION.html) - a detailed introduction about all the parameters
* [Cite](https://nextdenovo.readthedocs.io/en/latest/QSTART.html#cite) - if you get a good assembly with NextDenovo, please cite it

## Star

You can track updates by tab the `Star` button on the upper-right corner at this page.

## More

The complete user documentation is available [here](https://nextdenovo.readthedocs.io/en/latest/).
