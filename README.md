# NextDenovo
Fast and accurate *de novo* assembler for third generation sequencing (TGS) long reads

> **NextDenovo** includes two parts:

> * **NextCorrect**  (Binary released)  
>   This package aims to correct the TGS long reads with approximately 15% sequencing errors. The design philosophy is from [Falcon](https://github.com/PacificBiosciences/FALCON). After correction, for nanopore data, the identity is about 97-98%, most of remaining errors are located in the homopolymer or tandem repeat regions.
>
> * **NextGraph**   (ongoing)    
>      This package aims to construct the string graph with corrected data. 

* **INSTALL**  
`cd NextDenovo && make`

* **UNINSTALL**  
`cd NextDenovo && make clean`

* **RUN**  
`nextDenovo test_data/run.cfg -l log.txt`

* **INPUT** 

<pre>
	[General]                # global options
	job_type = sge           # [sge or local]. (default: sge)
	job_prefix = nextDenovo  # prefix tag for jobs. (default: nextDenovo)
	task = all               # task need to run [all, correct or graph]. (default: all)
	rewrite = no             # overwrite existed directory [yes, no]. (default: no)
	rerun = 3                # re-run unfinished jobs untill finished or reached ${rerun} loops, 0=no. (default: 3)
	input_type = raw         # input reads type [raw, corrected]. (default: raw)
	input_fofn = input.fofn  # input file, one line one file. (<b>required</b>)
	workdir = 01.workdir     # work directory. (default: ./)
	usetempdir = /tmp/test   # temporary directory in compute nodes to avoid high IO wait. (default: no)
	sge_options = -l vf={vf} -q all.q -pe smp {cpu}
	                         # options used in SGE.

	[correct_option]         # options using only in corrected step.
	read_cuoff = 1k          # filter reads with length < read_cuoff. (default: 1k)
	seed_cutoff = 25k        # minimum seed length. (<b>required</b>)
	seed_cutfiles = 10       # split seed reads into ${seed_cutfiles} subfiles. (default: ${pa_correction})
	blocksize = 10g          # block size for parallel running. (default: 10g)
	pa_raw_align = 10        # number of alignment tasks used to run in parallel. (default: 10)
	pa_correction = 15       # number of corrected tasks used to run in parallel. (default: 15)
	sort_options = -m 5g -t 2 -k 50   
	                         # sort_overlap options, see below. (defalut: -m 40G -t 8 -k 40)
	minimap2_options = -x ava-ont -t 10   
	                         # minimap2 options, used to set PacBio/Nanopore read overlap. (<b>required</b>)
	correction_options = -p 10            
	                         # nextCorrector options, see below. (default: -p 10)
</pre>


<pre>
	[sort_options]
	-t                       # number of threads to use. [8]
	-m                       # set maximum available buffer size, larger buffer size will accelerate sort process, suffix K/M/G. [40G]
	-k                       # maximum depth of each overlap, larger depth will produce more accurate and more corrected data with slower speed. [40]

	[correction_options]
	-p , --process           # set the number of processes used for correcting. (default: 10)
	-s, --split              # split the corrected seed with un-corrected regions. (default: False)
	-fast                    # 0.5-1 times faster mode with a little lower accuracy. (default: False)
	-max_lq_length           # maximum length of a continuous low quality region in a corrected seed, larger max_lq_length will produce
	                           more corrected data with lower accuracy. (default: auto [pb/1k, ont/10k])
</pre>

* **OUTPUT**    
cns.fasta with fasta format, the fasta header includes primary seqID, length, approximate corrected bases percentage (%). The two flanking un-corrected sequences are trimed.

* **HELP**   
Please raise an issue at the issue page.

* **Contact information**    
For additional help, please send an email to huj_at_grandomics_dot_com.

* **COPYRIGHT**    
NextDenovo is freely available for academic use and other non-commercial use. For commercial use, please contact [NextOmics](https://www.nextomics.cn/en/).

* **STAR**  
