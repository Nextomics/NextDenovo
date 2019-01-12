# NextDenovo
Fast and accurate *de novo* assembler for third generation sequencing (TGS) long reads

> **NextDenovo** includes two parts:

> * **NextCorrect**  (Binary released)   
>   This package aims to correct the TGS long reads with approximately 15% sequencing errors. The design philosophy is from [Falcon](https://github.com/PacificBiosciences/FALCON). After correction, for nanopore data, the identity is about 97-98%, most of remaining errors are located in the homopolymer or tandem repeat regions.
>
> * **NextGraphy**   (ongoing)    
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
	job_queue = all.q        # bind job to queue(s). (default: all.q)
	task = all               # task need to run [all, correct or graph]. (default: all)
	rewrite = no             # overwrite existed directory [yes, no]. (default: no)
	rerun = 3                # re-run unfinished jobs untill finished or reached ${rerun} loops, 0=no. (default: 3)
	input_type = raw         # input reads type [raw, corrected]. (default: raw)
	input_fofn = input.fofn  # input file, one line one file. (<b>required</b>)
	workdir = 01.workdir     # work directory. (default: ./)
	usetempdir = /tmp        # temporary directory in compute nodes to avoid high IO wait. (default: no)

	[correct_option]         # options using only in corrected step.
	read_cuoff = 1k          # filted reads with length < read_cuoff. (default: 1k)
	seed_cutoff = 25k        # Minimum seed length. (<b>required</b>)
	seed_cutfiles = 10       # split seed reads into ${seed_cutfiles} subfiles. (default: ${pa_correction})
	blocksize = 10g          # block size for parallel running. (default: 10g)
	pa_raw_align = 10        # number of processes used for aligning. (default: 10)
	pa_correction = 15       # number of processes used for correcting. (default: 15)
	minimap2_options = -x ava-ont -t 10   # minimap2 options, used to set PacBio/Nanopore read overlap (<b>required</b>)
	sge_correction_options = -pe smp 15   # request slot range for parallel jobs. (default: -pe smp 15)
	correction_options = -p 10            # nextCorrector options, see below. (default: -p 10)
</pre>

<pre>
	[correction_options]
	-p , --process              set the number of processes used for correcting. (default: 10)
	-s, --split                 split the corrected seed with un-corrected regions. (default: False)
	-fast                       0.5-1 times faster mode with a little less accuracy. (default: False)
	-max_cov_aln                maximum depth of the seed read, only use reads up to the MAX_COV_ALN average depth ranked by alignment length. (default: 130)
	-max_lq_length              maximum length of a continuous low quality region in a corrected seed. (default: auto [pb/1k, ont/10k])
	-min_cov_seed               minimum depth requirement of the seed read. (default: 10)
	-min_len_seed               minimum length of a seed read. (default: 10000)
	-min_len_aln                minimum length of a alignment used to correct. (default: 500)
	-min_cov_base               minimum depth to correct a raw base. (default: 4)
	-min_error_corrected_ratio  minimum corrected ratio of a corrected seed. (default: 0.8)
</pre>

* **OUTPUT**    
cns.fasta with fasta format, the fasta header includes primary seqID, length, corrected bases percentage (%). The two flanking un-corrected sequences are trimed.

* **HELP**   
Please raise an issue at the issue page.

* **Contact information**
For additional help, please send an email to huj_at_grandomics_dot_com.
