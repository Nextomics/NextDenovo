# NextDenovo Parameter Reference

NextDenovo requires at least one read file (option: input_fofn) as input, it works with gzip'd FASTA and FASTQ formats and uses a config file to pass options, see [here](./run.cfg) for an example.

* **INPUT**    
    - reads files list (one file one line)  
    `ls reads1.fasta reads2.fastq reads3.fasta.gz reads4.fastq.gz ... > input.fofn`
    - config file   
    A config file is a text file that contains a set of parameters (key=value pairs) to set runtime parameters for NextDenovo (see OPTION section).

* **OUTPUT**    
	- `workdir/03.ctg_graph/nd.asm.fasta`, contigs with fasta format, the fasta header includes ID, type, length, node count, a consecutive lowercase region in the sequence implies a weak connection, and a low quality base is marked with a single lowercase base.
	- `workdir/03.ctg_graph/nd.asm.fasta.stat`, some basic statistical information (N10-N90, Total size et al.).  
	
* **OPTION** 
<pre>
	[General]                # global options
	job_type = sge           # [local, sge, pbs...]. (default: sge)
	job_prefix = nextDenovo  # prefix tag for jobs. (default: nextDenovo)
	task = all               # task need to run [all, correct or assemble]. (default: all)
	rewrite = no             # overwrite existed directory [yes, no]. (default: no)
	deltmp = yes             # delete intermediate results. (default: yes)
	rerun = 3                # re-run unfinished jobs untill finished or reached ${rerun} loops, 0=no. (default: 3)
	parallel_jobs = 10       # number of tasks used to run in parallel. (default: 10)
	input_type = raw         # input reads type [raw, corrected]. (default: raw)
	input_fofn = input.fofn  # input file, one line one file. (<b>required</b>)
	workdir = 01.workdir     # work directory. (default: ./)
	usetempdir = /tmp/test   # temporary directory in compute nodes to avoid high IO wait. (default: no)
	nodelist = avanode.list.fofn
	                         # a list of hostnames of available nodes, one node one line, used with usetempdir for non-sge job_type.
	cluster_options = auto
	                         # a template to define the resource requirements for each job, which will pass to <a href="https://github.com/pygridtools/drmaa-python/wiki/FAQ">DRMAA</a> as the nativeSpecification field.

	[correct_option]         # options using only in corrected step.
	read_cutoff = 1k         # filter reads with length < read_cutoff. (default: 1k)
	seed_cutoff = 25k        # minimum seed length. (<b>required</b>)
	seed_cutfiles = 10       # split seed reads into ${seed_cutfiles} subfiles. (default: ${pa_correction})
	blocksize = 10g          # block size for parallel running. (default: 10g)
	pa_correction = 15       # number of corrected tasks used to run in parallel, overwrite parallel_jobs only for this step. (default: 15)
	minimap2_options_raw = -x ava-ont -t 10   
	                         # minimap2 options, used to find overlaps between raw reads and set PacBio/Nanopore read overlap, see <a href="./UTILITY.md/#minimap2-nd">here</a> for details. (<b>required</b>)
	sort_options = -m 5g -t 2 -k 50   
	                         # sort options, see <a href="./UTILITY.md/#ovl_sort">here</a> for details.  
	correction_options = -p 10            
	                         # -p, --process, set the number of processes used for correcting. (default: 10)
	                         # -b, --blacklist, disable the filter step and increase more corrected data.
	                         # -s, --split, split the corrected seed with un-corrected regions. (default: False)
	                         # -fast, 0.5-1 times faster mode with a little lower accuracy. (default: False)
	                         # -dbuf, Disable caching 2bit files and reduce ~TOTAL_INPUT_BASES/4 bytes of memory usage. (default:False)
	                         # -max_lq_length, maximum length of a continuous low quality region in a corrected seed, larger max_lq_length will produce more corrected data with lower accuracy. (default: auto [pb/1k, ont/10k])

	[assemble_option]
	minimap2_options_cns = -x ava-ont -t 8 -k17 -w17 
	                         # minimap2 options, used to find overlaps between corrected reads. (default: -k17 -w17)
	nextgraph_options = -a 1 # nextgraph options, see <a href="./UTILITY.md/#nextgraph">here</a> for details.
</pre>
