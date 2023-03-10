#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import argparse
from multiprocessing import Pool
from kit import *
from ctypes import *

log = plog()

class HelpFormatter(
		argparse.RawDescriptionHelpFormatter,
		argparse.ArgumentDefaultsHelpFormatter):
	pass

class ctg_cns_cfg (Structure):
	pass

class consensus_trimed (Structure):
	_fields_ = [
		("len", c_uint),
		("identity", c_float),
		("seq", c_char_p)
	]

class consensus_trimed_data (Structure):
	_fields_ = [
		("data", POINTER(consensus_trimed)),
		("i_m", c_int),
	]

class ref_qv (Structure):
	pass

class ref_ (Structure):
	_fields_ = [
		("n", c_char_p),
		("s", POINTER(c_uint32)),
		("qv", POINTER(ref_qv)),
		("qv_l", c_uint32),
		("length", c_uint32),
	]

class refs_ (Structure):
	_fields_ = [
		("ref", POINTER(ref_)),
		("i", c_uint32),
		("i_m", c_uint32),
	]

global P, CFG, REFS
P = CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "ctg_cns.so")
P.read_ref.argtypes = [c_char_p, POINTER(c_char_p), c_int]
P.read_ref.restype = POINTER(refs_)
P.refs_destroy.argtypes = [POINTER(refs_)]

P.ctg_cns_init.argtypes = [c_int, c_int, c_int, c_float, c_float, c_float]
P.ctg_cns_init.restype = POINTER(ctg_cns_cfg)
P.ctg_cns_destroy.argtypes = [POINTER(ctg_cns_cfg)]

P.ctg_cns_core.argtypes = [POINTER(ctg_cns_cfg), POINTER(ref_), c_char_p]
P.ctg_cns_core.restype = POINTER(consensus_trimed_data)
P.free_consensus_trimed_data.argtypes = [POINTER(consensus_trimed_data)]

def set_window_process(args):
	import psutil
	p = args.process
	w = args.window
	max_mem = psutil.virtual_memory().available/1536
	max_cpu = psutil.cpu_count()
	if args.process > max_cpu:
		args.process = max_cpu
	if args.window < 5000000 or args.process * args.window > max_mem:
		args.window = 5000000

	process = int(max_mem/args.window)
	if args.process > process:
		args.process = process
	# args.process = int(max_mem/args.window)
	# if args.process > max_cpu:
	# 	args.process = max_cpu
	# 	args.window = int(max_mem/args.process)

	if p != args.process or w != args.window:
		log.warning(
			'Adjust -p from %d to %d, -w from %d to %d, logical CPUs:%d, available RAM:~%dG, ' 
			% (p, args.process, w, args.window, max_cpu, max_mem/1024/682) + 
			'use -a to disable automatic adjustment.')

def iter_uncorrected_seqs(uncorrected_seqs, bam_list):
	i = REFS.contents.i
	for i in range(i):
		if string_at(REFS.contents.ref[i].n) in uncorrected_seqs:
			yield (i, bam_list)

def read_uncorrected_seqs(infile, index, corrected_seqs):
	uncorrected_seqs = set()
	if args.block_index != 'all':
		with open(infile) as IN:
			for line in IN:
				lines = line.strip().split()
				if lines:
					if lines[0] not in corrected_seqs and lines[1] == index:
						uncorrected_seqs.add(str2byte(lines[0]))
	else:
		with open(infile) as IN:
			for line in IN:
				if line.startswith('>'):
					seq_name = line.strip().split()[0][1:]
					if seq_name not in corrected_seqs:
						uncorrected_seqs.add(str2byte(seq_name))
	return uncorrected_seqs

def read_corrected_seqs(infile, corrected_seqs):
	last_seq = ''
	cur_seq_offset = last_seq_position = 0
	with open(infile) as IN:
		for line in IN:
			if line.startswith('>'):
				lines = line.split()[0].split('_s')
				last_seq = seq_name = lines[0][1:]
				if (len(lines) == 1 or lines[1] == '0'):
					# last_seq_position = len(line)
					last_seq_position += cur_seq_offset
					cur_seq_offset = len(line)
					log.warning('Skip corrected seq: ' + seq_name)
				else:
					cur_seq_offset += len(line)
				corrected_seqs.add(seq_name)
			else:
				cur_seq_offset += len(line)

	if last_seq:
		corrected_seqs.remove(last_seq)
	return last_seq_position

def worker(args):
	n, bam_list = args
	c_seq = P.ctg_cns_core(CFG, REFS.contents.ref[n], bam_list)
	c_seq_data = []
	for i in range(c_seq.contents.i_m):
		seq = byte2str(string_at(c_seq.contents.data[i].seq))
		seq_len = c_seq.contents.data[i].len
		seq_name = byte2str(string_at(REFS.contents.ref[n].n))
		if c_seq.contents.i_m != 1:
			seq_name += '_s%d' % (i)
		c_seq_data.append([seq_name, seq, seq_len])
	P.free_consensus_trimed_data(c_seq)
	return c_seq_data

def start():
	log.info(
		'Start a corrected worker in %d from parent %d' %
		(os.getpid(), os.getppid()))

def main(args):
	log.info('Corrected step options:')
	log.info(args)

	OUT = sys.stdout
	corrected_seqs = set()
	last_seq_position = 0
	if args.out != 'stdout':
		if os.path.exists(args.out):
			last_seq_position = read_corrected_seqs(args.out, corrected_seqs)
			OUT = open(args.out, 'r+')
			# OUT.seek(-1 * last_seq_position, 2)
			OUT.seek(last_seq_position, os.SEEK_SET)
		else:
			OUT = open(args.out, 'w')
	
	blockfile = args.block
	if args.block_index == 'all' or not args.block:
		args.block_index = 'all'
		blockfile = args.genome
	uncorrected_seqs = read_uncorrected_seqs(blockfile, args.block_index, corrected_seqs)
	
	global CFG, REFS
	uncorrected_seqs_len = len(uncorrected_seqs)
	c_uncorrected_seqs = (c_char_p * uncorrected_seqs_len)()
	c_uncorrected_seqs[:] = list(uncorrected_seqs)
	REFS = P.read_ref(str2byte(args.genome), c_uncorrected_seqs, uncorrected_seqs_len);
	del c_uncorrected_seqs

	if args.auto:
		set_window_process(args)

	CFG = P.ctg_cns_init(args.window, args.read_type, args.split, args.alignment_identity_ratio, 
		args.alignment_score_ratio, args.alignment_score_ratio)
	pool = Pool(args.process, initializer=start)
	for seq_data in pool.imap_unordered(worker, 
			iter_uncorrected_seqs(uncorrected_seqs, str2byte(args.bam_list)), chunksize=1):
		for seq_name, seq, seq_len in seq_data:
			if args.uppercase:
				seq = seq.upper()
			if seq_len > 10:
				print('>%s %d\n%s' % (seq_name, seq_len, seq), file=OUT)
			else:
				log.error('Failed to correct sequence: %s' % seq_name)
				sys.exit(1)

	pool.close()
	pool.join()
	if args.out != 'stdout':
		OUT.close()
	P.ctg_cns_destroy(CFG)
	P.refs_destroy(REFS)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		formatter_class = HelpFormatter,
		description = '''
%(prog)s:
	correct structural & base errors in the genome with long reads using multi-processor.

exmples:
	%(prog)s -g genome.fa -l lgs.sort.bam.list -r ont -p 10
'''
	)
	parser.add_argument('-g', '--genome', metavar = 'FILE', required=True, type=str,
		help='genome file, the reference of bam alignments.')
	parser.add_argument('-l', '--bam_list', metavar = 'FILE', required=True, type=str,
		help='sorted bam file list of long reads, one file one line, require index file.')
	parser.add_argument('-r', '--read_type', required=True, type=str.lower, choices=['clr', 'clr-rs', 'hifi', 'ont'],
		help='reads type, clr=PacBio continuous long read from sequel, clr-rs=PacBio continuous long read from RS, \
		hifi=PacBio highly accurate long reads, ont=NanoPore 1D reads')
	parser.add_argument('-b', '--block', metavar = 'FILE', type=str,
		help='genome block file, each line includes [seq_id, index].')
	parser.add_argument('-i', '--block_index', type=str, default = 'all',
		help='index of seqs need to be corrected in genome block file.')
	parser.add_argument('-o', '--out', metavar = 'FILE', default='stdout',
		help='output file, corrected seqs in output file will be skipped.')
	parser.add_argument('-p', '--process', metavar = 'INT', type=int, default=10,
		help='number of processes used for correcting.')
	parser.add_argument('-u', '--uppercase', action='store_true', default=False,
		help='output uppercase sequences.')
	parser.add_argument('-w', '--window', metavar = 'STR', type=str, default="5M",
		help='size of window (>=5M) to split super-long contigs, shorter size requires less memory and more CPU time.')
	parser.add_argument('-a', '--auto', action='store_false', default=True,
		help='automatically adjust window size (-w) and processes (-p).')
	parser.add_argument('-sp', '--split', action='store_false', default=True,
		help='split the corrected contig with un-corrected regions.')
	parser.add_argument('-id', '--alignment_identity_ratio', metavar = 'FLOAT', type=float, default=0.8,
		help='split the corrected contig if alignment_identity/median_alignment_identity < $identity_ratio, co-use with --split.')
	parser.add_argument('-as', '--alignment_score_ratio', metavar = 'FLOAT', type=float, default=0.8,
		help='split the corrected contig if alignment_score/max_alignment_score < $alignment_score_ratio, co-use with --split.')
	args, unknown = parser.parse_known_args()
	args.window = parse_num_unit(args.window)
	args.split = 1 if args.split else 0
	if args.read_type == 'ont':
		args.read_type = 1
	elif args.read_type == 'clr':
		args.read_type = 2
	elif args.read_type == 'hifi':
		args.read_type = 3
	elif args.read_type == 'clr-rs':
		args.read_type = 4
	else:
		log.error('Unrecognized option: -r %s' % args.read_type)
		sys.exit(1)
	main(args)
