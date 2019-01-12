#!/usr/bin/env python

import os
import sys
import re
import argparse
from multiprocessing import Pool
from kit import *
from ctypes import *

log = plog()

class HelpFormatter(
		argparse.RawDescriptionHelpFormatter,
		argparse.ArgumentDefaultsHelpFormatter):
	pass

class consensus_trimed_data (Structure):
	_fields_ = [
		("len", c_uint),
		("identity", c_float),
		("seq", c_char_p)
	]

class index_tag (Structure):
	_fields_ = [
		("f2bit", POINTER(c_void_p)),
		("offset", c_uint64),
	]

class ovl (Structure):
	_fields_ = [
		("ovl", c_void_p),
		("f2bits", POINTER(c_void_p)),
		("bases_arr", POINTER(c_char_p)),
		("decode_tbl", POINTER(c_uint32)),
		("handle_index", c_int),
		("indexs", POINTER(index_tag)),
	]

class ids (Structure):
	_fields_ = [
		("prev_qname", c_uint32),
		("prev_tname", c_uint32),
	]

CNS = CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "nextCorrector.so")
CNS.nextCorrect.argtypes = [POINTER(c_char_p), POINTER(c_uint), POINTER(c_uint), c_uint, c_uint, \
	c_uint, c_uint, c_uint, c_uint, c_float, c_uint, c_uint]
CNS.nextCorrect.restype = POINTER(consensus_trimed_data)
CNS.free_consensus_trimed.argtypes = [POINTER(consensus_trimed_data)]

OVL = CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "ovlSeq.so")
OVL.init_ovls.argtypes = [c_char_p, c_char_p]
OVL.init_ovls.restype = POINTER(ovl)
OVL.destory_ovls.argtypes = [POINTER(ovl)]
OVL.decode_ovl.argtypes = [POINTER(ovl), POINTER(ids), POINTER(c_uint32)]
OVL.decode_ovl.restype = c_int
OVL.bit2seq.argtypes = [POINTER(ovl), POINTER(c_uint32)]
OVL.bit2seq.restype = c_void_p
OVL.destory_seq.argtypes = [c_void_p]

def correct(seqs, aln_start, aln_end, count, max_aln_length, min_len_aln, max_cov_aln, min_cov_base, \
		max_lq_length, min_error_corrected_ratio, split, fast):

	c_seqs = (c_char_p * count)()
	c_seqs[:] = seqs
	c_aln_start = (c_uint * count)(*aln_start)
	c_aln_end = (c_uint * count)(*aln_end)

	c_consensus = CNS.nextCorrect(c_seqs, c_aln_start, c_aln_end, count, max_aln_length, min_len_aln, max_cov_aln, \
		min_cov_base, max_lq_length, min_error_corrected_ratio, split, fast)
	lens = c_consensus.contents.len
	identity = c_consensus.contents.identity
	sequence = string_at(c_consensus.contents.seq)
	
	CNS.free_consensus_trimed(c_consensus)
	del c_seqs
	del c_aln_start
	del c_aln_end
	return (lens, identity, sequence)

def read_seq_data(args, corrected_seeds):

	count = total_length = seed_length = max_aln_length = 0
	used_reads = []
	seed_name = ''
	seqs = []
	aln_start = []
	aln_end = []

	ovls = OVL.init_ovls(args.idxs, args.ovl)

	ids_ = ids(0, 0)
	tmp_arr = (c_uint32 * 7)()
	last_seed = -1

	while OVL.decode_ovl(ovls, byref(ids_), tmp_arr) >= 0:
		t_name, _, t_s, t_e, q_name, q_s, q_e = tmp_arr
		if seed_name == '+' or (last_seed != -1 and t_name != last_seed):
			if total_length / seed_length >= args.min_cov_seed:
				# print seed_name
				# log.info(seed_name)
				yield (seed_name, count, seqs, aln_start, aln_end, max_aln_length, args)
			count = 0
			seed_name = ''
			seqs = []
			aln_start = []
			aln_end = []
			total_length = 0
			seed_length = 0
			max_aln_length = 0
			used_reads = []
		if not str(seed_name):
			seed_length = t_e + 1
			total_length += seed_length
			max_aln_length = seed_length
			seed_name = t_name if seed_length >= args.min_len_seed and t_name not in corrected_seeds else '+'
		if t_e - t_s < args.min_len_aln or total_length / \
				seed_length > args.max_cov_aln * 1.5 or q_name in used_reads or seed_name == '+':
			continue
		else:
			seq = OVL.bit2seq(ovls, tmp_arr)
			used_reads.append(q_name)
			aln_start.append(t_s)
			aln_end.append(t_e)
			seqs.append(string_at(seq))
			total_length += t_e - t_s + 1
			count += 1
			OVL.destory_seq(seq)
			_ = t_e - t_s + q_e - q_s + 2
			if  _ > max_aln_length and t_name != q_name:
				max_aln_length = _
		last_seed = t_name
	if seed_length and total_length / seed_length >= args.min_cov_seed:
		yield (seed_name, count, seqs, aln_start, aln_end, max_aln_length, args)
	del ids_
	del tmp_arr
	OVL.destory_ovls(ovls)

def read_corrected_seeds(infile, corrected_seeds):
	last_seed = ''
	last_seed_position = 0
	with open(infile) as IN:
		for line in IN:
			if line.startswith('>'):
				last_seed_position = len(line)
				last_seed = seed_name = re.split(r'[_,\s]+', line.strip())[0][1:]
				log.warning('Skip corrected seed: ' + seed_name)
				corrected_seeds.append(int(seed_name))
			else:
				last_seed_position += len(line)

	if last_seed:
		corrected_seeds.remove(int(last_seed))
	return last_seed_position

def worker(args):
	seed_name, count, seqs, aln_start, aln_end, max_aln_length, args = args

	args.split = 1 if args.split else 0
	args.fast = 1 if args.fast else 0
	if args.max_lq_length < 5:
		args.max_lq_length = 5
		log.warning('Reset max_lq_length to %d.' % args.max_lq_length)

	lens, identity, sequence = correct(seqs, aln_start, aln_end, count, max_aln_length, args.min_len_aln, \
			args.max_cov_aln, args.min_cov_base, args.max_lq_length, args.min_error_corrected_ratio, args.split, args.fast)
	return (seed_name, lens, identity, sequence)

def start():
	log.info(
		'Start a correctted worker in %d from parent %d' %
		(os.getpid(), os.getppid()))

def main(args):
	log.info('Corrected step options:')
	log.info(args)
	corrected_region = re.compile(r"[ACGT]+")

	OUT = sys.stdout
	corrected_seeds = []
	last_seed_position = 0
	if args.out != 'stdout':
		if os.path.exists(args.out):
			last_seed_position = read_corrected_seeds(args.out, corrected_seeds)
			OUT = open(args.out, 'r+')
			OUT.seek(-1 * last_seed_position, 2)
		else:
			OUT = open(args.out, 'w')

	pool = Pool(args.process, initializer=start)
	for seed_name, lens, identity, seq in pool.imap_unordered(
			worker, read_seq_data(args, corrected_seeds), chunksize=1):
		if lens >= args.min_len_seed and identity >= args.min_error_corrected_ratio:
			if args.split:
				corrected_regions = corrected_region.findall(seq)
				for i in range(len(corrected_regions)):
					lens = len(corrected_regions[i])
					identity = 1
					if lens >= args.min_len_seed:
						print >>OUT, '>%s_%d %d %f\n%s' % (
							seed_name, i + 1, lens, identity, corrected_regions[i])
			else:
				print >>OUT, '>%s %d %f\n%s' % (seed_name, lens, identity, seq.lstrip('atgc'))
		elif lens == 3:
			log.warning('No enough memory and fail to correct seed: ' + str(seed_name))

	pool.close()
	pool.join()
	if args.out != 'stdout':
		OUT.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		formatter_class=HelpFormatter,
		description='''
%(prog)s:
	correct seed reads with multi-reads-block using multi-processor.
exmples:
	%(prog)s -f input.idxs -i sorted.ovl -p 24 > seed.cns.fa
'''
	)
	parser.add_argument('-f', '--idxs', metavar = 'FILE', required=True,
		help='set the reads index file.')
	parser.add_argument('-i', '--ovl', metavar = 'FILE', required=True,
		help='set the reads overlap file.')
	parser.add_argument('-o', '--out', metavar = 'FILE', default='stdout',
		help='set the output file, corrected seeds in output file will be skipped.')
	parser.add_argument('-p', '--process', type=int, default=10,
		help='set the number of processes used for correcting.')
	parser.add_argument('-s', '--split', action='store_true', default=False,
		help='split the corrected seed with un-corrected regions.')
	parser.add_argument('-fast', action='store_true', default=False,
		help='0.5-1 times faster mode with a little less accuracy.')
	parser.add_argument('-max_cov_aln', type=int, default=130,
		help='maximum depth of the seed read, only use reads up to the MAX_COV_ALN average depth ranked by alignment length.')
	parser.add_argument('-max_lq_length', type=int, default=10000,
		help='maximum length of a continuous low quality region in a corrected seed.')
	parser.add_argument('-min_cov_seed', type=int, default=10,
		help='minimum depth requirement of the seed read.')
	parser.add_argument('-min_len_seed', type=str, default=10000,
		help='minimum length of a seed read.')
	parser.add_argument('-min_len_aln', type=str, default=500,
		help='minimum length of a alignment used to correct.')
	parser.add_argument('-min_cov_base', type=int, default=4,
		help='minimum depth to correct a raw base.')
	parser.add_argument('-min_error_corrected_ratio', type=float, default=0.8, 
		help='minimum corrected ratio of a corrected seed.')
	args = parser.parse_args()
	args.max_lq_length = parse_num_unit(args.max_lq_length)
	args.min_len_seed = parse_num_unit(args.min_len_seed)
	args.min_len_aln = parse_num_unit(args.min_len_aln)
	main(args)
