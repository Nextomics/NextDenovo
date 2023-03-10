#!/usr/bin/env python
from __future__ import print_function

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

class buffer_t (Structure):
	pass

class sbbuf(Structure):
	pass

class ovl (Structure):
	_fields_ = [
		("ovl", c_void_p),
		("f2bits", POINTER(c_void_p)),
		("f2bits_", POINTER(POINTER(c_uint32))),
		("decode_tbl", POINTER(c_uint32)),
		("handle_index", c_int),
		("indexs", POINTER(index_tag)),
		("ovlbuf", POINTER(buffer_t)),
		("bitbuf", POINTER(sbbuf)),
	]

class ids (Structure):
	_fields_ = [
		("prev_qname", c_uint32),
		("prev_tname", c_uint32),
	]

CNS = CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "nextcorrect.so")
CNS.nextCorrect.argtypes = [POINTER(c_char_p), POINTER(c_uint), POINTER(c_uint), c_uint, c_uint, \
	c_uint, c_uint, c_uint, c_uint, c_float, c_uint, c_uint, c_int]
CNS.nextCorrect.restype = POINTER(consensus_trimed_data)
CNS.free_consensus_trimed.argtypes = [POINTER(consensus_trimed_data)]

OVL = CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "ovlseq.so")
OVL.init_ovls.argtypes = [c_char_p, c_char_p, c_int]
OVL.init_ovls.restype = POINTER(ovl)
OVL.destory_ovls.argtypes = [POINTER(ovl)]
OVL.decode_ovl.argtypes = [c_void_p, POINTER(c_uint32), POINTER(ids), POINTER(c_uint32), POINTER(buffer_t), c_int]
OVL.decode_ovl.restype = c_int
OVL.getseq.argtypes = [POINTER(ovl), POINTER(c_uint32)]
OVL.getseq.restype = c_void_p
global OVLDB

def correct(seqs, aln_start, aln_end, count, max_aln_length, min_len_aln, max_cov_aln, min_cov_base, \
		max_lq_length, min_error_corrected_ratio, split, fast, read_type):

	c_seqs = (c_char_p * count)()
	c_seqs[:] = seqs
	c_aln_start = (c_uint * count)(*aln_start)
	c_aln_end = (c_uint * count)(*aln_end)

	c_consensus = CNS.nextCorrect(c_seqs, c_aln_start, c_aln_end, count, max_aln_length, min_len_aln, max_cov_aln, \
		min_cov_base, max_lq_length, min_error_corrected_ratio, split, fast, read_type)
	identity = c_consensus.contents.identity
	sequence = byte2str(string_at(c_consensus.contents.seq))
	lens = c_consensus.contents.len
	
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

	ids_ = ids(0, 0)
	tmp_arr = (c_uint32 * 8)()
	last_seed = -1

	while OVL.decode_ovl(OVLDB.contents.ovl, OVLDB.contents.decode_tbl, byref(ids_), tmp_arr, OVLDB.contents.ovlbuf, 8) >= 0:
		t_name, _, t_s, t_e, q_name, q_s, q_e, match = tmp_arr
		if seed_name == '+' or (last_seed != -1 and t_name != last_seed):
			if total_length / seed_length >= args.min_cov_seed:
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
			total_length = 0
			max_aln_length = seed_length
			seed_name = t_name if seed_length >= args.min_len_seed and t_name not in corrected_seeds else '+'
		if t_e - t_s < args.min_len_aln or total_length / \
				seed_length > args.max_cov_aln * 1.5 or q_name in used_reads or seed_name == '+':
			continue
		else:
			# seq = OVL.getseq(OVLDB, tmp_arr) # here is too slow if args.process > ~7
			# seqs.append(string_at(seq))
			seqs.append(list(tmp_arr))
			used_reads.append(q_name)
			aln_start.append(t_s)
			aln_end.append(t_e)
			total_length += t_e - t_s + 1
			count += 1
			_ = t_e - t_s + q_e - q_s + 2
			if  _ > max_aln_length and t_name != q_name:
				max_aln_length = _
		last_seed = t_name
	if seed_length and total_length / seed_length >= args.min_cov_seed:
		yield (seed_name, count, seqs, aln_start, aln_end, max_aln_length, args)
	del ids_
	del tmp_arr

def	read_blacklist(infile, corrected_seeds):
	i = 0
	if os.path.exists(str(infile)): 
		with open(infile) as IN:
			for line in IN:
				seed_name = line.strip().split()[0]
				# log.warning('Skip seed in blacklist: ' + seed_name)
				corrected_seeds.add(int(seed_name))
				i += 1
		log.warning('Skip %d seeds in blacklist.', i)

def read_corrected_seeds(IN, corrected_seeds):
	i = 0
	offset = -1
	seed_name = 0
	last_line = ''
	last_validy_lines = [0, 0, -1]
	line = IN.readline()
	while line:
		offset = IN.tell()
		lines = line.strip().split()
		seed_name = re.split(r'[_,\s]+', lines[0])[0]
		corrected_seeds.add(int(seed_name))
		if last_line:
			last_validy_lines = last_line.strip().split()
		last_line = line
		i += 1
		line = IN.readline()

	if offset != -1:
		IN.seek(offset - len(last_line), 0)
		IN.truncate()
		corrected_seeds.remove(int(seed_name))
		i -= 1
	if i:
		log.warning('Skip %d corrected seeds.', i)
	return int(last_validy_lines[1]) + int(last_validy_lines[2]) + 1

def worker(args):
	seed_name, count, seqs_, aln_start, aln_end, max_aln_length, args = args

	args.split = 1 if args.split else 0
	args.fast = 1 if args.fast else 0
	args.max_lq_length = min(int(aln_end[0]/2), args.max_lq_length)
	seqs = []
	log.debug("init seed: %s" % seed_name)
	for tmp_arr_ in seqs_:
		tmp_arr = (c_uint32 * 8)(*tmp_arr_)
		seqs.append(string_at(OVL.getseq(OVLDB, tmp_arr)))
	log.debug("ovl2block done: %s" % seed_name)
	lens, identity, sequence = correct(seqs, aln_start, aln_end, count, max_aln_length, args.min_len_aln, \
			args.max_cov_aln, args.min_cov_base, args.max_lq_length, args.min_error_corrected_ratio, args.split, \
			args.fast, args.read_type)
	log.debug("cns done: %s" % seed_name)
	return (seed_name, lens, identity, sequence)

def start():
	log.info(
		'Start a cns worker in %d from parent %d' %
		(os.getpid(), os.getppid()))

def main(args):
	log.info('Corrected step options:')
	log.info(args)
	if args.debug:
		log.setLevel(10)
	corrected_region = re.compile(r"[ACGT]+")

	OUT = sys.stdout
	IDX = None
	corrected_seeds = set()
	read_blacklist(args.blacklist, corrected_seeds)
	if args.out != 'stdout':
		if os.path.exists(args.out):
			IDX = open(args.out + '.idx', 'r+')
			last_seed_position = read_corrected_seeds(IDX, corrected_seeds)
			OUT = open(args.out, 'r+')
			OUT.seek(last_seed_position, 0)
			OUT.truncate()
		else:
			OUT = open(args.out, 'w')
			IDX = open(args.out + '.idx', 'w')

	global OVLDB
	OVLDB = OVL.init_ovls(str2byte(args.idxs), str2byte(args.ovl), 0 if args.dbuf else 1)

	fail_seed = 0
	pool = Pool(args.process, initializer=start)
	for seed_name, lens, identity, seq in pool.imap_unordered(
			worker, read_seq_data(args, corrected_seeds), chunksize = 
				min(10, int(len(corrected_seeds)/100/args.process) + 1)):
		if lens >= args.min_len_seed and identity >= args.min_error_corrected_ratio:
			if args.split:
				corrected_regions = corrected_region.findall(seq)
				for i in range(len(corrected_regions)):
					lens = len(corrected_regions[i])
					identity = 1
					if lens >= args.min_len_seed:
						print('>%s_%d %d %f\n%s' % (
							seed_name, i + 1, lens, identity, corrected_regions[i]), file=OUT)
						if IDX:
							cur_offset = OUT.tell()
							print('%s_%d\t%d\t%d' % (seed_name, i+1, cur_offset - lens - 1, lens), file=IDX)
			else:
				print('>%s %d %f\n%s' % (seed_name, lens, identity, seq), file=OUT)
				if IDX:
					cur_offset = OUT.tell()
					print('%d\t%d\t%d' % (seed_name, cur_offset - lens - 1, lens), file=IDX)

		else:
			if lens == 3:
				fail_seed += 1
				log.warning('No enough memory and fail to correct seed: ' + str(seed_name))
			elif IDX:
				print('%d\t%d\t%d' % (seed_name, 0, 0), file=IDX)
		# IDX.flush()

	pool.close()
	pool.join()
	OVL.destory_ovls(OVLDB)
	if args.out != 'stdout':
		OUT.close()
		IDX.close()
	if fail_seed > 5:
		log.error('A total of %d seeds are failed to be corrected, increass RAM or decrease "-p" and rerun.' % fail_seed)
		sys.exit(1)

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
	parser.add_argument('-r', '--read_type', required=True, type=str.lower, choices=['clr', 'hifi', 'ont'],
		help='reads type, clr=PacBio continuous long read, hifi=PacBio highly accurate long reads, ont=NanoPore 1D reads')
	parser.add_argument('-b', '--blacklist', action='store_false', default=True,
		help='skip the seeds in blacklist file.')
	parser.add_argument('-o', '--out', metavar = 'FILE', default='stdout',
		help='set the output file, corrected seeds in output file will be skipped.')
	parser.add_argument('-p', '--process', type=int, default=10,
		help='set the number of processes used for correcting.')
	parser.add_argument('-s', '--split', action='store_true', default=False,
		help='split the corrected seed with un-corrected regions.')
	parser.add_argument('-dbuf', action='store_true', default=False,
		help='Disable caching 2bit files and reduce ~TOTAL_INPUT_BASES/4 bytes of memory usage.')
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
	parser.add_argument('-debug', action='store_true', default=False,
		help='debug mode')
	args = parser.parse_args()
	args.max_lq_length = parse_num_unit(args.max_lq_length)
	args.min_len_seed = parse_num_unit(args.min_len_seed)
	args.min_len_aln = parse_num_unit(args.min_len_aln)
	if args.read_type == 'ont':
		args.read_type = 1
	elif args.read_type == 'clr':
		args.read_type = 2
	elif args.read_type == 'hifi':
		args.read_type = 3
	else:
		log.error('Unrecognized option: -r %s' % args.read_type)
		sys.exit(1)
	if args.blacklist:
		args.blacklist = args.ovl + '.bl'
	main(args)
