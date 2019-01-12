#!/usr/bin/env python

import os
import time
import signal
import logging

__all__ = ['plog', 'pmkdir', 'write2file', 'parse_num_unit', 'parse_options_value']

class TimedOutExc(Exception):
	pass

class deco(object):
	@staticmethod
	def _deadline(timeout):
		def _decorate(f):
			def handler(signum, frame):
				raise TimedOutExc
			def new_f(*args):
				signal.signal(signal.SIGALRM, handler)
				signal.alarm(timeout)
				try:
					return f(*args)
				except TimedOutExc:
					return 0
				finally:
					signal.alarm(0)
			return new_f
		return _decorate

def plog(path = False):

	formatter = logging.Formatter('[%(levelname)s] %(asctime)s %(message)s')
	log = logging.getLogger()
	if log.handlers:
		return log

	log.setLevel(logging.DEBUG)
	console_handler = logging.StreamHandler()
	console_handler.setFormatter(formatter)
	log.addHandler(console_handler)

	if path:
		fileHandler = logging.FileHandler(path, mode='w')
		fileHandler.setFormatter(formatter)
		log.addHandler(fileHandler)
	return log

def pmkdir(path):
	if not os.path.exists(path):
		os.makedirs(path)
		return True
	else:
		return False

def ptime(init_time = 0):
	if init_time:
		return time.strftime("%H:%M:%S",time.gmtime((time.time() - init_time)))
	else:
		return time.strftime("%Y-%m-%d %H:%M:%S",time.gmtime((time.time())))

def write2file(content, path):
	with open(path, 'w') as OUT:
		print >>OUT, content

def parse_options_value(content, option, last = True, index = False):
	'''-p1 [4] -p2 [smp [5]] -p3 [ad] -d4 [f] '''
	contents = content.strip().split()
	contents_len = len(contents)
	for i in range(contents_len):
		if contents[i] == option:
			i += 1
			break
	j = i
	while (j <= contents_len):
		if j == contents_len or contents[j].startswith('-'):
			if index:
				return j - 1 if last else (i, j)
			else:
				return contents[j - 1] if last else ' '.join(contents[i:j])
		j += 1

def parse_num_unit(content):
	'''2Gb 2kb 2.3 kb 3.5g'''
	def parse_unit(unit):
		if unit[0] in ['k', 'K']:
			return 1000
		elif unit[0] in ['m', 'M']:
			return 1000000
		elif unit[0] in ['g', 'G']:
			return 1000000000

	if str(content)[-1].isdigit():
		return int(content)
	value, unit = 1, 1
	contents = str(content).strip().split()
	if len(contents) != 1:
		value = float(contents[0])
		unit = parse_unit(contents[1])
	else:
		if contents[0][-2].isdigit():
			value = float(contents[0][:-1])
			unit = parse_unit(contents[0][-1])
		else:
			value = float(contents[0][:-2])
			unit = parse_unit(contents[0][-2:])
	return int(value * unit + .499)