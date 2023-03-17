#!/usr/bin/env python
from __future__ import print_function

import os, re
import time
import signal
import logging

import sys
if sys.version_info[0] == 2:
	PYTHON_VERSION = 2
elif sys.version_info[0] == 3:
	PYTHON_VERSION = 3
else:
	raise Exception("Unknown Python version")

__all__ = ['plog', 'pmkdir', 'write2file', 'parse_num_unit', 'parse_options_value', 
			'getver', 'pypath', 'cal_n50_info', 'which', 'str2byte', 'byte2str', 
			'calgs', 'cal_minlen_from_idx']

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

class ExitOnCritical(logging.StreamHandler):
	def emit(self, record):
		if isinstance(record.msg, list):
			record.msg = ' '.join(record.msg)
		elif isinstance(record.msg, dict):
			record.msg ="\n" +  "\n".join(("%-30s%s" % (str(k).strip() + ':', str(v).strip()) for k, v in \
				sorted(record.msg.items(), key = lambda x: len(str(x[0]) + str(x[1])))))
		elif hasattr(record.msg, '__dict__'):
			record.msg ="\n" + "\n".join(("%-30s%s" % (str(k).strip() + ':', str(v).strip()) for k, v in \
				sorted(vars(record.msg).items(), key = lambda x: len(str(x[0]) + str(x[1])))))

		if record.levelno >= logging.ERROR:
			msg = record.msg
			record.msg = '\033[35m%s\033[0m' % msg
			super(ExitOnCritical, self).emit(record)
			record.msg = msg
		else:
			super(ExitOnCritical, self).emit(record)
		if record.levelno >= logging.CRITICAL:
			need_emit = False
			for handler in logging.getLogger().handlers:
				if need_emit:
					handler.emit(record)
				elif handler == self:
					need_emit = True
			raise Exception(record.msg)

def plog(path=None):
	formatter = logging.Formatter('[%(process)d %(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
	log = logging.getLogger()

	if not log.handlers:
		log.setLevel(logging.INFO)
		# console_handler = logging.StreamHandler()
		console_handler = ExitOnCritical()
		console_handler.setFormatter(formatter)
		log.addHandler(console_handler)

	if path:
		has_path_logger = False
		for logger in log.handlers:
			if isinstance(logger, logging.FileHandler) and logger.baseFilename == path:
				has_path_logger = True
				break

		if not has_path_logger:
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

def popen(path, *args, **kwargs):
	pmkdir(os.path.dirname(path))
	return open(path, *args, **kwargs)

def ptime(init_time = 0):
	if init_time:
		return time.strftime("%H:%M:%S",time.gmtime((time.time() - init_time)))
	else:
		return time.strftime("%Y-%m-%d %H:%M:%S",time.gmtime((time.time())))

def write2file(content, path, mode = 'w'):
	with open(path, mode) as OUT:
		print(content, file=OUT)

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
		if j == contents_len or contents[j].startswith('-') or contents[j].startswith('>') or contents[j].startswith('1>'):
			if index:
				return j - 1 if last else (i, j)
			else:
				return contents[j - 1] if last else ' '.join(contents[i:j])
		j += 1

def parse_num_unit(content, base=1000):
	'''2Gb 2kb 2.3 kb 3.5g'''
	def parse_unit(unit, base=1000):
		if unit[0] in ['k', 'K']:
			return base
		elif unit[0] in ['m', 'M']:
			return base * base
		elif unit[0] in ['g', 'G']:
			return base * base * base

	if str(content)[-1].isdigit():
		return int(content)
	value, unit = 1, 1
	contents = str(content).strip().split()
	if len(contents) != 1:
		value = float(contents[0])
		unit = parse_unit(contents[1], base)
	else:
		if contents[0][-2].isdigit():
			value = float(contents[0][:-1])
			unit = parse_unit(contents[0][-1], base)
		else:
			value = float(contents[0][:-2])
			unit = parse_unit(contents[0][-2:], base)
	return int(value * unit + .499)

def latestver(url):
	if PYTHON_VERSION == 3:
		from urllib.request import urlopen
	else:
		from urllib2 import urlopen

	try:
		if 'releases' in url:
			g = re.search(r'download/(.*?)/', urlopen(url, timeout=1).read().decode("utf-8"))
		else:
			g = re.search(r'([\S]+)', urlopen(url, timeout=1).read())
		if g:
			return g.group(1).lower().strip('v')
		else:
			return 'Unknown'
	except Exception:
		return 'Unknown'

def getver(path):
	ver = 'Unknown'
	readme = path + '/VERSION'
	if os.path.exists(readme):
		with open(readme) as f:
			ver = f.read().lower().strip().strip('v')

	latest = latestver('https://api.github.com/repos/Nextomics/NextDenovo/releases/latest')
	if latest != 'Unknown' and ver != latest:
		print(('\033[35m Please update to the latest version: %s, current version: %s \033[0m') % (latest, ver))
	return ver

def pypath():
	return sys.executable

def cal_n50_info(stat, outfile = None):
	stat.sort(key = int, reverse = True)
	gs = sum(stat)
	k = 1
	l = j = 0
	out = "%-5s %20s %20s\n" % ("Type", "Length (bp)", "Count (#)")
	for i in stat:
		l += i
		j += 1
		while l >= gs * 0.1 * k and k < 10:
			out += "N%d0 %20d%20d\n" % (k, i, j)
			k += 1
	out += '\n'
	out += "%-5s %18d%20s\n" % ("Min.", stat[-1], '-')
	out += "%-5s %18d%20s\n" % ("Max.", stat[0], '-')
	out += "%-5s %18d%20s\n" % ("Ave.", gs/j, '-')
	out += "%-5s %18d%20d\n" % ("Total", gs, j)
	if outfile:
		write2file(out, outfile)
	return out

def which(program):
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file
	return None

def str2byte(string, ignore = True):
	if ignore:
		try:
			string = string.encode("UTF-8")
		except Exception:
			pass
		finally:
			return string
	else:
		return string.encode("UTF-8")

def byte2str(byte, ignore = True):
	if ignore:
		try:
			byte = byte.decode("UTF-8")
		except Exception:
			pass
		finally:
			return byte
	else:
		return byte.decode("UTF-8")

def calgs(infile):
	import ctypes
	KIT = ctypes.CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "ckit.so")
	KIT.calgs.argtypes = [ctypes.c_char_p]
	KIT.calgs.restype = ctypes.c_uint64
	return KIT.calgs(str2byte(infile))

def cal_minlen_from_idx(fofn, file_count, total_len):
	import ctypes
	KIT = ctypes.CDLL(os.path.dirname(os.path.realpath(__file__)) + '/' + "ckit.so")
	KIT.cal_minlen_from_idx.argtypes = [ctypes.POINTER(ctypes.c_char_p), ctypes.c_int, ctypes.c_uint64]
	KIT.cal_minlen_from_idx.restype = ctypes.c_uint32
	c_flies = (ctypes.c_char_p * file_count)(*[str2byte(infile) for infile in fofn])
	minlen = KIT.cal_minlen_from_idx(c_flies, file_count, total_len)
	del c_flies
	return minlen
