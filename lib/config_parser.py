"""Configuration file parser."""
import re, os, sys
from kit import *

__all__ = ["ConfigParser"]

log = plog()

class Config(dict):

	def __repr__(self):
		ret = "\n".join(("%-30s%s" % (str(k).strip() + ':', str(v).strip()) for k, v in \
			sorted(self.items(), key = lambda x: len(str(x[0]) + str(x[1]))) if v and not str(k).startswith('_')))
		return ret

class ConfigParser:
	def __init__(self, cfgfile = None):
		self.cfg = self._defaultcfg()
		if cfgfile:
			self.cfgdir = os.path.dirname(os.path.abspath(cfgfile))
			self._read(cfgfile)
			self._check()

	def update(self, total_base, seed_cutoff, seed_depth):

		def __cal_blocksize(b, p, s):
			bs = 0
			i = 2
			while not bs:# or bs >= 16000000000: #16G
				i += 1
				while i * p <= s * (s + 1) / 2:
					i += 1
				bs = float(b) * s / (i * p - s * (s + 1) / 2)
			return int(bs + 1000000) if bs > 10000000 else 10000000

		gs = parse_num_unit(self.cfg['genome_size'])
		total_depth = int(total_base / gs)
		self.cfg['seed_cutoff'] = str(seed_cutoff)
		self.cfg['seed_depth'] = str(seed_depth)
		self.cfg['blocksize'] = __cal_blocksize(total_base - seed_depth * gs, 
				int(self.cfg['parallel_jobs']), int(self.cfg['seed_cutfiles']))
		
		self._set_minlen()
		self.cfg['sort_options'] += ' -k %d' % ((total_depth - 2) if total_depth <= 30 else min(total_depth - 5, 40))

		if self.cfg['read_type'] == 'hifi':
			self.cfg['minimap2_options_raw'] += ' -f %d' % (seed_depth * 20)

	def _set_minlen(self):
		if self.cfg['read_type'] in ['clr', 'clr-rs', 'ont']:
			if '-k' not in self.cfg['minimap2_options_cns'].split():
				self.cfg['minimap2_options_cns'] += ' -k 17'
			if '-w' not in self.cfg['minimap2_options_cns'].split():
				self.cfg['minimap2_options_cns'] += ' -w 17'
			if '--minlen' not in self.cfg['minimap2_options_cns']:
				self.cfg['minimap2_options_cns'] += ' --minlen %d' % (min(2000, int(self.cfg['seed_cutoff'])/10))
			if '--maxhan1' not in self.cfg['minimap2_options_cns']:
				self.cfg['minimap2_options_cns'] += ' --maxhan1 %d' % (min(5000, int(self.cfg['seed_cutoff'])/2))
			if '-min_len_seed' not in self.cfg['correction_options']:
				self.cfg['correction_options'] += ' -min_len_seed %d' % (int(self.cfg['seed_cutoff'])/2)
		else:
			if '--minide' not in self.cfg['minimap2_options_cns']:
				self.cfg['minimap2_options_cns'] += ' --minide 0.1'
			if '--maxhan1' not in self.cfg['minimap2_options_cns']:
				self.cfg['minimap2_options_cns'] += ' --maxhan1 %d' % (min(1000, int(self.cfg['seed_cutoff'])/2))

	def _defaultcfg(self):
		cfg = Config()
		cfg['job_type'] = 'sge'
		cfg['job_prefix'] = 'nextDenovo'
		cfg['task'] = 'all'
		cfg['rewrite'] = 1
		cfg['deltmp'] = 1
		cfg['rerun'] = 3
		cfg['parallel_jobs'] = '10'
		cfg['workdir'] = os.getcwd()
		cfg['input_type'] = 'raw'
		cfg['use_drmaa'] = False
		cfg['submit'] = None
		cfg['kill'] = None
		cfg['check_alive'] = None
		cfg['job_id_regex'] = None

		cfg['read_cutoff'] = '1k'
		cfg['seed_depth'] = '45'
		cfg['blocksize'] = '10g'
		cfg['pa_correction'] = '3'
		cfg['nodelist'] = ''
		cfg['cluster_options'] = ''
		cfg['_sge_queue'] = ''
		cfg['seed_cutfiles'] = cfg['pa_correction']
		cfg['correction_options'] = '-p 10'
		cfg['sort_options'] = '-m 40g -t 8 -k 40'
		cfg['minimap2_options_raw'] = ''

		cfg['_random_round_with_less_accuracy'] = 0
		cfg['minimap2_options_cns'] = ''

		cfg['minimap2_options_map'] = ''
		cfg['ctg_cns_options'] = ''
		cfg['nextgraph_options'] = ''
		cfg['_nextgraph_out_format'] = 'fasta'
		return cfg

	def _read(self, cfgfile):
		with open(cfgfile) as IN:
			for line in IN:
				line = line.strip()
				if not line or line[0].startswith('#'):
					continue
				group = re.search(r'([^;\s]+)\s*[=:]\s*([^;#\n]+)(\s*|#.*)$', line) # a option = value1 value2 # annotation
				if group and group.groups()[1].strip():
					self.cfg[group.groups()[0]] = group.groups()[1].strip()

	def _check(self):
		self.cfg['seed_cutfiles'] = str(max(int(self.cfg['pa_correction']), int(self.cfg['seed_cutfiles'])))
		self.cfg['input_fofn'] = self.cfg['input_fofn'] if self.cfg['input_fofn'].startswith('/') else self.cfgdir + '/' + self.cfg['input_fofn']
		self.cfg['workdir'] = self.cfg['workdir'] if self.cfg['workdir'].startswith('/') else self.cfgdir + '/' + self.cfg['workdir']
		self.cfg['raw_aligndir'] = self.cfg['workdir'] + '/01.raw_align'
		self.cfg['cns_aligndir'] = self.cfg['workdir'] + '/02.cns_align'
		self.cfg['ctg_graphdir'] = self.cfg['workdir'] + '/03.ctg_graph'

		if 'usetempdir' in self.cfg:
			if str(self.cfg['usetempdir']).lower() in ['no', '0', 'false']:
				del self.cfg['usetempdir']
			elif self.cfg['job_type'].lower() == 'local':
				log.error('Error, usetempdir cannot be used with local.')
				sys.exit(1)
			else:
				if not self.cfg['usetempdir'].startswith('/'):
					log.error('Error, usetempdir must be absolute path')
					sys.exit(1)
				if self.cfg['job_type'] != 'sge' and not os.path.exists(self.cfg['nodelist']):
					log.error('Error, usetempdir must be used with nodelist for non-sge job_type.')
					sys.exit(1)
				self.cfg['correction_options'] += ' -dbuf '
		elif '-dbuf' in self.cfg['correction_options']:
			log.warning('using \'-dbuf\' parameter but not setting \'usetempdir\' parameter usually' + \
				' significantly reduces NextDenovo\'s performance.')

		if 'input_fofn' not in self.cfg or not os.path.exists(self.cfg['input_fofn']):
			log.error('Error, can not find input_fofn')
			sys.exit(1)

		if self.cfg['task'] not in ['all', 'correct', 'assemble']:
			log.error('Error, task only accept: all|correct|assemble')
			sys.exit(1)

		if self.cfg['input_type'] not in ['raw', 'corrected']:
			log.error('Error, input_type only accept: raw|corrected')
			sys.exit(1)

		if 'read_type' not in self.cfg:
			log.error('Error, can not find read_type')
			sys.exit(1)
		elif self.cfg['read_type'] not in ['clr', 'clr-rs', 'ont', 'hifi', 'ccs']:
			log.error('Error, read_type only accept: clr|clr-rs|ont|hifi|ccs')
			sys.exit(1)
		if self.cfg['read_type'] == 'ccs':
			self.cfg['read_type'] = 'hifi'
		if self.cfg['read_type'] == 'hifi':
			# self.cfg['input_type'] = 'corrected'
			if self.cfg['seed_depth'] == '45':
				self.cfg['seed_depth'] = '40' #TODO check 40?
			if '-R' not in self.cfg['nextgraph_options'] and '--max_ide_ratio' not in self.cfg['nextgraph_options']:
				self.cfg['nextgraph_options'] += ' -R 0.7'
			if '-x' in self.cfg['minimap2_options_cns'] and \
					parse_options_value(self.cfg['minimap2_options_cns'], '-x') != 'ava-hifi':
				log.error('Error, Hifi/CCS reads must use "-x ava-hifi" in minimap2_options_cns.')
				sys.exit(1)
			if '-sp' not in self.cfg['ctg_cns_options']:
				self.cfg['ctg_cns_options'] += ' -sp'
			if self.cfg['input_type'] == 'corrected':
				log.warning('HiFi data still need to do self-correction,' + \
					' so, it is not recommended to set input_type as corrected, continue anyway...')
			# if 'genome_size' not in self.cfg:
			# 	log.error('Error, genome_size must be set for Hifi/CCS reads.')
			# 	sys.exit(1)

		if 'seed_cutoff' not in self.cfg or str(self.cfg['seed_cutoff']).startswith('-'):
			self.cfg['seed_cutoff'] = 0
		else:
			self.cfg['seed_cutoff'] = str(parse_num_unit(self.cfg['seed_cutoff']))

		if int(self.cfg['seed_cutoff']) <= 0 and 'genome_size' not in self.cfg:
			log.error('Error, genome_size or seed_cutoff must be set.')
			sys.exit(1)
		elif int(self.cfg['seed_cutoff']) > 0:
			self._set_minlen()

		minimap2_threads1 = minimap2_threads2 = 8
		if 'minimap2_options_raw' in self.cfg and '-t' in self.cfg['minimap2_options_raw']:
			minimap2_threads1 = int(parse_options_value(self.cfg['minimap2_options_raw'], '-t'))
		else:
			self.cfg['minimap2_options_raw'] += ' -t 8'
		if 'minimap2_options_cns' in self.cfg and '-t' in self.cfg['minimap2_options_cns']:
			minimap2_threads2 = int(parse_options_value(self.cfg['minimap2_options_cns'], '-t'))
		else:
			self.cfg['minimap2_options_cns'] += ' -t 8'
		self.cfg['_minimap2_threads'] = (minimap2_threads1, minimap2_threads2)

		minimap2_x = ''
		if '-x' in self.cfg['minimap2_options_raw']:
			minimap2_x = parse_options_value(self.cfg['minimap2_options_raw'], '-x')
		elif '-x' in self.cfg['minimap2_options_cns']:
			minimap2_x = parse_options_value(self.cfg['minimap2_options_cns'], '-x')
		elif self.cfg['read_type'].startswith('clr'):
			minimap2_x = 'ava-pb'
		elif self.cfg['read_type'] == 'ont':
			minimap2_x = 'ava-ont'
		else:
			minimap2_x = 'ava-hifi'
					
		if '-x' not in self.cfg['minimap2_options_raw']:
			self.cfg['minimap2_options_raw'] += ' -x ' + minimap2_x
		if '-x' not in self.cfg['minimap2_options_cns']:
			self.cfg['minimap2_options_cns'] += ' -x ' + minimap2_x

		if '-max_lq_length' not in self.cfg['correction_options']:
			if 'ont' in minimap2_x:
				self.cfg['correction_options'] += ' -max_lq_length 10000'
			else:
				self.cfg['correction_options'] += ' -max_lq_length 1000'
		if '-r' not in self.cfg['correction_options']:
			self.cfg['correction_options'] += ' -r %s' % self.cfg['read_type']

		if 'pb' in minimap2_x:
			self.cfg['minimap2_options_map'] += ' -x map-pb'
		elif 'ont' in minimap2_x:
			self.cfg['minimap2_options_map'] += ' -x map-ont'
		else:
			self.cfg['minimap2_options_map'] += ' -x asm20' #check asm 10? asm5?


		self.cfg['_cns_threads'] = int(parse_options_value(self.cfg['correction_options'], '-p'))
		self.cfg['_map_threads'] = int(parse_options_value(self.cfg['minimap2_options_map'], '-t')) \
			if '-t' in self.cfg['minimap2_options_map'] else self.cfg['_cns_threads']
		if '-p' not in self.cfg['ctg_cns_options']:
			self.cfg['ctg_cns_options'] += ' -p ' + str(self.cfg['_cns_threads'])


		if '-t' not in self.cfg['sort_options']:
			self.cfg['sort_options'] += ' -t ' + str(self.cfg['_cns_threads'])
		self.cfg['_sort_threads'] = int(parse_options_value(self.cfg['sort_options'], '-t'))
		if '-m' not in self.cfg['sort_options']:
			self.cfg['sort_options'] += ' -m %dg' % self.cfg['_sort_threads']
		self.cfg['_sort_mem'] = parse_options_value(self.cfg['sort_options'], '-m')

		if str(self.cfg['rewrite']).lower() in ['no', '0', 'false']:
			self.cfg['rewrite'] = 0
		else:
			self.cfg['rewrite'] = 1
			log.warning('Re-write workdir')

		if str(self.cfg['deltmp']).lower() in ['no', '0', 'false']:
			self.cfg['deltmp'] = 0
		else:
			self.cfg['deltmp'] = 1
		
		if str(self.cfg['rerun']).lower() in ['no', '0', 'false']:
			self.cfg['rerun'] = 0
		else:
			self.cfg['rerun'] = min(int(self.cfg['rerun']), 10)

		if self.cfg['input_type'] == 'corrected':
			if self.cfg['task'] != 'assemble':
				log.warning('Change task "%s" to "assemble", becasue the input_type is "%s"', self.cfg['task'], self.cfg['input_type'])
				self.cfg['task'] = 'assemble'
		else:
			if self.cfg['task'] == 'assemble':
				log.warning('Change task "%s" to "all", becasue the input_type is "%s"', self.cfg['task'], self.cfg['input_type'])
				self.cfg['task'] = 'all'

		if int(self.cfg['pa_correction']) > int(self.cfg['parallel_jobs']):
			self.cfg['pa_correction'] = self.cfg['parallel_jobs']

		if str(self.cfg['use_drmaa']).lower() in ['no', '0', 'false']:
			self.cfg['use_drmaa'] = None
		for opt in ['submit', 'kill', 'check_alive', 'job_id_regex']:
			if self.cfg[opt] and self.cfg[opt].lower() == 'auto':
				self.cfg[opt] = None
		
		if '-a' in self.cfg['nextgraph_options']:
			format_value = parse_options_value(self.cfg['nextgraph_options'], '-a')
			if format_value == '0':
				self.cfg['_nextgraph_out_format'] = 'non'
			elif format_value == '1':
				self.cfg['_nextgraph_out_format'] = 'fasta'
			elif format_value == '2':
				self.cfg['_nextgraph_out_format'] = 'graphml'
			elif format_value == '3':
				self.cfg['_nextgraph_out_format'] = 'gfa'
			elif format_value == '4':
				self.cfg['_nextgraph_out_format'] = 'path'
			else:
				log.error('Error, unaccepted option value -a %s in nextgraph_options.' % format_value)
				sys.exit(1)

