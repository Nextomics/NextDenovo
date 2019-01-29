#!/usr/bin/env python

import sys
import os
import re
import time
import psutil
import signal
from kit import *

__all__ = ['Task', 'Run']

log = plog()

class Task(object):

	def __init__(self, path, group = 1, max_subtask = 300, prefix = 'job', bash = '/bin/bash', convertpath = True): # need add output input
		self.tasks = []
		self.subtasks = []
		self.path = [path]
		self.group = group
		self.max_subtask = max_subtask
		self.prefix = prefix
		self.bash = bash
		self.run = ''
		self._check(convertpath)
		self._subtasks()

	def add_task(self, path, prefix = 0, group = 0, max_subtask = 0, bash = 0, convertpath = True):
		self.tasks = []
		self.path.append(path)
		if group:
			self.group = group
		if max_subtask:
			self.max_subtask = max_subtask
		if prefix:
			self.prefix = prefix
		if bash:
			self.bash = bash
		self._check(convertpath)
		self._subtasks()

	def check(self, path = False):
		if path:
			return True if os.path.exists(path + '.done') else False
		else:
			return True if os.path.exists(self.path[-1] + '.done') else False

	def set_task_done(self):
		for path in self.path:
			os.system('touch ' + path + '.done')

	def _check(self, convertpath = True):
		with open(self.path[-1]) as IN:
			for line in IN:
				line = line.strip()
				if not line or line[0].startswith('#'):
					continue
				lines = line.split()
				if convertpath:
					for i in range(len(lines)):
						if not re.search(r'\/',lines[i]):
							if os.path.exists(lines[i]) or lines[i - 1] == '>' or lines[i - 1] == '1>' or lines[i - 1] == '2>':
								lines[i] = os.path.abspath(lines[i])
							elif lines[i].startswith('>') and len(lines[i]) > 1 :
								lines[i] = lines[i][0] + ' ' + os.path.abspath(lines[i][1:])
							elif lines[i].startswith(('1>','2>')) and len(lines[i]) > 2:
								lines[i] = lines[i][:2] + ' ' + os.path.abspath(lines[i][2:])
							elif re.search(r'=',lines[i]):
								paramters = lines[i].split("=")
								lines[i] = paramters[0] + '=' + os.path.abspath(paramters[1]) if \
								paramters[1] and len(paramters) == 2 and os.path.exists(paramters[1]) else lines[i]
						elif re.search(r'=',lines[i]):
							paramters = lines[i].split("=")
							lines[i] = paramters[0] + '=' + os.path.abspath(paramters[1]) if \
								paramters[1] and len(paramters) == 2 else lines[i]
						else:
							lines[i] = os.path.abspath(lines[i])
				self.tasks.append(" ".join(lines))

	def _subtasks(self):
		tasks = []
		if isinstance(self.group, int):
			task = []
			for i in range(len(self.tasks)):
				if task and i % self.group == 0:
					tasks.append(task)
					task = []
				task.append(self.tasks[i])
			tasks.append(task)
		elif isinstance(self.group, list):
			j = 0
			last_i = 0
			task = []
			for i in range(len(self.tasks)):
				if i - last_i == self.group[j]:
					tasks.append(task)
					task = []
					last_i = i
					j += 1
				task.append(self.tasks[i])
			tasks.append(task)
		elif isinstance(self.group, str):
			task = []
			for i in range(len(self.tasks)):
				if task and i == self.group:
					tasks.append(task)
					task = []
				task.append(self.tasks[i])
			tasks.append(task)

		else:
			log.error('incorrect allocated task group')

		self.tasks = tasks
		log.info('analysis tasks done')

	def set_subtasks(self, job_prefix = 'work', maxfile = 300):
		d = os.path.abspath(self.path[-1]) + '.work/'
		subtask_file = job_prefix + '.sh'
		task_count = len(self.tasks)
		split = 1 if task_count > maxfile else 0
		for i in range(task_count):
			if split:
				subtask_dir = self.prefix + '{:0>{}}/'.format(split/maxfile, len(str(task_count/maxfile))) + self.prefix + '{:0>{}}'.format(i, len(str(task_count)))
				split += 1
			else:
				subtask_dir = self.prefix + '{:0>{}}'.format(i, len(str(task_count)))

			subtask_finish_lable = d + subtask_dir + '/' + subtask_file + '.done'
			if not os.path.exists(subtask_finish_lable):
				pmkdir(d + subtask_dir)
				subtask = '#!' + self.bash + '\n' + \
						'set -xve' + '\n' + \
						'hostname' + '\n' + \
						'cd ' + d + subtask_dir + '\n'
				for task in self.tasks[i]:
					subtask += 'time ' + task + '\n'
				subtask += 'touch ' + d + subtask_dir + '/' + subtask_file + '.done\n'
				with open(d + subtask_dir + '/' + subtask_file, 'w') as OUT:
					print >>OUT, subtask
			self.subtasks.append(d + subtask_dir + '/' + subtask_file)

	def set_run(self, max_pa_jobs = 5, bash = '/bin/bash', job_type = 'sge', interval = 30, cpu = 1, vf = '', sge_options = ''):
		self.run = Run(self.subtasks, max_pa_jobs, bash, job_type, interval, cpu , vf, sge_options)

class Run(object):
	"""docstring for Run"""
	RUNNINGTASK = {'sge':[], 'local':[]}

	def __init__(self, tasks, max_pa_jobs, bash, job_type, interval, cpu, vf, sge_options):
		self.tasks = tasks
		self.unfinished_tasks = []
		self.job_type = job_type
		self.max_pa_jobs = int(max_pa_jobs)
		self.interval = int(interval)
		self.cpu = str(cpu)
		self.vf = str(vf) if vf else self.cpu + 'G'
		self.bash = str(bash)
		self.sge_options = str(sge_options)
		self.option = self._getoption()
		self.check()
		
	def start(self):
		log.info('total jobs: ' + str(len(self.unfinished_tasks)))
		if self.job_type == 'sge':
			self._sge()
		else:
			self._local()

	def rerun(self):
		# Here we can analysis log file to check the reason of unfinished jobs, and then decided how to rerun unfinished jobs.
		if self.job_type == 'sge':
			vfs = re.split(r'(\d+)', self.vf)
			self.vf = str(int(int(vfs[-2]) * 1.5)) + vfs[-1] #most unfinished jobs were caused by Memory in SGE system. TODO: should avoid larger the total memory of computer-node
			self.option = self._getoption()
		# self.check()
		self.start()

	def check(self):
		self.unfinished_tasks = []
		for task in self.tasks:
			if not os.path.exists(task + '.done'):
				self.unfinished_tasks.append(task)
		return False if self.unfinished_tasks else True

	@classmethod
	def kill(cls, signum, frame):
		log.warning('Accept killed signal and kill all running jobs, please waite...')
		for pid in Run.RUNNINGTASK['sge']:
			os.popen('qdel ' + pid)
		for pid in Run.RUNNINGTASK['local']:
			os.kill(pid, signal.SIGKILL)
		Run.RUNNINGTASK = {'sge':[], 'local':[]}
		log.warning('Kill running jobs done')
		sys.exit(1)

	def _sge(self):
		j = 0
		runnedJop = []
		while j < len(self.unfinished_tasks):
			task = self.unfinished_tasks[j]
			if j < self.max_pa_jobs or len(self._check_running(runnedJop)) <= self.max_pa_jobs:
				qsubCmd = self.option + ' -o ' + task + '.o ' + '-e ' + task + '.e ' + task
				jobId = os.popen(qsubCmd).read().strip().split()[2]
				runnedJop.append(jobId)
				Run.RUNNINGTASK['sge'].append(jobId)
				log.info('Throw jobID:[' + jobId + '] jobCmd:['  + qsubCmd + '] in the qsub_cycle.')
				j += 1
			else:
				time.sleep(self.interval)
		else:
			while (1):
				if len(self._check_running(runnedJop)) != 0:
					time.sleep(self.interval)
				else:
					break

	def _check_running(self, runnedJop):
		runningJop = []
		user = os.popen('whoami').read().strip()
		qstatContent = os.popen('qstat -u ' + user).read().strip().split('\n')
		for i in qstatContent:
			i = i.strip().split()
			if i and i[0] in runnedJop and i[3] == user:
				if i[4] in ['qw', 'r', 't', 'hqw']:
					runningJop.append(i[0])
				elif i[4] == 'Eqw' or (len(i) > 7 and i[4] != 'qw' and i[7].split('@')[1] in self._get_died_nodes):
					log.error('Died job: ' + i[0] + ', ' + ' '.join(i))
					log.error('Delete job: ' + i[0])
					os.popen('qdel ' + i[0])
		return runningJop

	@property
	def _get_died_nodes(self):
		died_nodes = []
		for node in os.popen("qhost").read().split('\n'):
			nodes = node.strip().split()
			if len(nodes) > 9 and '-' in nodes[6:9]:
				died_nodes.append(nodes)
		return died_nodes

	def _local(self):
		subids = []
		for i in range(len(self.unfinished_tasks)):
			newpid = os.fork() 
			if newpid == 0:
				os.popen('sh ' + self.unfinished_tasks[i] + ' > ' + self.unfinished_tasks[i] + '.o ' + '2> ' + self.unfinished_tasks[i] + '.e ')
				sys.exit(0)
			else:
				subids.append(newpid)
				Run.RUNNINGTASK['local'].append(newpid)
				log.info('Throw jobID:[' + str(newpid) + '] jobCmd:['  + self.unfinished_tasks[i] + '] in the local_cycle.')
				if i >= self.max_pa_jobs - 1:
					os.wait()
			time.sleep(0.5)

		while subids:
			subid = subids.pop()
			if psutil.pid_exists(subid) and psutil.Process(subid).ppid() == os.getpid():
				os.waitpid(subid, 0)
				time.sleep(0.5)
			Run.RUNNINGTASK['local'].remove(subid)

	def _getoption(self):
		cmd = self.sge_options if 'qsub' in self.sge_options else 'qsub ' + self.sge_options
		
		if '-S' not in self.sge_options:
			cmd += ' -S {bash}'
		if self.vf and ('virtual_free' not in self.sge_options and 'vf' not in self.sge_options):
			cmd += ' -l virtual_free={vf}'
		if '-pe' not in self.sge_options:
			cmd + ' -pe smp {cpu}'
		return cmd.format(cpu=self.cpu, vf=self.vf, bash=self.bash)
