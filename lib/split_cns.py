#!/usr/bin/env python
from __future__ import print_function

import os, sys, re
import gzip
import argparse
from kit import *

log = plog()

class HelpFormatter(
        argparse.RawDescriptionHelpFormatter,
        argparse.ArgumentDefaultsHelpFormatter):
    pass

class File(object):

    def __init__(self, path):
        self.file = path
        self.name = os.path.basename(path)
        self.dir = os.path.dirname(path)
        self.type = self._set_type()
        self.handle = None
        self._check_exist()

    def _set_type(self):
        if self.name.endswith(('.gz', '.z', '.gzip', '.GZ', '.Z', '.GZIP')):
            return 'gz'
        else:
            return 'txt'

    @property
    def spath(self):
        names = self.file.split('.')
        if self.type == 'gz':
            names = names[:-1]
        if names[-1].lower() in ['fa', 'fasta', 'fq', 'fastq']:
            names = names[:-1]
        return '.'.join(names)

    def open_file(self):
        self.handle = open(self.file) if self.type == 'txt' else gzip.GzipFile(self.file)

    def close_file(self):
        self.handle.close()

    def _check_exist(self):
        if not os.path.isfile(self.file):
            log.error('Error, ' + self.file + ' is not exists.')
            sys.exit(1)

class Fastx(File):
    
    def __init__(self, path, fofn = False):
        self.fofn = fofn
        super(Fastx, self).__init__(path)

    def _read(self):
        name = seq = ''
        self.open_file()
        for line in self.handle:
            if self.type == 'gz':
                line = byte2str(line)
            line = line.strip()
            if line.startswith('>') or line.startswith('@'):
                if name and seq:
                    yield (name, seq)
                name = line.split()[0][1:]
                seq = ''
            elif line == '+':
                next(self.handle) #skip quality score line
                continue
            else:
                seq += line
        yield (name, seq)
        self.close_file()

    def read(self):
        if not self.fofn:
            for r in self._read():
                yield r
        else:
            self.open_file()
            for line in self.handle:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if self.dir and not line.startswith('/'):
                    line = self.dir + '/' + line
                for r in Fastx(line)._read():
                    yield r
            self.close_file()

    def cutf(self, count, pdir, rn = 0, ml = 0, index = True):
        fa_files = []
        idx_files = []
        path = pdir + '/cns{:0>{}}.fasta'
        for i in range(count):
            ppath = path.format(i, len(str(count)))
            fa_files.append(open(ppath, 'w'))
            if index:
                idx_files.append(open(ppath + '.idx', 'w'))
        i = 0
        t = 0
        for name, seq in self.read():
            i = i + 1 if i + 1 < count else 0
            lens = len(seq)
            if lens < ml:
                continue
            if rn:
                t += 1
                print('>%d %d %f pid=%s\n%s' % (t, lens, 1, name, seq), file=fa_files[i])
                name = str(t)
            else:
                print('>%s %d %f\n%s' % (name, lens, 1, seq), file=fa_files[i])
            if index:
                cur_offset = fa_files[i].tell()
                print('%s\t%d\t%d' % (name, cur_offset - lens - 1, lens), file=idx_files[i])

        for i in range(count):
            fa_files[i].close()
            if index:
                idx_files[i].close()
        log.info('split %d reads to %d subfiles', t, count)
        return fa_files

def main(args):
    log.info('Split step options:')
    log.info(args)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    f = Fastx(args.fofn, fofn = True)
    f.cutf(args.count, rn = args.rename, ml = args.min_len, pdir = args.outdir, index = args.index)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
        description='''
%(prog)s:
    split reads with a specified number of subfile in total.
exmples:
    %(prog)s -f input.fofn -c 24 -o out_dir
'''
    )
    parser.add_argument('-f', '--fofn', metavar='FILE', required=True,
        help='set the reads file.')
    parser.add_argument('-o', '--outdir', type=str, default='./',
        help='set the output directory.')
    parser.add_argument('-c', '--count', type=int, default=10,
        help='set the number of subfile in total.')
    parser.add_argument('-l', '--min_len', type=str, default='1kb',
        help='minimum length of a read.')
    parser.add_argument('-i', '--index', action='store_false', default=True,
        help='produce the index file.')
    parser.add_argument('-r', '--rename', action='store_false', default=True,
        help='rename seq ID with integer.')
    args = parser.parse_args()
    args.min_len = parse_num_unit(args.min_len)
    main(args)
