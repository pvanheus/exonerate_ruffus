#!/usr/bin/env python

import sys
import os
import os.path
import argparse
import glob
from subprocess import Popen, PIPE
import shlex

from Bio import SeqIO
from ruffus import *
from genblastA_to_gff3 import genblastA_process
# import Environment Modules stuff
MODULES_KEY = 'MODULESHOME'
if MODULES_KEY in os.environ:
	modules_init = os.path.join(os.environ[MODULES_KEY], 'init/python.py')
	execfile(modules_init)
	# need this for faToTwoBit
	module('load', 'blat/default')

parser = argparse.ArgumentParser(description='Use Ruffus to process .out files from genblastA')
parser.add_argument('--input_pattern', '-I', default='.out')
parser.add_argument('--working_directory', '-W', default='.')
parser.add_argument('--num_threads', '-N', type=int, default=1)
parser.add_argument('genome_filename', help='FASTA format genome file, must end in .fa or .fasta')
parser.add_argument('query_filename', help='FASTA format file of query sequences, must end in .fa or .fasta')
args = parser.parse_args()

FASTA_RE=r'\.(fa|fasta)$'

os.chdir(args.working_directory)
starting_files = glob.glob('*' + args.input_pattern)

def safe_open(filename, mode='r'):
	try:
		file_obj = open(filename, mode)
	except IOError as e:
		sys.stderr.write('Failed to open {}: {}'.format(filename, str(e)))
		sys.exit(1)
	return file_obj

@transform(starting_files, 
			suffix(args.input_pattern),
			'.genblastA.gff3')
def genblastA_to_gff3(input_file, output_file):
	in_file = safe_open(input_file)
	out_file = safe_open(output_file, 'w')
	genblastA_process(in_file, out_file, min_perc_coverage=80.0)

@transform(args.genome_filename,
	       regex(FASTA_RE),
	       '.2bit')
def make_twobit(input_file, output_file):
	cmdline = 'faToTwoBit -noMask {} {}'.format(input_file, output_file)
	cmdline_list = shlex.split(cmdline)
	proc = Popen(cmdline_list, stdout=PIPE, stderr=PIPE)
	(output, error) = proc.communicate()
	retcode = proc.wait()
	if (retcode != 0):
		sys.stderr.write('Failed to run faToTwoBit (cmdline: {}): output:\n{}\nerror:{}\n'.format(cmdline, output, error))

@transform(args.query_filename,
	       regex(FASTA_RE),
		   '.idx')
def make_index(input_file, output_file):
	SeqIO.index_db(output_file, input_file, 'fasta')

pipeline_run(multiprocess=args.num_threads)
