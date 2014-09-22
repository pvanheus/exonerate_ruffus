#!/usr/bin/env python

import sys
import os
import os.path
import argparse
import glob
from subprocess import Popen, PIPE
import shlex
import re
import logging

from Bio import SeqIO
from ruffus import *
from ruffus.drmaa_wrapper import run_job, error_drmaa_job

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
parser.add_argument('--num_jobs', '-J', type=int, default=1, help='Number of exonerate jobs to run')
parser.add_argument('--run_local', '-L', action='store_true', default=False, help='Force to run jobs locally, not on cluster')
parser.add_argument('--scripts_dir', '-S', default='/cip0/research/pvh/software/sanbi-scripts', help='Location of pvh sanbi-scripts')
parser.add_argument('--modules_home', '-M', default='/cip0/software/x86_64/modules/Modules/3.2.9', help='Value of Environment Modules MODULESHOME env variable')
parser.add_argument('--query_type', default='dna', choices=['dna','protein'], help='exonerate query type (dna or protein)')
parser.add_argument('--debug', action='store_true', default=False, help='Set log level to DEBUG')
parser.add_argument('--printout', action='store_true', default=False, help='Print what will be run, but do not run anything')
parser.add_argument('--load_database', action='store_true', default=False, help='Load results into MySQL database')
parser.add_argument('--db_prefix', default='', help='Prefix to prepend to database names (e.g. project name)')
parser.add_argument('--db_hostname', '-H', default='localhost', help='Hostname of database to use for bp_seqfeature_load')
parser.add_argument('--db_username', '-U', required=True, help='Username to log into database for bp_seqfeature_load')
parser.add_argument('--db_password', '-P', required=True, help='Password to log into database for bp_seqfeature_load')
parser.add_argument('--exonerate_upstream', default=3000, type=int)
parser.add_argument('--exonerate_downstream', default=3000, type=int)
parser.add_argument('--extra_exonerate_args', default=None, type=str)
parser.add_argument('genome_filename', help='FASTA format genome file, must end in .fa or .fasta')
parser.add_argument('query_filename', help='FASTA format file of query sequences, must end in .fa or .fasta')
args = parser.parse_args()

if args.debug:
	logging.basicConfig(level=logging.DEBUG)
else:
	logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('genblastA_process_ruffus')

sys.path.append(args.scripts_dir) # this is where scripts like genblastA_to_gff are stored
from genblastA_to_gff.genblastA_to_gff3 import genblastA_process

drmaa_session = None
if not args.run_local:
	try:
		import drmaa
		drmaa_session = drmaa.Session()
		drmaa_session.initialize()
	except ImportError:
		sys.stderr.write('Failed to import drmaa module, must run locally')
		args.run_local = True


FASTA_RE=r'\.(fa|fasta)$'

os.chdir(args.working_directory)
starting_files = glob.glob('*' + args.input_pattern)

def safe_open(filename, mode='r'):
	try:
		file_obj = open(filename, mode)
	except IOError as e:
		logger.error('Failed to open {}: {}'.format(filename, str(e)))
		sys.exit(1)
	return file_obj

def merge_gff(infiles, output_file):
	out_file = safe_open(output_file,'w')
	out_file.write('##gff-version 3\n')
	for input_filename in infiles:
		in_file = safe_open(input_filename)
		first_line = in_file.readline() # discard 
		if first_line == '':
			# we've got an empty file here
			continue
		assert first_line.startswith('##gff-version'), "Invalid GFF header in {}\n".format(input_filename)
		line = in_file.readline()
		while line.startswith('#'):
			line = in_file.readline()
		# we've now skipped over all the header lines in the input GFF3, so line holds a non-header line
		out_file.write(line)
		out_file.write(in_file.read()) # write the rest of the input file to the output
		in_file.close()
	out_file.close()

def bp_seqfeature_load(input_file, output_file, dbname, hostname, username, password):
	seqfeature_load_cmdline = 'bp_seqfeature_load -f -d "dbi:mysql:database={};host={}" -c -i -u {} -p {} {}'.format(
		                       dbname, hostname, username, password, input_file)
	seqfeature_cmd = shlex.split(seqfeature_load_cmdline)
	proc = Popen(seqfeature_cmd, stdout=PIPE, stderr=PIPE)
	(output, error) = proc.communicate()
	retcode = proc.wait()
	if retcode != 0:
		raise Exception('bp_seqfeature_load failed: error: {} output: {}'.format(error, output))
	else:
		out_file = safe_open(output_file, 'w')
		out_file.write('output:\n' + output + '\nerror:\n' + error)
		out_file.close()

@transform(starting_files, 
			suffix(args.input_pattern),
			'.genblastA.gff3')
def genblastA_to_gff3(input_file, output_file):
	in_file = safe_open(input_file)
	out_file = safe_open(output_file, 'w')
	genblastA_process(in_file, out_file, min_perc_coverage=80.0)

@merge(genblastA_to_gff3, 'genblastA.all.gff3')
def merge_genblastA_gff3(infiles, output_file):
	merge_gff(infiles, output_file)

@transform(merge_genblastA_gff3,
	       suffix('.gff3'),
	       '.log', args.db_prefix+'genblastA', args.db_hostname,
	       args.db_username, args.db_password)
def load_genblastA_db(input_file, output_file, dbname, hostname, username, password):
	bp_seqfeature_load(input_file, output_file, dbname, hostname, username, password)

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
		raise Exception('Failed to run faToTwoBit (cmdline: {}): output:\n{}\nerror:{}\n'.format(cmdline, output, error))

@transform(args.query_filename,
	       regex(FASTA_RE),
		   '.idx')
def make_index(input_file, output_file):
	SeqIO.index_db(output_file, input_file, 'fasta')

FASTA_RE_COMPILED = re.compile(FASTA_RE)
PATH_val = os.environ['PATH'] + ':' + os.path.join(args.scripts_dir,'run_est_mapping')
PYTHONPATH_val = os.path.join(args.scripts_dir, 'lib')
@follows(make_index, make_twobit)
@transform(genblastA_to_gff3,
	       suffix('.genblastA.gff3'),
	       '.exonerate.gff3',
	       args.genome_filename, args.query_filename)
def run_exonerate(input_file, output_file, genome_filename, query_filename):
	twobit_filename = FASTA_RE_COMPILED.sub('.2bit', genome_filename)
	job_name = input_file.replace('.genblastA.gff3', '.sge')
	job = 'run_est_mapping.py --query_type {}'.format(args.query_type)
	job += ' --upstream {} --downstream {} --mapper exonerate --save_mapper_output --augustus_hints'.format(
		   args.exonerate_upstream, args.exonerate_downstream)
	if args.extra_exonerate_args:
		job += '--extra_mapper_args "{}"'.format(args.extra_exonerate_args)
	job += ' {} {} {} {}'.format(query_filename, input_file, twobit_filename, output_file)
	job_queue = 'all.q'
	job_env = dict(PATH=PATH_val, PYTHONPATH=PYTHONPATH_val)
	if not args.run_local:
		job_env['MODULESHOME'] = args.modules_home
	run_job(job, job_name=job_name, job_other_options='-q {}'.format(job_queue),
		    job_environment=job_env, drmaa_session=drmaa_session, 
		    working_directory=args.working_directory,
		    run_locally=args.run_local, logger=logger)

@merge(run_exonerate, 'exonerate.all.gff3')
def merge_exonerate_gff3(infiles, output_file):
	merge_gff(infiles, output_file)

@transform(merge_exonerate_gff3,
	       suffix('.gff3'),
	       '.log', args.db_prefix+'exonerate', args.db_hostname,
	       args.db_username, args.db_password)
def load_exonerate_db(input_file, output_file, dbname, hostname, username, password):
	bp_seqfeature_load(input_file, output_file, dbname, hostname, username, password)

if args.printout:
	pipeline_printout()
	pipeline_printout_graph('exonerate_ruffus.jpg', output_format='jpg', pipeline_name='Exonerate')
else:
	pipeline_run([merge_genblastA_gff3], multiprocess=args.num_threads)
	if args.load_database:
		pipeline_run([load_genblastA_db])
	# can't use multithread with a DRMAA task, causes script to hang
	pipeline_run([merge_exonerate_gff3], multithread=args.num_jobs)
	if args.load_database:
		pipeline_run([load_exonerate_db])

if not args.run_local:
	drmaa_session.exit()
